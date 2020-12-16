//! SNP/InDel support stats based on CIGAR and MD tag.
//!
//! ## Getting started
//!
//! ### Install
//!
//! Clone repo and install with cargo.
//!
//! ```shell
//! cd validate-variant-cli
//! cargo install --path .
//! ```
//!
//! ### Run
//!
//! ```shell
//! vav \
//!     --bam tests/W080446T.many_variants.bam \
//!     --var "2:29474101C>A"
//! ```
extern crate bam;
#[macro_use]
extern crate clap;
extern crate env_logger;
extern crate log;
extern crate once_cell;
extern crate regex;
extern crate serde;
extern crate serde_json;

use std::collections::HashMap;
use std::error::Error;

use bam::bam_reader::{ModificationTime, Region};
use bam::header::Header as BamHeader;
use bam::record::AlignmentEntry;
use bam::record::Record as BamRecord;
use bam::IndexedReader as BamReader;
use clap::Clap;
use once_cell::sync::OnceCell;
use serde::Serialize;

mod error;
mod seq;
mod variant;

use crate::error::opterr;
use crate::seq::{Base, Ordering};
use crate::variant::Variant;

static MAPQ: OnceCell<u8> = OnceCell::new();
static MARGIN: OnceCell<u32> = OnceCell::new();

trait MakeRegion {
    fn make_region(&self, header: &BamHeader) -> Result<Region, Box<dyn Error>>;
}

impl MakeRegion for Variant {
    fn make_region(&self, header: &BamHeader) -> Result<Region, Box<dyn Error>> {
        let rid = header.reference_id(self.chrom()).ok_or_else(opterr)?;
        Ok(Region::new(rid, self.pos(), self.end()))
    }
}

#[derive(Serialize, Default)]
struct Summary {
    /// Ref support.
    reference: u32,
    /// Alt proper support.
    proper: u32,
    /// Alt support in margin.
    margin: u32,
    /// Alt support of low mapq.
    lowq: u32,
    /// Alt support of excessive support.
    /// Example, expecting chr1:12345A>C, got chr1:12345AT>CG, chr1:12345AT>C, etc.
    excessive: u32,
    /// Other alleles support.
    alleles: u32,
    /// Unknown support or exception.
    unknown: u32,
}

impl Summary {
    fn total_count(&self) -> u32 {
        self.reference
            + self.proper
            + self.margin
            + self.lowq
            + self.excessive
            + self.alleles
            + self.unknown
    }

    fn alt_count(&self) -> u32 {
        self.proper + self.margin + self.lowq + self.excessive
    }

    fn alt_freq(&self) -> f32 {
        let v = self.alt_count() as f32 / self.total_count() as f32;
        (v * 10000.0).round() / 10000.0
    }

    fn proper_freq(&self) -> f32 {
        let v = self.proper as f32 / self.total_count() as f32;
        (v * 10000.0).round() / 10000.0
    }

    fn margin_freq(&self) -> f32 {
        let v = self.margin as f32 / self.total_count() as f32;
        (v * 10000.0).round() / 10000.0
    }

    fn lowq_freq(&self) -> f32 {
        let v = self.lowq as f32 / self.total_count() as f32;
        (v * 10000.0).round() / 10000.0
    }

    fn ref_count(&self) -> u32 {
        self.reference
    }

    fn ref_freq(&self) -> f32 {
        let v = self.ref_count() as f32 / self.total_count() as f32;
        (v * 10000.0).round() / 10000.0
    }

    /// Validate record supportion for variant.
    ///
    /// ## Examples
    ///
    /// ```rust
    /// use bam::record::Record as BamRecord;
    ///
    /// use crate::variant::Variant;
    /// use crate::Summary;
    ///
    /// let var = Variant::try_parse("chr1:123456AT>-")?;
    /// let record = BamRecord::new();
    /// let sum = Summary::default();
    /// sum.valdiate(&record, &var)?;
    /// ```
    ///
    /// ## Warn
    ///
    /// Crate `bam` bam reader reading alignemnt with 0-based position, while variant is 1-based.
    /// So alignment `+1` or variant `-1` is necessary in some places.
    ///
    fn validate(&mut self, record: &BamRecord, var: &Variant) -> Result<(), Box<dyn Error>> {
        // Unmapped read
        if (!record.flag().is_mapped())
            || (record.start() + 1) as u32 > var.end()
            || (record.calculate_end() as u32) < var.pos()
        {
            return Ok(());
        }
        // Record ref
        let mut rref: Vec<Base> = Vec::with_capacity(var.refs().len());
        // Record alt
        let mut ralt: Vec<Base> = Vec::with_capacity(var.alts().len());
        // Front margin and end margin
        let mut front = 0;
        let mut end = 0;
        let mut iter = if let Ok(v) = record.alignment_entries() {
            v.skip_while(|i| {
                front += 1;
                i.ref_pos() < Some(var.pos() - 1)
            })
        } else {
            self.unknown += 1;
            return Ok(());
        };
        let mut next: Option<AlignmentEntry> = if let Some(v) = iter.next() {
            Some(v)
        } else {
            return Ok(());
        };

        let mut preskip = true;

        while let Some(curr) = next {
            next = iter.next();
            if let Some(ref v) = curr.record_pos() {
                end = *v;
            };

            if preskip && var.is_abbr_deletion() {
                preskip = false;
                log::info!("Skipping first base due to variant deletion format like `1:12345C>-`");
                continue;
            };

            if curr.is_insertion() {
                ralt.push(Base::from_byte(curr.record_nt().ok_or_else(opterr)?)?)
            } else if curr.is_deletion() {
                rref.push(Base::from_byte(curr.ref_nt().ok_or_else(opterr)?)?)
            } else {
                ralt.push(Base::from_byte(curr.record_nt().ok_or_else(opterr)?)?);
                rref.push(Base::from_byte(curr.ref_nt().ok_or_else(opterr)?)?)
            };

            if let Some(ref v) = next {
                if !v.is_seq_match() {
                    continue;
                }
            };

            if rref.len() >= var.refs().len() || ralt.len() >= var.alts().len() {
                break;
            }
        }
        end = record.aligned_query_end() - end;

        match (var.ref_cmp(&rref), var.alt_cmp(&ralt), rref == ralt) {
            // Record ref does not accord with variant ref.
            (Ordering::Nul, _, _) => {
                log::error!(
                    "Bam record `{}` ref {:?} does not accord with variant ref {:?}.",
                    String::from_utf8_lossy(&record.name().to_vec()),
                    rref,
                    var.refs()
                );
                self.alleles += 1;
                Ok(())
            }
            // Fully supported Alt
            (Ordering::Equ, Ordering::Equ, _) => {
                log::debug!(
                    "Fully supported alt by record `{}`",
                    String::from_utf8_lossy(&record.name().to_vec())
                );
                if Some(&record.mapq()) < MAPQ.get() {
                    self.lowq += 1;
                } else if Some(&front) < MARGIN.get() || Some(&end) < MARGIN.get() {
                    self.margin += 1;
                } else {
                    self.proper += 1;
                }
                Ok(())
            }
            // Fully supported Ref
            (Ordering::Equ, _, true) => {
                log::debug!(
                    "Fully supported ref by record `{}`",
                    String::from_utf8_lossy(&record.name().to_vec())
                );
                self.reference += 1;
                Ok(())
            }
            // Excessively supported ref
            // FIXME: Extra base considered the same with genome reference
            (Ordering::Sub, _, true) => {
                log::debug!(
                    "Excessively supported ref by record `{}`",
                    String::from_utf8_lossy(&record.name().to_vec())
                );
                self.reference += 1;
                Ok(())
            }
            // Partially supported Ref
            (_, _, true) => {
                log::debug!(
                    "Partially supported ref by record `{}`",
                    String::from_utf8_lossy(&record.name().to_vec())
                );
                self.reference += 1;
                Ok(())
            }
            // Partially supported Alt
            (Ordering::Sub, Ordering::Equ, false) => {
                log::debug!(
                    "Partially supported alt by record `{}`",
                    String::from_utf8_lossy(&record.name().to_vec())
                );
                self.excessive += 1;
                Ok(())
            }
            // Excessively supported Alt
            (_, Ordering::Sub, false) => {
                log::debug!(
                    "Excessively supported alt by record `{}`",
                    String::from_utf8_lossy(&record.name().to_vec())
                );
                self.excessive += 1;
                Ok(())
            }
            // Partially supported Alt
            (_, Ordering::Sup, false) => {
                log::debug!(
                    "Partially supported alt (interpreted as other allele) by record `{}`",
                    String::from_utf8_lossy(&record.name().to_vec())
                );
                self.alleles += 1;
                Ok(())
            }
            _ => {
                log::debug!(
                    "Other allele by record `{}`",
                    String::from_utf8_lossy(&record.name().to_vec())
                );
                self.alleles += 1;
                Ok(())
            }
        }
    }
}

#[derive(Clap)]
#[clap(name = crate_name!(), version = crate_version!(), author = crate_authors!(), about = crate_description!())]
struct Opts {
    #[clap(
        long,
        number_of_values = 1,
        about = "Input genome variant, e.g. 'chr1:12345AT>-'."
    )]
    var: Vec<String>,
    #[clap(long, default_value = "30", about = "Minimum read mapping quality.")]
    mapq: u8,
    #[clap(
        long,
        default_value = "10",
        about = "Minimum margin base distance for alt support. Margin stands for read start/end, softclip start/end etc."
    )]
    margin: u32,
    #[clap(short, long, about = "Print verbose info.")]
    verbose: bool,
    #[clap(about = "Input bam file.")]
    bam: String,
}

fn main() -> Result<(), Box<dyn Error>> {
    let mut opts = Opts::parse();
    MAPQ.set(opts.mapq).map_err(|_| opterr())?;
    MARGIN.set(opts.margin).map_err(|_| opterr())?;

    env_logger::Builder::new()
        .filter_level(if opts.verbose {
            log::LevelFilter::Debug
        } else {
            log::LevelFilter::Info
        })
        .init();

    log::warn!("Reading bam file {}.", &opts.bam);
    let mut sam = BamReader::build()
        .modification_time(ModificationTime::warn(|e| eprintln!("{}", e)))
        .from_path(&opts.bam)?;

    let mut varsum: HashMap<String, Summary> = HashMap::new();
    while let Some(each) = opts.var.pop() {
        if varsum.contains_key(&each) {
            continue;
        };
        let variant = Variant::try_parse(&each)?;
        let mut sum = Summary::default();
        log::warn!("Variant {} Parsed as {:?}", &each, variant);

        log::warn!("Fetching variant adjcent reads.");
        let reg = variant.make_region(sam.header())?;
        for i in sam.fetch(&reg)? {
            let record = i?;
            if ((record.start() + 1) as u32 > variant.pos())
                || ((record.calculate_end() as u32) < variant.end())
            {
                break;
            };
            match sum.validate(&record, &variant) {
                Ok(_) => {}
                Err(e) => {
                    log::error!("{}", e)
                }
            }
        }

        log::warn!(
            "Variant {} total {}; Ref {}({}); Proper alt {}({}); Margin alt {}({}); Lowq alt {}({})",
            &each,
            sum.total_count(),
            sum.reference,
            sum.ref_freq(),
            sum.proper,
            sum.proper_freq(),
            sum.margin,
            sum.margin_freq(),
            sum.lowq,
            sum.lowq_freq(),
        );
        varsum.insert(each, sum);
    }

    if varsum.len() == 1usize {
        println!(
            "{}",
            serde_json::to_string_pretty(&varsum.values().next().ok_or_else(opterr)?)?
        );
    } else {
        println!("{}", serde_json::to_string_pretty(&varsum)?);
    }
    Ok(())
}
