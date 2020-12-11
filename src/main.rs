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
extern crate clap;
extern crate env_logger;
#[macro_use]
extern crate lazy_static;
extern crate log;
extern crate regex;
extern crate serde;
extern crate serde_json;

use std::error::Error;

use bam::bam_reader::{ModificationTime, Region};
use bam::header::Header as BamHeader;
use bam::record::AlignmentEntry;
use bam::record::Record as BamRecord;
use bam::IndexedReader as BamReader;
use clap::Clap;
use serde::Serialize;

mod error;
mod seq;
mod variant;

use crate::error::err;
use crate::seq::{Base, Ordering};
use crate::variant::Variant;

trait MakeRegion {
    fn make_region(&self, header: &BamHeader) -> Result<Region, Box<dyn Error>>;
}

impl MakeRegion for Variant {
    fn make_region(&self, header: &BamHeader) -> Result<Region, Box<dyn Error>> {
        let rid = header.reference_id(self.chrom()).ok_or_else(err)?;
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
        self.alt_count() as f32 / self.total_count() as f32
    }

    fn proper_freq(&self) -> f32 {
        self.proper as f32 / self.total_count() as f32
    }

    fn ref_count(&self) -> u32 {
        self.reference
    }

    fn ref_freq(&self) -> f32 {
        self.ref_count() as f32 / self.total_count() as f32
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
        let mut iter = if let Ok(v) = record.alignment_entries() {
            v.skip_while(|i| i.ref_pos() < Some(var.pos() - 1))
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

            if preskip && var.is_abbr_deletion() {
                preskip = false;
                log::warn!("Skipping first base due to variant deletion format like `1:12345C>-`");
                continue;
            };

            if curr.is_insertion() {
                ralt.push(Base::from_byte(curr.record_nt().ok_or_else(err)?)?)
            } else if curr.is_deletion() {
                rref.push(Base::from_byte(curr.ref_nt().ok_or_else(err)?)?)
            } else {
                ralt.push(Base::from_byte(curr.record_nt().ok_or_else(err)?)?);
                rref.push(Base::from_byte(curr.ref_nt().ok_or_else(err)?)?)
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

        match (var.ref_cmp(&rref), var.alt_cmp(&ralt), rref == ralt) {
            // Record ref does not accord with variant ref.
            (Ordering::Nul, _, _) => {
                log::error!(
                    "Bam record `{:?}` ref {:?} does not accord with variant ref {:?}.",
                    String::from_utf8(record.name().to_vec()),
                    rref,
                    var.refs()
                );
                self.alleles += 1;
                Ok(())
            }
            // Fully supported Alt
            (Ordering::Equ, Ordering::Equ, _) => {
                self.proper += 1;
                Ok(())
            }
            // Fully supported Ref
            (Ordering::Equ, _, true) => {
                self.reference += 1;
                Ok(())
            }
            // Excessively supported ref
            // FIXME: Extra base considered the same with genome reference
            (Ordering::Sub, _, true) => {
                self.reference += 1;
                Ok(())
            }
            // Partially supported Ref
            (_, _, true) => {
                self.reference += 1;
                Ok(())
            }
            // Partially supported Alt
            (Ordering::Sub, Ordering::Equ, false) => {
                self.excessive += 1;
                Ok(())
            }
            // Excessively supported Alt
            (_, Ordering::Sub, false) => {
                self.excessive += 1;
                Ok(())
            }
            // Partially supported Alt
            (_, Ordering::Sup, false) => {
                self.alleles += 1;
                Ok(())
            }
            _ => {
                self.alleles += 1;
                Ok(())
            }
        }
    }
}

#[derive(Clap)]
struct Opts {
    #[clap(long, about = "Input bam file.")]
    bam: String,
    #[clap(long, about = "Input genome variant, e.g. 'chr1:12345AT>-'.")]
    var: String,
    #[clap(long, default_value = "0", about = "Minimum read mapping quality.")]
    mapq: u8,
    #[clap(
        long,
        default_value = "10",
        about = "Minimum margin base distance for alt support. Margin stands for read start/end, softclip start/end etc."
    )]
    margin: u32,
}

fn main() -> Result<(), Box<dyn Error>> {
    let opts = Opts::parse();

    env_logger::Builder::new()
        .filter_level(log::LevelFilter::Info)
        .init();

    log::warn!("Reading bam file {}.", &opts.bam);
    let mut sam = BamReader::build()
        .modification_time(ModificationTime::warn(|e| eprintln!("{}", e)))
        .from_path(&opts.bam)?;

    let mut sum = Summary::default();
    let var = Variant::try_parse(&opts.var)?;
    log::warn!("Parsing variant as {:?}", var);

    log::warn!("Fetching variant adjcent reads.");
    let reg = var.make_region(sam.header())?;
    for i in sam.fetch(&reg)? {
        let record = i?;
        if ((record.start() + 1) as u32 > var.pos())
            || ((record.calculate_end() as u32) < var.end())
        {
            break;
        };
        match sum.validate(&record, &var) {
            Ok(_) => {}
            Err(e) => {
                log::error!("{}", e)
            }
        }
    }

    println!("{}", serde_json::to_string_pretty(&sum)?);
    Ok(())
}
