use std::error::Error;

use bam::record::AlignmentEntry;
use bam::record::Record as BamRecord;

use crate::error::err;
use crate::seq::{Base, Ordering, Support};
use crate::variant::Variant;

/// Validate alignment supportion for variant.
pub trait VariantValidate {
    fn validate(&self, v: &Variant) -> Result<Support, Box<dyn Error>>;
}

impl VariantValidate for BamRecord {
    /// Validate record supportion for variant.
    ///
    /// ## Examples
    ///
    /// ```rust
    /// use bam::record::Record as BamRecord;
    ///
    /// use crate::variant::Variant;
    ///
    /// let var = Variant::try_parse("chr1:123456AT>-")?;
    /// let record = BamRecord::new();
    ///
    /// // should raise error
    /// record.validate(var)?;
    /// ```
    ///
    /// ## Warn
    ///
    /// Crate `bam` bam reader reading alignemnt with 0-based position, while variant is 1-based.
    /// So alignment `+1` or variant `-1` is necessary in some places.
    ///
    fn validate(&self, var: &Variant) -> Result<Support, Box<dyn Error>> {
        // Unmapped read
        if (!self.flag().is_mapped())
            || (self.start() + 1) as u32 > var.end()
            || (self.calculate_end() as u32) < var.pos()
        {
            return Ok(Support::Nul);
        }
        // Record ref
        let mut rref: Vec<Base> = Vec::with_capacity(var.refs().len());
        // Record alt
        let mut ralt: Vec<Base> = Vec::with_capacity(var.alts().len());
        let mut iter = if let Ok(mut v) = self.alignment_entries() {
            v.skip_while(|i| i.ref_pos() < Some(var.pos() - 1))
        } else {
            return Ok(Support::Unk);
        };
        let mut next: Option<AlignmentEntry> = if let Some(v) = iter.next() {
            Some(v)
        } else {
            return Ok(Support::Nul);
        };

        let mut preskip = true;

        loop {
            let curr = match next {
                Some(v) => v,
                None => break,
            };
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

            match next {
                Some(ref v) => {
                    if !v.is_seq_match() {
                        continue;
                    }
                }
                _ => {}
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
                    String::from_utf8(self.name().to_vec()),
                    rref,
                    var.refs()
                );
                Ok(Support::Oth)
            }
            // Fully supported Alt
            (Ordering::Equ, Ordering::Equ, _) => Ok(Support::Alt),
            // Fully supported Ref
            (Ordering::Equ, _, true) => Ok(Support::Ref),
            // Excessively supported ref
            // FIXME: Extra base considered the same with genome reference
            (Ordering::Sub, _, true) => Ok(Support::Ree),
            // Partially supported Ref
            (_, _, true) => Ok(Support::Rep),
            // Partially supported Alt
            (Ordering::Sub, Ordering::Equ, false) => Ok(Support::Ale),
            // Excessively supported Alt
            (_, Ordering::Sub, false) => Ok(Support::Ale),
            // Partially supported Alt
            (_, Ordering::Sup, false) => Ok(Support::Alp),
            _ => Ok(Support::Oth),
        }
    }
}
