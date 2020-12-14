use std::fmt;

use once_cell::sync::Lazy;
use regex::Regex;

use crate::error::err;
use crate::seq::{Base, Ordering};

static VAREX: Lazy<Regex> = Lazy::new(|| {
    Regex::new(r"(?i)^(?P<chrom>(?:chr|)[\w\.-]+):(?P<pos>\d+)(?P<refs>(?:[ATCGN]+|-))>(?P<alts>(?:[ATCGN]+|-))$").unwrap()
});

#[derive(PartialEq, Debug)]
pub struct Variant {
    chrom: String,
    pos: u32,
    refs: Vec<Base>,
    alts: Vec<Base>,
}

impl Variant {
    /// # Parse variant from a string.
    ///
    /// ## Format
    ///
    /// > chr1:12345AT>G
    ///
    /// ## Examples
    ///
    /// ```rust
    /// Variant::try_parse("1:12345AT>GC")?;
    /// ```
    pub fn try_parse(v: &str) -> Result<Self, Box<dyn std::error::Error>> {
        if let Some(c) = VAREX.captures(v) {
            Ok(Self {
                chrom: String::from(c.name("chrom").ok_or_else(err)?.as_str()),
                pos: c.name("pos").ok_or_else(err)?.as_str().parse()?,
                refs: Base::try_parse(c.name("refs").ok_or_else(err)?.as_str())?,
                alts: Base::try_parse(c.name("alts").ok_or_else(err)?.as_str())?,
            })
        } else {
            Err(Box::new(err()))
        }
    }

    /// Reference chromosome.
    pub fn chrom(&self) -> &String {
        &self.chrom
    }

    /// Reference start position.
    pub fn pos(&self) -> u32 {
        self.pos
    }

    /// Reference end
    ///
    /// ## Warn
    ///
    /// Maximum variant length: u32::Max.
    pub fn end(&self) -> u32 {
        if !self.refs.is_empty() {
            self.pos + self.refs.len() as u32 - 1
        } else {
            self.pos + 1
        }
    }

    /// Reference sequence.
    pub fn refs(&self) -> &Vec<Base> {
        &self.refs
    }

    /// Alternative sequence.
    pub fn alts(&self) -> &Vec<Base> {
        &self.alts
    }

    /// Reference sequence as string
    pub fn ref_str(&self) -> String {
        if !self.refs.is_empty() {
            self.refs
                .iter()
                .map(|i| i.stringify())
                .collect::<Vec<String>>()
                .concat()
        } else {
            String::from("-")
        }
    }

    /// Alternative sequence as string
    pub fn alt_str(&self) -> String {
        if !self.alts.is_empty() {
            self.alts
                .iter()
                .map(|i| i.stringify())
                .collect::<Vec<String>>()
                .concat()
        } else {
            String::from("-")
        }
    }

    pub fn is_abbr_deletion(&self) -> bool {
        self.alts.is_empty()
    }

    pub fn ref_cmp(&self, v: &[Base]) -> Ordering {
        if self.refs == v {
            Ordering::Equ
        } else if self.refs.is_empty() || v.is_empty() {
            Ordering::Emp
        } else if self.refs.starts_with(v) {
            Ordering::Sup
        } else if v.starts_with(&self.refs) {
            Ordering::Sub
        } else {
            Ordering::Nul
        }
    }

    pub fn alt_cmp(&self, v: &[Base]) -> Ordering {
        if self.alts == v {
            Ordering::Equ
        } else if self.alts.is_empty() || v.is_empty() {
            Ordering::Emp
        } else if self.alts.starts_with(v) {
            Ordering::Sup
        } else if v.starts_with(&self.alts) {
            Ordering::Sub
        } else {
            Ordering::Nul
        }
    }
}

impl fmt::Display for Variant {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}:{}{}>{}",
            self.chrom,
            self.pos,
            self.ref_str(),
            self.alt_str()
        )
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_try_parse() {
        assert_eq!(
            Variant::try_parse("chr1:12345AT>GC").unwrap(),
            Variant {
                chrom: String::from("chr1"),
                pos: 12345,
                refs: vec![Base::A, Base::T],
                alts: vec![Base::G, Base::C],
            }
        )
    }
}
