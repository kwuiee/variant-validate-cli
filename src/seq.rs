use std::io::Error as IOError;

use crate::error::err;

/// CIGAR Operations.
///
/// ## Warn
///
/// Not exactly the same with SAM specifications.
#[derive(PartialEq, Debug)]
pub enum CigarOp {
    M,
    I,
    D,
    X,
    S,
}

/// Base Types.
///
/// ## Warn
///
/// N for sequencer unknown base.
#[derive(PartialEq, Debug)]
pub enum Base {
    A,
    T,
    C,
    G,
    N,
}

impl<'a> Base {
    /// Parse base sequence from a str.
    ///
    /// ## Examples
    ///
    /// ```rust
    /// use crate::seq::Base;
    ///
    /// Base::try_parse("ATCC")?;
    /// ```
    /// ## Note
    ///
    /// `-` stands for null, otherwise a sequence of ATCGN is required.
    pub fn try_parse(v: &'a str) -> Result<Vec<Self>, IOError> {
        match v {
            "-" => Ok(vec![]),
            _ => {
                let mut r: Vec<Self> = Vec::with_capacity(v.len());
                v.bytes().try_for_each(|i| {
                    r.push(match i {
                        b'A' | b'a' => Base::A,
                        b'T' | b't' => Base::T,
                        b'C' | b'c' => Base::C,
                        b'G' | b'g' => Base::G,
                        b'N' | b'n' => Base::N,
                        _ => {
                            return Err(err(&format!("Error parsing `{}` as as Base sequence.", v)))
                        }
                    });
                    Ok(())
                })?;
                Ok(r)
            }
        }
    }

    pub fn stringify(&self) -> String {
        match self {
            Self::A => String::from("A"),
            Self::T => String::from("T"),
            Self::C => String::from("C"),
            Self::G => String::from("G"),
            Self::N => String::from("N"),
        }
    }

    pub fn from_byte(v: u8) -> Result<Self, IOError> {
        match v {
            b'A' | b'a' => Ok(Base::A),
            b'T' | b't' => Ok(Base::T),
            b'C' | b'c' => Ok(Base::C),
            b'G' | b'g' => Ok(Base::G),
            b'N' | b'n' => Ok(Base::N),
            _ => Err(err(&format!("Error parsing `{}` as valid Base", v))),
        }
    }
}

/// Position Query Base.
///
/// Including ref base, alt base, cigar, reference position and query position.
#[derive(PartialEq, Debug)]
pub struct QueryBase {
    // ref base
    r: Option<Base>,
    // alt base
    a: Option<Base>,
    cigar: CigarOp,
    refpos: u32,
    querypos: u32,
}

/// Sequence support enum.
///
/// ## Notes
///
/// - Ref: Expected reference sequence.
/// - rRef: Read/Record actual reference sequence.
/// - Alt: Expected alt.
/// - rAlt: Read/record actual alt sequence.
/// - End: Reach the end of read/record, read/record may not be long enough to cover Ref/Alt.
///
/// ## Implementation
///
/// Program shall be implemented in following orders:
///
/// ### Nul
///
/// > rRef != Ref || not rRef.startswith(Ref) || not Ref.startswith(rRef)
///
/// - 异常情况.
/// - 来源: 1. 参考基因组不同, 2. 给错ref序列, 3. 给错变异起始位点
/// - 条件: rRef不等于Ref 或 rRef不是Ref的子序列 或 Ref不是rRef的子序列
///
/// ### Alt
///
/// > (Ref == rRef && Alt == rAlt)
///
/// - 完全支持变异.
/// - 条件: Ref等于rRef 且 Alt等于rAlt
///
/// ### Ale
///
/// > rAlt.startswith(Alt)
///
/// - 结余支持变异.
/// - 条件: Alt等于rAlt 或 Alt是rAlt的子序列
///
/// ### Alp
///
/// > Alt.startswith(rAlt)
///
/// - 部分支持变异.
/// - 条件: rAlt等于Alt 或 rAlt是Alt的子序列
///
/// ### Oth:
///
/// > rRef != rAlt
///
/// - 其它支持.
/// - 条件: rRef 不等于 rAlt
///
/// ### Ref
///
/// > Ref == rRef
///
/// - 完全支持参考.
/// - 条件: Ref, rRef, rAlt三者相等
///
/// ### Ree
///
/// > rRef.startswith(Ref)
///
/// - 结余支持变异.
/// - 条件: Ref是rRef的开头子序列
///
/// ### Rep
///
/// > Ref.startswith(rRef)
///
/// - 部分支持变异.
/// - 条件: rRef是Ref的子序列
///
/// ### Oth
///
/// - 其它
///
/// ### Unk
///
/// Unknown
///
/// ### Nul
///
/// Support nothing or exceptions.
#[derive(PartialEq)]
pub enum Support {
    /// Ref fully, ref is fully supported
    Ref,
    /// partial ref, ref is partially supported, e.g. expecting "ATC" found "AT"
    Rep,
    /// Excessive ref, ref is excessively supported, e.g. expecting "ATC" found "ATCG"
    Ree,
    /// alt, alt is fully supported
    Alt,
    /// partial alt, alt is partially supported, e.g. expecting "ATC" found "AT"
    Alp,
    /// Excessive alt, alt is excessively supported, e.g. expecting "ATC" found "ATCG"
    Ale,
    /// Other alleles, a Rep or Alp may be considered as Oth?
    Oth,
    /// Unknown, e.g. when read MD tag missing
    Unk,
    /// Exception or support Null
    Nul,
}

impl Support {
    pub fn is_ref(&self) -> bool {
        matches!(self, Self::Ref)
    }

    pub fn may_ref(&self) -> bool {
        matches!(self, Self::Rep | Self::Ree)
    }

    pub fn any_ref(&self) -> bool {
        self.is_ref() | self.may_ref()
    }

    pub fn is_alt(&self) -> bool {
        matches!(self, Self::Alt)
    }

    pub fn may_alt(&self) -> bool {
        matches!(self, Self::Alp | Self::Ale)
    }

    pub fn any_alt(&self) -> bool {
        self.is_alt() | self.may_alt()
    }

    pub fn is_oth(&self) -> bool {
        matches!(self, Self::Oth)
    }

    pub fn is_nul(&self) -> bool {
        matches!(self, Self::Nul)
    }
}

/// Alignment sequence cmp.
#[derive(PartialEq, Debug)]
pub enum Ordering {
    /// One out of two sequence is empty.
    Emp,
    /// One is superset of the other.
    Sup,
    /// One equals the other.
    Equ,
    /// One is subset of the other.
    Sub,
    /// None above stands.
    Nul,
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn base_try_parse() {
        assert_eq!(
            Base::try_parse("ATC").unwrap(),
            vec![Base::A, Base::T, Base::C]
        )
    }
}
