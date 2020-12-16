# Variant Validate

Validate a SNP/InDel from bam using CIGAR and MD tag.

## Getting Started

```shell
$ vav --help
vav 0.1.1
slyo <sean.lyo@outlook.com>
Validate a SNP/InDel from bam using CIGAR and MD tag.

USAGE:
    vav [FLAGS] [OPTIONS] <bam>

ARGS:
    <bam>    Input bam file.

FLAGS:
    -h, --help       Prints help information
    -v, --verbose    Print verbose info.
    -V, --version    Prints version information

OPTIONS:
        --mapq <mapq>        Minimum read mapping quality. [default: 30]
        --margin <margin>    Minimum margin base distance for alt support. Margin stands for read
                             start/end, softclip start/end etc. [default: 10]
        --var <var>...       Input genome variant, e.g. 'chr1:12345AT>-'.
```

## Examples

With multiple variants.

```shell
$ vav tests/many_variants.bam --var '1:156843458A>G' --var "2:29474101C>A"
[2020-12-16T09:46:56Z WARN  vav] Reading bam file tests/many_variants.bam.
[2020-12-16T09:46:56Z WARN  vav] Variant 2:29474101C>A Parsed as Variant { chrom: "2", pos: 29474101, refs: [C], alts: [A] }
[2020-12-16T09:46:56Z WARN  vav] Fetching variant adjcent reads.
[2020-12-16T09:46:56Z WARN  vav] Variant 2:29474101C>A total 7724; Ref 7638(0.9889); Proper alt 21(0.0027); Margin alt 17(0.0022); Lowq alt 0(0)
[2020-12-16T09:46:56Z WARN  vav] Variant 1:156843458A>G Parsed as Variant { chrom: "1", pos: 156843458, refs: [A], alts: [G] }
[2020-12-16T09:46:56Z WARN  vav] Fetching variant adjcent reads.
[2020-12-16T09:46:56Z WARN  vav] Variant 1:156843458A>G total 6465; Ref 6389(0.9882); Proper alt 25(0.0039); Margin alt 4(0.0006); Lowq alt 0(0)
{
  "1:156843458A>G": {
    "reference": 6389,
    "proper": 25,
    "margin": 4,
    "lowq": 0,
    "excessive": 0,
    "alleles": 47,
    "unknown": 0
  },
  "2:29474101C>A": {
    "reference": 7638,
    "proper": 21,
    "margin": 17,
    "lowq": 0,
    "excessive": 0,
    "alleles": 48,
    "unknown": 0
  }
}
```

With single variant.

```shell
vav tests/many_variants.bam --var "2:29474101C>A"
[2020-12-16T09:51:07Z WARN  vav] Reading bam file tests/many_variants.bam.
[2020-12-16T09:51:07Z WARN  vav] Variant 2:29474101C>A Parsed as Variant { chrom: "2", pos: 29474101, refs: [C], alts: [A] }
[2020-12-16T09:51:07Z WARN  vav] Fetching variant adjcent reads.
[2020-12-16T09:51:07Z WARN  vav] Variant 2:29474101C>A total 7724; Ref 7638(0.9889); Proper alt 21(0.0027); Margin alt 17(0.0022); Lowq alt 0(0)
{
  "reference": 7638,
  "proper": 21,
  "margin": 17,
  "lowq": 0,
  "excessive": 0,
  "alleles": 48,
  "unknown": 0
}
```
