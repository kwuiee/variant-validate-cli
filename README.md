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

```shell
$ vav tests/many_variants.bam --var '1:156843458A>G' --var "2:29474101C>A"
```
