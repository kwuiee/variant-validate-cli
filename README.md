# Variant Validate

Validate a SNP/InDel from bam using CIGAR and MD tag.

## Getting Started

```shell
vav --help
vav

USAGE:
    vav [FLAGS] [OPTIONS] --bam <bam> --var <var>

FLAGS:
    -h, --help       Prints help information
    -v               Print verbose info.
    -V, --version    Prints version information

OPTIONS:
        --bam <bam>          Input bam file.
        --mapq <mapq>        Minimum read mapping quality. [default: 30]
        --margin <margin>    Minimum margin base distance for alt support. Margin stands for read
                             start/end, softclip start/end etc. [default: 10]
        --var <var>          Input genome variant, e.g. 'chr1:12345AT>-'.
```

## Examples

```shell
vav \
	--bam  tests/many_variants.bam \
    --var "2:29474101C>A"
```
