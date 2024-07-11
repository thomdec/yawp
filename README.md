# YAWP

Yet Another Window Pi (part of the YAYAS familly of programs -- Yet Another Yet Another Software)

## Installation

```
git clone https://github.com/thomdec/yawp
```

### Dependencies:

- numpy
- sckit-allel
- tabix
- bedtools and bcftools (only when using `--mask`)

`yawp.py` requires `func.py` in the same directory 

## How to run YAWP?

YAWP requires three files to run: 

- A VCF file
- A mutlicallable bed file
- A population file

By default, YAWP compute global estimates. YAWP will run in windows if a windows size size is supplied with `-w`.  To obtain chromosome-wide estimates the window size must be set to a length larger than all of the chromosomes. 

Optional parameters: 

- Mask bed file (eg 4-fold degenerate sites)

### Examples

```bash
# global
python yawp.py -v in.vcf.gz -b multicallable.bed.gz -p populations.csv -o output_prefix

# global with mask
python yawp.py -v in.vcf.gz -b multicallable.bed.gz -p populations.csv -m mask.bed -o output_prefix

# windows
python yawp.py -v in.vcf.gz -b multicallable.bed.gz -p populations.csv -w 100_000 -o output_prefix
```

## What is a multicallabe bed file?

A multicallable bed file array is composed of the chromosome id, the start and ends of intervals (zero-based half-open intervals, following BED conventions), and the callability of each sample in the interval. Sample that are callable within the interval are recorded with $1$, sample with missing data in the interval are recorded with $0$. The order of the samples is the same as in the VCF. 

```
# example of a multi-callable bed array
# the order of samples is in the callable matrix is the same as in the VCF file, for ease of computation
contig1 0   100 1   0   1
contig1 100 150 1   1   1
```

## How to create the multicallable bed array?

For gIMble aficionado, the `ForgeMulticallable.py` will create a multicallable bed file for the gimbleprep filtered VCF using the temporary directory of gimbleprep (kept with the `-k` flag).

```bash
# using the gimble output
python ForgeMulticallable.py -g tmp_gimble_directory
```

[I will write a `yawpprep` module ASAP]

## The output 

YAWP outputs two files:

- a .pi.tsv with: chromosome, start, end, pop, pi, n_diff, n_comp
- a .dxy.tsv with: chromosome, start, end, pop1, pop2, dxy, n_diff, n_comp

Plotting in R:

```R
library(tidyverse)

df <- read_tsv("output.pi.tsv") %>%
	mutate(mid = (start + end) / 2)

ggplot(df, aes(mid, pi, col = pop)) +
  geom_line() +
  facet_wrap(~chromosome, scales = "free_x")
```

One can compute pairwise $F_{ST}$ from the outputs: 

$F_{ST} = (d_{xy} - \hat\pi) / (d_{xy} + \hat\pi)$, (where $\hat\pi$ is the mean value of $\pi$ for both populations)


Note that to obtain a global estimate, one can always run YAWP in windows (any size) and divide the sum of pairwise difference by the sum of pairwise comparison. 


## Why Yet Another Window Program (YAWP) for computing common popgen statistics?

State of the art scripts to compute $\pi$ and $d_{xy}$ already exist (eg. [Pixy](https://github.com/ksamuk/pixy) and [Simon Martin's `popgenWindows.py`](https://github.com/simonhmartin/genomics_general)). Why YAWP then? 

- Both pixy and `popgenWindows.py` require VCFs with invariant sites. I find that including invariant sites defeats the very purpose of the VCF (**Variant** Call Format)! Yet, encoding invariant sites is important to account for missing data [see](https://doi.org/10.1111/1755-0998.13326). YAWP uses a bed file to encode callable variants less redundantly. 

For now, YAWP is just a set of scripts. However, I will try to make it into a package asap! 


## Acknowledgement

Reading through pixy and gIMble repository helped me a lot in this project! I am grateful for the work of coders involved in both packages
