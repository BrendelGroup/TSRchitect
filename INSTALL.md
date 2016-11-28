# __TSRchitect__ Installation and Setup

## Preliminary Steps
__TSRchitect__ is package for the [R software environment for statistical
computing](https://www.r-project.org/).
We assume that you have a current version of [R](https://www.r-project.org/)s
 installed and are familiar with its use.
Several other packages are prerequisite for __TSRchitect__, and these packages
are available through [_Bioconductor_](http://bioconductor.org/).
We assume that you have a current version of [R](https://www.r-project.org/)
installed and are familiar with its use.
Typical preliminary steps to install or updated these packages are as follows
(note: for system-wide installation, you would need to invoke
[R](https://www.r-project.org/) as root):

```bash
R
>source("http://bioconductor.org/biocLite.R")
>biocLite(c("GenomeInfoDb", "GenomicAlignments", "GenomicRanges", "Rsamtools", "BiocParallel", "ggbio", "Rsubread", "gtools"))
q()
```
## Obtaining TSRchitect
Eventually, __TSRchitect__ will also be available as a
[_Bioconductor_](http://bioconductor.org/) package.
For now (and in general, for the latest development version), you can install
__TSRchitect__ as follows (again, use root privileges on the install command for
system-wide installation):

```bash
git clone https://github.com/BrendelGroup/TSRchitect
R CMD INSTALL TSRchitect
```

To check on successful installation, load the package in your
[R](https://www.r-project.org/) console:

```bash
R
>library(TSRchitect)
q()
```
Either fix problems according to displayed error messages, or go on to the next
step.

## Finally

proceed to the [HOWTO](./HOWTO.md) document, and we'll tell you how to execute
sample workflows (or, equally easy, your very own data analyses).