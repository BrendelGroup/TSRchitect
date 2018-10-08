# __TSRchitect__ Installation and Setup

## Preliminary Steps
__TSRchitect__ is a package for the [R software environment for statistical
computing](https://www.r-project.org/).
We assume that you have a current version of [R](https://www.r-project.org/)
 installed and are familiar with its use.
Several other packages are prerequisite for __TSRchitect__, and these packages
are available through
[CRAN](https://cran.r-project.org/) or
[Bioconductor](http://bioconductor.org/).
Typical preliminary steps to install or update these packages are as follows
(note: for system-wide installation, you would need to invoke
[R](https://www.r-project.org/) as root):

```{r eval=FALSE}
#installing CRAN packages
install.packages(c("gtools","knitr"))

```{r eval=FALSE}                              
if (!requireNamespace("BiocManager", quietly=TRUE))
   install.packages("BiocManager")
#installing Bioconductor packages
BiocManager::install(c("AnnotationHub", "BiocGenerics", "BiocParallel",
"ENCODExplorer", "GenomicAlignments", "GenomeInfoDb",
"GenomicRanges", "IRanges", "Rsamtools", "rtracklayer",
"S4Vectors", "SummarizedExperiment"))
```

## Obtaining TSRchitect
__TSRchitect__ is available as a
[Bioconductor](http://bioconductor.org/) package and thus can be installed in
the same way as the prerequisite packages:

```{r eval=FALSE}
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("TSRchitect")
```

Optionally, you can install __TSRchitect__ directly from our group's GitHub
repository as follows (again, use root privileges on the install command for
system-wide installation):

```{bash eval=FALSE}
git clone https://github.com/BrendelGroup/TSRchitect
R CMD INSTALL TSRchitect
```

Note that the `master` branch should be a stable version, possibly somewhat
ahead of the Biocondcutor version in case of minor edits or addtions.
You may also see a `devel` branch, which would contain recently developed
addtions that are not yet fully tested.
To clone and install just the 'devel' branch, please do the following:
```{bash eval=FALSE}
git clone https://github.com/BrendelGroup/TSRchitect --branch devel --single-branch
R CMD INSTALL TSRchitect
```
In either case, please use the _GitHub Issues_ button to report problems or
suggest novel features.


## Installation as a singularity container

Assuming a current version of _singularity_ is installed on your system (if not,
it should be, and it's easy to do so), you can get the TSRchitect container from
[Singularity Hub](https://www.singularity-hub.org/collections/1204) as follows:

```bash
singularity pull --name tsr.simg shub://BrendelGroup/TSRchitect
```

For a gentle introduction to singularity, see our group
[handbook article](https://github.com/BrendelGroup/bghandbook/blob/master/doc/06.2-Howto-Singularity-run.md).


## Finally

Proceed to the [TSRchitectUsersGuide](./inst/doc/TSRchitectUsersGuide.Rmd) for
examples of how the use TSRchitect functions.
If you are interested in developing and running entire workflows for
transcription start site profiling, please see our paper or simple worked
examples in the [HOWTO](./demo/HOWTO.md) document.
