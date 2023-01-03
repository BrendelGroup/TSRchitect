# TSRchitect HOWTO - embedding TSRchitect in workflows

## Preparation

At this stage, you should have completed the TSRchitect installation steps
documented in the [INSTALL](../INSTALL.md) document; we'll assume that you have
downloaded the `tsr.sif` singularity container.
We explain how you can use _bash_ and _Rscript_ scripts to exectute the
examples in the [TSRchitectUsersGuide](../inst/doc/TSRchitectUsersGuide.Rmd).


## Example 1
This example only requires the _xdemo1_ script you will find in this directory.
Take a look at the script: it's simply a compilation of the instructions from
the relevant section in the
[TSRchitectUsersGuide](../inst/doc/TSRchitectUsersGuide.Rmd), to be executed in
your shell.
The following command may need to be adjusted to your needs.
Here we assume that the the _tsr.sif_ container is located in the parent of your
current working directory.

```
singularity exec -e -B${PWD}/.. ../tsr.sif  ./xdemo1
```

## Example 2
As with the previous example, first we must download the raw data. In this case we have only a single alignment file to retrieve, which is found here: https://oregonstate.app.box.com/s/3lb3spmqbiuofhbubovc1z8bfthmxh6f .
Once download of the file `peat.sorted.bam` is complete, please move it to subdirectory "PEATbam/".
Then execution of _xdemo2_ should work (see __Example 1__ for setup details):

```
singularity exec -e -B /projects/vbrendel/TMP /DATA/GROUP/prj/SINGULARITY/tsr.simg  ./xdemo2
```
## Example 3

```
singularity exec -e -B /projects/vbrendel/TMP /DATA/GROUP/prj/SINGULARITY/tsr.simg  ./xdemo3

