# TSRchitect HOWTO - embedding TSRchitect in workflows

## Preparation

At this stage, you should have completed the TSRchitect installation steps
documented in the [INSTALL](../INSTALL.md) document; we'll assume that you have
downloaded the `tsr.simg` singularity container.
We explain how you can use _bash_ and _Rscript_ scripts to exectute the
examples in the [TSRchitectUsersGuide](../inst/doc/TSRchitectUsersGuide.Rmd).


## Example 1
This example only requires the _xdemo1_ you will find in this directory.
Take a look at the script: it's simply a compilation of the instructions from
the relevant section in the
[TSRchitectUsersGuide](../inst/doc/TSRchitectUsersGuide.Rmd), to be executed in
your shell.
The following command needs to be adjusted to your needs.
Here we assume that the working directory is somewhere under the bound directory
_/projects/vbrendel/TMP_ and that the _tsr.simg_ container is located in
_/DATA/GROUP/prj/SINGULARITY_.

```
singularity exec -e -B /projects/vbrendel/TMP /DATA/GROUP/prj/SINGULARITY/tsr.simg  ./xdemo1
```

## Example 2
As with the previous example, we first must download the raw data. In this case we have only a single alignment file to retrieve, which is found here: https://oregonstate.app.box.com/s/3lb3spmqbiuofhbubovc1z8bfthmxh6f
Once download of the file `peat.sorted.bam` is complete, please move it to subdirectory "PEATbam/".
Then execution of _xdemo2_ should work (see __Example 1__ for setup details):

```
singularity exec -e -B /projects/vbrendel/TMP /DATA/GROUP/prj/SINGULARITY/tsr.simg  ./xdemo2
```
## Example 3

```
singularity exec -e -B /projects/vbrendel/TMP /DATA/GROUP/prj/SINGULARITY/tsr.simg  ./xdemo3

