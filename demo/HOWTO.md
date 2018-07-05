# TSRchitect HOWTO - embedding TSRchitect in workflows

## Preparation

At this stage, you should have completed the TSRchitect installation steps
documented in the [INSTALL](../INSTALL.md) document; we'll assume that you have
downloaded the `tsr.simg` singularity container.
We explain how you can use _bash_ and _Rscript_ scripts to exectute the
examples in the [TSRchitectUsersGuide](../inst/doc/TSRchitectUsersGuide.Rmd).


## Example 1
```
singularity exec -e -B /projects/vbrendel/TMP /DATA/GROUP/prj/SINGULARITY/tsr.simg  ./xdoit1
```

