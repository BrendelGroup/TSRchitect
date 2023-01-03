# TSRchitect: Promoter identification from diverse types of large-scale TSS profiling data

The TSRchitect repository encompasses an [R](https://www.r-project.org/)
package developed in the [Brendel Group](http://brendelgroup.org/) for analyses
of transcription start site data.
The code conforms to our [RAMOSE](https://brendelgroup.github.io/)
philosophy: it generates __reproducible__, __accurate__, and __meaningful__
results; it is __open__ (source) and designed to be __scalable__ and
__easy__ to use.


## Quick Start

Input to TSRchitect will be transcription profiling read alignment data in `bam`
or `bed` format as well as the appropriate genome annotation (if available).
Output consists of predicted Transcription Start Sites (TSS) and Transcription
Start Regions (TSR) as well as statistics summarizing the distribution and
characteristics of identified TSSs and TSRs.

The simplest way to get going is to use the TSRchitect
[Singularity](https://apptainer.org/) container, e.g. as follows:

```bash
git clone https://github.com/BrendelGroup/TSRchitect
cd TSRchitect/demo
wget https://BrendelGroup.org/SingularityHub/tsr.sif
alias rws="singularity exec -e -B{PWD}/..  tsr.sif"
rws R
```

In the above example, you clone this repository into your current directory,
go into the TSRchitect/demo directory that has been created, download the TSRchitect
apptainer, define the bash alias _rws_ ("run with singularity"), and check that
everything works by launching an [R](https://www.r-project.org/) console from within
the container.

Of course this assumes that you have [Apptainer/Singularity](https://apptainer.org/) installed on your system.
Check whether there is package built for your system.
Otherwise, follow the instructions to [install Singularity from source code](https://apptainer.org/user-docs/master/quick_start.html#quick-installation-steps).

The advantage of this approach is that the TSRchitect library and all its prerequisites
are available within the container, so that there is no further installation necessary
on your part to follow our examples and run your own analyses.
For example, in that console, you should see

```
R version 4.2.2 Patched (2022-11-10 r83330) -- "Innocent and Trusting"
...
> packageVersion("TSRchitect")
[1] '1.17.3'
> packageVersion("GenomicRanges")
[1] '1.50.2'
>
```

## Realistic Start

Please find detailed installation instructions and options in the
[INSTALL](./INSTALL.md) document.
Once all preparatory steps are taken care of, see the [HOWTO](./demo/HOWTO.md)
document for examples of how to load data into TSRchitect and predict and
characterize promoters.


## FAQ and References

Please see
[V. Brendel's TSRchitect FAQ](https://github.com/vpbrendel/TSRchitect/wiki/FAQ)
for usage examples and suggestions.

If you find _TSRchitect_ useful, you may cite:

Raborn RT, Sridharan K, Brendel VP (2017)
_TSRchitect: Promoter identification from large-scale TSS profiling data._
doi: 10.18129/B9.bioc.TSRchitect, https://doi.org/doi:10.18129/B9.bioc.TSRchitect. 


## Contact

Please direct all comments and suggestions to
[Volker Brendel](<mailto:vbrendel@indiana.edu>)
at [Indiana University](http://brendelgroup.org/) and
[Taylor Raborn](<mailto:rtraborn@asu.edu>) at his current address.
