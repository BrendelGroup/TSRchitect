# TSRchitect: Promoter identifcation from diverse types of large-scale TSS profiling data

The TSRchitect repository encompasses an [R](https://www.r-project.org/)
package developed in the [Brendel Group](http://brendelgroup.org/) for analyses
of transcription start site data.
The code conforms to our [RAMOSE](https://brendelgroup.github.io/)
philosophy: it generates __reproducible__, __accurate__, and __meaningful__
results; it is __open__ (source) and designed to be __scalable__ and
__easy__ to use.


## Quick Start [![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/1204)

Input to TSRchitect will be transcription profiling read alignment data in `bam`
or `bed` format as well as the appropriate genome annotation (if
available).
Output consists of predicted Transcription Start Sites (TSS) and Transcription
Start Regions (TSR) as well as statistics summarizing the distribution and
characteristics of identified TSSs and TSRs.

All the TSRchitect dependencies are encapsulated in a
[Singularity](https://www.sylabs.io/docs/) container available from
[Singularity Hub](https://singularity-hub.org/).
Thus, once you know what you are doing, execution could be as simple as

```
singularity pull --name tsr.simg shub://BrendelGroup/TSRchitect
singularity exec tsr.simg R
```

which will bring up an [R](https://www.r-project.org/) console with the
TSRchitect library and all its prerequisites available.
For example, in that console, you should see

```
R version 3.5.3 (2019-03-11) -- "Great Truth"
...
> packageVersion("TSRchitect")
[1] '1.8.9'
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

Our own publications will be linked here in due course.

## Contact

Please direct all comments and suggestions to
[Volker Brendel](<mailto:vbrendel@indiana.edu>)
at [Indiana University](http://brendelgroup.org/) and
[Taylor Raborn](<mailto:rtraborn@asu.edu>) at his current address.
