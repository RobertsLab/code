`code/apptainer_definition_files`

### Apptainer (Singularity) Definition Files

This is a repository for all of [the Roberts Lab](https://robertslab.github.io/resources/) definition files used for creating [Apptainer](https://apptainer.org/) images.

---

- [`bedtools-2.31.1.def`](https://github.com/RobertsLab/code/blob/master/apptainer_definition_files/bedtools-2.31.1.def): Creates image with [bedtools v2.31.1](https://bedtools.readthedocs.io/en/latest/) and installs in system `$PATH`. Bedtools can be called by entering `bedtools`.

- [`isoform-expression-hisat2-stringtie.def`](https://github.com/RobertsLab/code/blob/master/apptainer_definition_files/isoform-expression-hisat2-stringtie.def): Creates an image geared toward isoform expression analysis with the following software:

  - [`hisat2`](https://daehwankimlab.github.io/hisat2/)

  - [`gffcompare`](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml)

  - [`samtools`](https://www.htslib.org/doc/samtools.html)

  - [`stringtie`](https://ccb.jhu.edu/software/stringtie/)

- [`ubuntu-22.04-base.def`](https://github.com/RobertsLab/code/blob/master/apptainer_definition_files/ubuntu-22.04-base.def): The "base" image upon which all other images are built. Runs Ubuntu 22.04 and contains the necessary programs/libraries needed for the building/installation of other software.

