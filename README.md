# ChIPnCHOP 

![header_chipnchop](https://github.com/jvazpa/chipnchop/blob/main/format/header.png)

[![Project Status: Concept – Minimal or no implementation has been done yet, or the repository is only intended to be a limited example, demo, or proof-of-concept.](https://www.repostatus.org/badges/latest/concept.svg)](https://www.repostatus.org/#concept)

    <a href="https://discord.gg/CZQhuCYK">
        <img src="https://img.shields.io/discord/308323056592486420?logo=discord"

### ✍️ Authors

Ana González Toro (anagtoro7@gmail.com), Joaquín Tamargo Azpilicueta (joatamazp@alum.us.es), José Vázquez Pacheco (josvazpac@gmail.com)

### 🧩 About & Usage

ChIPnCHOP is a Chromatin Immuno-Precipitation analyisis pipeline for histones, transcription factors and other proteins that interact with DNA in *Arabidopsis thaliana* col-0 developed by Biochemistry Degree students at the University of Seville.

Once you have downloaded the repository, using it is easy as pie 🍰. First, you must specify the path of the parameters that are needed in the test_params.txt file, located in the test folder. Then, in the chipnchop folder you must call chipnchop by typing: `bash chipnchop <FULL/PATH/to/params/test_params.txt> [options]`. Remember, in `<FULL/PATH/to/params/test_params.txt>`, you must CLARIFY THE COMPLETE DIRECTORY TO PARAMETERS FILE. 

With regard to other [options], you must write `-TF` if the analyzed proteins are transcription factors or interact with the chromatine in a similar way, or `-HI` if the proteins analyzed are histones or other specific epigenetic marks that affect large portions of DNA, in contrast to transcription factors. By default, chipnchop will assume you are working with a transcription factor (-TF) unless you state otherwise by using option -HI. Basically, changing this option includes slight changes in some of the functions that are used so that they are optimized for the case study.

There are two examples available that you can check before you perform your own analysis. Those can be found at "test" file, within this repository. Depending on your interests, you could either try the analysis on transcription factors or epigenetic marks. 

### 🔗 Requiered Dependencies and Specifications

For running this pipeline, you have to check that the following software are installed in your system:
* [R](https://www.r-project.org/)
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [SAMTOOLS](https://sourceforge.net/projects/samtools/files/samtools/)
* [MACS2](https://github.com/macs3-project/MACS)
* [Homer](http://homer.ucsd.edu/homer/download.html)

At R, you must have downloaded these following packages previously:
* [BiocManager](https://cran.r-project.org/web/packages/BiocManager/vignettes/BiocManager.html)
* [ChIPseeker](https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html)
* [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html)
* [TxDb.Athaliana.BioMart.plantsmart28](https://bioconductor.org/packages/release/data/annotation/html/TxDb.Athaliana.BioMart.plantsmart28.html)
* [org.At.tair.db](https://bioconductor.org/packages/release/data/annotation/html/org.At.tair.db.html)
* [DOSE](https://bioconductor.org/packages/release/bioc/html/DOSE.html)
* [enrichplot](https://bioconductor.org/packages/release/bioc/html/enrichplot.html)

### 💻 Installation

First of all, you must open the terminal and get into the folder where you want to download this repository. In order to get it, you must write in the terminal "git clone https://github.com/username/chipnchop.git". You can find this link in the CODE green button above. Then, you must specify your github account and password and the download will begin. That's all!

### 🎯 Troubleshooting

#### * grep: ../../parameters/parameter_file.txt: No such file or directory

This error message comes up when the parameter file is not correctly specified or, most likely, if the path is not FULLY specified. In other words, you must type the whole path, and not going back and forth with double dots (..). 

##### * ChIPnCHOP doesn't read properly the parameters file

Once you have checked you have specified the correct path to the parameter file, you must ascertain that there are one space between the colon and the parameter. For instance, "working_directory:/home/mickeymouse/tmp/" won't be read properly, whereas "working_directory: /home/anagontor1/tmp/" will.

### 📍 Status

ChIPnCHOP has been tested on MacOS Catalina, MacOS Big Sur and Ubuntu 20.04. 

### 🗺 Roadmap

* Make multiple sampling parallelization available by using Sun Grid Engine (SGE),  Simple Linux Utility for Resource Management (Slurm) or similar. Yet there were a first version in which we included SGE parallelization, it has been tested with obsolete versions of the software. Thus, it has to be re-tested so that it can be used that way.
* Make the the pipeline executable for working with pair-end samples.
* Make it possible to work with other ecotype or species different from Arabidopsis thaliana col-0.

### License

This software is licensed 



