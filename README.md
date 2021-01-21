# ChIPnCHOP 

![header_chipnchop](https://github.com/jvazpa/chipnchop/blob/main/format/header.png)

A Chromatin Immuno-Precipitation analyisis pipeline for histones, transcription factors and other proteins that interact with DNA developed by Biochemistry Degree students at the University of Seville.

### ✍️ Authors

Ana González Toro, Joaquín Tamargo Azpilicueta (joatamazp@alum.us.es), José Vázquez Pacheco (josvazpac@gmail.com)

### About & Usage

When you have downloaded the repository, you must specify the path of the parameters that are needed in the test_params.txt file, located in the test folder. Then, in the chipnchop folder you must write in the terminal "bash chipnchop <FULL/PATH/to/params/test_params.txt> [options]". As for <FULL/PATH/to/params/test_params.txt>, you must CLARIFY THE COMPLETE DIRECTORY TO PARAMETERS FILE. With regards to [options], you must write -TF if the analyzed proteins are transcription factors or interact with the chromatine in a similar way, or -HI if the proteins analyzed are histones or other specific epigenetic marks that affect large portions of DNA, in contrast to transcription factors.

### Requiered Dependencies and Specifications

For running this pipeline, you have to check that the following softwares are installed in your system:
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
[Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
[SAMTOOLS](https://sourceforge.net/projects/samtools/files/samtools/)
[MACS2](https://github.com/macs3-project/MACS)
[R](https://www.r-project.org/)
[BiocManager](https://cran.r-project.org/web/packages/BiocManager/vignettes/BiocManager.html)
[ChIPseeker](https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html)
[clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html)
[TxDb.Athaliana.BioMart.plantsmart28](https://bioconductor.org/packages/release/data/annotation/html/TxDb.Athaliana.BioMart.plantsmart28.html)
[org.At.tair.db](https://bioconductor.org/packages/release/data/annotation/html/org.At.tair.db.html)
[DOSE](https://bioconductor.org/packages/release/bioc/html/DOSE.html)
[enrichplot](https://bioconductor.org/packages/release/bioc/html/enrichplot.html)
[Homer](http://homer.ucsd.edu/homer/download.html)

### Installation

First of all, you must open the terminal and enter the folder where you want to download this repository. In order to get it, you must write in the terminal "git clone https://github.com/jvazpa/chipnchop.git". You can find this link in the CODE section of this repository. Then, you must specify your github account and password and the download will begin.

### Troubleshooting



##### ChIPnCHOP doesn't read properly the parameters file

Once you have checked you have specified the correct path to the parameter file, you must ascertain that there are one space between the colon and the parameter. For instance, "working_directory:/home/mickeymouse/tmp/" won't be read properly, whereas "working_directory: /home/anagontor1/tmp/" will.

### Status

ChIPnCHOP has been tested on MacOS Catalina, MacOS Big Sur and Ubuntu 16.04. If you find a bug, please do report it on the GitHub issue tracker. 

### Roadmap

* Make multiple sampling parallelization available.  




