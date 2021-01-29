**Analysis of epigenetic marks with chipnchop**

**1. Introduction**

Transcriptional repressive mark H3K27me3 (histone H3 lysine 27 trimethylation) is essential for inheritance of the silencing memory from mother to daughter cells. In this file, we included genome and annotation files as well as 2 ChIP-seq samples corresponding to the CHIP and INPUT. 

**2. Usage**

Usage is really simple. You just have to call chipnchop function (‘bash chipnchop’) plus the COMPLETE PATH to the parameters file. Apart from that, you will need to write -HI in order to activate the option that will allow you to work with theses epigenetic marks. 

After checking all parameters are properly read, just wait… Depending on how many samples you have, as well as the weight of genome and annotation files, it should last longer or be faster.  

**3. Usage example**

We suggest using GNU screen as you are not yet allowed to parallelize the process.

After downloading (and opening a screen), go to chipnchop directory and type:

`bash chipnchop  full/installation/path/to/test/epigenetic\_mark/test\_params.txt -HI	`

Check if directories and parameters are correctly read.

After that, chill. This is going to take a bit. This example has been adapted, so it shouldn’t take more than 10 min. 
