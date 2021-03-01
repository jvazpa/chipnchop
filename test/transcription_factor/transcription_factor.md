**Analysis of transcription factors with chipnchop**

**1. Introduction**

PRR5 is a Pseudo-Response Regulator whose mutation affects various circadian-associated biological events such as: flowering time in the long-day photoperiod conditions, red light sensitivity of seedlings during early photomorphogenesis, and acts as transcriptional repressor. regulating the hypocotyl growth under photoperiodic conditions.


**2. Usage**

Usage is really simple. You just have to call chipnchop function (‘bash chipnchop’) plus the COMPLETE PATH to the parameters file. Apart from that, you will need to write -TF in order to activate the option that will allow you to work with these transcription factors. 

After checking all parameters are properly read, just wait… Depending on how many samples you have, as well as the weight of genome and annotation files, it should last longer or be faster.  

**3. Usage example**

We suggest using GNU screen as you are not yet allowed to parallelize the process.

After downloading (and opening a screen), go to chipnchop directory and type:

`bash chipnchop  full/installation/path/to/test/transcription_factor/test_params.txt -TF	`

Check if directories and parameters are correctly read.

After that, chill. This is going to take a bit. This example has been adapted, so it shouldn’t take more than 10 min. 
