<p align="center">
  <a href="http://www.sinotechgenomics.com">
    <img height="70" src="http://www.sinotechgenomics.com/Upload/0/WebsiteLogo/WebsiteLogo_20170620092534731.png">
  </a>
  <h1 align="center">BreakID</h1>
</p>


BreakID is an genomic breakpoint prediction method that can discover and genotype translocations, tandem duplications, inversions and translocations at single-nucleotide resolution in short-read massively parallel sequencing data. It uses discordant read pairs and split-reads to sensitively and accurately delineate genomic rearrangements. 

Prerequisits for BreakID
----------------
nib files
you can prepare them by the following steps:

`mkdir nib_files`

`cd nib_files`

`rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/ .`

`gunzip *.gz`

`wget https://users.soe.ucsc.edu/~kent/src/blatSrc35.zip`

`unzip blatSrc35.zip`

Installing BreakID
----------------

The easiest way to get BreakID is to download a statically linked binary from the [BreakID github release page](https://github.com/SinOncology/BreakID/releases/). You can also build BreakID from source. 

`git clone https://github.com/SinOncology/BreakID.git`

`cd BreakID/`

`sh install.sh`

Running BreakID
--------------------------
BreakID needs a sorted, indexed and duplicate marked bam file for every input sample. An indexed reference genome is required to identify split-reads. The output is in plain txt format. 

`bin/BreakID -i input_bam -o out_prefix  -n /path/to/nib/directory`

FAQ
---
* BreakID is running too slowly what can I do?
You should exclude telomere and centromere regions and also all unplaced contigs. BreakID ships with such an exclude list for human and mouse samples. In addition, you can filter input reads more stringently using -q 20 and -s 15. You can also deactivate the small InDel calling using -n if you are only interested in large SVs (>=300bp).


* What pre-processing of bam files is required?    
Bam file must be produced by the short reads mapping algorithm that can mark alignments of Split Reads with SA tag. We suggest users bwa's mem mode for mapping. Bam files need to be sorted, indexed and ideally duplicate marked. If multiple libraries are present for a single sample these need to be merged in a single bam file with unique ReadGroup tags.

Citation
--------

Linfang Jin, Dai Heng, Bingding Huang.  
[BreakID: genomics breakpoints identification to detect gene fusion events using discordant pairs and split read.](http://bioinformatics.oxfordjournals.org/)  
Bioinformatics 2018 xx: xxxx-xxxx.


License
-------
BreakID is distributed under the GPLv3. Consult the accompanying [LICENSE](https://github.com/fionakim/BreakID/LICENSE) file for more details.
