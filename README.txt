This directory contains the Feb. 2009 assembly of the human genome (hg19,
GRCh37 Genome Reference Consortium Human Reference 37 (GCA_000001405.1))
in one gzip-compressed FASTA file per chromosome.

The Feb. 2009 human reference sequence (GRCh37) was produced by the
Genome Reference Consortium:
	http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/

Note on chrM:
Since the release of the UCSC hg19 assembly, the Homo sapiens mitochondrion
sequence (represented as "chrM" in the Genome Browser) has been replaced in
GenBank with the record NC_012920. We have not replaced the original sequence,
NC_001807, in the hg19 Genome Browser. We plan to use the Revised Cambridge
Reference Sequence (rCRS, http://mitomap.org/bin/view.pl/MITOMAP/HumanMitoSeq)
in the next human assembly release.

Files included in this directory:

  - chr*.fa.gz: compressed FASTA sequence of each chromosome.

    Repeats from RepeatMasker and Tandem Repeats Finder (with period
    of 12 or less) are shown in lower case; non-repeating sequence is
    shown in upper case.

    RepeatMasker was run with the -s (sensitive) setting.
    Using: Jan 29 2009 (open-3-2-7) version of RepeatMasker and
    RELEASE 20090120 of library RepeatMaskerLib.embl

------------------------------------------------------------------
See also: Wellcome Trust Sanger Institute MHC Haplotype Project
for additional information on the chr6 alternate haplotype assemblies:
	http://www.sanger.ac.uk/HGP/Chr6/MHC/

The chr*_random sequences are unplaced sequence on those reference
chromosomes.

The chrUn_* sequences are unlocalized sequences where the corresponding
reference chromosome has not been determined.

See also: NCBI discussion of assembly procedures:
    http://www.ncbi.nlm.nih.gov/genome/assembly/assembly.shtml 

------------------------------------------------------------------
If you plan to download a large file or multiple files from this 
directory, we recommend that you use ftp rather than downloading the 
files via our website. To do so, ftp to hgdownload.cse.ucsc.edu, then 
go to the directory goldenPath/hg19/chromosomes. To download multiple 
files, use the "mget" command:

    mget <filename1> <filename2> ...
    - or -
    mget -a (to download all the files in the directory)

Alternate methods to ftp access.
    
Using an rsync command to download the entire directory:
    rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/ .
For a single file, e.g. chrM.fa.gz
    rsync -avzP 
        rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chrM.fa.gz .
    
Or with wget, all files:
    wget --timestamping 
        'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/*'
With wget, a single file:
    wget --timestamping 
        'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chrM.fa.gz' 
        -O chrM.fa.gz
    
To uncompress the fa.gz files:
    gunzip <file>.fa.gz

All the files in this directory are freely available for public use.
