# CHANGELOG

## 01/11/2022
* Added bowtie2 component, changes made to bin/ docs/ and pipeline templates
* Changes to run on Rhel8 to config file and shell scripts in bin/
* Removed expression plot component
* Replace deseq with deseq2
* Edits to the create scripts to accomodate changes. 
* Edits to prok config file to add components
* Additional bug fixes as reported on hisat component

## 03/29/18
* Moved samtools bin from /usr/src/samtools to /opt/packages/samtools to keep all package locations consistent
* Moved fastqc bin from /usr/bin/ to /opt/packages/fastqc since FastQC will be installed by source instead of by apt-get install

## 03/26/18
* Fixing issues where software.config was pointing to internal paths instead of the paths Docker needs

## 10/12/18
* Added BDBAG component
* changed create_euk script to make bag
* changed config file and added xml files

## 09/07/18
* Deseq2 replaces Deseq
* edited htseq shell script

## 9/21/17
* Various changes with Bedtools scripts to make use of modern versions of Samtools and Bedtools
* Fixed bug that doubled read counts for single-stranded reads in align_hisat_(split_)stats.pl

## 4/28/17
* Updated 'samtools_file_convert' to be compatible with modern versions of SAMtools

## 4/13/17
* Changed various options in pipeline configuration scripts
* Fixed various bugs

## 3/20/17
* Uncommented out code in align_tophat_stats.pl that dealt with single-end reads

## 2/16/17
* Initial commit of package-rnaseq

