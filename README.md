# fortuna

fortuna is an alignment-free method to quantify anambiguous building blocks of transcripts and to detect and annotate novel splicing
events. It proceeds in the steps shown in the figure below.

<p align="center">
<img src="https://user-images.githubusercontent.com/37735817/149656945-bba990cf-d608-4901-9fc5-a3adfe568fd0.jpg" width=45% height=50%>
</p>

fortuna starts (A) by "guessing" novel transcripts based on annotated splice sites. It then (B) creates a set of sequence fragments of annotated and guessed novel transcripts that contain all possible combinations of unspliced exonic segments. From this set of fragments we build a kallisto [2] index (C) and use it to efficiently pseudoalign reads to fragments (D), which yields counts of the most elementary splicing units, the **signature counts**. Optionally, fortuna can further incorporate novel splice sites (e.g. segment s<sub>2</sub>) identified by **any** spliced aligner from reads that remained unmapped in step (D). Computed counts can be directly used for alternative splicing analysis or added up to larger units such as those used by DEXSeq [3] (E1). In addition, fortuna annotates all novel events (E2) based on precise definitions of event types.

#### GENCODE fortuna index

If you want to use fortuna with human transcripts annoated in GENCODE (https://www.gencodegenes.org/human/), release 42 (GRCh38.p13, access date 17.12.2022), you can download here the fortuna index for different read lengths Other species will follow.
- [*H. sapiens*, read length 75bp](https://fizika-my.sharepoint.com/:u:/g/personal/lborozan_unios_hr/ERy8k_wj6DBGj03jVFphQ14BKjdGtogzU_43M7Wqm6CyHQ?e=zWYYBA&download=1),
- [*H. sapiens*, read length 100bp](https://fizika-my.sharepoint.com/:u:/g/personal/lborozan_unios_hr/EfRsiEzTOyJOvDoqflZkfCABTPLiIIKxMqUX6QvZk2jdpA?e=jwe3id&download=1),
- [*H. sapiens*, read length 150bp](https://fizika-my.sharepoint.com/:u:/g/personal/lborozan_unios_hr/ESzjW3JAQUNHrsuMN1Vum2kBLaKx7-E_r0tcd16FUAsH9g?e=q6hCne&download=1).

## Dependencies

Fortuna requires C++11 to run and CMake version >=2.8.12 to compile. 

Also, fortuna requires
* zlib library (included in most UNIX distributions)
* autoconf 2.69 (included in most UNIX distributions, Mac OS users may have to run ```brew install autoconf@2.69```)

Additional dependencies are included with fortuna, including
* kallisto (https://github.com/pachterlab/kallisto)
* htslib (https://github.com/samtools/htslib)


## Installation

Mac OS users may have to install autoconf 2.69 first (see Dependencies):

```
brew install autoconf@2.69
PATH="/opt/homebrew/opt/autoconf@2.69/bin:$PATH“
```

After cloning the repository, run the following commands in the fortuna folder. 

``` 
mkdir build 
cd build 
cmake .. 
make 
```

## Pre-compiled fortuna binaries
In addition to compiling from source, pre-compiled fortuna binaries can be downloaded for the x86_64 architecture running [Linux](https://github.com/canzarlab/fortuna/files/10284967/fortuna.zip) or [Mac OS](https://github.com/canzarlab/fortuna/files/10295151/fortuna_mac_x86.zip) or for [Mac OS on the Arm64 architecture](https://github.com/canzarlab/fortuna/files/10289050/fortuna_mac.zip) (M1/M2).


## Docker

Alternatively, fortuna can be run in a docker container. The image can be downloaded [here](https://fizika-my.sharepoint.com/:u:/g/personal/lborozan_unios_hr/EZapKOTGS05AhDK7821YRW4BcBx_nV6Yk6IQXz3YiPU_Lg?e=5s3nLR&download=1). It has been built on Linux and was tested on Linux and Mac OS, but not on Windows.

Then fortuna can be run on a small example dataset (see below) using the following commands:

``` 
# Download and load docker image
wget -O fortuna1-00.tar "https://fizika-my.sharepoint.com/:u:/g/personal/lborozan_unios_hr/EZapKOTGS05AhDK7821YRW4BcBx_nV6Yk6IQXz3YiPU_Lg?e=5s3nLR&download=1"
docker load -i fortuna1-00.tar

# Download sample
wget -O sample.zip "https://fizika-my.sharepoint.com/:u:/g/personal/lborozan_unios_hr/EaprBlTdJU5Ekt6-EKRkyr8B1KFRbbptoCDLjrFMgv-pqQ?e=vmPJVY&download=1"
unzip sample.zip -d sample

# Run docker image with the sample
docker run -it -v "$PWD"/sample:/fortuna/sample fortuna:1.00

# Once inside the docker image, run the following commands
cd sample
bash sample.sh
``` 
## Example dataset
 
An example dataset can be downloaded [here](https://fizika-my.sharepoint.com/:u:/g/personal/lborozan_unios_hr/EaprBlTdJU5Ekt6-EKRkyr8B1KFRbbptoCDLjrFMgv-pqQ?e=vmPJVY&download=1). It contains approximately 240,000 reads taken from sample 4 (patient id 11352.p1) in [4] that map to gene STAT1. fortuna will map over 99% of those reads and find slightly over 1000 reads supporting a novel exon skipping at chr2:190984393-190986854. After extracting the data to the root fortuna folder, use the following commands.
 
``` 
cd sample
bash sample.sh
```

Two tab-separated output files will be created after a couple of seconds - cnt.tsv and alt.tsv. 

File cnt.tsv contains signature counts in the format documented below. An example line from the file follows.
```
223832|223833|223836|   2274
``` 
The above example implies that a signature comprised of subexons 223832, 223833, and 223836 has a count of 2274. The IDs of the subexons are taken from the input GTF file sample.rg.gtf (field NodeId). Note that the GTF file has been pre-processed using script ``` bin/processGTF.sh ```.

File alt.tsv contains novel splicing information. To obtain the information regarding the previously mentioned novel exon skipping, use the following command. Mac OS does not support ```grep -P```. Omit ```-P``` and, depending on your implementation of grep, you might have to replace ```\t``` by ```^V<tab>``` (press ```ctrl+V``` then press ```tab```).
```
grep -P "chr2\t190984393\t190986854\t" alt.tsv
```
The result should be:
```
chr2 190984393 190986854 NM_007315.3,NM_139266.2,XM_006712718.1,XM_017004783.2,XR_001738914.2,XR_001738915.2 Locus_18066 ES 1256
```
which means that on chromosome 2 there exists a novel intron between coordinates 190984393 and 190986854 on gene Locus_18066 that is a novel exon skipping with regards to transcripts NM_007315.3, NM_139266.2, XM_006712718.1, XM_017004783.2, XR_001738914.2, and XR_001738915.2 which is supported by 1256 reads. Just like subexon IDs in the cnt.tsv file, chromosome, gene (gene_id field), and transcript IDs (transcript_id field) are taken from the input GTF file.

Below is an example sashimi plot of the novel exon skipping we generated using IGV.
<img width="853" alt="sashimi" src="https://user-images.githubusercontent.com/37735817/208976495-17478028-b299-4f54-a6b2-1169f33e034b.png">


## Usage

Below we describe how to run the individual steps shown in the top figure. For a complete list of options, run fortuna without any parameters.

 ``` ./fortuna ```

### Preprocessing

A transcript annotation in GTF format needs to be preprocessed by computing minimal disjoint segments bounded by splice sites (segments s<sub>1</sub>, s<sub>2</sub>, s<sub>3</sub>, s<sub>4</sub> in the figure).  Running bash script ``` bin/processGTF.sh ``` with the path to the annotation

``` bash processGTF.sh /path/to/the/annotation.gtf ```

will generate a segmented annotation in ``` /path/to/the/annotation.rg.gtf ```.


### Indexing

Steps B and C are executed using fortuna's "--index" command. Given a genome sequence in FASTA format, a preprocessed transcript annotation (GTF),
and a target read length, it builds an index of annotated and novel transcript fragments that can be used for efficient quantification using kallisto's pseudoalignments (next step). The 
following options are available, mandatory arguments are marked by (*). 

* ```-ingtf <FILE>``` input preprocessed GTF file (*)
* ```-infa <FILE>``` input chromosome FASTA file (*)
* ```-outfa <FILE>``` outputs intermediate index fasta (*)
* ```-ind <FILE>``` outputs the final index file for quantification (*)
* ```-rl <INT>``` target read length (*) 
* ```-tmp <FOLD>``` temporary folder (defaults to the current folder)
* ```-Mc <INT>``` maximum number of fragments per gene (default 25000)
* ```-lgo <STR>``` large gene optimization strategy (default "11")

If input reads are of variable length (e.g. after trimming), we recommend setting the target read length (parameter -rl) to the length of the longest read in the sample.

The last parameter "-lgo" controls the extention of the catalog of known trascripts by novel "guessed" ones (step A). It takes as argument a string of up to two numbers from the set {0, 2, 1} where 0 is the least restrictive option, 2 is the middle ground and 1 is the most restrictive option. The less restrictive optimizations are, the more additional isoforms will be generated. Option 0 denotes the most comprehensive extension. It creates fragments by combining known donor and acceptor sites of the same gene (set <img src="https://render.githubusercontent.com/render/math?math=T^{ap}_{g}"> in the manuscript).  Option 2 only allows to combine splice sites that lie within the boundaries of a known transcript (set <img src="https://render.githubusercontent.com/render/math?math=T^{as}_{g}">). Finally, option 1 does not create any novel fragments but only uses annotated transcripts. If the provided argument consists of 2 numbers, fragments are created initially using the strategy specified by the first number, gradually reducing to the strategy specified by the second if -Mc is exceeded.  

The following options further restrict the set of generated fragments: 

* ```-mil <INT>``` minimum generated intron length (default unrestricted)
* ```-Mil <INT>``` maximum generated intron length (default unrestricted)
* ```-mel <INT>``` minimum generated exon length (default unrestricted)
* ```-Mel <INT>``` maximum generated exon length (default unrestricted) 
* ```-exs <INT>``` maximum number of skipped exons (default unrestricted)
* ```-Moh <INT>``` maximum overhang on fragment ends (default unrestricted)
* ```--pmrna``` generate fragments from known premrna (default off)
* ```--fse``` generate all singleton subexons despite the limits defined by -Mc (default off)

Optionally, all generated transcript fragments can be output in GTF format.

* ```-outgtf <FILE>``` outputs a gtf file describing fragments

Here is an example call to fortuna's indexing step: 

``` ./fortuna --index -ingtf seq/a.gtf -infa seq/a.fa -outfa ind/a.fa -ind ind/a -rl 75 -tmp tmp/ -lgo "22" -exs 8 --fse ```


### Quantification

During the quantification step ("--quant"), reads will be pseudoaligned to fragments using kallisto (step D). From pseudoalignments fortuna derives signature counts and annotates all novel events (step E2). Optionally, unaligned reads will be written to disk and can be realigned across novel splice sites in the refinement step (below).

* ```-rl <INT>``` input read length (*)
* ```-gtf <FILE>``` preprocessed GTF file (*)
* ```-infq <FILE>``` FASTQ file with reads (*)
* ```-ind <FILE>``` index computed in previous step (*)
* ```-fa <FILE>``` transcriptome FASTA file
* ```-outfq <FILE>``` output FASTQ file with misaligned reads for the refinement step
* ```-cnt <FILE>``` output file with signature counts
* ```-alt <FILE>``` output file with novel event annotation
* ```-bam <FILE>``` output pseudoalignment BAM file
* ```-miss <INT>``` maximum number of allowed mismatches during alignment; reads exceeding this number will be considered unaligned
* ```-thr <INT>``` number of threads
* ```--p ``` outputs an intermediate file which is used to quantify paired-end reads

As in the indexing step, we recommend setting -rl to the length of the longest read.

Signature counts are output in a tab separated file specified by --cnt.  It contains in each line the number of reads with a specific mapping signature. Mapping signatures are specified by a sequence of subexon IDs taken from the segmented GTF file, separated by symbol "|". In the following example, the alignment of 7 reads overlap subexons a,b,c, and d:

``` a|b|c|d| 7 ```

The tab separated event file contains in each line the information of a novel event:

 ``` chromosome left-flank right-flank transcripts event-type count gene-id ```

Each event has one of the following types: ES (exon skipping), IR (intron retention), AA (novel acceptor), AD (novel donor), AP (novel pair), IE (intron-in-exon) and XX (unknown). Field transcripts is a comma separated list of annotated transcript IDs in which the event was observed. Field count is the number of reads supporting the event. In case of ES, AA, AD, AP and IE, left and right flanks are the genomic coordinates of donor and acceptor sites of the novel intron, while in the case of IR and XX, 5' and 3' coordinates of the region of the read supporting the novel exonic segment are reported.

Here is an example call of the quantification step: 

 ``` ./fortuna --quant -rl 75 -gtf seq/a.gtf -infq dat/a.fq -ind ind/a -fa seq/a.fa -outfq dat/a.unaligned.fq -cnt res/cnt.tsv -alt res/alt.tsv -miss 3 -thr 4 ```


### Refinement

Optionally ("--refine"), fortuna incorporates novel splice sites identified by any spliced aligner (such as STAR) from previously unmapped reads.

* ```-rl <INT>``` input read length (*)
* ```-gtf <FILE>``` preprocessed GTF file (*)
* ```-bam <FILE>``` genome BAM file with realigned reads (*)
* ```-outcnt <FILE>``` output file with refined signature counts 
* ```-outalt <FILE>``` output file with novel event annotation
* ```-ref <FILE>``` output refinement reference file
* ```-incnt <FILE>``` input file with signature counts
* ```-inalt <FILE>``` input file with novel event annotation
 
Input files specified by "-incnt" and "-inalt", if specified, will be augmented by the new splice sites taken from the BAM file and stored to files specified by "-outcnt" and "--outalt". The refinement reference file is a tab separated file which contains all newly created or refined subexons and introns. They are listed
with their IDs and genomic coordinates (chromosome, start, end). We recommend setting -rl to the length of the longest read.

Here is an example call to the refinement step:

 ``` ./fortuna --refine -rl 75 -gtf seq/a.gtf -bam res/star.bam -incnt res/cnt.tsv -outcnt res/cnt.refined.tsv -inalt res/alt.tsv -outalt res/alt.refined.tsv -ref res/ref.txt ```
 
 
### Paired-end reads

fortuna can use two intermediate single-ended outputs obtained by ```--quant``` and merge them into a paired-end count file. Command ```--pairs```, which does that, has following inputs:

* ```-lf <FILE>``` input file with left reads (*)
* ```-rf <FILE>``` output file with subexon counts (*)
* ```-out <FILE>``` output file with subexon counts (*)
* ```-lp <STRING>``` left file read suffix
* ```-rp <STRING>``` right file read suffix

Sometimes, left and right read names have specific suffixes at their respective ends, e.g. r1/1 (left) and r1/2 (right) could represent a read pair. Optional arguments ```--lp``` and ```--rp``` are used to filter out postfixes ("/1" and "/2" in the example).

Here is an example call:
 
 ``` ./fortuna --pairs -lf res/left.cnt -rf res/right.cnt -out res/paired.cnt ```


### Transform

fortuna also allows to transform (```--trans```) signature counts to subexon counts as used by DEXSeq for alternative downstream analysis. It requires two parameters:
* ```-incnt <FILE>``` input file with signature counts (*)
* ```-outcnt <FILE>``` output file with subexon counts (*)

Here is an example conversion of counts:
 
 ``` ./fortuna --trans -incnt res/a.cnt -outcnt res/a.subexon.cnt ```

 
## References
 [1] A. Dobin et al. Star: ultrafast universal rna-seq aligner. Bioinformatics, 29(1):15–21, 2013.
 
 [2] N. L. Bray, H. Pimentel, P. Melsted, and L. Pachter. Near-optimal probabilistic rna-seq quantification. Nature Biotechnology, 34:525–527, 2016.
 
 [3] S. Anders, A. Reyes, and W. Huber. Detecting differential usage of exons from rna-seq data. Genome Research, 22:2008–2017, 2012.
 
 [4] K. Jaganathan et al. Predicting splicing from primary sequence with deep learning. Cell, 176(3):535–548, 2019.
 

## Developer
* Luka Borozan (lborozan@mathos.hr)
