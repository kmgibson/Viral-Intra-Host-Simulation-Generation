# Scripts to generate simulated intra-host NGS reads.

### These were used for the haplotype comparison project.

See article on [Infection, Genetics and Evolution](https://doi.org/10.1016/j.meegid.2020.104277).

<br>

_**We have plans to convert all these scripts into a snakemake workflow. For now, we have just documented the scripts and structure of it all.**_

General workflow is:

1. Create directory structure with parameter file and start sequence (MRCA).
2. Simulate intra-host populations.
3. Simulate reads from the population.

--- 
 
*Authors of scripts: Keylie M. Gibson and Matthew L. Bendall <br>*
*Copyright (C) 2019 Keylie M. Gibson and Matthew L. Bendall*

---
---

<br>

##########################################################################################
## Step 1: Create directory structure
##########################################################################################

```
#!/usr/bin/env bash

## Mutation Rate
for m in "1e-3" "1e-4" "1e-5" "1e-6" "1e-7" "1e-8" "3e-4" "3e-5" "3e-8" "5e-3" "5e-4" "5e-5" "5e-6" "5e-7" "5e-8"; do
    ## Sample size
    for i in "100" "200" "400" "600" "650" "700" "900" "1100" "1300" "1500" "1600" "1800" "2000"; do
        ## Effective Population Size
        for j in "500" "1000" "2500" "5000" "7500" "10000"; do 
            mkdir -p u${m}_s${i}_Ne${j}
            ## Copy the MRCA into each parameter directory
            cp seqGMRCA u${m}_s${i}_Ne${j}/seqGMRCA
            ## Copy the parameter config file into each directory and change the needed parameters
            cat paramtemplate.txt | sed "s/XXMXX/$m/" | sed "s/XXSSXX/$i/" | sed "s/XXNEXX/$j/" > u${m}_s${i}_Ne${j}/parameters
        done
    done
done

## Create a list that contains all the directories (i.e., parameter names)
ls -d u* > list.txt
```

<br>
The directory structure looks like such:

```
./viral_simulation
	|
	├── scr
	|	└── CoalEvol7.3.5
	|
	├── list.txt
	|
	├── paramtemplate.txt
	|
	├── seqGMRCA
	|
	├── count_haps.py
	|
	├── calc_avg_hap.py
	|
	├── u1e-3_s100_Ne1000
	|	├── parameters
	|	└── seqGMRCA
	|
	├── u1e-3_s100_Ne10000
	|	├── parameters
	|	└── seqGMRCA
	|
	├── u1e-3_s100_Ne10000
	|	├── parameters
	|	└── seqGMRCA
	|
	├── u1e-3_s100_Ne2500
	|	├── parameters
	|	└── seqGMRCA
	|
	|
```



##########################################################################################
## Step 2: Run CoalEvol
##########################################################################################

### Running CoalEvol sequence simulator on each parameter for 5 replicates. 


```
!/usr/bin/env bash

for name in $(cat list.txt)
    cd $name
    for count in `seq 1 5`; do
		echo ${count}
		../../src/CoalEvol7.3.5
		Newname="replicate_${count}"
		echo $Newname
		mkdir -p $Newname
		mv Results/ $Newname/
    done
done

```

<br>
The directory structure will look like this after:

```
./viral_simulation
	|
	├── scr
	|	└── CoalEvol7.3.5
	|
	├── list.txt
	|
	├── paramtemplate.txt
	|
	├── seqGMRCA
	|
	├── count_haps.py
	|
	├── calc_avg_hap.py
	|
	├── u1e-3_s100_Ne1000
	|	├── parameters
	|	├── seqGMRCA
	|	├── replicate_1
	|	|	└── Results
	|	├── replicate_2
	|	|	└── Results
	|	├── replicate_3
	|	|	└── Results
	|	├── replicate_4
	|	|	└── Results
	|	└── replicate_5
	|		└── Results
	|
	|

```



##########################################################################################
## Step 3: Count number of haplotypes
##########################################################################################

We need to make sure every sequence has it's own directory to store the read files with the fasta file.

```
!/usr/bin/env bash

##--- Making New Directories to complete ART in
## Parameter
for par in $name; do #parameter
	    ## Sequence
	    for seq in $name/Results/seq*; do
           fseqname=${seq##*/} #seq.fas
           seqname=${fseqname%.fas} #seq
           ## Make directory with seqname
           mkdir -p $name/$seqname
           ## Moves sequence into new directory
           mv $seq $name/$seqname/
	    done
    done
done
```
<br>

File name = `count_haps.py`

#### This script counts the number of haplotypes in each sequence file.

```
!/usr/bin/env python2

from Bio import SeqIO
from collections import Counter
import argparse

# run like:
# python count_haps.py -input [input_fasta_file] -output [output_file_name.fasta] > count_haps.txt
# you need a count_haps.txt file for each sequence directory in each replicate in each parameter directory

parser = argparse.ArgumentParser()
parser.add_argument("-input", help="input file")
parser.add_argument("-output", help="output file")
args = parser.parse_args()


counts = Counter()
for s in SeqIO.parse(args.input, "fasta"):
    counts[str(s.seq)] += 1

num = 1

with open(args.output, "w") as outh:
    seq = args.input.split(".")[0]
    for sequence, count in counts.most_common():
        print "%s\t%d" % (sequence, count)
        if count > 1:
            print >> outh, ">duplicate_from_args.input_%s" % seq, num
            print >> outh, "%s" % sequence
            num += 1

```
<br>
The directory structure will look like this after:

```
./viral_simulation
	|
	├── scr
	|	└── CoalEvol7.3.5
	|
	├── list.txt
	|
	├── paramtemplate.txt
	|
	├── seqGMRCA
	|
	├── count_haps.py
	|
	├── calc_avg_hap.py
	|
	├── u1e-3_s100_Ne1000
	|	├── parameters
	|	├── seqGMRCA
	|	├── replicate_1
	|	|	├── Results
	|	|	├── sequences00001
	|	|	|	├── sequences00001.fas
	|	|	|	└── count_haps.py
	|	|	├── sequences00002
	|	|	|	├── sequences00002.fas
	|	|	|	└── count_haps.py
	|	|	├── .....
	|	|	|
	
```
            
##########################################################################################
## Step 4: Calculate average number of haplotypes per parameter
##########################################################################################

#### Create `hap_summary.txt` for each replicate in each parameter.

```
!/usr/bin/env bash

for p in u*; do
    for rep in $p/replicate_?; do
        for s in $rep/sequences*; do
            hapcount=$(cat $s/count_haps.txt | wc -l)
            echo -e "${p}\t$(basename $s)\t${hapcount}" >> $rep/hap_summary.txt
        done
    done
done

```

<br>

#### Count the average haplotype count for each parameter set.

File = `calc_avg_hap.py`

```
!/usr/bin/env python2

from Bio import SeqIO
from collections import Counter
from glob import glob

s2000 = glob('u*s2000*/Results/hap_summary.txt')

outh = open('avg_num_haps_summary.txt', 'w')

print >>outh, "Mutation_rate\tEffective_pop_size\tAvg_num_haps"

for p in s2000:
    mu = p.split('/')[0].split('_')[0]
    Ne = p.split('/')[0].split('_')[2]
    num_seq_counts = Counter()
    num_haps = []
    for s in open(p, 'rU'):
        seq = s.split('\t')[1]
        hap = int(s.split('\t')[2])
        num_seq_counts[seq] += 1
        num_haps.append(hap)
    if len(num_seq_counts) == len(num_haps):
        hap_sum = sum(num_haps)
        avg_hap = int(hap_sum/(len(num_haps)))
        print >>outh, "%s\t%s\t%s" % (mu, Ne, avg_hap)
    else:
        print "ERROR: %s" % p
```

<br>
The directory structure will look like this after:

```
./viral_simulation
	|
	├── scr
	|	└── CoalEvol7.3.5
	|
	├── list.txt
	|
	├── paramtemplate.txt
	|
	├── seqGMRCA
	|
	├── count_haps.py
	|
	├── calc_avg_hap.py
	|
	├── count_haps.txt
	|
	├── u1e-3_s100_Ne1000
	|	├── parameters
	|	├── seqGMRCA
	|	├── replicate_1
	|	|	├── hap_summary.txt
	|	|	├── Results
	|	|	├── sequences00001
	|	|	|	├── sequences00001.fas
	|	|	|	└── count_haps.py
	|	|	├── sequences00002
	|	|	|	├── sequences00002.fas
	|	|	|	└── count_haps.py
	|	|	├── .....
	|	|	|
	
```



##########################################################################################
## Step 5: Runs ART read simulator on the output files created from the CoalEvol.sh script IN PARALLEL

##########################################################################################

`$name` is a line from `list.txt`. We run this in array format in SLURM.

```
!/usr/bin/env bash

module load art # loads the art module
module load parallel # loads parallelization module
module load samtools/1.3.1 # loads samtools


##--- Running ART in Parallel and compressing files
function runart () {
    local f=$1
    local nseq=$(wc -l < $f)
    head -n $(( ($nseq - 1) * 2 )) $f > $(dirname $f)/tmp.fas
    ####### Run art. Adjust parameters for ART here.
    art_illumina -ss MSv1 -sam -i $(dirname $f)/tmp.fas -ef -rs 04142017 -p  -l 150 -c 100 -m 215 -s 120 -o $(dirname $f)/read
    #######
    du -csh $(dirname $f)/*
    rm -f $(dirname $f)/tmp.fas && echo "removed $(dirname $f)/tmp.fas"
    rm -f $(dirname $f)/read1.aln && echo "removed $(dirname $f)/read1.aln"
    rm -f $(dirname $f)/read2.aln && echo "removed $(dirname $f)/read2.aln"
    echo "compressing read.sam" &&\
      samtools view -b $(dirname $f)/read.sam > $(dirname $f)/read.bam &&\
      rm -f $(dirname $f)/read.sam && echo "removed $(dirname $f)/read.sam"
    echo "compressing read_errFree.sam" &&\
      samtools view -b $(dirname $f)/read_errFree.sam > $(dirname $f)/read_errFree.bam &&\
      rm -f $(dirname $f)/read_errFree.sam && echo "removed $(dirname $f)/read_errFree.sam"
    echo "compressing $(dirname $f)/read1.fq" &&\
        gzip -f $(dirname $f)/read1.fq
    echo "compressing $(dirname $f)/read2.fq" &&\
        gzip -f $(dirname $f)/read2.fq
    du -csh $(dirname $f)/*
}
export -f runart

##--- Runs are in parallel on all the parameters
for d in $name/replicate_1/sequences*/sequences*.fas; do echo $d done | parallel -j $(nproc) runart
```

<br>
The directory structure will look like this after:

```
./viral_simulation
	|
	├── scr
	|	└── CoalEvol7.3.5
	|
	├── list.txt
	|
	├── paramtemplate.txt
	|
	├── seqGMRCA
	|
	├── count_haps.py
	|
	├── calc_avg_hap.py
	|
	├── count_haps.txt
	|
	├── u1e-3_s100_Ne1000
	|	├── parameters
	|	├── seqGMRCA
	|	├── replicate_1
	|	|	├── hap_summary.txt
	|	|	├── Results
	|	|	├── sequences00001
	|	|	|	├── sequences00001.fas
	|	|	|	├── count_haps.py
	|	|	|	├── read_errFree.bam
	|	|	|	├── read.bam
	|	|	|	├── read1.fq.gz
	|	|	|	└── read2.fq.gz
	|	|	├── sequences00002
	|	|	|	├── sequences00002.fas
	|	|	|	├── count_haps.py
	|	|	|	├── read_errFree.bam
	|	|	|	├── read.bam
	|	|	|	├── read1.fq.gz
	|	|	|	└── read2.fq.gz
	|	|	├── .....
	|	|	|
	
```



##########################################################################################
## Step 6: Run the simulated reads through HAPHPIPE.
##########################################################################################

For more information regarding HAPHPIPE, see [github](https://github.com/gwcbi/haphpipe).

Each pair of reads for each sample can be run through one of the two pipelines available within HAPHPIPE. We used the reference-based pipeline, `haphpipe_assemble_02`. Following the HAPHPIPE pipeline, we put the each sample through the various haplotype reconstruction tools. 


`$name` is a line from `list.txt`. We run this in array format in SLURM. Here is an example of a single sample running through `haphpipe_assemble_02`. 

```
!/usr/bin/env bash

$ haphpipe_assemble_02 -h
USAGE:
haphpipe_assemble_02 [read1] [read2] [amplicons_fasta] [samp_id] <outdir>

----- HAPHPIPE assembly pipeline 02 -----

This pipeline implements amplicon assembly using a reference-based approach.
Reads are error-corrected and aligned to provided amplicon reference with up to
five refinement steps.

Input:
read1:             Fastq file for read 1. May be compressed (.gz)
read2:             Fastq file for read 2. May be compressed (.gz)
amplicons_fasta:   Amplicon reference sequence (fasta)
samp_id:           Sample ID
outdir:            Output directory (default is sample_dir/haphpipe_assemble_02)

# This is an example with the parameter set used in the structure examples.
$ haphpipe_assemble_02 \
 u1e-3_s100_Ne1000/replicate_1/sequences00001/read1.fq.gz \
 u1e-3_s100_Ne1000/replicate_1/sequences00001/read2.fq.gz \
 seqGMRCA \
 u1e-3_s100_Ne1000_sequences00001

```

You will likely only care about the reads, either corrected, trimmed or raw paired reads, `final.fna`, `final.bam` for inputs into the haplotype reconstruction programs.

<br>
The directory structure will look like this after:

```
./viral_simulation
	|
	├── scr
	|	└── CoalEvol7.3.5
	|
	├── list.txt
	|
	├── paramtemplate.txt
	|
	├── seqGMRCA
	|
	├── count_haps.py
	|
	├── calc_avg_hap.py
	|
	├── count_haps.txt
	|
	├── u1e-3_s100_Ne1000
	|	├── parameters
	|	├── seqGMRCA
	|	├── replicate_1
	|	|	├── hap_summary.txt
	|	|	├── Results
	|	|	├── sequences00001
	|	|	|	├── sequences00001.fas
	|	|	|	├── count_haps.py
	|	|	|	├── read_errFree.bam
	|	|	|	├── read.bam
	|	|	|	├── read1.fq.gz
	|	|	|	├── read2.fq.gz
	|	|	|	└── haphpipe_assemble_02
	|	|	|		├── corrected_1.fastq
	|	|	|		├── corrected_2.fastq
	|	|	|		├── corrected_U.fastq
	|	|	|		├── final.bam
	|	|	|		├── final.bam.bai
	|	|	|		├── final_bt2.out
	|	|	|		├── final.fna
	|	|	|		├── final.vcf.gz
	|	|	|		├── final.vcf.gz.tbi
	|	|	|		├── haphpipe.out
	|	|	|		├── refined.01.fna
	|	|	|		├── refined.02.fna
	|	|	|		├── refined.03.fna
	|	|	|		├── refined.04.fna
	|	|	|		├── refined_bt2.01.out
	|	|	|		├── refined_bt2.02.out
	|	|	|		├── refined_bt2.03.out
	|	|	|		├── refined_bt2.04.out
	|	|	|		├── refined_bt2.out
	|	|	|		├── refined.fna
	|	|	|		├── refined_summary.out
	|	|	|		├── trimmed_1.fastq
	|	|	|		├── trimmed_2.fastq
	|	|	|		├── trimmed_U.fastq
	|	|	|		└── trimmomatic_summary.out
	|	|	├── .....
	|	|	|
	
```





