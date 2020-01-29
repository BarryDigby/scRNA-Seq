# Kallisto | Bustools: An Unbridled Rage

This tutorial is designed to guide the user through a `kallisto | bustools` pipeline to generate the 3 files required for downstream analysis: `matrix.mtx`, `genes.tsv` and `barcodes.txt`. As of writing this tutorial (28th January 2020) the Pachter lab is pushing it's new tool `kb-python` which is essentially a wrapper script that executes `kallisto | bustools` for the user with some python scripts under the hood for parsing output into the correct format. That being said, the [bustools github page](https://github.com/BUStools/BUS_notebooks_python/tree/master/dataset-notebooks) is home to numerous `ipynb` files for conducting scRNA-Seq analysis using `kallisto | bustools` which are inconsistent and result in incomplete output files. Perhaps more aggregious is the [getting started](https://www.kallistobus.tools/getting_started) bustools tutorial (directed from the [bustools website](https://bustools.github.io/) which leaves the user with unannotated gene ID outputs and no counsel on how to rectify the problem. 

Firstly, I will illustrate why the tutorials are lacking before providing a solution to the problems identified. 

*** 

## The "getting started tutorial". 
My direcory has been set up with all the the prerequisite files:

```
barry@NUIG:~/scRNA$ tree */
data/
├── SRR8599150_S1_L001_R1_001.fastq.gz
└── SRR8599150_S1_L001_R2_001.fastq.gz
files/
├── 10xv2_whitelist.txt
└── transcripts_to_genes.txt
reference/
└── Mus_musculus.GRCm38.cdna.all.fa

0 directories, 5 files
```

To run this tutorial, I will be using a singularity container. To build the container, download the `scRNA.def` file available in the **<> code** section of this repository and run `sudo singularity build scRNA.simg scRNA.def`. (install singularity at [this link](https://singularity.lbl.gov/install-linux). 

#### 2. Build an Index
```
kallisto index -i reference/Mus_musculus.GRCm38.cdna.all.fa.idx -k 31 reference/Mus_musculus.GRCm38.cdna.all.fa

[build] loading fasta file reference/Mus_musculus.GRCm38.cdna.all.fa
[build] k-mer length: 31
[build] warning: clipped off poly-A tail (longer than 10)
        from 641 target sequences
[build] warning: replaced 3 non-ACGUT characters in the input sequence
        with pseudorandom nucleotides
[build] counting k-mers ... done.
[build] building target de Bruijn graph ...  done 
[build] creating equivalence classes ...  done
[build] target de Bruijn graph has 734746 contigs and contains 100614952 k-mers 
```

Check: 
```
barry@NUIG:~/scRNA$ tree *
data
├── SRR8599150_S1_L001_R1_001.fastq.gz
└── SRR8599150_S1_L001_R2_001.fastq.gz
files
├── 10xv2_whitelist.txt
└── transcripts_to_genes.txt
reference
├── Mus_musculus.GRCm38.cdna.all.fa
└── Mus_musculus.GRCm38.cdna.all.fa.idx
```

#### 3. Run Kallisto
```
kallisto bus -i reference/Mus_musculus.GRCm38.cdna.all.fa.idx -o bus_output/ -x 10xv2 -t 8 data/SRR8599150_S1_L001_R1_001.fastq.gz data/SRR8599150_S1_L001_R2_001.fastq.gz

[index] k-mer length: 31
[index] number of targets: 118,489
[index] number of k-mers: 100,614,952
[index] number of equivalence classes: 433,624
[quant] will process sample 1: data/SRR8599150_S1_L001_R1_001.fastq.gz
                               data/SRR8599150_S1_L001_R2_001.fastq.gz
[quant] finding pseudoalignments for the reads ... done
[quant] processed 8,860,361 reads, 3,431,849 reads pseudoaligned
```
check
```
barry@NUIG:~/scRNA$ tree *
bus_output
├── matrix.ec
├── output.bus
├── run_info.json
└── transcripts.txt
data
├── SRR8599150_S1_L001_R1_001.fastq.gz
└── SRR8599150_S1_L001_R2_001.fastq.gz
files
├── 10xv2_whitelist.txt
└── transcripts_to_genes.txt
reference
├── Mus_musculus.GRCm38.cdna.all.fa
└── Mus_musculus.GRCm38.cdna.all.fa.idx
```

#### Run Bustools
```
bustools correct -w ../files/10xv2_whitelist.txt -p output.bus | bustools sort -T tmp/ -t8 -p - | bustools count -o genecount/genes -g ../files/transcripts_to_genes.txt -e matrix.ec -t transcripts.txt --genecounts -

Found 737280 barcodes in the whitelist
Number of hamming dist 1 barcodes = 20550336
Processed 3431849 bus records
In whitelist = 3281671
Corrected = 36927
Uncorrected = 113251
Read in 3318598 BUS records
```
check
```
barry@NUIG:~/scRNA$ tree *
bus_output
├── genecount
│   ├── genes
│   ├── genes.barcodes.txt
│   ├── genes.genes.txt
│   └── genes.mtx
├── matrix.ec
├── output.bus
├── run_info.json
├── tmp
└── transcripts.txt
data
├── SRR8599150_S1_L001_R1_001.fastq.gz
└── SRR8599150_S1_L001_R2_001.fastq.gz
files
├── 10xv2_whitelist.txt
└── transcripts_to_genes.txt
reference
├── Mus_musculus.GRCm38.cdna.all.fa
└── Mus_musculus.GRCm38.cdna.all.fa.idx
```

## The "genes.genes.txt" problem
Nice! We have followed the tutorial to a tee and all of the output files seem kosher. However the main problem (and the theme of this tutorial) is that the `genes.gene.txt` file has only got ENG gene ID's present in the file:
```
barry@NUIG:~/scRNA$ head bus_output/genecount/genes.genes.txt 
ENSMUSG00000037736.18
ENSMUSG00000029804.17
ENSMUSG00000056124.5
ENSMUSG00000052516.19
ENSMUSG00000031511.15
ENSMUSG00000037126.16
ENSMUSG00000034164.17
ENSMUSG00000044034.11
ENSMUSG00000025269.16
ENSMUSG00000030446.17
```
This file is used as the annotation file for the `anndata` object in `scanpy`. A key step to scRNA-Seq analysis is to filter out all of the mitochondrial genes in your cells. This is easily achieved in python when you have gene symbols, simply use the python `startswith("MT")` to flag each mitochondrial gene. The task is made far mroe complicated (and the output much less interpretable) if only ENS gene ID's are used. The same is true for marker gene analysis. 

So, after your first tutorial using `bustools`, you are stuck with output that is not properly annotated. This cannot be fixed by using `biomaRt` to annotate the file, as non unique gene symbols truncate the file. Attemting to load the truncated file into `scanpy` will result in an error, it must maintain the same length as bustools output originally defined. 

> As an introduction tutorial, it leaves far too many things unexplained and results in unannotated gene ID's. 

## The "transcripts_to_gene.txt" problem 
Intuitively, one can see that the `transcripts_to_genes.txt` file containing ENST ID's, ENSG ID's and gene symbols is how the output files are annotated - by using ENST as a key for ENSG + gene symbol values. Lets look at some of the online suggestions on how to generate one. 

Using the python code at this [notebook](https://github.com/BUStools/BUS_notebooks_python/blob/master/dataset-notebooks/10x_pbmc_1k_v3chem_python/10x_pbmc_1k_v3chem.ipynb) it instructs you to implement the script on a GTF file. After copying python code cells `ln[9]` and `ln[10]` into a script (and editing the with.open() line to newest GTF version .99) the output does not contain  gene symbols. 
```
barry@NUIG:~/scRNA$ python py.py 
Created transcript_to_gene.tsv file
barry@NUIG:~/scRNA$ head transcript_to_gene.tsv 
ENST00000651793.1	ENSG00000286217.2
ENST00000651769.1	ENSG00000100142.15
ENST00000527779.1	ENSG00000137692.12
ENST00000338863.12	ENSG00000007866.21
ENST00000535093.1	ENSG00000134759.14
ENST00000646354.1	ENSG00000185532.19
ENST00000346219.7	ENSG00000008226.20
ENST00000570899.1	ENSG00000040633.13
ENST00000557761.1	ENSG00000258690.1
ENST00000625998.2	ENSG00000120685.20
```

Further down the analysis, you will see the `scanpy` plots have gene symbols. This is because they (yet again pulling a file out of thin air) use this block of code `ln[29]` to annotate the object:
```
gene_names = {}
with open('../index/mart_export_human.txt') as f:
    f.readline()
    for line in f:
        g,t,gn = line.split()
        gene_names[g] = gn
 ```
 There is no mention of where they sourced, or how they generated the `mart_export_human.txt` file. 

## Issues thus far
 * Scripts provided do not make a comprehensive `transcripts_to_genes.txt` file. 
 * If kallisto aligns to the cDNA index, why does it need ncRNA information from the GTF?
 * `gene.genes.txt` file only has ENSG ID
 * The getting started tutorial and the `ipynb` use   `bustools count` and `bustools text` interchangeably. 
 * One tutorial uses `bustools correct`, another does not. 
 
## The solution:
To run a pipeline with successfully annotated output for scanpy analysis, you will need `fastq files`, `reference cDNA transcriptome` and a `10xv(X)_whitelist`. 

#### Index the transcriptome
This step of the analysis remains unchanged. Proceed to index your `.cDNA.all.fa` file. **N.B** this analysis uses Ensembl reference files. 

#### Generating a transcripts_to_gene.txt file
To generate the file, do not use the Ensembl GTF file. Instead, use the ensembl transcriptome file used as input in the indexing step. Use this bash command to generate a correct `transcriptrs_to_genes.txt` file:

``` bash
cat Homo_sapiens.GRCh38.cdna.all.fa | awk '{if($1~/>/)print $1"\t"$4"\t"$7}' > t2g.txt;
sed -i 's/>//g' t2g.txt; sed -i 's/gene://g' t2g.txt; sed -i 's/gene_symbol://g' t2g.txt
```

```
bdigby@compute05:/data/MA5112/Practicals/scRNA-Seq$ head Assets/t2g.txt 
ENST00000434970.2	ENSG00000237235.2	TRDD2
ENST00000415118.1	ENSG00000223997.1	TRDD1
ENST00000448914.1	ENSG00000228985.1	TRDD3
ENST00000631435.1	ENSG00000282253.1	TRBD1
ENST00000632684.1	ENSG00000282431.1	TRBD1
ENST00000390583.1	ENSG00000211923.1	IGHD3-10
ENST00000431440.2	ENSG00000232543.2	IGHD4-11
ENST00000632524.1	ENSG00000282455.1	IGHD7-27
ENST00000633009.1	ENSG00000282323.1	IGHD1-26
ENST00000634070.1	ENSG00000282724.1	IGHD6-25
```

#### Run kallisto bus
This step of the analysis remains unchanged. 

#### Run bustools
After trial and error and the help of Sarah Ennis, `bustools text` has emerged as the command to use for annotation of the output `gene.genes.tsv` file:
``` nextflow
bustools correct -w $whitelist -o ${bus}/output.corrected.bus ${bus}/output.bus
sort_file="${bus}/output.corrected.bus"
bustools sort -T tmp/ -t 8 -o ${bus}/output.sort.bus $sort_file
bustools text -o output.sort.txt ${bus}/output.sort.bus
```
`${bus}` refers to the directory where the output of `kallisto bus` resides. 
 
Run the following python script to format your `bustools text` output. Change the paths in the file to represent your `transcripts_to_gene.txt` file generated above, and to the `matrix.ec` file output by `bustools sort`. 
``` python 
#!/usr/bin/env python

gene_min = 200
gene_max = 10000

#setup working directory
import os
os.chdir("/data/MA5112/Practicals/scRNA-Seq/")

from subprocess import call
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys, collections

tr2g = {}
trlist = []
with open("Assets/t2g.txt") as f:
    for line in f:
        l = line.split()
        tr2g[l[0]] = l[1]
        trlist.append(l[0])

genes = list(set(tr2g[t] for t in tr2g))

# load equivalence classes
ecs = {}
with open('Analysis_output/sorted/bus_output/matrix.ec') as f:
    for line in f:
        l = line.split()
        ec = int(l[0])
        trs = [int(x) for x in l[1].split(',')]
        ecs[ec] = trs
        
def ec2g(ec):
    if ec in ecs:
        return list(set(tr2g[trlist[t]] for t in ecs[ec]))        
    else:
        return []

cell_gene = collections.defaultdict(lambda: collections.defaultdict(float))
pbar=None
pumi=None
with open('Analysis_output/counts/output.sort.txt') as f:
    gs = set()
    for line in f:
        l = line.split()
        barcode,umi,ec,count = line.split()
        ec = int(ec)
        
        if barcode == pbar:
            # same barcode
            if umi == pumi:
                # same UMI, let's update with intersection of genelist
                gl = ec2g(ec)
                gs.intersection_update(gl)
            else:
                # new UMI, process the previous gene set
                for g in gs:
                    cell_gene[barcode][g] += 1.0/len(gs)
                # record new umi, reset gene set
                pumi = umi
                gs = set(ec2g(ec))
        else:
            # work with previous gene list
            for g in gs:
                cell_gene[pbar][g] += 1.0/len(gs)
            
            if sum(cell_gene[pbar][g] for g in cell_gene[pbar]) < 10:
                del cell_gene[pbar]
            
            pbar = barcode
            pumi = umi
            
            gs = set(ec2g(ec))

    for g in gs:
        cell_gene[pbar][g] += 1.0/len(gs)
        
    if sum(cell_gene[pbar][g] for g in cell_gene[pbar]) < 10:
        del cell_gene[pbar]

barcode_hist = collections.defaultdict(int)
for barcode in cell_gene:
    cg = cell_gene[barcode]
    s = len([cg[g] for g in cg])
    barcode_hist[barcode] += s

#Output a gene count histogram
bcv = [x for b,x in barcode_hist.items() if x > gene_min and x < gene_max]
plt.switch_backend('agg')
fig = plt.figure()
ax = fig.add_subplot(111)
ax.hist(bcv,bins=100)
ax.set_title("Histogram")
plt.xlabel("number of genes detected")
plt.ylabel("number of barcodes")
fig.savefig('scanpy/gene_hist.png')

outfile = 'scanpy/matrix.mtx'

gene_to_id = dict((g,i+1) for i,g in enumerate(genes))
barcodes_to_use = [b for b,x in barcode_hist.items() if x > gene_min and x < gene_max]

num_entries = 0
for barcode in barcodes_to_use:
    num_entries += len([x for x in cell_gene[barcode].values() if x>0])

with open(outfile, 'w') as of:
    of.write('%%MatrixMarket matrix coordinate real general\n%\n')
    #number of genes
    of.write("%d %d %d\n"%(len(genes), len(barcodes_to_use), round(num_entries)))
    bcid = 0
    for barcode in barcodes_to_use:
        bcid += 1
        cg = cell_gene[barcode]
        gl = [(gene_to_id[g],cg[g]) for g in cg if cg[g] > 0]
        gl.sort()
        for x in gl:
            of.write("%d %d %f\n"%(x[0],bcid,x[1]))

gene_names = {}
with open("Assets/t2g.txt") as f:
    f.readline()
    for line in f:
        t,g,gn = line.split()
        gene_names[g] = gn

id_to_genes = dict((i,g) for (g,i) in gene_to_id.items())
gl = []
for i in range(1,len(genes)+1):
    g = id_to_genes[i]
    gid = g
#    gid = g[:g.find('.')]
    if gid in gene_names:
        gn = gene_names[gid]
    else:
        gn = ''
    gl.append((g,gn))

with open('scanpy/genes.tsv','w') as of:
    for g,gn in gl:
        of.write("%s\t%s\n"%(g,gn))
        
with open('scanpy/barcodes.tsv','w') as of:
    of.write('\n'.join(x + '' for x in barcodes_to_use))
    of.write('\n')
```

#### Output:
```
bdigby@lugh:/data/MA5112/Practicals/scRNA-Seq$ head scanpy/genes.tsv 
ENSG00000234100.4	OR14J1
ENSG00000282035.1	AC090958.4
ENSG00000288143.1	ABCD1P2
ENSG00000224557.7	HLA-DPB2
ENSG00000144010.9	TRIM43B
ENSG00000065320.9	NTN1
ENSG00000123600.19	METTL8
ENSG00000105520.10	PLPPR2
ENSG00000231022.1	RPS3AP9
ENSG00000229103.1	WASF5P
``` 

The `genes.tsv` file is now fully annotated and ready to be loaded into `scanpy` with `matrix.mtx` and `barcodes.txt` for single-cell analysis. 

