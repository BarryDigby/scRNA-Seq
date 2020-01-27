#!/usr/bin/env nextflow

// V2 chemistry
// params.chemistry = '10xv2'
// params.barcode_whitelist = "/data/MA5112/Practicals/scRNA-Seq/Assets/10xv2_whitelist.txt"

// V3 chemistry
params.chemistry = '10xv3'
params.barcode_whitelist = "/data/MA5112/Practicals/scRNA-Seq/Assets/10xv3_whitelist.txt"

Channel
        .fromPath(params.barcode_whitelist)
        .set{barcode_whitelist_bustools }

// Read in reference fasta file

params.ref_fa = "/data/MA5112/Practicals/scRNA-Seq/Reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
ref_fasta = files( params.ref_fa )

params.ref_gtf = "/data/MA5112/Practicals/scRNA-Seq/Reference/Homo_sapiens.GRCh38.99.gtf"
Channel 
	.fromPath( params.ref_gtf )
	.into{ ref_gtf; ref_gtf1 }

params.ref_gff3 = "/data/MA5112/Practicals/scRNA-Seq/Reference/Homo_sapiens.GRCh38.99.gff3"
ref_gff3 = files( params.ref_gff3 )

// Read in fastq files separately

params.read1 = "/data/MA5112/Practicals/scRNA-Seq/Data/5k_pbmc_protein_v3_fastqs/5k_pbmc_protein_v3_gex_fastqs/*_R1*"
Channel
	.fromPath( params.read1 )
	.set { read1_ch }

params.read2 = "/data/MA5112/Practicals/scRNA-Seq/Data/5k_pbmc_protein_v3_fastqs/5k_pbmc_protein_v3_gex_fastqs/*_R2*"
Channel
	.fromPath( params.read2 ) 
	.set { read2_ch }




process make_transcriptome{
	publishDir "Reference/", mode:'copy'

	input:
	file ref_fasta
	file gff from ref_gff3

	output:
	file "gffread.transcriptome.fa" into extracted_transcriptome

	script:
	"""
	gffread -F -T $gff -w "gffread.transcriptome.fa" -g $ref_fasta
	"""
}
	

process index{
	publishDir "Reference/", mode:'copy'

	input:
	file transcriptome from extracted_transcriptome

	output:
	file "Homo_sapiens.cDNA.idx" into indexed_transcriptome

	script:
	"""
	kallisto index -i Homo_sapiens.cDNA.idx ${transcriptome}
	"""
}


process tx2gene{
	echo true
	publishDir "Assets/", mode:'copy'

	input:
	file gtf from ref_gtf1

	output:
	file "transcripts_to_genes.txt" into tx2gene

	script:
	"""
	pypath="/data/MA5112/Practicals/scRNA-Seq/tx2gene.py"
	cat $gtf | python \$pypath > transcripts_to_genes.txt
	"""
} 


combined = read1_ch.merge(read2_ch)
combined_flat = combined.flatten().collect()



process kallisto{
	echo true
	publishDir "Analysis_output/raw", mode:'copy'

	input:
	file(index) from indexed_transcriptome.collect()
	file(reads) from combined_flat
	
	output:
	file "bus_output" into (kallisto_bus_sort, kallisto_bus_sort1)
	file "kallisto.log" into output
	script:
	"""
	kallisto bus \
	-i $index \
	-o bus_output/ \
	-x ${params.chemistry} \
	-t 8 \
	$reads | tee kallisto.log
	"""
}


process fix{
	publishDir "Analysis_output/raw", mode:'copy'

	input:
	file bus from kallisto_bus_sort

	output:
	file "transcripts2.txt" into tx

	script:
	"""
	cut -d':' -f2 ${bus}/transcripts.txt > transcripts2.txt
	"""
}



process bustools_correct_sort{
	publishDir "Analysis_output/sorted", mode:'copy'

	input:
	file bus from kallisto_bus_sort1
	file whitelist from barcode_whitelist_bustools

	output:
	file bus into (kallisto_corrected, kallisto_metrics)

	script:
	correct = "bustools correct -w $whitelist -o ${bus}/output.corrected.bus ${bus}/output.bus" 
	sort_file = "${bus}/output.corrected.bus"
	"""
	$correct
	mkdir -p tmp
	bustools sort \
	-T tmp/ \
	-t 8 \
	-o ${bus}/output.corrected.sort.bus \
	$sort_file
	"""
} 



process bustools_count{
	publishDir "Analysis_output/counts", mode:'copy'

	input:
	file bus from kallisto_corrected
	file t2g from tx2gene.collect()
	file t from tx
	
	output:
	file "${bus}_genecount"

	script:
	"""
	mkdir -p ${bus}_genecount
	
	bustools count \
	-o ${bus}_genecount/gene \
	-g $t2g \
	-e ${bus}/matrix.ec \
	-t $t \
	--genecounts \
	${bus}/output.corrected.sort.bus
	"""
}

process metrics{
	publishDir "Analysis_output/metrics", mode:'copy'

	input:
	file bus from kallisto_metrics

	output:
	file "${bus}.json"

	script:
	"""
	bustools inspect \
	-o ${bus}.json \
	${bus}/output.corrected.sort.bus
	"""
}
