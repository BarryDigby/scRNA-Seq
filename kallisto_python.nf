#!/usr/bin/env nextflow

// V3 chemistry
params.chemistry = '10xv3'
params.barcode_whitelist = "/data/MA5112/Practicals/scRNA-Seq/Assets/10xv3_whitelist.txt"

Channel
        .fromPath(params.barcode_whitelist)
        .set{barcode_whitelist_bustools}

params.ref_fa = "/data/MA5112/Practicals/scRNA-Seq/Reference/Homo_sapiens.GRCh38.cdna.all.fa"
Channel
	.fromPath( params.ref_fa )
	.into{ ref_fa; ref_fa1 } 

params.read1 = "/data/MA5112/Practicals/scRNA-Seq/Data/5k_pbmc_protein_v3_fastqs/5k_pbmc_protein_v3_gex_fastqs/*_R1*"
Channel
	.fromPath( params.read1 )
	.set { read1_ch }

params.read2 = "/data/MA5112/Practicals/scRNA-Seq/Data/5k_pbmc_protein_v3_fastqs/5k_pbmc_protein_v3_gex_fastqs/*_R2*"
Channel
	.fromPath( params.read2 ) 
	.set { read2_ch }




process index{
	publishDir "Reference/", mode:'copy'

	input:
	file cDNA_fasta from ref_fa

	output:
	file "Homo_sapiens.cDNA.idx" into indexed_transcriptome

	script:
	"""
	kallisto index -i Homo_sapiens.cDNA.idx ${cDNA_fasta}
	"""
}


process tx2gene{
	echo true
	publishDir "Assets/", mode:'copy'

	input:
	file "Homo_sapiens.GRCh38.cdna.all.fa" from ref_fa1

	output:
	file "t2g.txt" into tx2gene

	shell:
	$/
	cat Homo_sapiens.GRCh38.cdna.all.fa | awk '{if($1~/>/)print $1"\t"$4"\t"$7}' > t2g.txt;
	sed -i 's/>//g' t2g.txt; sed -i 's/gene://g' t2g.txt; sed -i 's/gene_symbol://g' t2g.txt
	/$
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
	file "bus_output" into kallisto_bus_sort
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



process bustools_correct_sort{
        publishDir "Analysis_output/sorted", mode:'copy'

        input:
        file bus from kallisto_bus_sort
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
        -o ${bus}/output.sort.bus \
        $sort_file
        """
}
 


process bustools_count{
	publishDir "Analysis_output/counts", mode:'copy'

	input:
	file bus from kallisto_corrected
	file t2g from tx2gene.collect()

	output:
	file "output.sort.txt" into out

	script:
	"""
	bustools text \
	-o output.sort.txt \
	${bus}/output.sort.bus
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
	${bus}/output.sort.bus
	"""
}
