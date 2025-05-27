#!usr/bin/env nextflow
nextflow.enable.dsl=2

workflow {
	def read_pairs= Channel.fromFilePairs(params.reads, checkIfExists:true)
	read_pairs.view()
	fastqc(read_pairs)
	star_index_out = build_star_index(tuple("hg38", file(params.genome_fasta), file(params.gtf)))
	salmon_index_out = build_salmon_index(tuple("hg38", file(params.transcriptome_fasta)))
	align = align_star(read_pairs, star_index_out)
	counts = feature_counts(align)
	tpms = salmon_quant(read_pairs, salmon_index_out)
}

process fastqc {
	publishDir "${params.outdir}/fastqc"
	tag "$sra_id"

	input:
	tuple val(sra_id), path(reads)

	output:
	path("*.html")
	path("*.zip")

	script:
	"""
	fastqc -t 8 ${reads[0]} ${reads[1]} -o .
	"""
}

// Genome alignment with STAR
process align_star {
	publishDir "$params.outdir/align"
	tag "$sra_id"

	input:
	tuple val(sra_id), path(reads)
	path(index_dir)

	output:
	tuple val(sra_id), path("*.bam")

	script:
	"""
	STAR --runThreadN 16 \\
	--genomeDir ${index_dir} \\
	--readFilesIn ${reads[0]} ${reads[1]} \\
	--readFilesCommand zcat \\
	--outFileNamePrefix ${sra_id}_ \\
	--outSAMtype BAM SortedByCoordinate \\
	--twopassMode Basic
	"""
}

// FeatureCounts read summarization
process feature_counts {
	tag "$sra_id"

	input:
	tuple val(sra_id), path(bam)

	output:
	path("*.txt")

	script:
	"""
	featureCounts -T 16 \\
	-a ${params.gtf} \\
	-o ${sra_id}_counts.txt \\
	-g gene_id \\
	-t exon -p -B -C \\
	${bam}
	"""
}

// TPM quantification using Salmon
process salmon_quant {
	tag "$sra_id"

	input:
	tuple val(sra_id), path(reads)
	path(index_dir)

	output:
	path("${sra_id}_salmon/quant.sf")

	script:
	"""
	salmon quant -i ${index_dir} -l A \\
	-1 ${reads[0]} -2 ${reads[1]} \\
	-p 8 --validateMappings --gcBias \\
	-o ${sra_id}_salmon
	"""
}

process build_star_index {
    publishDir "${params.outdir}/star_index"
    tag "$prefix"

    input:
    tuple val(prefix), path(genome_fasta), path(gtf_file)

    output:
    path("star_index", emit: index_dir)

    script:
    """
    mkdir -p star_index
    STAR --runThreadN 16 \\
         --runMode genomeGenerate \\
         --genomeDir star_index \\
         --genomeFastaFiles ${genome_fasta} \\
         --sjdbGTFfile ${gtf_file} \\
         --sjdbOverhang 48
    """
}

process build_salmon_index {
    publishDir "${params.outdir}"
    tag "$prefix"

    input:
    tuple val(prefix), path(transcriptome_fasta)

    output:
    path("salmon_index", emit: index_dir)

    script:
    """
    mkdir -p salmon_index
    salmon index -t ${transcriptome_fasta} -i salmon_index
    """
}
