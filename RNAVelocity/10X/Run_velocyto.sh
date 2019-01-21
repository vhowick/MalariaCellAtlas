bash

CELLRANGER_OUTPUT=/lustre/scratch117/cellgen/team218/TA/Malaria_RNAVelocity/10XOutput/cellranger210_count_24726_4472STDY7251827_Pberghei_at10_20170115
CELLRANGER_OUTPUT=/lustre/scratch117/cellgen/team218/TA/Malaria_RNAVelocity/10XOutput/cellranger210_count_24726_4472STDY7251827_Pknowlesi_at10_20170115/
CELLRANGER_OUTPUT=/lustre/scratch117/cellgen/team218/TA/Malaria_RNAVelocity/10XOutput/cellranger210_count_24612_4472STDY7230041_3D7_Jan16v3-ensembl_37_transcriptome
BARCODE_FILE=$CELLRANGER_OUTPUT/filtered_gene_bc_matrices/Pberghei_at10_20170115/barcodes.tsv
BARCODE_FILE=$CELLRANGER_OUTPUT/filtered_gene_bc_matrices/Pknowlesi_at10_20170115/barcodes.tsv
BARCODE_FILE=$CELLRANGER_OUTPUT/filtered_gene_bc_matrices/3D7_Jan16v3/barcodes.tsv
#CELLRANGER_OUTPUT=/lustre/scratch117/cellgen/team218/TA/course/RNATranscriptomics/test001
#BARCODE_FILE=$CELLRANGER_OUTPUT/outs/filtered_gene_bc_matrices/mm10/barcodes.tsv
#GTF=/lustre/scratch117/cellgen/team218/TA/my_software/refdata-cellranger-mm10-1.2.0/genes/genes.gtf
GTF=/lustre/scratch117/cellgen/team218/TA/genomebuilding/Pberghei_forvelocyto10X.gtf
GTF=/lustre/scratch117/cellgen/team218/TA/genomebuilding/Pk_forvelocyto10X.gtf
GTF=/lustre/scratch117/cellgen/team218/TA/genomebuilding/Plasmodium_falciparum_forvelocyto10X.gtf

# echo "Pipestance completed successfully!" >> $CELLRANGER_OUTPUT/_log

conda activate base
# general running
#samtools sort -t CB -O BAM -o cellsorted_possorted_genome_bam.bam possorted_genome_bam.bam
#bsub -R"select[mem>10000] rusage[mem=10000]" -M10000 -o sort.%J.out -e sort.%J.err samtools sort -t CB -O BAM -o cellsorted_possorted_genome_bam.bam possorted_genome_bam.bam


bsub -R"select[mem>60000] rusage[mem=60000] span[hosts=1]" -M60000 -n 4 \
	-q normal -o velocyto10x.%J.out -e velocyto10x.%J.err \
	velocyto run --samtools-memory 1024 --samtools-threads 2 \
	--umi-extension chr -o velocyto_output \
	-b $BARCODE_FILE \
	$CELLRANGER_OUTPUT/possorted_genome_bam.bam \
	$GTF
