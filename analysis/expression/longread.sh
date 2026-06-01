
##MINIMAP2
# paths to reference genome
REF_GENOME="/scratch/bioschool/phd/blz227562/hg38/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
THREADS=12 # Change this to match your server's available cores

# Loop through all fastq or fastq.gz files in the current directory
for fastq_file in *.fastq.gz *.fastq; do

    # Check if files exist to avoid running on empty extensions
    [ -e "$fastq_file" ] || continue

    # Extract the base name of the file (strips out .fastq or .fastq.gz)
    base_name=$(basename "$fastq_file" | sed -E 's/\.fastq(\.gz)?$//')

    echo "========================================================"
    echo "Processing: $fastq_file -> Output prefix: $base_name"
    echo "========================================================"

    # Run minimap2 with IsoQuant-required MD tags, sort, and save as BAM
    minimap2 -ax splice:hq -uf -Y -G 200000 --MD -t $THREADS "$REF_GENOME" "$fastq_file" | \
    samtools view -Su - | \
    samtools sort -@ $THREADS -o "${base_name}.sorted.bam" -

    # Index the sorted BAM file immediately
    echo "Indexing: ${base_name}.sorted.bam"
    samtools index "${base_name}.sorted.bam"

done

echo "All files processed, sorted, and indexed successfully!"

##ISOQUANT
/home/bioschool/phd/blz227562/.conda/envs/longread_env/bin/isoquant.py --bam SGNex_A549_directcDNA_replicate1_run3.sorted.bam SGNex_A549_directcDNA_replicate2_run1.sorted.bam SGNex_A549_directcDNA_replicate3_run2.sorted.bam SGNex_A549_directcDNA_replicate5_run3.sorted.bam SGNex_H9_directcDNA_replicate1_run2.sorted.bam SGNex_H9_directcDNA_replicate2_run3.sorted.bam SGNex_H9_directcDNA_replicate3_run3.sorted.bam SGNex_H9_directcDNA_replicate4_run3.sorted.bam SGNex_Hct116_directcDNA_replicate1_run4.sorted.bam SGNex_Hct116_directcDNA_replicate3_run2.sorted.bam SGNex_Hct116_directcDNA_replicate4_run1.sorted.bam SGNex_Hct116_directcDNA_replicate5_run1.sorted.bam SGNex_HepG2_directcDNA_replicate1_run1.sorted.bam SGNex_HepG2_directcDNA_replicate1_run2.sorted.bam SGNex_HepG2_directcDNA_replicate4_run1.sorted.bam SGNex_HepG2_directcDNA_replicate4_run2.sorted.bam SGNex_HepG2_directcDNA_replicate5_run3.sorted.bam SGNex_HEYA8_directcDNA_replicate1_run3.sorted.bam SGNex_HEYA8_directcDNA_replicate2_run3.sorted.bam SGNex_HEYA8_directcDNA_replicate3_run2.sorted.bam SGNex_K562_directcDNA_replicate1_run2.sorted.bam SGNex_K562_directcDNA_replicate2_run1.sorted.bam SGNex_K562_directcDNA_replicate3_run1.sorted.bam SGNex_K562_directcDNA_replicate4_run2.sorted.bam SGNex_MCF7_directcDNA_replicate1_run2.sorted.bam SGNex_MCF7_directcDNA_replicate3_run3.sorted.bam SGNex_MCF7_directcDNA_replicate4_run2.sorted.bam --output isoquant_bam --data_type nanopore --report_novel_unspliced true --threads 24 --reference /scratch/bioschool/phd/blz227562/hg38/GCF_000001405.40_GRCh38.p14_genomic.fna.gz --genedb /scratch/bioschool/phd/blz227562/hg38/GCF_000001405.40_GRCh38.p14_genomic.gtf --complete_genedb --polya_trimmed none --sqanti_output


