COUNTING BIFUNCTIONAL GENES

REQUIRES: R, PYTHON, conda environment from nm\_nr.yaml, additional python/R libraries are explicitly mentioned in the scripts

1. Count the noncoding transcripts for humans and other organisms, separates the noncoding transcripts for all genes and bifunctional genes: count_noncoding_transcripts.py
Input data was rna.fna files from the RefSeq Assembly of the respective organisms

Output files in /data/counts/

* 18122023\_T2T-CHM13v2.0\_genewise\_counts.csv
* 18122023\_T2T-CHM13v2.0\_gene\_accessions.csv
* 23122023\_GRCh38.p14\_genewise\_counts.csv
* 23122023\_GRCh38.p14\_gene\_accessions.csv
* 23112024\_GRCh38.p14\_genewise\_status.csv
* 22032025\_rhizopus\_transcripts\_accessions.csv
* 22032025\_rhizopus\_transcripts\_counts.csv
* 22032025\_rhizopus\_transcripts\_headers.csv
* 22032025\_rice\_transcripts\_accessions.csv
* 22032025\_rice\_transcripts\_counts.csv
* 22032025\_rice\_transcripts\_headers.csv
* 22032025\_spombe\_transcripts\_accessions.csv
* 22032025\_spombe\_transcripts\_counts.csv
* 22032025\_spombe\_transcripts\_headers.csv
* 22032025\_wheat\_transcripts\_accessions.csv
* 22032025\_wheat\_transcripts\_counts.csv
* 22032025\_wheat\_transcripts\_headers.csv
* 05032024\_chimpanzee\_transcripts\_accessions.csv"
* 05032024\_chimpanzee\_transcripts\_counts.csv"
* 05032024\_chimpanzee\_transcripts\_headers.csv"
* 05032024\_mouse\_transcripts\_accessions.csv"
* 05032024\_mouse\_transcripts\_counts.csv"
* 05032024\_mouse\_transcripts\_headers.csv"
* 05032024\_zebrafish\_transcripts\_accessions.csv"
* 05032024\_zebrafish\_transcripts\_counts.csv"
* 05032024\_zebrafish\_transcripts\_headers.csv"
* 18022025\_celegans\_transcripts\_accessions.csv"
* 18022025\_celegans\_transcripts\_counts.csv"
* 18022025\_celegans\_transcripts\_headers.csv"
* 18022025\_chicken\_transcripts\_accessions.csv"
* 18022025\_chicken\_transcripts\_counts.csv"
* 18022025\_chicken\_transcripts\_headers.csv"
* 18022025\_fruitfly\_transcripts\_accessions.csv"
* 18022025\_fruitfly\_transcripts\_counts.csv"
* 18022025\_fruitfly\_transcripts\_headers.csv"
* 18022025\_rat\_transcripts\_accessions.csv"
* 18022025\_rat\_transcripts\_counts.csv"
* 18022025\_rat\_transcripts\_headers.csv"
* 18022025\_xenopus\_transcripts\_accessions.csv"
* 18022025\_xenopus\_transcripts\_counts.csv"
* 18022025\_xenopus\_transcripts\_headers.csv"
* 18022025\_yeast\_transcripts\_accessions.csv"
* 18022025\_yeast\_transcripts\_counts.csv"
* 18022025\_yeast\_transcripts\_headers.csv"
* 20042025\_fruitfly\_transcripts\_accessions.csv"
* 20042025\_fruitfly\_transcripts\_counts.csv"
* 20042025\_fruitfly\_transcripts\_headers.csv"
* 22032025\_afumi\_transcripts\_accessions.csv"
* 22032025\_afumi\_transcripts\_counts.csv"
* 22032025\_afumi\_transcripts\_headers.csv"
* 22032025\_aniger\_transcripts\_accessions.csv"
* 22032025\_aniger\_transcripts\_counts.csv"
* 22032025\_aniger\_transcripts\_headers.csv"
* 22032025\_arabidopsis\_transcripts\_accessions.csv"
* 22032025\_arabidopsis\_transcripts\_counts.csv"
* 22032025\_arabidopsis\_transcripts\_headers.csv"
* 22032025\_candida\_transcripts\_accessions.csv"
* 22032025\_candida\_transcripts\_counts.csv"
* 22032025\_candida\_transcripts\_headers.csv"
2. Assign gene status to all genes in other organisms: count\_gene\_status\_other\_species.R
This file took outputs from previous script and also was used to plot Figures 2A and 2B in the manuscript.

Output files in /data/counts/other species/

3. Conservation of bifunctional genes and the status of their orthologs: conservation of species.R
This file took outputs from previous script and was also used  Figures 2C and 2D in the manuscript.
Output files in /data/counts/
* 2025-04-29\_counts\_status\_Across\_species.csv

=============================================================================================================================================
Comparison of GRCh38.p14 and T2T

1. Load both genes with annotations and compare gene biotype annotations \& plot them: FigS2.R
No output files except for plots.



=============================================================================================================================================
ENRICHMENT ANALYSIS FOR BIFUNCTIONAL, CODING \& NONCODING GENES

REQUIRES: R, clusterProfiler, ggplot2 alongwith other libraries mentioned in the script

1. All enrichment/over-representation analysis was done by the script:Fig3.R

Outputs are saved in /enrichment/
Files are also supplied as Data S5, S6 and S7 accompanying the manuscript.

* results\_GO\_BP.csv
* results\_GO\_MF.csv
* results\_GO\_CC.csv

=============================================================================================================================================
EXPRESSION ANALYSIS FROM GTEX: NON-CODING JUNCTIONS OF BIFUNCTIONAL GENES

REQUIRES: R, dplyr, ggplot2, additional tools and software mentioned in exp\_env.yaml.
Please create and activate conda environment before trying to replicate the analysis.
Outputs are saved by tissue types.

1. Annotate the junctions as coding, non-coding or shared: get\_unique\_junctions.R

Outputs in /junctions\_nnotated/

2. Map data from recount3 to annotated junctions: recount\_junction\_data\_retrieval\_mapping.R

Outputs in /data/GTEx exon junctions/

3. Calculate coding and noncoding read fractions for bifunctional genes: exon\_junction\_analysis.R

Outputs in /data/GTEx exon junctions/tissue\_medians \& Outputs in /data/GTEx exon junctions/stats/

4. Plot the outputs, calculate Tau specificity score and global non-coding fraction medians: Fig4.R

Oututs in /data/GTEx exon junctions/

=============================================================================================================================================


EXPRESSION OF SELECT NONCODING TRANSCRIPT VARIANTS THROUGH REAL-TIME QPCR

Requires: R, ggplot2

Input files are in /data/expression

1. Plot expression with error bars for n=3: Fig5.R

=============================================================================================================================================

ANALYSIS \& PREDICTION OF ORFs

REQUIRES: R, PYTHON, blast+, conda environment from nm\_nr.yaml, additional python/R libraries are explicitly mentioned in the scripts

1. ORFs were predicted with orfipy from noncoding transcript fasta files: https://github.com/urmi-21/orfipy

Output files in /data/ORFs/predicted/

* 03052024\_pep\_ATG\_bifunc\_genes
* 03052024\_pep\_CTG\_bifunc\_genes
* 03052024\_pep\_GTG\_bifunc\_genes
* 03052024\_pep\_TTG\_bifunc\_genes
* 03052024\_rna\_ATG\_bifunc\_genes.fa
* 03052024\_rna\_CTG\_bifunc\_genes.fa
* 03052024\_rna\_GTG\_bifunc\_genes.fa
* 03052024\_rna\_TTG\_bifunc\_genes.fa



2. Formatted ORFs for downstream analysis: format\_ORFs.py

Outputs in /data/ORFs/predicted/

* 03052024\_formatted\_pep\_ATG\_bifunc\_genes.csv
* 03052024\_formatted\_pep\_CTG\_bifunc\_genes.csv
* 03052024\_formatted\_pep\_GTG\_bifunc\_genes.csv
* 03052024\_formatted\_pep\_TTG\_bifunc\_genes.csv
3. Get unique ORFs from the predicted peptide sequence: get\_unique\_ORFs.R

Outputs in /data/ORFs/unqiue\_ORFs/

* 2024-05-05\_fasta\_uniqueORFs\_bifunc\_ATG.txt
* 2024-05-05\_fasta\_uniqueORFs\_bifunc\_CTG.txt
* 2024-05-05\_fasta\_uniqueORFs\_bifunc\_GTG.txt
* 2024-05-05\_fasta\_uniqueORFs\_bifunc\_TTG.txt
4. blastp of predicted peptides: unique\_ORFs\_blast.sh
* blastdb was created using the peptides/proteins from the human RefSeq assembly protein file

Outputs in /data/ORFs/blastp/05052024\_csv/

* 05052024\_blastp-fast\_hybridGenes\_uniqueORFs\_GTG.csv
* 05052024\_blastp-fast\_hybridGenes\_uniqueORFs\_TTG.csv
* 05052024\_blastp-fast\_hybridGenes\_uniqueORFs\_ATG.csv
* 05052024\_blastp-fast\_hybridGenes\_uniqueORFs\_CTG.csv
5. Format blastp hits and match NMD status as in NCBI to predicted unique ORFs: match\_blastp\_hits\_downstream\_analysis

Outputs in data/ORFs/unique\_ORFs/

* 2024-05-21\_blastpHybridGenes\_ATG\_uniqueORFsMatched.tsv
* 2024-05-21\_blastpHybridGenes\_CTG\_uniqueORFsMatched.tsv
* 2024-05-21\_blastpHybridGenes\_GTG\_uniqueORFsMatched.tsv
* 2024-05-21\_blastpHybridGenes\_TTG\_uniqueORFsMatched.tsv

Outputs in data/ORFs/uniqueORFs\_nmdMatch/

* 2024-05-15\_bifunc\_ATG\_uniqueORFs\_nmd.csv
* 2024-05-16\_bifunc\_CTG\_uniqueORFs\_nmd.csv
* 2024-05-16\_bifunc\_GTG\_uniqueORFs\_nmd.csv
* 2024-05-16\_bifunc\_TTG\_uniqueORFs\_nmd.csv
* 2024-05-21\_unique0RFs\_bifunc\_ATG\_nmd\_blast.csv
* 2024-05-21\_unique0RFs\_bifunc\_CTG\_nmd\_blast.csv
* 2024-05-21\_unique0RFs\_bifunc\_GTG\_nmd\_blast.csv
* 2024-05-21\_unique0RFs\_bifunc\_TTG\_nmd\_blast.csv
6. Calculate Kozak Score from noncoding transcript sequences: implemented from https://github.com/Agleason1/TIS-Predictor

Outputs in /kss/

* 14062024\_kss-len\_formatted\_pep\_ATG\_bifunc\_genes.csv
* 14062024\_kss-len\_formatted\_pep\_CTG\_bifunc\_genes.csv
* 14062024\_kss-len\_formatted\_pep\_GTG\_bifunc\_genes.csv
* 14062024\_kss-len\_formatted\_pep\_TTG\_bifunc\_genes.csv
7. Combine and predict NMD, blast and KSS outputs: mRNA\_KSS\_match.R
This script also takes predicted and formatted outputs from all coding transcripts of bifunctional genes to put a flag if a noncoding ORF is the same as the coding transcript ORF. This has been done using code similar to what is described above for noncoding transcripts.

Outputs in data/ORFs/uniqueORFs\_nmdMatch/

* 2024-07-17\_unique0RFs\_bifunc\_CTG\_nmd\_blast\_mRNAs\_KSS.csv
* 2024-07-17\_unique0RFs\_bifunc\_GTG\_nmd\_blast\_mRNAs\_KSS.csv
* 2024-07-17\_unique0RFs\_bifunc\_TTG\_nmd\_blast\_mRNAs\_KSS.csv
* 2024-07-17\_unique0RFs\_bifunc\_ATG\_nmd\_blast\_mRNAs\_KSS.csv
8. These files were used to plot: Fig7.R

=============================================================================================================================================


EXPRESSION ANALYSIS FOR LONGREAD SEQUENCING DATA FROM SGNEX

REQUIRES: AWS s3 CLI, minimap2, samtools, isoquant and tools mentioned in longread\_env.yaml

1. Downloaded fastq files for direct-cDNA sequencies of all cell lines and replicates mentioned in SGNEx at https://registry.opendata.aws/sgnex/
Sample details are at https://github.com/GoekeLab/sg-nex-data/blob/master/docs/samples.tsv
2. Align the reads to GRCh38.p14 using minimap2 and quantify transcripts using the isoQuant tool: longread.sh
Aligned and sorted bam files have not been uploaded. IsoQuant outputs are in /data/longread/



=============================================================================================================================================

Matching microRNA targets from microCLIP and TargetScan 

1. Download Target Scan default predictions and microCLIP predictions and lift coordinates to hg38. Downloaded files from the microCLIP paper (https://www.nature.com/articles/s41467-018-06046-y) and TargetScan website (Downloaded files: Predicted\_Targets\_Context\_Scores.default\_predictions.txt, Predicted\_Target\_Locations.default\_predictions.hg19.bed, miR\_Family\_Info.txt)

Output files in /microRNA/hg38\_liftOver/

2\. Match predictions to transcripts from NCBI genome gtf:microCLIP\_matching\_all\_optimized.R

Output files in /microRNA/

3\. Get all possible microRNA-gene-transcript combinations and count sites common to both predictions:microRNA\_common\_targets.R

Output files in /microRNA/

=============================================================================================================================================



MISCELLANEOUS ANALYSIS

1. Transcript expression from CHESS3 was downloaded from: https://ccb.jhu.edu/chess/ (file was chess3.1.3.GRCh38.gtf.gz). The expression for noncoding transcripts from RefSeq was filtered and plotted.
Script: CHESS3\_analysis.R

2\. Circular RNA matching to bifunctional genes: circular RNA sequences were retrieved from circAtlas 3.0: circRNA\_matching.R



=============================================================================================================================================



OTHER NOTES

* If there's a file reading error, please manually look for the files in the data folder as there might be some conflict in moving data from Lab PC to HPC and back!
* For issues, please email or contact Janki at blz227562.iitd.ac.in or jankinsan@gmail.com

