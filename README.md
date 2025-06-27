This is a data analysis pipeline developed by Pr Vijaykrishna Dhanasekaran Lab to process sequenced NBD114.24 barcoded RSV samples to consensus sequences & phylogenetic trees through quality filtering, subtyping/alignment/consensus generation using custom IRMA module, and optional mafft alignment/iqtree tree construction.
Scripts in use: config_empty.sh, rsv_irma_labelled_xargs.sh, rsv_mafft_iqtree.sh

0) 	obtain fastq/bam files of samples  (Dorado Post Sequencing Analysis - Barcoding)
1)	conda install -y samtools filtlong irma mafft iqtree
2)	which samtools filtlong irma mafft iqtree; samtools filtlong irma mafft iqtree --version
4)	put RSV_ont folder in miniconda3/bin/IRMA_RES/modules/
5)	put config_empty.sh in same directory; change details & rename config_empty.shto config.sh
7)	download rsv_irma_labelled_qc.sh; create script.sh for cd & bash commands
8)	cd [sample directory]; bash [script directory]/rsv_irma_labelled_qc.sh
9) 	copy /irma_consensus/s into one folder; bash rsv_mafft_iqtree.sh

After finishing, dorado barcoding_bam is in data acquisition unit, delete the bams in transfer hard drive & put the processed unfiltered fastq.gz files into new folder /hide/ (so when re-run script, script runs IRMA on /qc_reads/ data instead of filtlong on these data)


------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------
# Activate miniconda3 (for some terminal)
```bash
source ~/miniconda3/bin/activate  
```
# Add timestamp to WSL
```bash
PS1="[\u@\h \D{%Y%m%d-%H:%M:%S}]\$"
```
# Use all cores (inconsistent!)
```bash
THREADS=$(nproc --all)
```
# Remove irma_consensus & irma_results directories
```bash
find . -type d \( -name "irma_consensus" -o -name "irma_results" \) -exec rm -rf {} \;
```
# Extract fasta from irma_consensus
```bash
find . -type d -name "irma_consensus" -exec cp -n -v {}/* . \;
```
# Combine fasta files
```bash
find . -type f \( -name "*.fasta" -o -name "*.fa" \) -exec cat {} \; > pooled.fasta
echo "" >> pooled.fasta
```
# Frequent Commands
```bash
cd /mnt/e/
cd /Volumes/Extreme\ SSD/
cd /mnt/c/Users/volca/OneDrive/Academia/HKU_Influenza/RSV_results/
bash /mnt/c/Users/volca/OneDrive/Academia/HKU_Influenza/RSV_results/rsv_irma_labelled_xargs.sh
bash /mnt/c/Users/volca/OneDrive/Academia/HKU_Influenza/RSV_results/rsv_mafft_iqtree.sh

: <<'DISABLE'
DISABLE
```


