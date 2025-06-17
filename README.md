This is a data analysis pipeline developed by Pr Vijaykrishna Dhanasekaran Lab to process sequenced NBD114.24 barcoded RSV samples to consensus sequences & phylogenetic trees through quality filtering, subtyping/alignment/consensus generation using custom IRMA module, and optional mafft alignment/iqtree tree construction.
Scripts in use: config_empty.sh, rsv_irma_labelled_xargs.sh, rsv_mafft_iqtree.sh

0)  obtain fastq/bam files of samples  (Dorado Post Sequencing Analysis - Barcoding)
1)	conda install -y samtools filtlong irma mafft iqtree
2)	which samtools filtlong irma mafft iqtree; samtools filtlong irma mafft iqtree --version
3)	rename reference fasta to consensus.fasta & replace the one in /RSV_ont/reference/
4)	put /RSV_ont/ folder in /miniconda3/bin/IRMA_RES/modules/
5)	put in same directory config_empty.sh, rsv_irma_labelled_xargs.sh
6)	change details & rename config_empty.shto config.sh
7)	use rsv_commands.txt for commands
8)	bash rsv_irma_labelled_xargs.sh
9)  bash rsv_mafft_iqtree.sh


------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------
# Activate miniconda3 (for some terminal)
source ~/miniconda3/bin/activate  

# Add timestamp to WSL
PS1="[\u@\h \D{%Y%m%d-%H:%M:%S}]\$"

# Use all cores (inconsistent!)
THREADS=$(nproc --all)

cd /mnt/e/
cd /mnt/c/Users/volca/OneDrive/Academia/HKU_Influenza/RSV_results/
bash /mnt/c/Users/volca/OneDrive/Academia/HKU_Influenza/RSV_results/rsv_irma_labelled_xargs.sh
bash /mnt/c/Users/volca/OneDrive/Academia/HKU_Influenza/RSV_results/rsv_mafft_iqtree.sh


cd /mnt/e/20250322_2D2H_1009/20250322_2H_bam
bash /mnt/c/Users/volca/OneDrive/Academia/HKU_Influenza/RSV_results/threefour.sh
cd /mnt/e/20250323_1H_10/20250323_1H_bam
bash /mnt/c/Users/volca/OneDrive/Academia/HKU_Influenza/RSV_results/threefour.sh
cd /mnt/e/20250328_3G_18/20250328_3G_bam
bash /mnt/c/Users/volca/OneDrive/Academia/HKU_Influenza/RSV_results/threefour.sh
