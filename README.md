This is a data analysis pipeline developed by Pr Vijaykrishna Dhanasekaran Lab to process sequenced NBD114.24 barcoded RSV samples to consensus sequences & phylogenetic trees through quality filtering, subtyping/alignment/consensus generation using custom IRMA module, and optional mafft alignment/iqtree tree construction.
Scripts in use: config_empty.sh, rsv_irma_labelled_xargs.sh, rsv_mafft_iqtree.sh

0)  obtain fastq/bam files of samples  (Dorado Post Sequencing Analysis - Barcoding)
1)	conda install -y samtools filtlong irma mafft iqtree
2)	which samtools filtlong irma mafft iqtree; samtools filtlong irma mafft iqtree --version
3)	cd to folder with sample files
4)	put RSV_ont folder in miniconda3/bin/IRMA_RES/modules/
5)	put in same directory config_empty.sh, rsv_irma_labelled_xargs.sh
6)	change details & rename config_empty.shto config.sh
7)	use rsv_commands.txt for commands
8)	bash rsv_irma_labelled_xargs.sh
9)  bash rsv_mafft_iqtree.sh


------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------

PS1="[\u@\h \D{%Y%m%d-%H:%M:%S}]\$ "

cd /mnt/e/
cd /mnt/c/Users/volca/OneDrive/Academia/HKU_Influenza/RSV_results/

bash /mnt/c/Users/volca/OneDrive/Academia/HKU_Influenza/RSV_results/rsv_irma_labelled_xargs.sh
bash /mnt/c/Users/volca/OneDrive/Academia/HKU_Influenza/RSV_results/threefour.sh

cd /mnt/c/Users/volca/OneDrive/Academia/HKU_Influenza/RSV_results/20250425_1B1E_2222/20250425_1B
bash /mnt/c/Users/volca/OneDrive/Academia/HKU_Influenza/RSV_results/threefour.sh
cd /mnt/c/Users/volca/OneDrive/Academia/HKU_Influenza/RSV_results/20250425_1B1E_2222/20250425_1E
bash /mnt/c/Users/volca/OneDrive/Academia/HKU_Influenza/RSV_results/threefour.sh
cd /mnt/c/Users/volca/OneDrive/Academia/HKU_Influenza/RSV_results/20250515_2A2C_2422/20250515_2A
bash /mnt/c/Users/volca/OneDrive/Academia/HKU_Influenza/RSV_results/threefour.sh
cd /mnt/c/Users/volca/OneDrive/Academia/HKU_Influenza/RSV_results/20250515_2A2C_2422/20250515_2C
bash /mnt/c/Users/volca/OneDrive/Academia/HKU_Influenza/RSV_results/threefour.sh
cd /mnt/c/Users/volca/OneDrive/Academia/HKU_Influenza/RSV_results/20250526_1B1C_2414/20250526_1B
bash /mnt/c/Users/volca/OneDrive/Academia/HKU_Influenza/RSV_results/threefour.sh
cd /mnt/c/Users/volca/OneDrive/Academia/HKU_Influenza/RSV_results/20250526_1B1C_2414/20250526_1C
bash /mnt/c/Users/volca/OneDrive/Academia/HKU_Influenza/RSV_results/threefour.sh


cd /mnt/c/Users/volca/OneDrive/Academia/HKU_Influenza/RSV_results/test_irma
bash /mnt/c/Users/volca/OneDrive/Academia/HKU_Influenza/RSV_results/rsv_mafft_iqtree.sh

cd /mnt/e/20250322_2D2H_1009/20250322_2H_bam
bash /mnt/c/Users/volca/OneDrive/Academia/HKU_Influenza/RSV_results/zeroone.sh
cd /mnt/e/20250323_1H_10/20250323_1H_bam
bash /mnt/c/Users/volca/OneDrive/Academia/HKU_Influenza/RSV_results/zeroone.sh
cd /mnt/e/20250328_3G_18/20250328_3G_bam
bash /mnt/c/Users/volca/OneDrive/Academia/HKU_Influenza/RSV_results/zeroone.sh
cd /mnt/e/20250322_2D2H_1009/20250322_2H_bam
bash /mnt/c/Users/volca/OneDrive/Academia/HKU_Influenza/RSV_results/two.sh
cd /mnt/e/20250323_1H_10/20250323_1H_bam
bash /mnt/c/Users/volca/OneDrive/Academia/HKU_Influenza/RSV_results/two.sh
cd /mnt/e/20250328_3G_18/20250328_3G_bam
bash /mnt/c/Users/volca/OneDrive/Academia/HKU_Influenza/RSV_results/two.sh
cd /mnt/e/20250322_2D2H_1009/20250322_2H_bam
bash /mnt/c/Users/volca/OneDrive/Academia/HKU_Influenza/RSV_results/threefour.sh
cd /mnt/e/20250323_1H_10/20250323_1H_bam
bash /mnt/c/Users/volca/OneDrive/Academia/HKU_Influenza/RSV_results/threefour.sh
cd /mnt/e/20250328_3G_18/20250328_3G_bam
bash /mnt/c/Users/volca/OneDrive/Academia/HKU_Influenza/RSV_results/threefour.sh
