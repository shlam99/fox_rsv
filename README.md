This is a data analysis pipeline developed by Pr Vijaykrishna Dhanasekaran Lab to process sequenced NBD114.24 barcoded RSV samples to consensus sequences & phylogenetic trees through quality filtering, subtyping/alignment/consensus generation using custom IRMA module, and optional mafft alignment/iqtree tree construction.
Scripts in use: rsv_commands.txt, config_empty.sh, rsv_irma_labelled_xargs.sh, rsv_mafft_iqtree.sh

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
