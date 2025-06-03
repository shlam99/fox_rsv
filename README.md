0)  obtain fastq/bam files of samples  (Dorado Post Sequencing Analysis - Barcoding)
1)	conda install -y samtools filtlong irma mafft iqtree
2)	which samtools filtlong irma mafft iqtree; samtools filtlong irma mafft iqtree --version
3)	cd to folder with sample files
4)	put RSV_ont folder in miniconda3/bin/IRMA_RES/modules/
5)	put in same directory config_empty.sh, rsv_irma_labelled_xargs.sh
6)	change details & rename config_empty.shto config.sh
7)	bash rsv_irma_labelled_xargs.sh
8)  bash rsv_mafft_iqtree.sh
