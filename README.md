0)  obtain fastq or bam files of samples  (Dorado Post Sequencing Analysis - Barcoding)
1)	conda install -y samtools filtlong irma mafft iqtree
2)	which samtools filtlong irma mafft iqtree; samtools filtlong irma mafft iqtree --version
3)	cd (folder with sample files)
4)	Put RSV_ont folder in miniconda3/bin/IRMA_RES/modules/
5)	(put in same directory) config.sh, rsv_irma_labelled.sh
6)	(change details) config.sh
7)	bash rsv_irma_labelled.sh
8)	(either) bash script_tree_single_run_rsv.sh (for single run tree)
9)  (or) bash script_tree_all_runs_rsv.sh (for multiple runs tree)
