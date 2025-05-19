0)  Obtain fastq/bam sample files (Dorado Post Sequencing Analysis)
1)	conda install -y filtlong minimap2 samtools mafft iqtree
2)	which filtlong minimap2 samtools mafft iqtree; filtlong minimap2 samtools mafft iqtree --version
3)	cd (folder with sample fastq files transferred from dorado barcoding analysis)
4)	Put RSV_ont folder in miniconda3/bin/IRMA_RES/modules/
5)	(put in same directory) config.sh, script_consensus_ont_rsv.sh, tree script, reference fasta
6)	(change parameters) config.sh
7)	bash script_consensus_ont_rsv.sh
8)	(either) bash script_tree_single_run_rsv.sh (for single run tree)
9)  (or) bash script_tree_all_runs_rsv.sh (for multiple runs tree)
