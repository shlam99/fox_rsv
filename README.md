1)	conda install -y filtlong minimap2 samtools mafft iqtree
2)	which filtlong minimap2 samtools mafft iqtree; filtlong minimap2 samtools mafft iqtree --version
3)	cd (folder with sample fastq files transferred from dorado barcoding analysis)
4)	(put in same directory) config.sh, script_consensus_ont_rsv.sh, tree script, reference fasta
5)	(change parameters) config.sh
6)	bash script_consensus_ont_rsv.sh
7)	(either) bash script_tree_single_run_rsv.sh (for single run tree)
8)	(or) bash script_tree_all_runs_rsv.sh (for multiple runs tree)

