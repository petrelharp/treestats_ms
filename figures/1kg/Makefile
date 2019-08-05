chr20.strictMask.fasta.gz :
	# from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/accessible_genome_masks/README.20120515_phase1_pop_gen_mask
	# N - the base is an N in the reference genome GRCh37
	# L - depth of coverage is much lower than average
	# H - depth of coverage is much higher than average
	# Z - too many reads with zero mapping quality overlap this position
	# Q - the average mapping quality at the position is too low
	# P - the base passed all filters
	# 0 - an overlapping base was never observed in aligned reads
	wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/accessible_genome_masks/chr20.strictMask.fasta.gz


%.strictMask.binary.gz : %.strictMask.fasta.gz
	zcat $< | tr 'NLHZQ0P' '0000001' | tail -n +2 | gzip -c > $@
