.PHONY: figs


figs : swept.1000.999.1e-09.diversity.pdf swept.10000.765.1e-09.diversity.pdf divergence_diagram.pdf branch_site_diagram.pdf 
	
1kgfigs : figures/tgp_geva_chr20.branch.1000000.diversity.pdf 1kg_chr20.branch.1000000.diversity.pdf
	$(MAKE) -C 1kg

### Selective sweep simulations

swept.1000.%.1e-09.diversity.pdf : swept.1000.%.trees
	python3 plot_sweep.py $< 1000 1e-9

swept.1000.%.trees :
	slim -d seed=$* -d N=1000 sweep.slim

swept.10000.%.1e-09.diversity.pdf : swept.10000.%.trees
	python3 plot_sweep.py $< 10000 1e-9

swept.10000.%.trees :
	slim -d seed=$* -d N=10000 sweep.slim

### Introgression simulations

intro/introgressed.1000.%.trees :
	slim -d seed=$* -d N=1000 introgression.slim

intro/introgressed.10000.%.trees :
	slim -d seed=$* -d N=10000 introgression.slim

### 1000 genomes whatnot

%.branch.1000000.diversity.pdf : %.trees 
	./plot_1kg_diversity.py $< 1kg/chr20.strictMask.binary.gz branch 1e6

%.branch.1000000.divergences.tsv : %.trees
	./calc_1kg_divergences.py $< 1kg/chr20.strictMask.binary.gz branch 1e6

%.branch.1000000.D.pdf %.branch.1000000.f4.pdf : %.branch.1000000.divergences.tsv 
	./plot_1kg_D_and_f4.py $*.trees 1kg/chr20.strictMask.binary.gz branch 1e6

#### rules from above

include ../Makefile


