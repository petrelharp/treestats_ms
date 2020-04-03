.PHONY: figs


figs : swept.1000.999.1e-09.diversity.pdf swept.10000.765.1e-09.diversity.pdf divergence_diagram.pdf branch_site_diagram.pdf divergence_diagram_simple.pdf allele_frequency_diagram.pdf adding_weights.pdf

1kg_figs : 1kg/relate_chr20.site.1000000.region.diversity.pdf 1kg/relate_chr20.branch.1000000.region.diversity.pdf

1kg/relate_chr20_site_div_branch.1000000.diversity.pdf : 1kg/relate_chr20.branch.1000000.diversity.tsv 1kg/relate_chr20.site.1000000.diversity.tsv
	python3 compare_relate_diversity.py

%.site.1000000.region.diversity.pdf :
	./plot_1kg_diversity.py $*.trees chr20.strictMask.binary.gz site 1000000

%.branch.1000000.region.diversity.pdf :
	./plot_1kg_diversity.py $*.trees chr20.strictMask.binary.gz branch 1000000

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


