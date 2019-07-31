.PHONY: figs


figs : swept.1000.999.1e-09.diversity.pdf swept.10000.765.1e-09.diversity.pdf divergence_diagram.pdf

swept.1000.%.1e-09.diversity.pdf : swept.1000.%.trees
	python3 plot_sweep.py $< 1000 1e-9

swept.1000.%.trees :
	slim -d seed=$* -d N=1000 sweep.slim

swept.10000.%.1e-09.diversity.pdf : swept.10000.%.trees
	python3 plot_sweep.py $< 10000 1e-9

swept.10000.%.trees :
	slim -d seed=$* -d N=10000 sweep.slim

include ../Makefile
