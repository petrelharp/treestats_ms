.PHONY: figs

figs : swept.1000.999.1e-09.diversity.pdf 

# ones that take a long time
bigfigs : swept.10000.765.1e-09.diversity.pdf 

%.1e-09.diversity.pdf : %.trees
	python3 plot_sweep.py $< 1000

swept.1000.%.trees :
	slim -d seed=$* -d N=1000 sweep.slim
