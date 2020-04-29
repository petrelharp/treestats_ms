.PHONY: all, figures

treestats_paper.pdf : references.bib macros.tex figures review-responses.tex

treestats_paper_only.pdf : treestats_paper.pdf
	pdfjam --outfile $@ $< 1-31

cover_letter.pdf : treestats_paper.pdf
	pdfjam --outfile $@ $< 32

review_responses.pdf : treestats_paper.pdf
	pdfjam --outfile $@ $< 33-

figures :
	$(MAKE) -C figures

eps_figs : figures/swept.1000.999.1e-09.diversity.eps figures/swept.10000.765.1e-09.diversity.eps figures/divergence_diagram.eps figures/branch_site_diagram.eps figures/divergence_diagram_simple.eps figures/allele_frequency_diagram.eps figures/adding_weights.eps figures/1kg/relate_chr20.site.1000000.region.diversity.eps figures/1kg/relate_chr20.branch.1000000.region.diversity.eps figures/intro/introgressed.1000.23.10000.1e-09.f4sd.eps figures/1kg/relate_chr20_site_div_branch.1000000.diversity.eps figures/tree_0_init.eps figures/tree_0_out_0.eps tskit_stat_benchmarks/benchmarks_without_copy_longer_genome.eps figures/intro/introgressed.1000.23.1000000.1e-09.f4.eps figures/1kg/relate_chr20_GBR.branch.1000000.site.ratio.eps

%.pdf : %.tex %.bbl
	while ( pdflatex $<;  grep -q "Rerun to get" $*.log ) do true ; done

%.aux : %.tex
	-pdflatex $<

%.bbl : %.aux
	bibtex $<

%.html : %.md
	Rscript -e "templater::render_template(md.file='$<', output='$@')"

%.svg : %.pdf
	inkscape $< --export-plain-svg=$@

%.png : %.pdf
	convert -density 300 $< -flatten $@

%.pdf : %.ink.svg
	inkscape $< --export-pdf=$@

%.eps : %.pdf
	inkscape --without-gui --export-eps=$@ $<

treestats_paper-diff%.tex : treestats_paper.tex review-responses.tex
	latexdiff-git -r $* $<

treestats_paper-diff-to-submission.tex : treestats_paper-diffde01bff8f43713f6038dc45ef44ab1a7c5540eaa.tex
	mv $< $@

LATEX_MACROS = macros.tex
PANDOC_OPTS = 
PANDOC_HTML_OPTS = -c resources/pandoc.css --mathjax=https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML
PANDOC_PDF_OPTS = 
ifeq ($(wildcard $(LATEX_MACROS)),)
	# LATEX_MACROS doesn't exist
else
	PANDOC_HTML_OPTS += -H <(echo '\['; cat $(LATEX_MACROS); echo '\]')
	PANDOC_PDF_OPTS += -H $(LATEX_MACROS)
endif

%.pdf : %.md
	pandoc $(PANDOC_OPTS) $(PANDOC_PDF_OPTS) -f markdown -o $@ $<
