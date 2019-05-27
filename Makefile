.PHONY: all, clean

all: treestats_paper.pdf

# problem_statement.pdf: problem_statement.tex
# 	latexmk -pdf problem_statement

treestats_paper.pdf : references.bib macros.tex

clean: 
	rm problem_statement-1_0.pdf problem_statement-1.asy problem_statement-1.pre problem_statement-1.tex problem_statement.aux problem_statement.fdb_latexmk problem_statement.fls problem_statement.log problem_statement.pdf problem_statement.pre

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
