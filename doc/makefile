# Makefile for docs

all: notes.pdf

notes.pdf: notes.tex notes.bbl notes.aux figures/coords.tikz
	lualatex notes.tex
	lualatex notes.tex

notes.aux: notes.tex
	lualatex notes.tex

notes.bbl: notes-ref.bib notes.aux
	bibtex notes

clean:
	rm -f *.log
	rm -f *.bcf
	rm -f *.blg
	rm -f *.aux
	rm -f *.bbl
	rm -f *.out
	rm -f *.xml
