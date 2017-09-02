# Makefile for docs

all:
	lualatex notes.tex
	biber notes
	lualatex notes.tex
	lualatex notes.tex

clean:
	rm -f *.log
	rm -f *.bcf
	rm -f *.blg
	rm -f *.aux
	rm -f *.bbl
	rm -f *.out
	rm -f *.xml
