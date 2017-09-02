# Makefile for docs

all:
	lualatex Notes-FE.tex
	biber Notes-FE
	lualatex Notes-FE.tex
	lualatex Notes-FE.tex

clean:
	rm -f *.log
	rm -f *.bcf
	rm -f *.blg
	rm -f *.aux
	rm -f *.bbl
	rm -f *.out
	rm -f *.xml
