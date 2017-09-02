# Makefile for docs

all:
	lualatex notes.tex
	biber notes
	lualatex notes.tex
	lualatex notes.tex

clean:
	rm -f *.log
