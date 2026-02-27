all: manuscript.pdf
manuscript.pdf: manuscript.md bibliography.bib
	pandoc -V geometry:margin=1.5in -H format.sty --latex-engine=xelatex --bibliography bibliography.bib --csl=nature.csl manuscript.md -o manuscript.pdf
