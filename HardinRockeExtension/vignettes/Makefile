all: pdf clean

pdf: HardinRocke.pdf

%.pdf:	%.tex
	pdflatex $* 
	bibtex $* 
	pdflatex $* 
	pdflatex $* 
	rm $<

clean:
	rm -f *.aux *.bbl *.blg *.log *.toc *.out 
