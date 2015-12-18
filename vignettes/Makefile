all: pdf clean

pdf: HardinRocke.pdf GreenMartinSimulation.pdf

%.pdf:	%.tex
	rm -f $@
	pdflatex $* 
	bibtex $* 
	pdflatex $* 
	pdflatex $* 
	rm $<

clean:
	#rm -f *.aux *.bbl *.blg *.log *.toc *.out 
