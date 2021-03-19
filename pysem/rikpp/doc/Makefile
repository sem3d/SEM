DOC = manual

.SUFFIXES:

all: $(DOC).pdf

%.pdf: %.ps
	ps2pdf $<

%.ps: %.dvi
	dvips $* -o 

$(DOC).dvi: *.tex *.bib
	latex $(DOC)

# selfdoc.tex: ../SRC/*.f90 
# 	cd .. ; csh ./selfdoc.csh ; cd MANUAL

*.bib:
	bibtex $(DOC) ; latex $(DOC)


clean:
	rm *.pdf *.aux *.log *.dvi \
	*.bbl *.blg *.nlo *.toc

biblio: 
	bibtex $(DOC)
recom:
	pdflatex $(DOC)
