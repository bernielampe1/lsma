.PHONY: all
all:
	pdflatex presentation.tex
	bibtex presentation
	pdflatex presentation.tex
	pdflatex presentation.tex

.PHONY: clean
clean:
	rm -f *.bbl *.log *.blg *.aux *.toc *.snm *.out *.nav

