.SUFFIXES: .gp .tex .eps .pdf .eps.gz

.eps.pdf:
		epstopdf --outfile $@ $<

.eps.gz.pdf:
		gzip -dc $< | epstopdf --filter > $@

.pdf.eps:
		pdftops -eps $< $@

all:longdust.pdf

longdust.pdf:longdust.tex longdust.bib len.pdf
		pdflatex longdust; bibtex longdust; pdflatex longdust; pdflatex longdust;

len.eps:plot-tr.gp eval.txt
		gnuplot plot-tr.gp

clean:
		rm -fr *.toc *.aux *.bbl *.blg *.idx *.log *.out *~ longdust.pdf
