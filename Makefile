all:
	make cleanpdf && make report.pdf && make slides.pdf

draft.pdf:
	cd draft && latexmk -halt-on-error -f -shell-escape -pdf -quiet draft.tex && rsync draft.pdf ../draft.pdf 

report.pdf:
	cd report && latexmk -halt-on-error -f -shell-escape -pdf -quiet report.tex && rsync report.pdf ../report.pdf 

slides.pdf:
	cd slides && latexmk -halt-on-error -f -shell-escape -pdf -quiet slides.tex && rsync slides.pdf ../slides.pdf 

install:
	pip install -e modules/
	cp scripts/croco-ipy-load.py $$VIRTUAL_ENV/bin/croco-ipy-load
	chmod +x $$VIRTUAL_ENV/bin/croco-ipy-load

cleanpdf:
	rm -f report.pdf slides.pdf report/report.pdf slides/slides.pdf

cleanaux:
	rm -f */*.aux */*.fdb_latexmk */*.fls */*.log */*.bak* */*.bbl */*.blg */*.out */*Notes.bib */*blx.bib */*.run.xml */*.toc */*.spl */*.nav */*.snm

clean:
	make cleanpdf && make cleanaux




