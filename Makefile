all:
	make cleanpdf && make report.pdf && make draft.pdf

draft.pdf:
	cd draft && latexmk -halt-on-error -f -shell-escape -pdf -quiet draft.tex && rsync draft.pdf ../draft.pdf 

report.pdf:
	cd report && latexmk -halt-on-error -f -shell-escape -pdf -quiet report.tex && rsync report.pdf ../report.pdf 

install:
	pip install -e modules/
	@if [ -n "$$VIRTUAL_ENV" ]; then \
		cp modules/croco-ipy-load.py $$VIRTUAL_ENV/bin/croco-ipy-load; \
		chmod +x $$VIRTUAL_ENV/bin/croco-ipy-load; \
	elif [ -n "$$CONDA_PREFIX" ]; then \
		cp modules/croco-ipy-load.py $$CONDA_PREFIX/bin/croco-ipy-load; \
		chmod +x $$CONDA_PREFIX/bin/croco-ipy-load; \
	else \
		echo "No virtual or Conda environment found."; exit 1; \
	fi

cleanpdf:
	rm -f report.pdf slides.pdf draft.pdf report/report.pdf slides/slides.pdf draft/draft.pdf

cleanaux:
	rm -f */*.aux */*.fdb_latexmk */*.fls */*.log */*.bak* */*.bbl */*.blg */*.out */*Notes.bib */*blx.bib */*.run.xml */*.toc */*.spl */*.nav */*.snm */*.lof */*.lot

clean:
	make cleanpdf && make cleanaux




