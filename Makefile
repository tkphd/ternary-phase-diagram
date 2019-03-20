# GNU Makefile for ternary-phase-diagram

all: ternary-diagram.png
.PHONY: all

ternary-diagram.png: phase-diagram.py paraboloids.c paraboloids.so
	python $<

paraboloids.c: free-energies.py
	python $<

paraboloids.so: paraboloids.c
	gcc -Wall -fPIC -shared $< -o $@

.PHONY: docs
docs: ternary-diagram.tex
	pdflatex -interaction=nonstopmode $<

.PHONY: clean
clean:
	rm -f paraboloids.h paraboloids.c paraboloids.so ternary-diagram.pdf
