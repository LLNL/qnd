# Minimal Makefile for Sphinx documentation
# (requires gmake)

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = qnd
SOURCEDIR     = .
BUILDDIR      = _build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# rst2pdf may require pip install rst2pdf (not Anaconda default package)
# rst2html may be the module name rst2html.py
RST2PDF=rst2pdf
RST2HTML=rst2html.py
.SUFFIXES: .rst

%.pdf: %.rst
	$(RST2PDF) -s pdfstyle.json --use-floating-images $< -o $@

%.html: %.rst
	$(RST2HTML) $< $@

.PHONY: help Makefile clean

clean::
	rm -rf $(BUILDDIR)/*

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
