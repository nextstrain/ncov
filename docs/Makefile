# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = src
BUILDDIR      = build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help help-docker Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.ONESHELL:
docker-html:
	set -euox
	docker build -t nextstrain-docs-builder --network=host .
	docker run -it --rm \
	--name=nextstrain-docs-builder-$(shell date +%s) \
	--init \
	--user=$(shell id -u):$(shell id -g) \
	--volume=$(shell pwd):/home/user/src \
	--workdir=/home/user/src \
	--env 'TERM=xterm-256colors' \
	nextstrain-docs-builder
