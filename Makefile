# =============================================================================
# Makefile — reproduction du projet "ML from Scratch in R"
# Usage : make [cible]. Adapter RSCRIPT/QUARTO au besoin.
# =============================================================================
RSCRIPT ?= Rscript
QUARTO  ?= quarto

.PHONY: all tests sims derivations rapport book check clean help pipeline docker memoire site

help:            ## Affiche les cibles disponibles
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
	awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-14s\033[0m %s\n", $$1, $$2}'

all: tests sims  ## Lance les tests puis toutes les simulations

tests:           ## Suite de tests, un processus par fichier (isolation, 622 vérifs)
	$(RSCRIPT) run_tests.R

pipeline:        ## Reproduction à dépendances suivies (targets::tar_make)
	$(RSCRIPT) -e "targets::tar_make()"

docker:          ## Construit l'image Docker reproductible
	docker build -t mlfromscratch .

# memoire.qmd est aussi un chapitre du livre : le rendre depuis la racine
# recompile tout le livre. Quarto remonte l'arborescence jusqu'au premier
# _quarto.yml, d'où le dossier jetable avec son propre _quarto.yml qui masque
# celui du projet et isole le manuscrit.
memoire:         ## Rend le mémoire SEUL en PDF (manuscrit autonome, prêt pour arXiv)
	@rm -rf .manuscrit && mkdir -p .manuscrit
	@cp rapport/memoire.qmd references.bib .manuscrit/
	@printf 'project:\n  type: default\n' > .manuscrit/_quarto.yml
	$(QUARTO) render .manuscrit/memoire.qmd --to pdf
	@mv .manuscrit/memoire.pdf rapport/memoire.pdf
	@rm -rf .manuscrit
	@echo "--> rapport/memoire.pdf"

slides:          ## Rend le diaporama de soutenance (revealjs, HTML)
	$(QUARTO) render rapport/slides.qmd

jury:            ## Rend la fiche des questions du jury (HTML)
	$(QUARTO) render rapport/questions_jury.qmd --to html

site:            ## Construit le site de documentation pkgdown
	$(RSCRIPT) -e "pkgdown::build_site()"

sims:            ## Toutes les études Monte Carlo
	$(RSCRIPT) run_all.R sims

derivations:     ## Rend les 48 dérivations Quarto (HTML)
	$(QUARTO) render derivations

rapport:         ## Rend le rapport de synthèse et l'annexe
	$(QUARTO) render rapport/rapport.qmd
	$(QUARTO) render rapport/appendix_synthese.qmd

book:            ## Compile le livre Quarto complet (dérivations + rapport)
	$(QUARTO) render

check:           ## R CMD check du package
	$(RSCRIPT) -e "if (requireNamespace('rcmdcheck', quietly=TRUE)) rcmdcheck::rcmdcheck() else message('installer rcmdcheck')"

clean:           ## Supprime les sorties générées
	rm -rf simulations/output _book *_files derivations/*.html rapport/*.html

.DEFAULT_GOAL := help
