# =============================================================================
# Dockerfile — environnement reproductible pour "ML from Scratch in R"
# -----------------------------------------------------------------------------
# Image de base rocker/verse : R épinglé + Quarto + LaTeX (pour le livre).
# Les versions EXACTES des packages sont restaurées depuis renv.lock, garantissant
# une reproduction bit-à-bit indépendante de la machine hôte.
#
#   docker build -t mlfromscratch .
#   docker run --rm mlfromscratch                 # lance la suite de tests isolée
#   docker run --rm mlfromscratch make book       # compile le livre Quarto
# =============================================================================
ARG R_VERSION=4.6.1
FROM rocker/verse:${R_VERSION}

# Bibliothèques système pour quelques dépendances compilées (GMM/GLPK, XML).
RUN apt-get update && apt-get install -y --no-install-recommends \
        libglpk-dev libxml2-dev libgsl-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /project

# Restauration des packages AUX VERSIONS de renv.lock (couche cachée si inchangé).
COPY renv.lock renv.lock
RUN R -q -e "install.packages('renv', repos='https://cloud.r-project.org'); renv::restore(prompt = FALSE)"

# Code du projet.
COPY . .

# Par défaut : suite de tests complète, un processus par fichier (isolation).
CMD ["Rscript", "run_tests.R"]
