# Charge toutes les implémentations de R/ avant d'exécuter les tests, afin de
# pouvoir tester sans installer le package. Détecte la racine du projet en
# remontant depuis le répertoire courant jusqu'au dossier contenant R/.
.find_project_root <- function(start = getwd()) {
  d <- normalizePath(start, winslash = "/", mustWork = FALSE)
  for (i in 1:6) {
    if (dir.exists(file.path(d, "R")) && file.exists(file.path(d, "DESCRIPTION"))) {
      return(d)
    }
    d <- dirname(d)
  }
  stop("Racine du projet introuvable (dossier R/ + DESCRIPTION).")
}

local({
  root <- .find_project_root()
  for (f in list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE)) {
    source(f, local = FALSE)
  }
})
