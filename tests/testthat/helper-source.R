# En développement, le paquet n'est pas installé : on charge les implémentations
# depuis R/ pour pouvoir tester sans installer.
#
# Sous R CMD check (et via run_tests.R, qui utilise pkgload::load_all), le paquet
# est déjà chargé : il ne faut alors RIEN sourcer. Le dossier R/ d'un paquet
# *installé* ne contient que la base lazy-load (.rdb/.rdx) et aucune source : la
# détection de racine s'y accrochait, ne chargeait donc aucune fonction, et les
# tests échouaient en CI sur « could not find function » alors qu'ils passaient
# en local. D'où la garde ci-dessous et l'exigence de vrais fichiers .R.
if (!"mlfromscratch" %in% loadedNamespaces()) {

  .find_project_root <- function(start = getwd()) {
    d <- normalizePath(start, winslash = "/", mustWork = FALSE)
    for (i in 1:6) {
      if (file.exists(file.path(d, "DESCRIPTION")) &&
          length(list.files(file.path(d, "R"), pattern = "\\.R$")) > 0) {
        return(d)
      }
      d <- dirname(d)
    }
    stop("Racine du projet introuvable (dossier R/ avec sources + DESCRIPTION).")
  }

  local({
    root <- .find_project_root()
    for (f in list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE)) {
      source(f, local = FALSE)
    }
  })
}
