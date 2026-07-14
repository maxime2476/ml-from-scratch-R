# Point d'entrée pour R CMD check / testthat::test_dir.
# Les fonctions sont chargées depuis R/ (package non installé pendant le
# développement) via un helper dans tests/testthat/helper-source.R.
#
# NB : pour exécuter la suite COMPLÈTE de façon fiable, préférer `Rscript
# run_tests.R` (racine du projet), qui isole chaque fichier de test dans un
# sous-processus R frais. Chargés ensemble dans un même processus, les nombreux
# packages de référence (DoubleML, grf, iml, hdm, gmm...) se marchent dessus
# (setGeneric S4 sur coef/mean, etc.) et corrompent l'environnement — chaque
# fichier passe pourtant isolément (337 tests, 0 échec).
library(testthat)
test_dir("testthat")
