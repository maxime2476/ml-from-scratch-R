# Point d'entrée pour R CMD check / testthat::test_dir.
# Les fonctions sont chargées depuis R/ (package non installé pendant le
# développement) via un helper dans tests/testthat/helper-source.R.
library(testthat)
test_dir("testthat")
