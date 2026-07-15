# Point d'entrée pour R CMD check : le paquet est installé, on l'attache et on
# teste ce qui est réellement livré (idiome standard testthat).
#
# NB : pour exécuter la suite COMPLÈTE en développement, préférer `Rscript
# run_tests.R` (racine du projet), qui isole chaque fichier de test dans un
# sous-processus R frais. Chargés ensemble dans un même processus, les nombreux
# packages de référence (DoubleML, grf, iml, hdm, gmm...) se marchent dessus
# (setGeneric S4 sur coef/mean, etc.) et corrompent l'environnement — chaque
# fichier passe pourtant isolément.
library(testthat)
library(mlfromscratch)

test_check("mlfromscratch")
