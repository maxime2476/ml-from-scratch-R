# Attachement des packages qui exigent d'être « attachés » (et non seulement
# appelés via `::`) pour que `testthat::test_dir` fonctionne en une passe.
# mclust : son dispatch interne (mclustBIC) est introuvable sans attachement.
if (requireNamespace("mclust", quietly = TRUE)) suppressMessages(library(mclust))
