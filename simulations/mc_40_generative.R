# =============================================================================
# Monte Carlo — Module 40 : chaque classifieur generatif a son regime optimal.
# (1) Covariances homogenes -> LDA (moins de variance) ; heterogenes -> QDA.
# (2) Features correlees -> Naive Bayes (hypothese d'independance) se degrade.
# =============================================================================

for (f in c("40_generative", "mc_tools")) source(file.path("R", paste0(f, ".R")))
suppressMessages(library(ggplot2))
set.seed(2026)
out_dir <- "simulations/output"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

acc <- function(fit, pred, tr, te) mean(pred(fit(tr$X, tr$y), te$X) == te$y)

## (1) LDA vs QDA : covariances homogenes vs heterogenes ----------------------
R <- 100; res <- matrix(0, 2, 2, dimnames = list(c("homogene","heterogene"), c("LDA","QDA")))
for (r in seq_len(R)) {
  n <- 60
  # homogene : meme covariance
  A <- matrix(c(1, 0.5, 0.5, 1), 2); mk <- function(m) matrix(rnorm(n*2), n, 2) %*% chol(A) + rep(m, each = n)
  Xh <- rbind(mk(c(0,0)), mk(c(1.5,1.5))); yh <- factor(rep(1:2, each = n))
  # heterogene : covariances opposees
  Xe <- rbind(matrix(rnorm(n*2),n,2)%*%diag(c(0.4,2.5)), matrix(rnorm(n*2),n,2)%*%diag(c(2.5,0.4))+1)
  ye <- factor(rep(1:2, each = n))
  teh <- list(X = rbind(mk(c(0,0)), mk(c(1.5,1.5))), y = yh)
  tee <- list(X = rbind(matrix(rnorm(n*2),n,2)%*%diag(c(0.4,2.5)), matrix(rnorm(n*2),n,2)%*%diag(c(2.5,0.4))+1), y = ye)
  res[1,1] <- res[1,1] + acc(lda_fit, lda_predict, list(X=Xh,y=yh), teh)
  res[1,2] <- res[1,2] + acc(qda_fit, qda_predict, list(X=Xh,y=yh), teh)
  res[2,1] <- res[2,1] + acc(lda_fit, lda_predict, list(X=Xe,y=ye), tee)
  res[2,2] <- res[2,2] + acc(qda_fit, qda_predict, list(X=Xe,y=ye), tee)
}
res <- res / R
cat("=== (1) Precision selon le regime de covariance ===\n\n")
cat(sprintf("%-12s %8s %8s\n", "regime", "LDA", "QDA"))
cat(sprintf("%-12s %8.3f %8.3f  <- LDA prefere (moins de variance)\n", "homogene", res[1,1], res[1,2]))
cat(sprintf("%-12s %8.3f %8.3f  <- QDA prefere (covariances propres)\n\n", "heterogene", res[2,1], res[2,2]))

## (2) Naive Bayes vs LDA sous correlation des features -----------------------
R <- 200; nb_indep <- nb_corr <- lda_corr <- 0
for (r in seq_len(R)) {
  n <- 100
  mkc <- function(rho, m) { A <- matrix(c(1,rho,rho,1),2); matrix(rnorm(n*2),n,2)%*%chol(A) + rep(m,each=n) }
  # correlees
  Xc <- rbind(mkc(0.9, c(0,0)), mkc(0.9, c(1.2,-1.2))); yc <- factor(rep(1:2, each=n))
  tec <- list(X = rbind(mkc(0.9,c(0,0)), mkc(0.9,c(1.2,-1.2))), y = yc)
  nb_corr <- nb_corr + acc(naive_bayes_fit, naive_bayes_predict, list(X=Xc,y=yc), tec)
  lda_corr <- lda_corr + acc(lda_fit, lda_predict, list(X=Xc,y=yc), tec)
}
cat("=== (2) Features CORRELEES (rho=0.9) : Naive Bayes vs LDA ===\n")
cat(sprintf("  Naive Bayes (suppose l'independance) : %.3f\n", nb_corr / R))
cat(sprintf("  LDA (modelise la covariance)         : %.3f\n", lda_corr / R))
cat("\n=> Chaque modele a son regime : LDA si covariances egales, QDA si elles\n")
cat("   different, et Naive Bayes paie l'hypothese d'independance quand elle est fausse.\n")

df <- data.frame(regime = rep(c("covariances homogenes","covariances heterogenes"), each=2),
                 modele = rep(c("LDA","QDA"), 2), acc = c(res[1,],res[2,]))
gg <- ggplot(df, aes(modele, acc, fill = modele)) + geom_col(show.legend=FALSE) +
  facet_wrap(~regime) + coord_cartesian(ylim=c(0.5,1)) +
  labs(title="LDA vs QDA : chacun son regime de covariance",
       x=NULL, y="precision test") + theme_minimal(base_size=12)
ggsave(file.path(out_dir, "mc_40_generative.png"), gg, width=8, height=5, dpi=120)
cat("\nGraphique -> mc_40_generative.png\n")
