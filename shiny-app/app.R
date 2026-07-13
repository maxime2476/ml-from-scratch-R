# =============================================================================
# Application Shiny — « ML from Scratch in R » : explorer les idées clés
# -----------------------------------------------------------------------------
# Trois onglets interactifs, tous branchés sur les implémentations from-scratch :
#   1. Biais-variance du ridge   (Module 4 : ridge_bias_var)
#   2. Chemins de régularisation (Module 4 : lasso_fit / ridge_fit)
#   3. Orthogonalisation (DML)   (Modules 1, 4 : partialling-out de Neyman)
#
# Lancement :  shiny::runApp("shiny-app")   (depuis la racine du package)
# =============================================================================

library(shiny)
local({
  root <- if (dir.exists("R")) "R" else file.path("..", "R")   # racine du package
  for (f in c("00_linalg", "01_ols", "04_regularisation"))
    source(file.path(root, paste0(f, ".R")), local = FALSE)
})

# --- Générateurs de données --------------------------------------------------
gen_sparse <- function(n, p, s, snr, seed = 1) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  beta <- c(rep(1, s), rep(0, p - s))
  sig <- sqrt(sum((X %*% beta)^2) / n) / snr
  y <- as.numeric(X %*% beta + sig * rnorm(n))
  list(X = X, y = y, beta = beta, sigma2 = sig^2)
}

# =============================================================================
# UI
# =============================================================================
ui <- fluidPage(
  titlePanel("ML from Scratch in R — laboratoire interactif"),
  tabsetPanel(

    # --- Onglet 1 : biais-variance ridge ---
    tabPanel("1. Biais-variance (ridge)",
      sidebarLayout(
        sidebarPanel(
          helpText("Décomposition exacte de l'EQM du ridge (Module 4, éq. 4.4) le",
                   "long de la pénalité \\(\\lambda\\). L'optimum arbitre biais² vs variance."),
          sliderInput("bv_n", "n (observations)", 30, 500, 100, step = 10),
          sliderInput("bv_p", "p (variables)", 2, 80, 20, step = 1),
          sliderInput("bv_snr", "rapport signal/bruit", 0.5, 5, 2, step = 0.1),
          helpText("La ligne verticale marque le \\(\\lambda\\) qui minimise l'EQM.")),
        mainPanel(plotOutput("bv_plot", height = "460px"), verbatimTextOutput("bv_txt")))),

    # --- Onglet 2 : chemins de régularisation ---
    tabPanel("2. Chemins de régularisation",
      sidebarLayout(
        sidebarPanel(
          helpText("Chemins des coefficients quand \\(\\lambda\\) varie. Le lasso",
                   "(Module 4) met des coefficients EXACTEMENT à zéro (sélection) ;",
                   "le ridge les rétrécit sans les annuler."),
          radioButtons("path_method", "Méthode", c("Lasso" = "lasso", "Ridge" = "ridge")),
          sliderInput("path_p", "p (variables)", 4, 40, 12, step = 1),
          sliderInput("path_s", "s (variables actives)", 1, 10, 4, step = 1),
          sliderInput("path_n", "n (observations)", 30, 400, 120, step = 10)),
        mainPanel(plotOutput("path_plot", height = "460px")))),

    # --- Onglet 3 : orthogonalisation DML ---
    tabPanel("3. Orthogonalisation (DML)",
      sidebarLayout(
        sidebarPanel(
          helpText("Modèle partiellement linéaire \\(Y=\\theta D+g(X)+\\varepsilon\\),",
                   "\\(D=m(X)+v\\). L'estimateur NAÏF (régresser Y sur D seul) est",
                   "biaisé par la confusion \\(g,m\\) ; l'estimateur ORTHOGONAL",
                   "(partialling-out de Neyman, Modules 1/16) l'élimine."),
          sliderInput("dml_conf", "force de confusion", 0, 3, 1.5, step = 0.1),
          sliderInput("dml_n", "n (observations)", 100, 1000, 400, step = 50),
          sliderInput("dml_p", "p (covariables)", 2, 20, 6, step = 1),
          actionButton("dml_go", "Rééchantillonner")),
        mainPanel(plotOutput("dml_plot", height = "440px"), verbatimTextOutput("dml_txt"))))
  )
)

# =============================================================================
# SERVER
# =============================================================================
server <- function(input, output, session) {

  # --- Onglet 1 ---
  output$bv_plot <- renderPlot({
    dat <- gen_sparse(input$bv_n, input$bv_p, s = max(1, input$bv_p %/% 2), snr = input$bv_snr)
    Xs <- scale(dat$X)
    lams <- 10^seq(-3, 4, length.out = 120)
    dec <- t(sapply(lams, function(l) unlist(ridge_bias_var(Xs, dat$beta, dat$sigma2, l)[c("bias2", "variance", "mse")])))
    lopt <- lams[which.min(dec[, "mse"])]
    matplot(log10(lams), dec, type = "l", lwd = 2, lty = 1,
            col = c("#d1495b", "#00798c", "black"),
            xlab = expression(log[10](lambda)), ylab = "composante de l'EQM",
            main = "Décomposition biais-variance du ridge")
    abline(v = log10(lopt), lty = 2, col = "grey50")
    legend("top", c("biais²", "variance", "EQM"), col = c("#d1495b", "#00798c", "black"),
           lwd = 2, bty = "n", horiz = TRUE)
  })
  output$bv_txt <- renderText({
    dat <- gen_sparse(input$bv_n, input$bv_p, s = max(1, input$bv_p %/% 2), snr = input$bv_snr)
    Xs <- scale(dat$X)
    lams <- 10^seq(-3, 4, length.out = 120)
    mse <- sapply(lams, function(l) ridge_bias_var(Xs, dat$beta, dat$sigma2, l)$mse)
    mse0 <- ridge_bias_var(Xs, dat$beta, dat$sigma2, 1e-8)$mse
    sprintf("EQM(lambda=0) = %.3f | EQM optimale = %.3f | gain = %.0f %%",
            mse0, min(mse), 100 * (1 - min(mse) / mse0))
  })

  # --- Onglet 2 ---
  output$path_plot <- renderPlot({
    dat <- gen_sparse(input$path_n, input$path_p, input$path_s, snr = 2)
    Xs <- scale(dat$X); yc <- dat$y - mean(dat$y)
    lams <- 10^seq(-2, 1.2, length.out = 40) * max(abs(crossprod(Xs, yc))) / nrow(Xs)
    B <- sapply(lams, function(l) {
      if (input$path_method == "lasso")
        lasso_fit(Xs, yc, lambda = l * nrow(Xs), standardize = FALSE, intercept = FALSE)$beta
      else ridge_fit(Xs, yc, lambda = l * nrow(Xs), standardize = FALSE, intercept = FALSE)$coefficients
    })
    cols <- ifelse(dat$beta != 0, "#d1495b", "grey70")
    matplot(log10(lams), t(B), type = "l", lty = 1, lwd = 1.5, col = cols,
            xlab = expression(log[10](lambda)), ylab = "coefficient",
            main = paste0("Chemin de régularisation — ", toupper(input$path_method)))
    abline(h = 0, col = "grey40", lty = 3)
    legend("topright", c("actif (vrai != 0)", "nul (vrai = 0)"),
           col = c("#d1495b", "grey70"), lwd = 2, bty = "n")
  })

  # --- Onglet 3 ---
  dml_data <- eventReactive(list(input$dml_go, input$dml_conf, input$dml_n, input$dml_p), {
    n <- input$dml_n; p <- input$dml_p; cf <- input$dml_conf
    X <- matrix(rnorm(n * p), n, p)
    gX <- cf * (X[, 1] + 0.5 * X[, 2]^2); mX <- cf * 0.7 * X[, 1]
    d <- mX + rnorm(n); y <- 1.0 * d + gX + rnorm(n)      # theta0 = 1
    list(X = X, y = y, d = d)
  }, ignoreNULL = FALSE)

  dml_fit <- reactive({
    dd <- dml_data(); X <- dd$X; y <- dd$y; d <- dd$d
    naive <- coef(lm(y ~ d))[2]
    # partialling-out par lasso (résidualisation de Y et D sur X)
    Xs <- scale(X)
    lam <- 0.05 * max(abs(crossprod(Xs, y - mean(y))))
    ry <- (y - mean(y)) - as.numeric(Xs %*% lasso_fit(Xs, y - mean(y), lambda = lam, standardize = FALSE, intercept = FALSE)$beta)
    rd <- (d - mean(d)) - as.numeric(Xs %*% lasso_fit(Xs, d - mean(d), lambda = lam, standardize = FALSE, intercept = FALSE)$beta)
    ortho <- sum(rd * ry) / sum(rd^2)
    list(naive = as.numeric(naive), ortho = ortho)
  })

  output$dml_plot <- renderPlot({
    f <- dml_fit()
    barplot(c(`naïf (Y~D)` = f$naive, `orthogonal (Neyman)` = f$ortho),
            col = c("#d1495b", "#00798c"), ylim = range(0, 1, f$naive, f$ortho) * 1.15,
            ylab = expression(hat(theta)), main = "Estimation de l'effet causal theta")
    abline(h = 1.0, lty = 2, lwd = 2); text(0.2, 1.05, "vrai theta = 1", pos = 4)
  })
  output$dml_txt <- renderText({
    f <- dml_fit()
    sprintf("naïf = %.3f (biais %+.3f) | orthogonal = %.3f (biais %+.3f)  [vrai theta = 1]",
            f$naive, f$naive - 1, f$ortho, f$ortho - 1)
  })
}

shinyApp(ui, server)
