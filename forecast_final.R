library(tidyverse)
library(lubridate)
library(xts)
library(rugarch)
library(gridExtra)
library(flextable)

# Décommenter si session vide (rendements_xts absent) :
brent_raw <- read_csv("DCOILBRENTEU.csv") %>%
  rename(Date = observation_date, Prix = DCOILBRENTEU) %>%
  filter(!is.na(Prix)) %>%
  arrange(Date) %>%
  mutate(Log_Prix  = log(Prix),
         Rendement = c(NA, diff(Log_Prix))) %>%
  filter(!is.na(Rendement))

rendements_xts <- xts(brent_raw$Rendement, order.by = brent_raw$Date)
colnames(rendements_xts) <- "Rendement"


# ===========================================================================
# 1. DÉFINITION DE LA GRILLE
# ===========================================================================

arma_grid <- expand.grid(
  arma_p = 0:2,
  arma_q = 0:2
)

garch_grid <- expand.grid(
  garch_p = 1:2,
  garch_q = 1:2
)

variance_models <- c("eGARCH")
distributions   <- c("norm", "std")

total <- nrow(arma_grid) * nrow(garch_grid) * length(distributions)
cat(sprintf("TOTAL modèles à estimer : %d\n\n", total))


# ===========================================================================
# 2. BOUCLE D'ESTIMATION
# ===========================================================================

resultats <- data.frame(
  Modele_Variance  = character(),
  ARMA             = character(),
  GARCH_ordre      = character(),
  Distribution     = character(),
  AIC              = numeric(),
  BIC              = numeric(),
  LogLik           = numeric(),
  Convergence      = logical(),
  stringsAsFactors = FALSE
)

compteur <- 0
n_echecs <- 0

cat("Estimation en cours...\n")
cat(rep("-", 60), "\n", sep = "")

for (vm in variance_models) {
  for (gi in 1:nrow(garch_grid)) {
    gp <- garch_grid$garch_p[gi]
    gq <- garch_grid$garch_q[gi]
    
    for (ai in 1:nrow(arma_grid)) {
      ap <- arma_grid$arma_p[ai]
      aq <- arma_grid$arma_q[ai]
      
      for (dist in distributions) {
        
        compteur     <- compteur + 1
        label_modele <- sprintf("%-7s | ARMA(%d,%d) + %s(%d,%d) | %-4s",
                                vm, ap, aq, vm, gp, gq, dist)
        
        if (compteur %% 20 == 0 || compteur == 1) {
          cat(sprintf("[%3d/%d] %s\n", compteur, total, label_modele))
        }
        
        spec <- tryCatch(
          ugarchspec(
            variance.model     = list(model = vm, garchOrder = c(gp, gq)),
            mean.model         = list(armaOrder = c(ap, aq), include.mean = TRUE),
            distribution.model = dist
          ),
          error = function(e) NULL
        )
        
        if (is.null(spec)) { n_echecs <- n_echecs + 1; next }
        
        fit <- tryCatch(
          ugarchfit(spec = spec, data = rendements_xts,
                    solver = "hybrid", silent = TRUE),
          error   = function(e) NULL,
          warning = function(w) suppressWarnings(
            ugarchfit(spec = spec, data = rendements_xts,
                      solver = "hybrid", silent = TRUE)
          )
        )
        
        converge <- !is.null(fit) && fit@fit$convergence == 0
        
        if (!converge) {
          n_echecs <- n_echecs + 1
          resultats <- rbind(resultats, data.frame(
            Modele_Variance  = vm,
            ARMA             = sprintf("(%d,%d)", ap, aq),
            GARCH_ordre      = sprintf("(%d,%d)", gp, gq),
            Distribution     = dist,
            AIC              = NA_real_,
            BIC              = NA_real_,
            LogLik           = NA_real_,
            Convergence      = FALSE,
            stringsAsFactors = FALSE
          ))
          next
        }
        
        ic     <- infocriteria(fit)
        aic_v  <- ic[1]
        bic_v  <- ic[2]
        loglik <- likelihood(fit)
        
        resultats <- rbind(resultats, data.frame(
          Modele_Variance  = vm,
          ARMA             = sprintf("(%d,%d)", ap, aq),
          GARCH_ordre      = sprintf("(%d,%d)", gp, gq),
          Distribution     = dist,
          AIC              = round(aic_v,  6),
          BIC              = round(bic_v,  6),
          LogLik           = round(loglik, 3),
          Convergence      = TRUE,
          stringsAsFactors = FALSE
        ))
        
      }
    }
  }
}

cat(rep("-", 60), "\n", sep = "")
cat(sprintf("Estimation terminée : %d modèles, %d échecs de convergence.\n\n",
            compteur, n_echecs))


# ===========================================================================
# 3. TABLEAUX COMPARATIFS — TOP 5 AIC ET TOP 5 BIC
# ===========================================================================

res_ok <- resultats %>%
  filter(Convergence == TRUE) %>%
  arrange(AIC)

cat(sprintf("Modèles convergents : %d / %d\n\n", nrow(res_ok), total))

print(
  res_ok %>% arrange(AIC) %>% head(5) %>%
    mutate(Rang = row_number()) %>%
    select(Rang, ARMA, GARCH_ordre, Distribution, AIC, BIC, LogLik),
  row.names = FALSE
)

print(
  res_ok %>% arrange(BIC) %>% head(5) %>%
    mutate(Rang = row_number()) %>%
    select(Rang, ARMA, GARCH_ordre, Distribution, AIC, BIC, LogLik),
  row.names = FALSE
)


# ===========================================================================
# 4. GRAPHIQUE — SCATTER AIC vs BIC
# ===========================================================================

meilleur_aic <- res_ok %>% arrange(AIC) %>% slice(1)
meilleur_bic <- res_ok %>% arrange(BIC) %>% slice(1)

label_aic <- "ARMA(2,2) eGARCH(2,1)"
label_bic <- "ARMA(0,0) eGARCH(1,1)"

res_ok_plot <- res_ok %>%
  mutate(highlight = case_when(
    AIC == meilleur_aic$AIC & BIC == meilleur_aic$BIC ~ "Meilleur AIC",
    AIC == meilleur_bic$AIC & BIC == meilleur_bic$BIC ~ "Meilleur BIC",
    TRUE ~ "Autre"
  ))

g_scatter <- ggplot(res_ok_plot,
                    aes(x = AIC, y = BIC,
                        colour = Distribution,
                        shape  = highlight,
                        size   = highlight)) +
  geom_point(alpha = 0.85) +
  scale_colour_manual(
    values = c("std"  = "#d7191c", "norm" = "#2c7bb6"),
    labels = c("std"  = "Student-t", "norm" = "Normale")
  ) +
  scale_shape_manual(
    values = c("Meilleur AIC" = 17, "Meilleur BIC" = 17, "Autre" = 16),
    guide  = "none"
  ) +
  scale_size_manual(
    values = c("Meilleur AIC" = 3.5, "Meilleur BIC" = 3.5, "Autre" = 2.5),
    guide  = "none"
  ) +
  annotate("text",
           x = meilleur_aic$AIC + 0.0002, y = meilleur_aic$BIC + 0.0005,
           label = paste0("Meilleur AIC ", label_aic),
           colour = "#1a9641", size = 3, hjust = 0) +
  annotate("text",
           x = meilleur_bic$AIC + 0.0002, y = meilleur_bic$BIC - 0.0005,
           label = paste0("Meilleur BIC ", label_bic),
           colour = "#d97b00", size = 3, hjust = 0, vjust = 1) +
  labs(title  = "AIC vs BIC - 72 modèles eGARCH estimés",
       x = "AIC", y = "BIC", colour = "Distribution") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "right")

print(g_scatter)

# =============================================================================
# PRÉVISION OUT-OF-SAMPLE — EXPANDING WINDOW
# Modèle : eGARCH(1,1) - ARMA(0,0) - Student-t
# Train  : mars 2016 -> mars 2025  (~2270 observations)
# Test   : mars 2025 -> mars 2026  (~253 observations)
# Proxy  : Volatilité réalisée 10 jours (VR_10j)
# =============================================================================


# ===========================================================================
# 5. COUPURE TRAIN / TEST
# ===========================================================================

date_coupure <- as.Date("2025-03-01")

r_train <- rendements_xts[index(rendements_xts) <  date_coupure]
r_test  <- rendements_xts[index(rendements_xts) >= date_coupure]

cat(sprintf("Train : %s -> %s  (%d observations)\n",
            format(start(r_train)), format(end(r_train)), nrow(r_train)))
cat(sprintf("Test  : %s -> %s  (%d observations)\n",
            format(start(r_test)),  format(end(r_test)),  nrow(r_test)))

n_test <- nrow(r_test)


# ===========================================================================
# 6. SPÉCIFICATION DU MODÈLE
# ===========================================================================

spec_ew <- ugarchspec(
  variance.model     = list(model = "eGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "std"
)


# ===========================================================================
# 7. BOUCLE EXPANDING WINDOW
# ===========================================================================

cat(sprintf("Expanding window - %d itérations...\n\n", n_test))

sigma_oos <- rep(NA_real_, n_test)
dates_oos <- index(r_test)
r_full    <- rendements_xts

idx_fin_train <- which(index(r_full) < date_coupure)
n_train_init  <- length(idx_fin_train)

t_start <- proc.time()

for (t in 1:n_test) {
  
  n_fenetre <- n_train_init + (t - 1)
  r_fenetre <- r_full[1:n_fenetre]
  
  fit_t <- tryCatch(
    ugarchfit(spec = spec_ew, data = r_fenetre,
              solver = "hybrid", silent = TRUE),
    error   = function(e) NULL,
    warning = function(w) suppressWarnings(
      ugarchfit(spec = spec_ew, data = r_fenetre,
                solver = "hybrid", silent = TRUE)
    )
  )
  
  if (!is.null(fit_t) && fit_t@fit$convergence == 0) {
    fc           <- ugarchforecast(fit_t, n.ahead = 1)
    sigma_oos[t] <- as.numeric(sigma(fc))
  }
  
  if (t %% 25 == 0 || t == 1 || t == n_test) {
    elapsed <- (proc.time() - t_start)["elapsed"]
    pct     <- t / n_test * 100
    restant <- if (t > 1) (elapsed / t) * (n_test - t) else NA
    cat(sprintf("  [%3d/%d]  %.0f%%  |  ecoule : %ds  |  restant : %ds\n",
                t, n_test, pct, round(elapsed),
                ifelse(is.na(restant), 0, round(restant))))
  }
}

elapsed_total <- (proc.time() - t_start)["elapsed"]
cat(sprintf("\nBoucle terminée en %.1f secondes.\n\n", elapsed_total))

n_na <- sum(is.na(sigma_oos))
if (n_na > 0) {
  cat(sprintf("%d prévision(s) NA -> interpolation linéaire.\n\n", n_na))
  sigma_oos <- zoo::na.approx(sigma_oos, na.rm = FALSE)
}


# ===========================================================================
# 8. PROXY DE VOLATILITÉ RÉALISÉE — VR_10j
# ===========================================================================

r_test_vec <- as.numeric(r_test)
r2_test    <- r_test_vec^2
absr_test  <- abs(r_test_vec)

# VR_10j = sqrt(mean(r²)) sur 10 jours glissants — même unité que sigma_t
vol_realisee_10j <- zoo::rollapply(
  r2_test, width = 10,
  FUN  = function(x) sqrt(mean(x)),
  fill = NA, align = "right"
)

df_oos <- data.frame(
  Date        = dates_oos,
  Sigma_prevu = sigma_oos,
  r2          = r2_test,
  AbsR        = absr_test,
  VR_10j      = vol_realisee_10j
)


# ===========================================================================
# 9. MÉTRIQUES D'ÉVALUATION
# ===========================================================================

sigma2_oos <- sigma_oos^2

idx_valide <- !is.na(sigma2_oos) & !is.na(r2_test)
mae_var    <- mean(abs(sigma2_oos[idx_valide] - r2_test[idx_valide]))
rmse_var   <- sqrt(mean((sigma2_oos[idx_valide] - r2_test[idx_valide])^2))

idx_10j  <- !is.na(sigma_oos) & !is.na(vol_realisee_10j)
cor_10j  <- cor(sigma_oos[idx_10j], vol_realisee_10j[idx_10j])
cor_brut <- cor(sigma_oos[idx_valide], absr_test[idx_valide])

# QLIKE — filtre r² < 5e percentile pour éviter divisions par ~0
seuil_qlike <- quantile(r2_test[r2_test > 0], probs = 0.05)
idx_qlike   <- idx_valide & r2_test > seuil_qlike
qlike <- mean(
  sigma2_oos[idx_qlike] / r2_test[idx_qlike] -
    log(sigma2_oos[idx_qlike] / r2_test[idx_qlike]) - 1
)

cat(sprintf("MAE   (sigma2 vs r2)       : %.8f\n", mae_var))
cat(sprintf("RMSE  (sigma2 vs r2)       : %.8f\n", rmse_var))
cat(sprintf("Corrélation sigma / VR_10j : %.4f\n", cor_10j))
cat(sprintf("Corrélation sigma / |r_t|  : %.4f\n", cor_brut))
cat(sprintf("QLIKE corrigé              : %.6f\n", qlike))


# ===========================================================================
# 10. GRAPHIQUES
# ===========================================================================

# ---- 10.1  sigma_t prévu vs VR_10j -----------------------------------------
g_fw1 <- ggplot(df_oos, aes(x = Date)) +
  geom_line(aes(y = VR_10j,      colour = "Vol. réalisée 10j"),
            linewidth = 0.6, alpha = 0.9) +
  geom_line(aes(y = Sigma_prevu, colour = "sigma_t prévu (eGARCH)"),
            linewidth = 0.8) +
  scale_colour_manual(values = c(
    "Vol. réalisée 10j"      = "#999999",
    "sigma_t prévu (eGARCH)" = "#d7191c"
  )) +
  labs(title    = "Figure 7: Prévision de volatilité - Expanding Window",
       subtitle = "eGARCH(1,1) | Période test : mars 2025 -> mars 2026",
       x = "Date", y = "Volatilité (sigma_t)", colour = NULL) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top")
print(g_fw1)

# ---- 10.2  Erreurs de prévision dans le temps --------------------------------
df_oos <- df_oos %>%
  mutate(Erreur = Sigma_prevu^2 - r2)

g_fw2 <- ggplot(df_oos, aes(x = Date, y = Erreur)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_line(colour = "#2c7bb6", linewidth = 0.5, alpha = 0.8) +
  geom_smooth(method = "loess", span = 0.3,
              colour = "#d7191c", se = FALSE, linewidth = 1) +
  labs(title    = "Figure 8: Erreurs de prévision (sigma2_prévu - r2_t)",
       subtitle = "Rouge = tendance loess  |  > 0 : sur-estimation  |  < 0 : sous-estimation",
       x = "Date", y = "Erreur de prévision") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))
print(g_fw2)

# ---- 10.3  Scatter sigma_t prévu vs VR_10j -----------------------------------
g_fw3 <- ggplot(df_oos %>% filter(!is.na(VR_10j)),
                aes(x = Sigma_prevu, y = VR_10j)) +
  geom_point(alpha = 0.5, size = 1.8, colour = "#2c7bb6") +
  geom_smooth(method = "lm", colour = "#d7191c", se = TRUE, linewidth = 1) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey40") +
  labs(title    = "Figure 9: sigma_t prévu vs Volatilité réalisée 10j",
       subtitle = sprintf("Corrélation = %.4f  |  Ligne pointillée = prévision parfaite",
                          cor_10j),
       x = "sigma_t prévu (eGARCH)", y = "Volatilité réalisée 10j") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))
print(g_fw3)

# ---- 10.4  Vue complète train + test -----------------------------------------
fit_full_train <- ugarchfit(spec = spec_ew, data = r_train,
                            solver = "hybrid", silent = TRUE)

sigma_full_df <- rbind(
  data.frame(Date    = as.Date(index(r_train)),
             Sigma   = as.numeric(sigma(fit_full_train)),
             Periode = "In-sample (train)"),
  data.frame(Date    = as.Date(dates_oos),
             Sigma   = sigma_oos,
             Periode = "Out-of-sample (test)")
)

g_fw4 <- ggplot(sigma_full_df, aes(x = Date, y = Sigma, colour = Periode)) +
  geom_line(linewidth = 0.5) +
  geom_vline(xintercept = date_coupure,
             linetype = "dashed", colour = "black", linewidth = 0.8) +
  annotate("text",
           x     = date_coupure,
           y     = max(sigma_full_df$Sigma, na.rm = TRUE) * 0.95,
           label = "  Coupure\n  mars 2025",
           hjust = 0, size = 3.5, colour = "black") +
  scale_colour_manual(values = c(
    "In-sample (train)"    = "#2c7bb6",
    "Out-of-sample (test)" = "#d7191c"
  )) +
  labs(title    = "Volatilité conditionnelle - Vue complète (2016-2026)",
       subtitle = "Bleu = in-sample  |  Rouge = out-of-sample (expanding window)",
       x = "Date", y = "sigma_t", colour = NULL) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top")
print(g_fw4)

