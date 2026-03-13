library(tidyverse)
library(lubridate)
library(xts)
library(zoo)
library(rugarch)
library(moments)
library(gridExtra)
library(tseries)
library(FinTS)
library(forecast)
library(patchwork)
library(flextable)


# ===========================================================================
# 1. CHARGEMENT ET NETTOYAGE DES DONNÉES
# ===========================================================================

brent_raw <- read_csv(
  "DCOILBRENTEU.csv",
  col_types = cols(
    DATE         = col_character(),
    DCOILBRENTEU = col_double()
  )
)

head(brent_raw, 5)

brent_raw <- brent_raw %>%
  rename(Date = observation_date,
         Prix = DCOILBRENTEU) %>%
  mutate(Date = ymd(Date))

brent_clean <- brent_raw %>%
  filter(!is.na(Prix)) %>%
  arrange(Date)


# ===========================================================================
# 2. CALCUL DES LOG-RENDEMENTS
# ===========================================================================

brent_clean <- brent_clean %>%
  mutate(
    Log_Prix  = log(Prix),
    Rendement = c(NA, diff(Log_Prix))
  ) %>%
  filter(!is.na(Rendement))

rendements_xts <- xts(brent_clean$Rendement, order.by = brent_clean$Date)
colnames(rendements_xts) <- "Rendement"

prix_xts <- xts(brent_clean$Prix, order.by = brent_clean$Date)
colnames(prix_xts) <- "Prix"


# ===========================================================================
# 3. ANALYSE EXPLORATOIRE (EDA)
# ===========================================================================

# ---- 3.1  Prix -----------------------------------------------
g1 <- ggplot(brent_clean, aes(x = Date, y = Prix)) +
  geom_line(colour = "#2c7bb6", linewidth = 0.6) +
  labs(title    = "Figure 1: Prix quotidiens du Brent",
       subtitle = "Source : FRED - DCOILBRENTEU",
       x = "Date", y = "Prix (USD/baril)") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

# ---- 3.2  Log-rendements -------------------------------------
g2 <- ggplot(brent_clean, aes(x = Date, y = Rendement)) +
  geom_line(colour = "#d7191c", linewidth = 0.4, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  labs(title = "Figure 2: Log-rendements journaliers du Brent",
       x = "Date", y = "Log-rendement") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

g1 + g2


# ---- 3.3  Histogramme + densité normale ----------------------
mu_r <- mean(brent_clean$Rendement)
sd_r <- sd(brent_clean$Rendement)

g3 <- ggplot(brent_clean, aes(x = Rendement)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 80, fill = "#fdae61", colour = "white", alpha = 0.85) +
  stat_function(fun  = dnorm, args = list(mean = mu_r, sd = sd_r),
                colour = "#2c7bb6", linewidth = 1.1, linetype = "dashed") +
  labs(title    = "Figure 3: Distribution des log-rendements du Brent",
       subtitle = "Courbe bleue = densité normale théorique",
       x = "Log-rendement", y = "Densité") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))
print(g3)

# ---- 3.3b  Zoom sur les queues (axe Y limité) ----------------
g3b <- ggplot(brent_clean, aes(x = Rendement)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 80, fill = "#fdae61", colour = "white", alpha = 0.85) +
  stat_function(
    fun  = dnorm,
    args = list(mean = mu_r, sd = sd_r),
    colour = "#2c7bb6", linewidth = 1.1, linetype = "dashed"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Figure 4: Zoom sur les queues",
       x = "Log-rendement", y = "Densité") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))
print(g3b)

g3 + g3b


# ---- 3.4  Statistiques descriptives --------------------------
tableau_stats <- data.frame(
  Indicateur = c("Moyenne", "Variance", "Ecart-type", "Skewness", "Kurtosis"),
  Valeur = c(
    sprintf("%+.6f", mu_r),
    sprintf("%.6f",  var(brent_clean$Rendement)),
    sprintf("%.6f",  sd_r),
    sprintf("%+.4f", skewness(brent_clean$Rendement)),
    sprintf("%.4f",  kurtosis(brent_clean$Rendement))
  )
)

mon_tableau <- flextable(tableau_stats) %>% set_caption(caption = "Tableau 1 : Statistiques descriptives") %>% 
  theme_booktabs() %>%
  autofit() %>%
  align(align = "center", part = "all") %>%
  bold(part = "header")

mon_tableau


# ---- 3.5  ACF / PACF des rendements et rendements2 -----------
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 3, 0))
acf(coredata(rendements_xts),
    lag.max = 30, main = "ACF des rendements", col = "#2c7bb6", lwd = 2)
pacf(coredata(rendements_xts),
     lag.max = 30, main = "PACF des rendements", col = "#d7191c", lwd = 2)
acf(coredata(rendements_xts)^2,
    lag.max = 30, main = "ACF des rendements2 (effet ARCH)", col = "#1a9641", lwd = 2)
pacf(coredata(rendements_xts)^2,
     lag.max = 30, main = "PACF des rendements2", col = "#fdae61", lwd = 2)

mtext("Figure 5 : Analyse de l'autocorrélation et des effets ARCH", 
      side = 3, outer = TRUE, cex = 1, font = 2)

par(mfrow = c(1, 1))





lb_test <- Box.test(coredata(rendements_xts)^2, lag = 10, type = "Ljung-Box")
print(lb_test)


# ---- 3.6  Tests de stationnarité (ADF + KPSS) ----------------
# ADF  - H0 : racine unitaire -> rejet si p < 0.05 -> stationnaire
# KPSS - H0 : stationnaire   -> non-rejet si p > 0.05 -> stationnaire

adf_prix <- adf.test(coredata(log(prix_xts)), alternative = "stationary")
print(adf_prix)

adf_rend <- adf.test(coredata(rendements_xts), alternative = "stationary")
print(adf_rend)

kpss_rend <- suppressWarnings(kpss.test(coredata(rendements_xts)))
print(kpss_rend)


# ---- 3.7  Test ARCH-LM (Engle, 1982) -------------------------
# H0 : pas d'effet ARCH -> rejet si p < 0.05 -> GARCH justifie

arch_test <- ArchTest(coredata(rendements_xts), lags = 10)
print(arch_test)

jarque.bera.test(coredata(rendements_xts))


# ===========================================================================
# 4. MODELISATION GARCH(1,1)
# ===========================================================================

# ---- 4.1  Specification ----------------------------------------------
spec_garch <- ugarchspec(
  variance.model     = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "std"
)

# ---- 4.2  Estimation -------------------------------------------------
fit_garch <- ugarchfit(spec = spec_garch, data = rendements_xts)
print(fit_garch)


# ---- 4.3  Coefficients -----------------------------------------------
coef_garch     <- coef(fit_garch)
alpha_g        <- coef_garch["alpha1"]
beta_g         <- coef_garch["beta1"]
demi_vie_garch <- log(0.5) / log(alpha_g + beta_g)

cat(sprintf("Persistance (alpha + beta) = %.4f\n", alpha_g + beta_g))
cat(sprintf("Demi-vie de la volatilite  = %.1f jours\n", demi_vie_garch))


# ---- 4.4  Volatilite conditionnelle ----------------------------------
sigma_garch <- sigma(fit_garch)

sigma_garch_df <- data.frame(
  Date   = index(sigma_garch),
  Sigma  = as.numeric(sigma_garch),
  Modele = "GARCH(1,1)"
)


# ---- 4.5  Diagnostics des residus ------------------------------------
residus_std_garch <- residuals(fit_garch, standardize = TRUE)

print(Box.test(coredata(residus_std_garch),   lag = 10, type = "Ljung-Box"))
print(Box.test(coredata(residus_std_garch)^2, lag = 10, type = "Ljung-Box"))


# ---- 4.6  Test de Jarque-Bera - residus GARCH(1,1) ------------------
# H0 : normalite des residus (skewness = 0, kurtosis = 3)

jb_garch <- jarque.bera.test(coredata(residus_std_garch))
print(jb_garch)

sk_garch <- skewness(coredata(residus_std_garch))
ku_garch <- kurtosis(coredata(residus_std_garch))
cat(sprintf("Skewness : %.4f | Kurtosis : %.4f | JB p-value : %.4e\n",
            sk_garch, ku_garch, jb_garch$p.value))


# ===========================================================================
# 5. MODELISATION EGARCH(1,1)
# ===========================================================================

# ---- 5.1  Specification ----------------------------------------------
spec_egarch <- ugarchspec(
  variance.model     = list(model = "eGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "std"
)

# ---- 5.2  Estimation -------------------------------------------------
fit_egarch <- ugarchfit(spec = spec_egarch, data = rendements_xts)
print(fit_egarch)


# ---- 5.3  Coefficients -----------------------------------------------
coef_egarch <- coef(fit_egarch)

if ("gamma1" %in% names(coef_egarch)) {
  cat(sprintf("Parametre de levier gamma = %.4f\n", coef_egarch["gamma1"]))
  if (coef_egarch["gamma1"] < 0) {
    cat("-> gamma < 0 : les chocs negatifs augmentent davantage la volatilite.\n")
  } else {
    cat("-> gamma > 0 : les chocs positifs augmentent davantage la volatilite.\n")
  }
}

beta_e          <- coef_egarch["beta1"]
demi_vie_egarch <- log(0.5) / log(abs(beta_e))
cat(sprintf("Persistance beta = %.4f | Demi-vie = %.1f jours\n", beta_e, demi_vie_egarch))


# ---- 5.4  Volatilite conditionnelle EGARCH ---------------------------
sigma_egarch <- sigma(fit_egarch)

sigma_egarch_df <- data.frame(
  Date   = index(sigma_egarch),
  Sigma  = as.numeric(sigma_egarch),
  Modele = "EGARCH(1,1)"
)


# ---- 5.5  Diagnostics EGARCH -----------------------------------------
residus_std_egarch <- residuals(fit_egarch, standardize = TRUE)


print(Box.test(coredata(residus_std_egarch),   lag = 10, type = "Ljung-Box"))
print(Box.test(coredata(residus_std_egarch)^2, lag = 10, type = "Ljung-Box"))


# ---- 5.6  Test de Jarque-Bera - residus EGARCH(1,1) -----------------
jb_egarch <- jarque.bera.test(coredata(residus_std_egarch))
print(jb_egarch)

sk_egarch <- skewness(coredata(residus_std_egarch))
ku_egarch <- kurtosis(coredata(residus_std_egarch))

cat(sprintf("GARCH  : skewness = %.4f | kurtosis = %.4f | JB p = %.4e\n",
            sk_garch, ku_garch, jb_garch$p.value))
cat(sprintf("EGARCH : skewness = %.4f | kurtosis = %.4f | JB p = %.4e\n",
            sk_egarch, ku_egarch, jb_egarch$p.value))


# ===========================================================================
# 6. COMPARAISON GARCH vs EGARCH
# ===========================================================================

# ---- 6.1  Criteres d'information -------------------------------------
aic_garch  <- infocriteria(fit_garch)[1]
bic_garch  <- infocriteria(fit_garch)[2]
aic_egarch <- infocriteria(fit_egarch)[1]
bic_egarch <- infocriteria(fit_egarch)[2]

cat(sprintf("              GARCH(1,1)   EGARCH(1,1)\n"))
cat(sprintf("AIC         : %+9.4f   %+9.4f\n", aic_garch,  aic_egarch))
cat(sprintf("BIC         : %+9.4f   %+9.4f\n", bic_garch,  bic_egarch))
cat(sprintf("Log-Lik     : %+9.2f   %+9.2f\n",
            likelihood(fit_garch), likelihood(fit_egarch)))


# ---- 6.2  Graphique comparatif ---------------------------------------
sigma_compare <- rbind(sigma_garch_df, sigma_egarch_df)

g6 <- ggplot(sigma_compare, aes(x = Date, y = Sigma, colour = Modele)) +
  geom_line(linewidth = 0.45, alpha = 0.85) +
  scale_colour_manual(values = c("GARCH(1,1)" = "#2c7bb6", "EGARCH(1,1)" = "#d7191c")) +
  labs(title    = "Figure 6: Superposition des volatilites conditionnelles",
       subtitle = "GARCH(1,1) vs EGARCH(1,1) - Brent (2016-2026)",
       x = "Date", y = "sigma_t", colour = "Modele") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top")
print(g6)


# ---- 6.3  Graphiques cote-a-cote -------------------------------------
p_garch  <- ggplot(sigma_garch_df,  aes(x = Date, y = Sigma)) +
  geom_line(colour = "#2c7bb6", linewidth = 0.5) +
  labs(title = "GARCH(1,1)", x = NULL, y = "sigma_t") +
  theme_minimal(base_size = 11)

p_egarch <- ggplot(sigma_egarch_df, aes(x = Date, y = Sigma)) +
  geom_line(colour = "#d7191c", linewidth = 0.5) +
  labs(title = "EGARCH(1,1)", x = NULL, y = "sigma_t") +
  theme_minimal(base_size = 11)

gridx <- grid.arrange(p_garch, p_egarch, ncol = 1,
                      top = "Volatilite conditionnelle : comparaison")


# ---- 6.4  Test de Diebold-Mariano ------------------------------------
# H0 : pouvoir predictif identique entre GARCH et EGARCH
# Erreur : e_t = r2_t - sigma2_t (variance realisee vs prevue)

r2_insample <- coredata(rendements_xts)^2
e_garch_dm  <- as.numeric(r2_insample) - as.numeric(sigma_garch)^2
e_egarch_dm <- as.numeric(r2_insample) - as.numeric(sigma_egarch)^2

dm_result <- dm.test(e_garch_dm, e_egarch_dm,
                     alternative = "two.sided",
                     h = 1, power = 2)
print(dm_result)

mse_garch  <- mean(e_garch_dm^2)
mse_egarch <- mean(e_egarch_dm^2)
cat(sprintf("MSE GARCH  : %.8f\n", mse_garch))
cat(sprintf("MSE EGARCH : %.8f\n", mse_egarch))
cat(sprintf("Ratio MSE  : %.4f  (< 1 -> EGARCH meilleur)\n", mse_egarch / mse_garch))

