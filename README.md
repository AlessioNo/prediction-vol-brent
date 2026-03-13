# Modélisation et Prévision de la Volatilité du Pétrole Brent
### Approche comparative GARCH(1,1) vs EGARCH(1,1) — 2016-2026

---

## Description

Ce projet modélise la dynamique de la volatilité des rendements journaliers du pétrole Brent
sur la période mars 2016 – mars 2026, en comparant un modèle GARCH(1,1) symétrique et
un modèle EGARCH(1,1) asymétrique, tous deux estimés sous distribution de Student-t.
L'évaluation out-of-sample repose sur une procédure d'expanding window (253 prévisions,
mars 2025 – mars 2026).

---

## Données

- **Source** : FRED – Federal Reserve Economic Data
- **Série** : DCOILBRENTEU (prix quotidien du Brent en USD/baril)
- **Période** : 2 mars 2016 → 2 mars 2026
- **Observations** : 2609 brutes → 70 NA supprimés → 2538 log-rendements

---

## Structure du projet
```
├── DCOILBRENTEU.csv                  # Données brutes FRED
├── brent_volatilite_garch.R          # Script principal
├── appendix_model_selection.R        # Grille 72 modèles eGARCH
├── brent_expanding_window.R          # Prévision out-of-sample
├── tableaux_resultats.R              # Tableaux flextable
└── brent_volatilite_memoire.Rmd      # Rapport PDF
```

---

## Scripts

### `brent_volatilite_garch.R` — Script principal
- Chargement et nettoyage des données
- Calcul des log-rendements
- Analyse exploratoire (prix, rendements, histogramme, ACF/PACF)
- Tests de stationnarité (ADF, KPSS)
- Test ARCH-LM (Engle, 1982)
- Estimation GARCH(1,1) et EGARCH(1,1) sous Student-t
- Diagnostics des résidus (Ljung-Box, Jarque-Bera, Sign Bias Test)
- Comparaison AIC/BIC/MSE/Diebold-Mariano

### `appendix_model_selection.R` — Sélection automatique
- Grille de 72 modèles eGARCH
- Variation des ordres ARMA (p,q ∈ {0,1,2}), GARCH ({(1,1),(1,2),(2,1),(2,2)})
- Distributions : normale vs Student-t
- Scatter AIC vs BIC — validation du choix ARMA(0,0) et Student-t

### `brent_expanding_window.R` — Évaluation out-of-sample
- Coupure : train (mars 2016 – février 2025) / test (mars 2025 – mars 2026)
- 253 ré-estimations eGARCH(1,1), horizon h = 1 jour
- Proxy volatilité réalisée : VR_10j = √(mean(r²)) sur 10 jours glissants
- Métriques : MAE, RMSE, Corrélation, QLIKE

### `tableaux_resultats.R` — Tableaux académiques
- Tableau GARCH(1,1) — coefficients, t-stats, p-values
- Tableau EGARCH(1,1) — idem + paramètre γ
- Tableau comparaison GARCH vs EGARCH
- Tableau métriques out-of-sample

---

## Résultats principaux

| Critère | GARCH(1,1) | EGARCH(1,1) | Verdict |
|---|---|---|---|
| AIC | -4.8657 | -4.8702 | EGARCH ✓ |
| BIC | -4.8542 | -4.8564 | EGARCH ✓ |
| Log-Vraisemblance | 6179.60 | 6186.31 | EGARCH ✓ |
| MSE | 0.00008308 | 0.00008206 | EGARCH ✓ |
| Persistance | 0.956 | 0.970 | — |
| Demi-vie | 15.5 jours | 22.7 jours | — |
| γ (asymétrie) | — | +0.171*** | — |

**Corrélation OOS (σ_t vs VR_10j) : 0.731**

---

## Packages R requis
```r
install.packages(c(
  "tidyverse", "lubridate", "xts", "zoo",
  "rugarch", "moments", "gridExtra", "tseries",
  "FinTS", "forecast", "patchwork", "flextable"
))
```

---

## Ordre d'exécution
```r
# 1. Script principal
source("brent_volatilite_garch.R")

# 2. Sélection de modèle (optionnel — long)
source("appendix_model_selection.R")

# 3. Expanding window (~ 3 minutes)
source("brent_expanding_window.R")

# 4. Tableaux
source("tableaux_resultats.R")

# 5. Générer le rapport PDF
rmarkdown::render("brent_volatilite_memoire.Rmd")
```

---

## Auteurs

**Nocera Alessio & Brochud Corentin**
Master Économétrie — 2025/2026
