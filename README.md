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
├── DCOILBRENTEU.csv       # Données brutes FRED
├── script_final.R         # Script principal
└── forecast_final.R       # Sélection de modèle + prévision out-of-sample
```

---

## Scripts

### `script_final.R` — Script principal
- Chargement et nettoyage des données
- Calcul des log-rendements
- Analyse exploratoire : prix (Fig.1), rendements (Fig.2), histogramme (Fig.3-4), ACF/PACF (Fig.5)
- Statistiques descriptives — tableau flextable
- Tests de stationnarité (ADF sur log-prix et rendements, KPSS sur rendements)
- Test ARCH-LM (Engle, 1982) + Ljung-Box sur rendements²
- Estimation GARCH(1,1) — Student-t, ARMA(0,0)
- Estimation EGARCH(1,1) — Student-t, ARMA(0,0)
- Diagnostics des résidus : Ljung-Box, Jarque-Bera
- Comparaison AIC / BIC / Log-vraisemblance / MSE / Diebold-Mariano
- Figure 6 : superposition volatilités conditionnelles GARCH vs EGARCH

### `forecast_final.R` — Sélection de modèle + prévision out-of-sample
**Partie 1 — Grille de sélection (72 modèles)**
- Variation des ordres ARMA (p,q ∈ {0,1,2}), GARCH ({(1,1),(1,2),(2,1),(2,2)})
- Distributions : normale vs Student-t
- Top 5 AIC et Top 5 BIC
- Figure scatter AIC vs BIC — validation ARMA(0,0) et Student-t

**Partie 2 — Expanding window**
- Coupure : train (mars 2016 – février 2025) / test (mars 2025 – mars 2026)
- 253 ré-estimations eGARCH(1,1), horizon h = 1 jour
- Proxy volatilité réalisée : VR_10j = √(mean(r²)) sur 10 jours glissants
- Métriques : MAE, RMSE, Corrélation σ/VR_10j, Corrélation σ/|r_t|, QLIKE
- Figure 7 : σ_t prévu vs VR_10j
- Figure 8 : erreurs de prévision + tendance loess
- Figure 9 : scatter σ_t vs VR_10j
- Vue complète 2016-2026 : in-sample + out-of-sample

---

## Résultats principaux

| Critère            | GARCH(1,1)  | EGARCH(1,1) | Verdict  |
|--------------------|-------------|-------------|----------|
| AIC                | -4.8657     | -4.8702     | EGARCH ✓ |
| BIC                | -4.8542     | -4.8564     | EGARCH ✓ |
| Log-Vraisemblance  | 6179.60     | 6186.31     | EGARCH ✓ |
| MSE                | 0.00008308  | 0.00008206  | EGARCH ✓ |
| Persistance        | 0.956       | 0.970       | —        |
| Demi-vie           | 15.5 jours  | 22.7 jours  | —        |
| γ (asymétrie)      | —           | +0.171***   | —        |

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
# 1. Script principal — données, modèles, diagnostics, comparaison
source("script_final.R")

# 2. Sélection de modèle + expanding window (~3 minutes)
source("forecast_final.R")
```

> ⚠️ `forecast_final.R` recharge les données automatiquement.
> Il peut être lancé indépendamment de `script_final.R`.

---

## Auteurs

**Nocera Alessio & Brochud Corentin**
Master Économétrie — 2025/2026
