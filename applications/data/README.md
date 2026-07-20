# Données LaLonde (NSW) — provenance

Deux jeux embarqués pour rendre le chapitre `applications/lalonde.qmd`
reproductible sans dépendance externe.

| Fichier | Lignes | Contenu |
|---|---|---|
| `lalonde_exp.rds` | 445 | Échantillon **expérimental** NSW (185 traités, 260 contrôles randomisés) — l'étalon-or. |
| `lalonde_psid.rds` | 2675 | Traités NSW + groupe de comparaison **non expérimental** PSID-1. |

**Variables** : `treat`, `age`, `education`, `black`, `hispanic`, `married`,
`nodegree`, `re74`, `re75`, `re78` (revenus réels 1974/75/78).

## Source

Données de LaLonde (1986), échantillon de Dehejia & Wahba (1999, 2002).
Extraites du package R **`qte`** (Callaway ; jeux `lalonde.exp` et
`lalonde.psid`), **archivé de CRAN** en 2026 — d'où leur inclusion directe ici.
Les données elles-mêmes sont publiques et largement diffusées (identiques à
`Matching::lalonde`, `causalsens`, etc.).

## Référence

- LaLonde, R. (1986). *Evaluating the Econometric Evaluations of Training
  Programs with Experimental Data.* American Economic Review, 76(4), 604-620.
- Dehejia, R. & Wahba, S. (1999, 2002). *Causal Effects in Nonexperimental
  Studies.*
