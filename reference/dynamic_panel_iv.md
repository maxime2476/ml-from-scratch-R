# Panel dynamique par variables instrumentales (Anderson-Hsiao / Arellano-Bond)

Dans \\y\_{it}=\rho y\_{i,t-1}+\alpha_i+\varepsilon\_{it}\\,
l'estimateur **within** (effets fixes) est BIAISE (\\y\_{i,t-1}\\
demeaне est correle a l'erreur demeanee : **biais de Nickell**). On
differencie pour eliminer \\\alpha_i\\ puis on instrumente \\\Delta
y\_{i,t-1}\\ par le **niveau retarde** \\y\_{i,t-2}\\ (valide car non
correle a \\\Delta\varepsilon\\).

## Usage

``` r
dynamic_panel_iv(data, id = "id", time = "time", y = "y")
```

## Arguments

- data:

  data.frame

- id, time, y:

  noms des colonnes.

## Value

liste : `rho` (coefficient autoregressif).
