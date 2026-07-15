# Estimateur within (effets fixes) d'un panel dynamique — pour illustrer le biais

Estimateur within (effets fixes) d'un panel dynamique — pour illustrer
le biais

## Usage

``` r
dynamic_panel_fe(data, id = "id", time = "time", y = "y")
```

## Arguments

- data, id, time, y:

  cf. `dynamic_panel_iv`.

## Value

le coefficient AR estime (biaise : Nickell).
