## Individual region counterfactuals

These require the corresponding [`regional_fit`](../regional_fit) to have been run, matching both `short_run` and `region`.

```r
orderly::orderly_run("regional_counterfactual",
                     parameters = list(region = "london", short_run = TRUE))
```
