## Combined counterfactuals

To run this task you must have compiled each of the 7 individual region counterfactuals (see [`regional_counterfactual`](../regional_counterfactual)). The `short_run` parameter can be used to select beween short (debug) and full runs.

```r
orderly::orderly_run("combined_counterfactual",
                     parameters = list(short_run = TRUE))
```
