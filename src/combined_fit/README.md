## Combined region fits

To run this task you must have compiled each of the 7 individual region fits (see [`regional_fit`](../regional_fit)). The `short_run` parameter can be used to select beween short (debug) and full runs.

```r
orderly::orderly_run("combined_fit",
                     parameters = list(short_run = TRUE))
```
