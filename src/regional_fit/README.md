## Individual region fits

These fits take quite a long time to run. If you just want to run representative samples in order to be able to run the combined task, you probably want to use the `short_run = TRUE` parameter. The region parameter must be given.

Particles will be run in parallel based on the environment variable `MC_CORES` (defaulting to 2). This only has an effect if you compiled sircovid with OpenMP support, and if your machine supports it. This should be straightforward on Windows and Linux, but may require some extra work ([see the R docs](https://mac.r-project.org/openmp/), though installation with brew may be easier)

```
orderly::orderly_run("regional_fit",
                     parameters = list(region = "london", short_run = TRUE))
```
