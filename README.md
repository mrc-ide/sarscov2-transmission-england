# sarscov2-transmission-england

This is an [orderly](https://www.vaccineimpact.org/orderly/) repository which contains the analysis to our preprint

> The 2020 SARS-CoV-2 epidemic in England: key epidemiological drivers and impact of interventions; Knock, Whittles, Lees, Perez-Guzman et al.

Running the complete analysis will take many hours, and you may want to run `regional_fit` in parallel on separate nodes of an HPC.  The tasks are parallelised using OpenMP and will consume as many cores as you have available on each node, respecting the environment variable `MC_CORES`.  See the `regional_fit` task for more information.  On 32 core Xeons (dual 16-core 2.6 GHz), each copy of the `regional_fit` task takes about 12 hours, and this must be run 7 times, one for each region.

## Requirements

The core requirement is our [sircovid](https://mrc-ide.github.io/sircovid/) package and its dependencies. Because that package is in constant development you will probably want to pin your versions of the software to the versions we used for preparation:

```r
remotes::install_github(c(
  "mrc-ide/dust@v0.8.19",
  "mrc-ide/mcstate@v0.5.11",
  "mrc-ide/sircovid@v0.10.23"))
```

However, you can always install the versions that we are using with

```r
drat:::add("ncov-ic")
install.packages("sircovid")
```

You will also need a recent [orderly](https://www.vaccineimpact.org/orderly/) which can be installed with

```r
drat:::add("vimc")
install.packages("orderly")
```

To install _all_ required packages, you can use remotes:

```r
remotes::install_deps()
```

## License

MIT © Imperial College of Science, Technology and Medicine
