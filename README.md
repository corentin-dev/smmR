# Semi Markov Model

[![PLMLab build status](https://plmlab.math.cnrs.fr/lmrs/statistique/smmR/badges/master/pipeline.svg)](https://plmlab.math.cnrs.fr/lmrs/statistique/smmR/-/pipelines) [![CRAN versions](https://www.r-pkg.org/badges/version/smmR)](https://CRAN.R-project.org/package=smmR) [![CRAN logs](https://cranlogs.r-pkg.org/badges/smmR)](https://CRAN.R-project.org/package=smmR) [![RDocumentation](https://api.rdocumentation.org/badges/version/smmR)](https://lmrs.pages.math.cnrs.fr/statistique/smmR/)

## Description

Performs parametric and non-parametric estimation and simulation for multi-state discrete-time semi-Markov processes.

For the parametric estimation, several discrete distributions are considered for the sojourn 
times:

- Uniform,
- Geometric,
- Poisson,
- Discrete Weibull,
- Negative Binomial.

The non-parametric estimation concerns the sojourn time distributions, where no assumptions are done on the shape of distributions. Moreover, the estimation can be done on the basis of one or several sample paths, with or without censoring at the beginning or/and at the end of the sample paths. 

Reliability indicators such as reliability, maintainability, availability, BMP-failure rate, RG-failure rate, mean time to failure and mean time to repair are available as well.

The implemented methods are described in:

* Barbu, V. S., & Limnios, N. (2009). Semi-Markov chains and hidden semi-Markov models toward applications: their use in reliability and DNA analysis (Vol. 191). _Springer Science & Business Media_. doi:10.1007/978-0-387-73173-5. ([Journal version](https://link.springer.com/book/10.1007/978-0-387-73173-5))

* Barbu, V., & Limnios, N. (2006). Empirical estimation for discrete-time semi-Markov processes with applications in reliability. _Nonparametric Statistics_, 18(7-8), 483-498. doi:10.1080/10485250701261913. ([Journal version](https://www.tandfonline.com/doi/pdf/10.1080/10485250701261913))

* Trevezas, S., & Limnios, N. (2011). Exact MLE and asymptotic properties for nonparametric semi-Markov models. _Journal of Nonparametric Statistics_, 23(3), 719-739. doi:10.1080/10485252.2011.555543. ([Journal version](https://www.tandfonline.com/doi/pdf/10.1080/10485252.2011.555543))

Estimation and simulation of discrete-time k-th order Markov chains are 
also considered.

## Install

* Install from CRAN:

```R
install.packages('smmR')
```

* Install latest development version from `git`:

```R
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_git(
  url = "https://plmlab.math.cnrs.fr/lmrs/statistique/smmR",
  dependencies = TRUE,
  build_vignettes = FALSE)
```

## Quickstart

```R
library(smmR)
```

To create an object `smmparametric`, we need several pieces of information:

- a state space (here `states`)
- an initial distribution (here `alpha`)
- a transition matrix (here `p`)
- a distribution matrix stating the laws (here `dist`)
- parameters for each state

```R
states <- c("a", "b", "c") # state space

alpha <- c(1, 0, 0) # initial distribution
```

We define the transition matrix `p` and the distribution matrix `dist` (see [here](https://lmrs.pages.math.cnrs.fr/statistique/smmR/reference/smmparametric.html?q=geom#arguments) for details):
```R
p <- matrix(data = c(0.0, 0.5, 0.5,
                     1.0, 0.0, 0.0,
                     1.0, 0.0, 0.0),
            nrow = 3, ncol = 3, byrow = TRUE)

distr <- matrix(c(NA, "geom", "geom", 
                  "geom", NA, NA,
                  "geom", NA, NA), 
                nrow = 3, ncol = 3, byrow = TRUE) # Distribution matrix
```

There is only one parameter for the geometrical law:
```
parameters <- array(c( NA, 0.8, 0.8,
                       0.3,  NA,  NA,
                       0.5,  NA,  NA,
                       NA,  NA,  NA,
                       NA,  NA,  NA,
                       NA,  NA,  NA), 
                    c(3, 3, 2))
```

Finally, we create our object `smmparametric`:

```R
mysmmparam <- smmparametric(states = states, init = alpha, ptrans = p, 
                         type.sojourn = "fij", distr = distr, param = parameters)
```

From this `smmparametric` object, we can simulate sequences:
```R
M <- 10000
seq <- simulate(object = mysmmparam, nsim = M)
```

And we can fit another model, here a `smmnonparametric` object:
```R
estimate <- fitsmm(sequences = seq, states = states, type.sojourn = "fij")
```

You can find an application to reliability [here](https://lmrs.pages.math.cnrs.fr/statistique/smmR/articles/Textile-Factory.html).

## Contributing

Contributions to this package are warmly welcome. Do not hesitate to open an issue to discuss new features. 

The official repository is at [PLMLab](https://plmlab.math.cnrs.fr/lmrs/statistique/smmR/). But to help with issues and contributions, a mirror has been setup at [Github](https://github.com/corentin-dev/smmR).

If you want to contribute to the code, you can fork the repository, make some changes and create a pull request to have them integrated into the package. You can use the `devtools::check()` function in order to verify that tests are still passing. See also the contributing guidelines.

If you encounter a problem, open a new issue. Try to be concise and explain what the problem is. If you have an example code that shows the error, it can be helpful.

## Acknowledgements

We acknowledge the project AStERiCs _Apprentissage Statistique à l'Echelle pour la Représentation et la Classification non-supervisées_ (RIN project funded by the Normandy Region), DAISI on Biomedical Data Classification (co-financed by the European Union with the European Regional Development Fund (ERDF) and by the Normandy Region).
