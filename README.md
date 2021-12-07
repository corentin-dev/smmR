# Semi Markov Model

[![PLMLab build status](https://plmlab.math.cnrs.fr/lmrs/statistique/smmR/badges/master/pipeline.svg)](https://plmlab.math.cnrs.fr/lmrs/statistique/smmR/-/pipelines) [![CRAN versions](https://www.r-pkg.org/badges/version/smmR)](https://CRAN.R-project.org/package=smmR) [![CRAN logs](https://cranlogs.r-pkg.org/badges/smmR)](https://CRAN.R-project.org/package=smmR) [![RDocumentation](https://api.rdocumentation.org/badges/version/smmR)](https://www.rdocumentation.org/packages/smmR)

## Contribute

The official repository is at [PLMLab](https://plmlab.math.cnrs.fr/lmrs/statistique/smmR/). But to help with issues and contributions, a mirror has been setup at [Github](https://github.com/corentin-dev/smmR).

## Install

```R
install.packages('smmR')
```

## Get started

First, load the needed libraries:
```r
library(smmR)
library(DiscreteWeibull)
```

Then, let us create a **smmparametric** object to represent the semi-Markov chain associated to the system:
```r
states <- c("1", "2", "3") # State space

alpha <- c(1, 0, 0) # Initial distribution

p <- matrix(data = c(0, 1, 0, 
                     0.95, 0, 0.05, 
                     1, 0, 0), nrow = 3, byrow = TRUE) # Transition matrix

distr <- matrix(c(NA, "geom", NA, 
                  "dweibull", NA, "dweibull", 
                  "dweibull", NA, NA), 
                nrow = 3, ncol = 3, byrow = TRUE) # Distribution matrix

param1 <- matrix(c(NA, 0.8, NA, 
                   0.3, NA, 0.5,
                   0.6, NA, NA), 
                 nrow = 3, ncol = 3, byrow = TRUE)

param2 <- matrix(c(NA, NA, NA, 
                   0.5, NA, 0.7,
                   0.9, NA, NA), 
                 nrow = 3, ncol = 3, byrow = TRUE)

parameters <- array(c(param1, param2), c(3, 3, 2))

factory <- smmparametric(states = states, init = alpha, ptrans = p, 
                         type.sojourn = "fij", distr = distr, param = parameters)
```

After that, we are able to simulate a sequence of sample sizes $`M = 10,000`$:

```r
M <- 10000
seq <- simulate(object = factory, nsim = M)
```

Thanks to the `smmR` package, we can estimate any semi-Markov model
with one or several discrete sequences. In our case, we are going to
introduce a **non-parametric estimation**:

```r
estimate <- fitsmm(sequences = seq, states = states, type.sojourn = "fij")
```

The estimate $`\hat{p}`$ of the transition matrix $`p`$ is:

```r
print(x = estimate$ptrans, digits = 2)
```

See the vignette `textile-factory` for more details.

## Build manually

### Get sources

You need to have git. Then:
```bash
git clone https://plmlab.math.cnrs.fr/lmrs/statistique/smmR.git
```

### Build

```bash
R CMD check --as-cran .
```

### Check
```bash
R CMD check --as-cran .
```
