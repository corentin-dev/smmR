# Function to check if an object is of class smm
is.smm <- function(x) {
  inherits(x, "smm")
}

# Method to get the sojourn time distribution f
.get.f <- function(x, ...) {
  UseMethod(".get.f", x)
}

# Method to get the semi-Markov kernel q
.get.q <- function(x, ...) {
  UseMethod(".get.q", x)
}

# Method to get the number of parameters
# (useful for the computation of criteria such as AIC and BIC)
.getKpar <- function(x) {
  UseMethod(".getKpar", x)
}
