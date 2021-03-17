# smmR 1.0.2

## Bug fixes

* Remove a bad argument check of the transition matrix `ptrans` in the S3 
class `mm`: when the order of the Markov model is `k = 1`, we checked if all 
diagonal elements of the transition matrix `ptrans` was 0 (transition to the 
same states were impossible).

## Minor Improvements

* Add more examples in the documentation.

# smmR 1.0.1

## Minor Improvements

* Added a `NEWS.md` file to track changes to the package.

* Remove `Depends` field in DESCRIPTION file.

## Bug fixes

* For initial distribution, transition matrices and conditional sojourn time 
distributions, we allow the fact that the sum of the probabilities may not be 
equal to 1 exactly due to round-off errors (we introduce `.Machine$double.eps`).