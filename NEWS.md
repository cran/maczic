# maczic 1.0.0 

* Deleted the extra 'The' in the last sentence from the DESCPITION file.
* Added an argument 'digits' to `plot_sensitivity()` and `maczic_power()`  to round the values to be returned.
* Added boot.ci.type="perc" to the return of `mediate_zi()` and `mediate_zi_vcoef()` so using `summary()`from the `mediation` package will print 'Nonparametric Bootstrap Confidence Intervals with the Percentile Method' in the summary if boot=TRUE.
* Modified `mediate_iv()` to get a bootstrap sample for an iteration until there are no warnings from the fitted mediator and the outcome model with the bootstrap sample.
* Added import(mathjaxr) to maczic.R to remove the note 'All declared Imports should be used' from 'R CMD check'.
* Fixed @docType package in maczic.R by replacing it with @aliases maczic-package

# maczic 0.2.0 (4/23/2021) 

* Added a `NEWS.md` file to track changes to the package.
* Added grant information to the DESCPITION file.
* Removed assign("last.warning", NULL, envir = baseenv()) from `maczic_power()` because it was flagged as an error in r-devel due to recent changes which lock the base environment, so that new bindings can no longer be added there.
* Corrected power (rejection) calculation in `maczic_power()` from p-value<0.05 to p-value<=0.05. 

# maczic 0.1.0 (11/04/2020)

* Initial release on CRAN.
