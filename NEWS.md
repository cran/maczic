# maczic 0.2.0 

* Added a `NEWS.md` file to track changes to the package.
* Added grant information to the DESCPITION file.
* Removed assign("last.warning", NULL, envir = baseenv()) from `maczic_power()` because it was flagged as an error in r-devel due to recent changes which lock the base environment, so that new bindings can no longer be added there.
* Corrected power (rejection) calculation in `maczic_power()` from p-value<0.05 to p-value<=0.05. 

# maczic 0.1.0 (11/04/2020)

* Initial release on CRAN.
