source("renv/activate.R")
#source("~/.Rprofile")

options(
	repos = c(
		PMPC = "https://packagemanager.posit.co/cran/latest",
		CRAN = "https://cran.r-project.org"
	)
)

options(
	warnPartialMatchArgs = FALSE,
	warnPartialMatchDollar = FALSE,
	warnPartialMatchAttr = FALSE
)

tv <- function(...) {
	targets::tar_visnetwork(targets_only = TRUE, ...)
}

td <- function() {
	targets::tar_make(use_crew = FALSE, callr_function = NULL, as_job = FALSE)
}

tm <- function(...) {
	targets::tar_make(...)
}

tc <- function(...) {
	targets::tar_make(use_crew = FALSE, ...)
}
