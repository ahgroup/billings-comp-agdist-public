###
# Fitting brms models to data
# Zane Billings
# 2024-01-29
# This script fits the finalized brms models to the data. There are four
# models, each with the same structure. The only difference is the outcome and
# the outcome distribution.
# titerincrease and prevactiter both use Gaussian likelihoods
# seroconversion and seroprotection both use Binomial likelihoods
# In contrast to the previous modeling approach, all strains are fitted within
# the same model, and we allow the random effects to handle differences across
# the strains.
# Models rerun on: 2024-10-14
###

# Setup ========================================================================
## Declare dependencies ####
suppressPackageStartupMessages({
	library(hgp)
	library(cmdstanr)
	library(brms)
	library(rstan)
})

## Load settings ####
# This script sets up variables for brms settings that are reused, so they can
# all be adjusted in the same place.
#source(here::here("R", "common-functions", "brms-settings.R"))
source(here::here("R", "Utils.R"))

get_brms_priors <- function(model = c("lmm", "gamm")) {
	model <- match.arg(model)

	if (model == "gamm") {
		out <- c(
			brms::prior(normal(0,1), class = "Intercept"),
			brms::prior(normal(0,1), class = "b"),
			brms::prior(student_t(3, 0, 1), class = "sd", lb = 0),
			brms::prior(student_t(3, 0, 1), class = "sigma", lb = 0),
			brms::prior(student_t(3, 0, 0.25), class = "sds", lb = 0)
		)
	} else if (model == "lmm") {
		out <- c(
			brms::prior(normal(0,1), class = "Intercept"),
			brms::prior(normal(0,1), class = "b"),
			brms::prior(student_t(3, 0, 1), class = "sd", lb = 0),
			brms::prior(student_t(3, 0, 1), class = "sigma", lb = 0)
		)
	} else {
		stop(paste0(
			"Model should be one of 'lmm' or 'gamm', not '", model, "'."
		))
	}

	return(out)
}

get_brms_formula <- function(model = c("lmm", "gamm")) {
	model <- match.arg(model)

	if (model == "gamm") {
		out <- brms::brmsformula(
			y | cens(c, y2) ~ 1 +
				birth_year_c + age_c + sex_i + race_i +
				log_pretiter + s(d_norm, k = 5, by = strain_type) +
				(1 | strain_type) +
				(1 | study) + (1 | subject_id) +
				(1 | strain_type:vaccine_name) + (1 | strain_type:strain_name)
		)
	} else if (model == "lmm") {
		out <- brms::brmsformula(
			y | cens(c, y2) ~ 1 +
				birth_year_c + age_c + sex_i + race_i +
				log_pretiter + d_norm + (1 + d_norm | strain_type) +
				(1 | study) + (1 | subject_id) +
				(1 | strain_type:vaccine_name) + (1 | strain_type:strain_name)
		)
	} else {
		stop(paste0(
			"Model should be one of 'lmm' or 'gamm', not '", model, "'."
		))
	}

	return(out)
}

nest_model_data_for_fitting <- function(model_data, additional_nesting_terms) {
	if(missing(additional_nesting_terms)) {
		additional_nesting_terms <- character()
	}

	nested_model_data <-
		model_data |>
		# Add the censoring variables for brms
		hgp::format_hai_data(post_titer = "log_posttiter") |>
		tidyr::nest(
			model_data = -c(metric, tidyselect::all_of(additional_nesting_terms))
		)
	return(nested_model_data)
}

define_model_metadata <- function(
		nested_model_data, seed_list,
		run_gamm = TRUE, run_lmm = TRUE
) {
	# prior_args <- c("no", "only")[c(run_posterior, run_prior)]
	model_forms <- c("gamm", "lmm")[c(run_gamm, run_lmm)]
	#models_to_run <- c("overall", "simple")[c(run_overall, run_simple)]

	out <-
		tidyr::expand_grid(
			model_forms,
			nested_model_data
		)

	if(missing(seed_list)) {seed_list <- sample(1:100000, nrow(out), TRUE)}

	out$fitting_seed <- seed_list[1:nrow(out)]

	return(out)
}

generate_brms_model_filenames <- function(nested_model_data, prefix = NULL) {
	n <- nested_model_data |>
		nrow()

	i <-
		seq_len(n) |>
		stringr::str_pad(width = nchar(n), side = "left", pad = "0")

	f <- here::here(
		"results", "large-files", "brms-fits", prefix,
		paste0("brms-model-", i, ".Rds")
	)

	return(f)
}

get_brms_settings <- function(testing = FALSE, verbose = TRUE) {
	if (isTRUE(testing)) {
		brms_settings1 <- list(
			iter = 20,
			warmup = 10,
			chains = 4,
			cores = 4,
			control = list(adapt_delta = 0.8, max_treedepth = 10)
		)
	} else {
		brms_settings1 <- list(
			iter = 750,
			warmup = 350,
			chains = 48,
			cores = 48,
			control = list(adapt_delta = 0.99, max_treedepth = 15)
		)
	}

	if (isTRUE(verbose)) {
		brms_settings2 <- list(
			silent = 2,
			refresh = 100
		)
	} else {
		brms_settings2 <- list(
			silent = 0,
			refresh = 0
		)
	}

	brms_settings <- c(
		brms_settings1,
		brms_settings2
	)

	return(brms_settings)
}

fit_brms_model <- function(
		data, model, sample_priors, seed, test = FALSE
	) {
	model_formula <- get_brms_formula(model)
	model_priors <- get_brms_priors(model)
	brms_settings <- get_brms_settings(test, FALSE)

	fitted_model <- brms::brm(
		formula = model_formula,
		data = data,
		prior = model_priors,
		family = gaussian(),
		backend = "cmdstanr",
		algorithm = "sampling",
		seed = seed,
		save_pars = brms::save_pars(all = TRUE),
		sample_prior = sample_priors,
		silent = brms_settings$silent,
		refresh = brms_settings$refresh,
		iter = brms_settings$iter,
		warmup = brms_settings$warmup,
		chains = brms_settings$chains,
		cores = brms_settings$cores,
		control = brms_settings$control
	)

	return(fitted_model)
}



# # END OF FILE ==================================================================
