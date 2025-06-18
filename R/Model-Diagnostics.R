###
# Model Diagnostics
# Zane
# 2024-01-06
# Check (and report) common diagnostics for the models to ensure the fit is OK
# to make conclusions about. Diagnostics we will include are:
# 1. R-hat (target < 1.01 for all parameters)
# 2. ESS (target bulk and tail > 100*n_chains for all parameters)
# 3. prior vs. posterior shrinkage (no real target, just need to examine)
# 4. Number of divergent transitions (target <1%)
# 5. Number of treedepth exceedences (target <1%)
# 6. E-BFMI (target > 0.3 for all chains)
# Suggested links for Andreas to read to understand what this code does:
# https://www.jstatsoft.org/article/view/v080i01
# https://journal.r-project.org/archive/2018/RJ-2018-017/index.html
# https://mc-stan.org/learn-stan/diagnostics-warnings.html
# https://mc-stan.org/docs/reference-manual/execution.html#model-diagnostics
# https://adv-r.hadley.nz/fp.html
# https://purrr.tidyverse.org/articles/base.html
# https://www.jumpingrivers.com/blog/new-features-r410-pipe-anonymous-functions/
# (section "shorthand syntax for anonymous functions")
# https://coolbutuseless.github.io/2019/03/13/anonymous-functions-in-r-part-1/
# I WOULD type all of this out in the comments but I don't have enough time
# Targetized 2025-04-15
###

suppressPackageStartupMessages({
	library(rstan)
	library(brms)
	library(posterior, include.only = NULL)
	library(future, include.only = NULL)
	library(future.apply, include.only = NULL)
})

get_model_diagnostics <- function(fitted_brms_model) {
	model_summary <- posterior::summarize_draws(fitted_brms_model)
	out <- tibble::tibble(
		"Pct. Divergences" = formatC(
			rstan:::get_num_divergent(fitted_brms_model$fit) /
				brms::ndraws(fitted_brms_model) * 100,
			digits = 1, format = "f"
		) |> paste0("%"),
		"Pct. TD Exceeded" = formatC(
			rstan:::get_num_max_treedepth(fitted_brms_model$fit) /
				brms::ndraws(fitted_brms_model) * 100,
			digits = 1, format = "f"
		) |> paste0("%"),
		"min E-BFMI" = rstan::get_bfmi(fitted_brms_model$fit) |> min(),
		"min ESS (tail)" = min(model_summary$ess_tail, na.rm = TRUE),
		"min ESS (bulk)" = min(model_summary$ess_bulk, na.rm = TRUE),
		"max R_hat" = max(model_summary$rhat, na.rm = TRUE)
	)
	return(out)
}

make_diagnostics_table <- function(cleaned_diagnostics_list, model_metadata) {
	cleaned_diagnostics <-
		cleaned_diagnostics_list |>
		dplyr::bind_rows() |>
		tibble::add_column(
			model = model_metadata$model_forms,
			metric = model_metadata$metric
		) |>
		dplyr::transmute(
			Model = factor(
				model,
				levels = c("gamm", "lmm"),
				labels = c("GAMM", "LMM")
			),
			`Pct. Divergences`,
			dplyr::across(tidyselect::contains("ESS"), \(x) sprintf("%.0f", x)),
			dplyr::across(c(`min E-BFMI`, `max R_hat`), \(x) sprintf("%.2f", x))
		)

	out <-
		cleaned_diagnostics |>
		flextable::flextable() |>
		#flextable::align(j = 3, align = "right") |>
		flextable::merge_v(j = 1) |>
		flextable::valign(j = 1, valign = "top") |>
		flextable::fix_border_issues() |>
		flextable::fontsize(size = 10, part = "all")

	return(out)
}
