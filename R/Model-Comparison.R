###
# Model validation and comparison
# Zane Billings
# 2025-02-28
# Check the ELPD of each fitted model with LOO-CV and compare models
# Targetized 2025-04-16
###

suppressPackageStartupMessages({
	library(brms)
	library(rstan)
	library(here, include.only = NULL)
})

correct_waic <- function(waic, nan_values, model_metadata) {
	# Get the current waic thing
	x <- waic
	# Get the pointwise WAIC contribution matrix and remove the NaN values
	pointwise_finite <- x$pointwise[-nan_values, ]
	# Record how many points we are using
	N <- nrow(pointwise_finite)
	d <- attr(x, "dims")[[1]]
	# Get the overall estimates, which are columnwise sums of the pointwise
	# contributions for each statistic
	estimates <- colSums(pointwise_finite)
	# And now compute the standard error columnwise of each statistic
	SEs <- apply(pointwise_finite, 2, \(x) sqrt(N * var(x)))
	# Format the output
	out <- x
	out$estimates <- matrix(
		c(estimates, SEs),
		ncol = 2,
		dimnames = list(
			rownames(x$estimates),
			colnames(x$estimates)
		)
	)
	out$pointwise <- pointwise_finite
	out$elpd_waic <- estimates[[1]]
	out$p_waic <- estimates[[2]]
	out$waic <- estimates[[3]]
	out$se_elpd_waic <- SEs[[1]]
	out$se_p_waic <- SEs[[2]]
	out$se_waic <- SEs[[3]]
	attr(out, "dims") = c(d, N)
	# attr(out, "model_name") = paste0(
	# 	as.character(model_metadata[[i, 1]]), "_",
	# 	model_metadata[[i, 2]]
	# )
	return(out)
}

find_ll_nan_values <- function(waic_list) {
	nan_points <- lapply(
		waic_list,
		\(x) x$pointwise[, 1] |> is.na() |> which()
	)
	nan_vals <- unlist(nan_points) |> unique()
	return(nan_vals)
}


calculate_elpd_from_ll <- function(ll, nan_values) {
	idx_to_use <- setdiff(1:ncol(ll), nan_values)
	llmat_c <- ll[,idx_to_use]
	out <- loo::loo(llmat_c)
	return(out)
}

make_ic_comparisons <- function(loo_list, model_metadata, ic) {
	comparison_groups <-
		model_metadata |>
		dplyr::select(metric) |>
		dplyr::distinct()

	n <- nrow(comparison_groups)

	model_names <- paste(
		model_metadata$model_forms,
		model_metadata$metric,
		sep = "_"
	)

	elpd_comps <- lapply(
		seq_len(n),
		function(i) {
			g <- comparison_groups[i, ]
			equal_idxs <- purrr::map(
				1:ncol(g),
				\(x) {
					col <- names(g)[[x]]
					idx <- which(model_metadata[[col]] == as.character(g[[1, x]]))
					return(idx)
				}
			)
			all_idxs <- unique(do.call(c, equal_idxs))
			mod_names <- model_names[all_idxs]
			mod_elpd <- rlang::set_names(loo_list[all_idxs], mod_names)
			# Note that compare_ic() is apparently deprecated, but the replacement
			# function loo_compare() doesn't actually work on objects of class waic
			# TODO submit a brms issue about this
			out <- suppressWarnings(brms::compare_ic(x = mod_elpd, ic = ic))
			return(out)
		}
	)

	names(elpd_comps) <- comparison_groups$metric

	return(elpd_comps)
}

# Function to format the elpd comparison nicely
format_elpd <- function(elpd) {
	est_gam <- elpd$gamm$estimates
	est_lm <- elpd$lmm$estimates

	k_gam <- elpd$gamm$diagnostics$pareto_k
	k_lm <- elpd$lmm$diagnostics$pareto_k

	neff_gam <- elpd$gamm$diagnostics$n_eff
	neff_lm <- elpd$lmm$diagnostics$n_eff

	reff_gam <- elpd$gamm$diagnostics$r_eff
	reff_lm <- elpd$lmm$diagnostics$r_eff

	diff_est <- elpd$ic_diffs__

	# Make a nice-ish table for supplement
	supplement_diagnostics <-
		dplyr::bind_rows(
			"GAMM" = est_gam[3, ],
			"LMM" = est_lm[3, ],
			.id = "Model"
		) |>
		dplyr::mutate(
			"LOOIC" = paste0(
				formatC(Estimate, digits = 1, format = "f"), " ± ",
				formatC(SE, digits = 1, format = "f")
			),
			.keep = "unused"
		) |>
		dplyr::mutate(
			"Max. Pareto k" = c(
				max(k_gam),
				max(k_lm)
			) |>
				formatC(digits = 2, format = "f"),
			"Min. N_eff" = c(
				min(neff_gam),
				min(neff_lm)
			) |>
				formatC(digits = 1, format = "f"),
			"Max. R_eff" = c(
				max(reff_gam),
				max(reff_lm)
			) |>
				formatC(digits = 2, format = "f")
		)

	# Difference table row for main text
	diff_tbl <-
		diff_est |>
		as.data.frame() |>
		tibble::remove_rownames() |>
		dplyr::mutate(
			rat = formatC(abs(LOOIC / SE), digits = 1, format = "f"),
			LOOIC = formatC(LOOIC, digits = 2, format = "f"),
			SE = formatC(SE, digits = 2, format = "f"),
			est = paste0(LOOIC, " (±", SE, ")")
		)

	return(list(
		"supp" = supplement_diagnostics,
		"diff" = diff_tbl
	))
}

make_elpd_comparison_table <- function(formatted_elpd_comps) {
	loo_diff_tbl_data <-
		formatted_elpd_comps |>
		purrr::map(\(x) x$diff) |>
		dplyr::bind_rows(.id = "Metric") |>
		dplyr::transmute(
			"Metric" = Metric,
			"ΔELPD (LMM - GAMM)" = gsub("-", "", est),
			"ΔELPD / SE" = rat
		) |>
		dplyr::filter(
			Metric %in% c("Cartographic", "Grantham", "p-Epitope", "Temporal")
		)

	loo_tbl_diff <-
		loo_diff_tbl_data |>
		dplyr::arrange(`Metric`) |>
		flextable::flextable() |>
		flextable::autofit()

	return(loo_tbl_diff)
}

clean_supp_elpd_data <- function(formatted_elpd_comps) {
	tbl_data <-
		formatted_elpd_comps |>
		purrr::map(\(x) x$supp)

	tbl_data_comb <-
		tbl_data |>
		dplyr::bind_rows(.id = "metric") |>
		dplyr::filter(
			metric %in% c("Cartographic", "Grantham", "p-Epitope", "Temporal")
		)

	return(tbl_data_comb)
}

make_all_model_elpd_table <- function(supp_elpd_data) {
	out <-
		supp_elpd_data |>
		dplyr::select(metric, Model, LOOIC) |>
		tidyr::pivot_wider(
			names_from = Model,
			values_from = LOOIC
		) |>
		`names<-`(c("Metric", "GAMM LOO-IC", "LMM LOO-IC"))

	out_tbl <-
		out |>
		dplyr::select(Metric, `LMM LOO-IC`) |>
		dplyr::mutate(
			Metric = factor(
				Metric,
				levels = c("Cartographic", "p-Epitope", "Grantham", "Temporal")
			)
		) |>
		dplyr::arrange(Metric) |>
		flextable::flextable() |>
		flextable::autofit()

	return(out_tbl)
}

make_elpd_diagnostic_table <- function(supp_elpd_data) {
	out <-
		supp_elpd_data |>
		dplyr::select(-LOOIC)

	out_tbl <-
		out |>
		dplyr::mutate(
			Metric = factor(
				metric,
				levels = c("Cartographic", "p-Epitope", "Grantham", "Temporal")
			),
			.keep = "unused",
			.before = dplyr::everything()
		) |>
		dplyr::arrange(Metric) |>
		flextable::flextable() |>
		flextable::merge_v(j = 1) |>
		flextable::valign(j = 1, valign = "top") |>
		flextable::fix_border_issues() |>
		flextable::autofit()

	return(out_tbl)
}
