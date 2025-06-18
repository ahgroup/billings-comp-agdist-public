###
# Getting results from models
# Zane
# 2024-01-06
# Now that the models are fitted, we want to get the posterior predictions of
# interest and save them as a data set.
# Targetized 2025-04-16
###

suppressPackageStartupMessages({
	library(here)
	library(brms)
	library(rstan)
	library(ggplot2)
	library(patchwork)
})

generate_mem_grid <- function(model, dose_included = FALSE) {
	if (isTRUE(dose_included)) {
		grid <- marginaleffects::datagrid(
			model = model,
			dose = c("SD", "HD"),
			d_norm = seq(0, 1, 0.01),
			strain_type = c("H1N1", "H3N2", "B-Vic", "B-Yam")
		)
	} else if (isFALSE(dose_included)) {
		grid <- marginaleffects::datagrid(
			model = model,
			d_norm = seq(0, 1, 0.01),
			strain_type = c("H1N1", "H3N2", "B-Vic", "B-Yam")
		)
	} else {
		stop(paste(
			"'dose_included' argument should be TRUE or FALSE, not",
			dose_included
		))
	}

	return(tibble::tibble(grid))
}

summarize_prediction_draws <- function(prediction_epreds, keep_vars) {
	if (missing(keep_vars)) {
		keep_vars <- c("strain_type", "d_norm", ".epred")
	}

	out <-
		prediction_epreds |>
		dplyr::ungroup() |>
		# Using any_of() here will group by dose if it exists (we are on a dose
		# model) but will not use it if it doesn't exist
		dplyr::select(dplyr::any_of(keep_vars)) |>
		dplyr::group_by(dplyr::across(-.epred)) |>
		dplyr::summarize(
			tidybayes::mean_hdci(.epred),
			.groups = "drop"
		)

	return(out)
}

# This function gets the correct epred draws from the model.
get_predictions <- function(
		model, prediction_grid, N, summarize = TRUE, rf = ~strain_type, ...
	) {
	if (missing(N)) {N <- brms::ndraws(model)}
	pred_draws <-
		model |>
		tidybayes::epred_draws(
			newdata = prediction_grid,
			re_formula = rf,
			ndraws = N
		) |>
		dplyr::mutate(
			ti = .epred - log_pretiter
		)
	if (isTRUE(summarize)) {
		out <- summarize_prediction_draws(pred_draws, ...)
	} else {
		out <- pred_draws
	}

	return(out)
}

create_model_prediction_plot <- function(epreds, model_metadata, file_path) {
	model_data <-
		model_metadata |>
		# Just select the first model since they all have the same data
		dplyr::filter(
			model_forms == "gamm",
			metric %in% c("Temporal", "Cartographic", "Grantham", "p-Epitope")
		) |>
		dplyr::mutate(
			metric = factor(
				metric,
				levels = c("Cartographic", "Grantham", "p-Epitope", "Temporal")
			)
		) |>
		dplyr::select(metric, model_data) |>
		tidyr::unnest(model_data) |>
		dplyr::mutate(
			y = posttiter
		)

	model_preds <-
		model_metadata |>
		dplyr::select(metric, model_type = model_forms) |>
		tibble::add_column(epred_draws = epreds) |>
		dplyr::filter(
			metric %in% c("Temporal", "Cartographic", "Grantham", "p-Epitope")
		) |>
		dplyr::mutate(
			metric = factor(
				metric,
				levels = c("Cartographic", "Grantham", "p-Epitope", "Temporal")
			)
		) |>
		tidyr::unnest(epred_draws) |>
		dplyr::mutate(
			dplyr::across(c(y, ymin, ymax), hgp::hai_to_natural_scale),
			model_type = factor(
				model_type,
				levels = c("gamm", "lmm"),
				labels = c("GAMM", "LMM")
			)
		)

	main_model_fig <-
		ggplot2::ggplot() +
		ggplot2::aes(x = d_norm, y = y) +
		ggplot2::geom_point(
			data = model_data,
			position = ggplot2::position_jitter(0.01, 0.15, 340),
			size = 1, alpha = 0.01,
			show.legend = TRUE
		) +
		ggplot2::geom_ribbon(
			data = model_preds,
			mapping = ggplot2::aes(ymin = ymin, ymax = ymax, fill = model_type),
			color = "transparent", alpha = 0.25,
			show.legend = FALSE
		) +
		ggplot2::geom_line(
			data = model_preds,
			mapping = ggplot2::aes(color = model_type, linetype = model_type),
			lwd = 2.5, alpha = 0.75
		) +
		ggplot2::facet_grid(strain_type ~ metric) +
		ggplot2::scale_y_continuous(
			trans = "log2",
			breaks = hgp::hai_to_natural_scale(seq(0, 12, 2)),
			labels = hgp::hai_to_natural_scale(seq(0, 12, 2)),
			minor_breaks = hgp::hai_to_natural_scale(seq(0, 12, 1))
		) +
		ggplot2::scale_color_brewer(
			palette = "Dark2",
			aesthetics = c("color", "fill"),
			name = "Model"
		) +
		ggplot2::scale_linetype_manual(
			name = "Model",
			values = 1:2,
			guide = ggplot2::guide_legend(override.aes = list(lwd = 1.25))
		) +
		ggplot2::coord_cartesian(
			ylim = hgp::hai_to_natural_scale(c(0, 10))
		) +
		ggplot2::labs(
			x = "Normalized antigenic distance",
			y = "Post-vaccination titer"
		) +
		hgp::theme_ms() +
		ggplot2::theme(
			legend.key.width = ggplot2::unit(0.5, "in"),
			legend.key.height = ggplot2::unit(0.5, "in")
		)

	ggplot2::ggsave(
		filename = file_path,
		plot = main_model_fig,
		width = 6.5 * 1.5,
		height = 8 * 1.5
	)

	return(file_path)
}

calculate_contrast_preds <- function(
		model_metadata, model_type_to_use, fitted_models, prediction_grids
	) {
	# First get all of the unique combinations of metrics chosen pairwise
	metric_combs <-
		model_metadata |>
		dplyr::filter(
			metric %in% c("Temporal", "Cartographic", "Grantham", "p-Epitope")
		) |>
		dplyr::pull(metric) |>
		unique() |>
		combn(m = 2)

	n_combs <- ncol(metric_combs)

	# Get only the models of a certain model type(s)
	selected_models <-
		model_metadata |>
		dplyr::mutate(
			model = fitted_models,
			grid = prediction_grids
		) |>
		dplyr::filter(model_forms %in% model_type_to_use)

	# Get the full posterior predictions, not just the summaries
	# We have to do this to construct a valid CI for the contrasts
	# Could be refactored so that the predicting/summarizing step are done
	# separately and then we don't have to recompute this
	# But it's fast enough that it doesn't matter
	full_draws <- purrr::pmap(
		selected_models,
		\(model, grid, ...) model |>
			brms::posterior_epred(
				newdata = grid,
				re_formula = ~strain_type
			)
	)

	# Next get the contrast predictions
	contrast_preds <- lapply(
		seq_len(n_combs),
		function(i) {
			this_comb <- metric_combs[, i]
			m1_preds <- full_draws[[which(selected_models$metric == this_comb[[1]])]]
			m2_preds <- full_draws[[which(selected_models$metric == this_comb[[2]])]]
			ct_preds <- m1_preds - m2_preds
			y <- apply(ct_preds, 2, tidybayes::mean_hdci) |>
				dplyr::bind_rows()
			grid <- selected_models$grid[[1]]
			contr <- dplyr::bind_cols(grid, y) |>
				dplyr::select(d_norm, strain_type, y, ymin, ymax) |>
				dplyr::mutate(
					m1 = this_comb[[1]],
					m2 = this_comb[[2]]
				)
		}
	)

	# Combine and assign nice labels
	contrast_preds_comb <- dplyr::bind_rows(contrast_preds) |>
		dplyr::mutate(contrast_label = paste0(m1, " / ", m2))

	return(contrast_preds_comb)
}

make_contrast_preds_plot <- function(contrast_preds, file_path) {
	contrast_preds_plot <-
		ggplot2::ggplot(contrast_preds) +
		ggplot2::aes(x = d_norm, y = y) +
		ggplot2::geom_ribbon(
			ggplot2::aes(ymin = ymin, ymax = ymax, fill = strain_type),
			color = "transparent", alpha = 0.2,
			show.legend = FALSE
		) +
		ggplot2::geom_line(
			ggplot2::aes(color = strain_type, linetype = strain_type),
			lwd = 1.5, alpha = 1
		) +
		ggplot2::geom_hline(
			yintercept = c(-2, 0, 2),
			linetype = 1,
			color = "black",
			alpha = 0.5
		) +
		ggplot2::facet_wrap(~contrast_label) +
		ggplot2::scale_y_continuous(
			breaks = seq(-6, 6, 2),
			minor_breaks = seq(-6, 6, 1),
			labels = MASS::fractions(2 ^ seq(-6, 6, 2))
		) +
		ggplot2::coord_cartesian(
			ylim = c(-5, 5)
		) +
		ggplot2::scale_linetype_manual(
			name = "Strain type",
			values = 2:5
		) +
		ggplot2::scale_color_viridis_d(
			aesthetics = c("color", "fill"),
			name = "Strain type"
		) +
		ggplot2::labs(
			x = "Normalized antigenic distance",
			y = "Fold change in post-vaccination titers between metrics"
		) +
		hgp::theme_ms() +
		ggplot2::theme(
			axis.text = ggplot2::element_text(size = 16),
			axis.title = ggplot2::element_text(size = 24),
			strip.text = ggplot2::element_text(size = 24),
			legend.key.width = ggplot2::unit(0.5, "in"),
			legend.key.height = ggplot2::unit(0.5, "in"),
			legend.text = ggplot2::element_text(size = 18)
		)

	ggplot2::ggsave(
		file_path,
		plot = contrast_preds_plot,
		width = 13,
		height = 10
	)

	return(file_path)
}

make_lmm_parameters_plot <- function(
		model_metadata, fitted_models, file_path
	) {
	require(patchwork, quietly = TRUE)

	# Add the models to the model metadata and filter for only the LMMs
	# Probably should refactor to have an "augmented model metadata" target
	# to avoid repeating this step. but it is fast and that takes time
	lmms <- model_metadata |>
		tibble::add_column(model = fitted_models) |>
		dplyr::filter(model_forms == "lmm")

	# Now get the parameter samples
	lm_parm_samples <- purrr::map(
		lmms$model,
		\(x) brms::as_draws_df(
			x,
			variable = c(
				"b_Intercept", "^b_.*d_norm$",
				"^r_strain_type\\["
			),
			regex = TRUE
		) |>
			tibble::as_tibble() |>
			# Pivot everything longer except the control variables
			tidyr::pivot_longer(
				cols = !c(.chain, .draw, .iteration, b_Intercept, b_d_norm)
			) |>
			# Separate out the main components of the name
			tidyr::separate_wider_regex(
				name,
				patterns = c(
					".*\\[",
					strain_type = "[^,]+",
					",",
					parameter = "[^\\]]+",
					"\\]"
				),
				too_few = "align_start"
			) |>
			# Now reconstruct the strain-type specific slopes and intercepts
			dplyr::mutate(
				is_intercept = as.numeric(parameter == "Intercept"),
				st_value = value + (is_intercept * b_Intercept) +
					((1 - is_intercept) * b_d_norm)
			) |>
			dplyr::select(
				.chain, .iteration, .draw, strain_type, parameter, value = st_value
			) |>
			# Cleanup
			dplyr::mutate(
				parameter = factor(
					parameter,
					levels = c("Intercept", "d_norm"),
					labels = c("Intercept", "Slope")
				),
				strain_type = clean_subtype(strain_type)
			)
	)

	# Combine the list of data frames and
	# Get the posterior summaries
	lm_parm_samples_combined <-
		lm_parm_samples |>
		rlang::set_names(lmms$metric) |>
		dplyr::bind_rows(.id = "metric") |>
		dplyr::filter(
			metric %in% c("Temporal", "Cartographic", "Grantham", "p-Epitope")
		) |>
		dplyr::mutate(
			metric = factor(
				metric,
				levels = c("Cartographic", "Grantham", "p-Epitope", "Temporal")
			) |>
				forcats::fct_rev()
		) |>
		dplyr::summarise(
			tidybayes::mean_hdci(value),
			.by = c(metric, parameter, strain_type)
		)

	# Finally make the plot
	# Easiest to do it in two steps and patchwork together
	lmm_intercept_plot <-
		lm_parm_samples_combined |>
		dplyr::filter(parameter == "Intercept") |>
		ggplot2::ggplot() +
		ggplot2::aes(
			x = y, xmin = ymin, xmax = ymax,
			y = metric
		) +
		ggplot2::geom_errorbar(width = 0.25, linewidth = 1.5) +
		ggplot2::geom_point(size = 5) +
		ggplot2::facet_grid(~ strain_type) +
		ggplot2::labs(
			x = "LMM Intercept",
			y = NULL
		) +
		hgp::theme_ms()

	lmm_slope_plot <-
		lm_parm_samples_combined |>
		dplyr::filter(parameter == "Slope") |>
		ggplot2::ggplot() +
		ggplot2::aes(
			x = y, xmin = ymin, xmax = ymax,
			y = metric
		) +
		ggplot2::geom_errorbar(width = 0.25, linewidth = 1.5) +
		ggplot2::geom_point(size = 5) +
		ggplot2::facet_grid(~ strain_type) +
		ggplot2::labs(
			x = "LMM Slope",
			y = NULL
		) +
		hgp::theme_ms()

	lmm_parm_plot <- lmm_intercept_plot / lmm_slope_plot &
		ggplot2::theme(
			axis.text = ggplot2::element_text(size = 16),
			axis.title = ggplot2::element_text(size = 24),
			strip.text = ggplot2::element_text(size = 24),
			legend.key.width = ggplot2::unit(0.5, "in"),
			legend.key.height = ggplot2::unit(0.5, "in"),
			legend.text = ggplot2::element_text(size = 18)
		)

	ggplot2::ggsave(
		file_path,
		plot = lmm_parm_plot,
		width = 13,
		height = 7
	)

	return(file_path)
}

# Andreas wants to do a sensitivity analysis where we look at each vaccine
# component separately instead of just the subtypes.
create_vaccine_specific_prediction_plot <- function(
		model_metadata, fitted_models, file_path
	) {
	aug_models <-
		model_metadata |>
		tibble::add_column(model = fitted_models)

	vac_component_preds <-
		aug_models |>
		dplyr::mutate(
			unique_vacs_per_strain_type = purrr::map(
				model_data,
				\(d) tidyr::expand(d, tidyr::nesting(strain_type, vaccine_name))
			),
			grid = purrr::map2(
				model, unique_vacs_per_strain_type,
				\(m, d) marginaleffects::datagrid(
					model = m,
					d_norm = seq(0, 1, 0.01),
					study = "UGA"
				) |>
					dplyr::select(-vaccine_name, -strain_type) |>
					tidyr::expand_grid(d)
			)
		)

	vac_component_preds$epred_draws <- purrr::pmap(
		vac_component_preds,
		\(model, grid, ...) get_predictions(
			model, grid,
			rf = ~strain_type + vaccine_name,
			keep_vars = c("strain_type", "vaccine_name", "d_norm", ".epred")
		)
	)

	# This part is copied from create_model_prediction_plot() since it's too
	# much work to edit that function
	model_data <-
		model_metadata |>
		# Just select the first model since they all have the same data
		dplyr::filter(
			model_forms == "gamm",
			metric %in% c("Temporal", "Cartographic", "Grantham", "p-Epitope")
		) |>
		dplyr::select(metric, model_data) |>
		tidyr::unnest(model_data) |>
		dplyr::mutate(
			y = posttiter,
			strain_type = clean_subtype(strain_type),
			vaclab = paste0(strain_type, "\n", vaccine_name)
		)

	model_preds <-
		vac_component_preds |>
		dplyr::filter(
			metric %in% c("Temporal", "Cartographic", "Grantham", "p-Epitope")
		) |>
		tidyr::unnest(epred_draws) |>
		dplyr::select(
			model_type = model_forms,
			metric, strain_type, vaccine_name, d_norm,
			y, ymin, ymax
		) |>
		dplyr::mutate(
			dplyr::across(c(y, ymin, ymax), hgp::hai_to_natural_scale),
			model_type = factor(
				model_type,
				levels = c("gamm", "lmm"),
				labels = c("GAMM", "LMM")
			),
			strain_type = clean_subtype(strain_type),
			vaclab = paste0(strain_type, "\n", vaccine_name)
		)

	main_model_fig <-
		ggplot2::ggplot() +
		ggplot2::aes(x = d_norm, y = y) +
		ggplot2::geom_point(
			data = model_data,
			position = ggplot2::position_jitter(0.01, 0.15, 340),
			size = 0.5, alpha = 0.01,
			show.legend = TRUE
		) +
		ggplot2::geom_ribbon(
			data = model_preds,
			mapping = ggplot2::aes(ymin = ymin, ymax = ymax, fill = model_type),
			color = "transparent", alpha = 0.25,
			show.legend = FALSE
		) +
		ggplot2::geom_line(
			data = model_preds,
			mapping = ggplot2::aes(color = model_type, linetype = model_type),
			lwd = 1.5, alpha = 0.75
		) +
		ggplot2::facet_grid(vaclab ~ metric) +
		ggplot2::scale_y_continuous(
			trans = "log2",
			breaks = hgp::hai_to_natural_scale(seq(0, 12, 2)),
			labels = hgp::hai_to_natural_scale(seq(0, 12, 2)),
			minor_breaks = hgp::hai_to_natural_scale(seq(0, 12, 1))
		) +
		ggplot2::scale_color_brewer(
			palette = "Dark2",
			aesthetics = c("color", "fill"),
			name = "Model"
		) +
		ggplot2::scale_linetype_manual(
			name = "Model",
			values = 1:2,
			guide = ggplot2::guide_legend(override.aes = list(lwd = 1.25))
		) +
		ggplot2::coord_cartesian(
			ylim = hgp::hai_to_natural_scale(c(0, 10))
		) +
		ggplot2::labs(
			x = "Normalized antigenic distance",
			y = "Post-vaccination titer"
		) +
		hgp::theme_ms() +
		ggplot2::theme(
			legend.key.width = ggplot2::unit(0.5, "in"),
			legend.key.height = ggplot2::unit(0.5, "in"),
			strip.text.y = ggplot2::element_text(size = 12)
		)

	ggplot2::ggsave(
		filename = file_path,
		plot = main_model_fig,
		width = 6.5 * 1.5,
		height = 7 * 1.5
	)

	return(file_path)
}

# Function for summarizing the fixed effects coefficients from the lienar model
summarize_fixed_effects <- function(model_metadata, model_list) {
	aug_models <-
		model_metadata |>
		tibble::add_column(model = model_list) |>
		dplyr::filter(model_forms == "lmm")

	all_params <- brms::variables(aug_models$model[[1]])
	params_of_interest <- all_params[startsWith(all_params, "b_")]
	params_of_interest <- params_of_interest[params_of_interest != "b_Intercept"]

	model_summaries <-
		purrr::map(
			aug_models$model,
			\(x) brms::posterior_summary(x, variable = params_of_interest) |>
				as.data.frame() |>
				tibble::rownames_to_column(var = "param") |>
				tibble::as_tibble()
		) |>
		rlang::set_names(aug_models$metric) |>
		dplyr::bind_rows(.id = "metric")

	return(model_summaries)
}

create_fixed_effects_summary_table <- function(fixed_effects_summaries) {
	fes_formatted <- fixed_effects_summaries |>
		dplyr::mutate(
			dplyr::across(
				c(Estimate, Q2.5, Q97.5),
				\(x) formatC(x, digits = 2, format = "f", flag = "0- ")
			),
			metric = factor(
				metric,
				levels = c("Cartographic", "Grantham", "p-Epitope", "Temporal")
			)
		) |>
		tidyr::drop_na(metric) |>
		dplyr::transmute(
			Metric = metric,
			Parameter = dplyr::case_match(
				param,
				"b_birth_year_c" ~ "Birth Year",
				"b_age_c" ~ "Age",
				"b_sex_i" ~ "Sex",
				"b_race_i" ~ "Race/Ethnicity",
				"b_log_pretiter" ~ "Log pre-vaccination HAI titer",
				"b_d_norm" ~ "Normalized antigenic distance"
			),
			Estimate = paste0(Estimate, " (", Q2.5, ",", Q97.5, ")")
		) |>
		dplyr::arrange(Metric)

	fes_wide <- tidyr::pivot_wider(
		fes_formatted,
		names_from = Parameter,
		values_from = Estimate
	)

	fes_tbl <-
		fes_wide |>
		flextable::flextable() |>
		flextable::footnote(
			i = 1, j = 4:5,
			value = flextable::as_paragraph(c(
				"Reference: Male (vs. female)",
				"Reference: Non-Hispanic white (vs. other)"
			)),
			ref_symbols = c("1", "2"),
			part = "header"
		) |>
		flextable::align(j = 2:7, align = "right") |>
		flextable::align(j = 2:7, align = "center", part = "header") |>
		flextable::autofit()

	return(fes_tbl)
}

# Function for variance decomposition for the linear model
decompose_model_variance <- function(this_model, this_model_form, ndraws) {
	# First get the list of variable names and pull out the RE variances
	all_params <- brms::variables(this_model)
	re_vars <- all_params[startsWith(all_params, "sd_")]

	# The vector of all variance component names has all of the RE variances
	# plus the residual variance of the Gaussian outcome
	var_params <- c(re_vars, "sigma")

	# First, to get the variance explained by the fixed effects we need to get
	# the variance of η = X * beta, i.e. the linear predictor. We get this as the
	# linpreds (same as epreds since the model is Gaussian) of the inputted
	# dataset. This part gets the ndraws x nobs matrix of samples of each
	# value of the linear predictor
	eta_matrix <- brms::posterior_linpred(this_model)

	# Get the posterior samples of the linear predictor variance by taking the
	# variance of each sample of the linear predictor
	var_eta <- apply(eta_matrix, 1, var)

	# Now get the posterior samples of all of the random effects variances
	# and also the residual variance
	var_re <- this_model |>
		brms::as_draws_df(variable = var_params) |>
		tibble::as_tibble()

	# Now get the total re variance
	total_re_var <- var_re |>
		dplyr::select(tidyselect::starts_with("sd_")) |>
		apply(MARGIN = 1, FUN = sum)

	# Make everything into a nice data frame
	var_comps <- tibble::add_column(
		var_re,
		"var_eta" = var_eta,
		"total_re_var" = total_re_var
	)

	var_comps_names <- c("var_eta", var_params, "total_re_var")

	# First get the denominator for all of the variance ratios
	total_variance <- var_comps |>
		dplyr::select(dplyr::all_of(var_comps_names), -total_re_var) |>
		rowSums()

	# Next construct the posterior samples for all of the individual variance
	# ratios by iterating over var_comps_names -- we have to loop over the list
	# of dataframes each time though.
	var_ratios <- purrr::map(
		1:length(var_comps_names),
		\(i) {
			this_vc <- var_comps_names[[i]]
			vc_list <- var_comps[[this_vc]] / total_variance
		}
	) |>
		rlang::set_names(paste0("VR_", var_comps_names))

	# Deal with the subtype variance, which has two RE components in the LMM
	# but only one in the GAMM
	# TODO modify this to automatically match any REs for the same grouping
	# structure
	if (this_model_form == "lmm") {
		# Manually add together the two variance components for subtype
		var_ratios$subtype <- var_ratios$VR_sd_strain_type__Intercept +
			var_ratios$VR_sd_strain_type__d_norm
		var_ratios$VR_sd_strain_type__Intercept <- NULL
		var_ratios$VR_sd_strain_type__d_norm <- NULL
	} else if (this_model_form == "gamm") {
		var_ratios$subtype <- var_ratios$VR_sd_strain_type__Intercept
		var_ratios$VR_sd_strain_type__Intercept <- NULL
	}

	var_ratios_summary <-
		var_ratios |>
		purrr::map(tidybayes::mean_hdci) |>
		dplyr::bind_rows(.id = "vc")

	return(var_ratios_summary)
}

bind_variance_decompositions <- function(model_metadata, variance_decompositions_list) {
	out <-
		model_metadata |>
		dplyr::select(-model_data, -fitting_seed) |>
		tibble::add_column(vd = variance_decompositions_list) |>
		tidyr::unnest(vd)

	return(out)
}

create_variance_ratio_table <- function(var_ratios_df) {
	browser()
	var_ratios_fmt <-
		var_ratios_df |>
		dplyr::filter(model_forms == "lmm") |>
		dplyr::mutate(
			vc = factor(
				vc,
				levels = c(
					"VR_sigma",
					"VR_var_eta",
					"VR_total_re_var",
					"subtype",
					"VR_sd_strain_type:strain_name__Intercept",
					"VR_sd_strain_type:vaccine_name__Intercept",
					"VR_sd_study__Intercept",
					"VR_sd_subject_id__Intercept"
				),
				labels = c(
					"Residual variance",
					"Fixed effects",
					"Total random effects",
					"Specific random effects_Subtype",
					"Specific random effects_Assay strain",
					"Specific random effects_Vaccine strain",
					"Specific random effects_Study site",
					"Specific random effects_Subject"
				)
			),
			dplyr::across(
				c(y, ymin, ymax),
				\(x) formatC(x * 100, digits = 0, format = "d")
			)
		) |>
		dplyr::transmute(
			Metric = metric,
			"Variance component" = vc,
			"Est. Contribution" = paste0(y, "% (", ymin, ", ", ymax, ")")
		) |>
		dplyr::arrange(Metric, `Variance component`) |>
		dplyr::filter(
			Metric %in% c("Cartographic", "Grantham", "p-Epitope", "Temporal")
		)

	var_ratio_wide <-
		var_ratios_fmt |>
		dplyr::mutate(
			.id = dplyr::row_number(),
			.by = c(Metric, `Variance component`)
		) |>
		tidyr::pivot_wider(
			names_from = `Variance component`,
			values_from = `Est. Contribution`
		) |>
		dplyr::select(-.id)

	var_ratios_tbl <-
		var_ratio_wide |>
		flextable::flextable() |>
		flextable::separate_header() |>
		flextable::merge_v(j = 1) |>
		flextable::valign(j = 1, valign = "top") |>
		flextable::vline(j = c(1, 4)) |>
		flextable::align(j = 2:8, align = "right") |>
		flextable::fontsize(size = 8, part = "all") |>
		flextable::fix_border_issues() |>
		fit_flextable_to_page()

	return(var_ratios_tbl)
}
