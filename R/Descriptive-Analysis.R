###
# Data summary tables and/or figures
# Zane Billings
# 2024-01-29
# Targetized 2025-04-02
# This script makes Table 1 for the paper, in addition to any plots and tables
# that summarize the data but not the model fits.
###

# ---- Setup ----
# Declare package dependencies
suppressPackageStartupMessages({
	library("ggplot2", include.only = NULL)
	library("Hmisc", include.only = NULL)
	library(patchwork)
	library(brms)
})

# Data description tables ####
## Counts by studies table ####
# Table should show the number of HAI samples (paired), individuals, and
# person-years contributed by each study. We'll use a function to make a
# sub-table for each of these.
create_count_subtable <- function(data, label) {
	out <- gtsummary::tbl_cross(
		data,
		study, season_short,
		label = list(study = label, season_short = "Season")
	) |>
		gtsummary::modify_column_indent(
			columns = label, rows = c(F, T, T, T, T)#, indent = 4L
		)
	return(out)
}

# First we need to create the sub-table showing the number of HAI counts.
# First we need to do a bit of filtering to make sure we don't double count
# any b-pre assays, or count anything multiple times because of the ag-dist
# multiple metrics.
get_distinct_measurements <- function(dat) {
	distinct_measurements <-
		dat |>
		dplyr::mutate(
			strain_type = hgp::replace_strain_names(
				strain_name,
				from = "short",
				to = "subtype"
			)
		) |>
		dplyr::distinct(
			study, subject_id, season_short, strain_type, strain_name,
			vaccine_name, pretiter, posttiter, fold_change
		)

	return(distinct_measurements)
}

# Next we'll count the number of person years contributed by each study in
# each season -- this treats each repeat by a known subject as a new
# individual.
get_distinct_personyears <- function(dat) {
	distinct_personyears <-
		dat |>
		dplyr::distinct(study, season_short, subject_id)

	return(distinct_personyears)
}

# Next we'll count the number of new unique participants who arrived at each
# study site in each season.
get_distinct_new_participants <- function(dat) {
	new_participants <-
		dat |>
		# Since the season variable is an ordered factor, we find the season an
		# individual first occurred in by taking the minimum of all of the season
		# observations for each individual.
		dplyr::distinct(subject_id, season_short, study) |>
		dplyr::group_by(subject_id, study) |>
		dplyr::mutate(first_time_point = min(season_short)) |>
		dplyr::ungroup() |>
		# Now we extract the records which represent an individual's first visit
		# and we count only those by strata.
		dplyr::filter(first_time_point == season_short)

	return(new_participants)
}

# Now we create the combined table
create_count_by_study_table <- function(
		distinct_measurements, distinct_personyears, distinct_new_participants,
		file_path
) {
	set_gtsummary_theme()

	hai_count_subtbl <-
		distinct_measurements |>
		create_count_subtable("Paired HAI assays")

	py_count_subtbl <-
		distinct_personyears |>
		create_count_subtable("Person years")

	ind_count_subtbl <-
		distinct_new_participants |>
		create_count_subtable("New participants")

	combined_counts_tbl <-
		list(hai_count_subtbl, py_count_subtbl, ind_count_subtbl) |>
		# Stack the three tables on top of each other since they have the
		# same columns
		gtsummary::tbl_stack() |>
		# Convert to a flextable object for Word compatibility
		gtsummary::as_flex_table() |>
		# Make some final edits to change the row headers that say "Total, n" to
		# just say "Overall"
		flextable::compose(
			i = c(5, 10, 15), j = 1,
			part = "body",
			value = flextable::as_paragraph("Overall")
		)

	save_file_default(
		combined_counts_tbl,
		file_path
	)

	return(file_path)
}

## Data counts ####
# First get counts of total measurements, number of people, and number of
# person years so we can refer to them inline in the MS
write_data_counts <- function(
		distinct_measurements, distinct_personyears, distinct_new_participants,
		file_path
	) {
	data_counts <- list(
		"measurements" = distinct_measurements |> nrow(),
		"individuals" = distinct_new_participants |> nrow(),
		"personyears" = distinct_personyears |> nrow()
	)

	save_file_default(
		data_counts,
		file_path
	)

	return(file_path)
}

# Get the summary statistics about the number of measurements added per person
# year so we can include it in-text.
calculate_measurements_per_person_year <- function(
		distinct_measurements, file_path
	) {
	avg_measurements_per_py <-
		distinct_measurements |>
		dplyr::group_by(study, subject_id, season_short) |>
		dplyr::count() |>
		dplyr::ungroup() |>
		dplyr::reframe(
			quibble(n)
		) |>
		tibble::deframe()

	save_file_default(
		avg_measurements_per_py,
		file_path
	)

	return(file_path)
}

## Demographics per person table ####
# Make a table that summarizes the features of the study population
make_demographics_table <- function(model_data, file_path) {
	set_gtsummary_theme()

	tbl_demographics_data <-
		model_data |>
		# First we need to calculate an age variable that is invariant over seasons,
		# so we'll use the age at first enrollment (aafe for short), which is just
		# the minimum age ever reported for a person.
		dplyr::group_by(subject_id, study) |>
		dplyr::mutate(
			aafe = min(age),
			race = forcats::fct_recode(
				race,
				"Black or African American" = "Black or American American"
			)
		) |>
		dplyr::ungroup() |>
		# Now we need to do some counting, so first we'll ensure that we get only
		# one record per unique HAI assay
		dplyr::transmute(
			subject_id, study, sex, race, aafe, birth_year,
			season_short, strain_name, vaccine_name
		) |>
		dplyr::distinct() |>
		# And then add a count variable of the number of HAIs per unique person
		dplyr::add_count(subject_id, study, name = "meas") |>
		# Now we'll count the number of person years per person, so we need to remove
		# some variables and get only distinct records per person-year.
		dplyr::select(-strain_name, -vaccine_name) |>
		dplyr::distinct() |>
		# Now count the number of records per person, which gives us the number of
		# person-years contributed by an individual.
		dplyr::add_count(subject_id, study, name = "py") |>
		# Now we filter down to only have one unique record per person with the
		# demographic variables that are season-invariant.
		dplyr::select(-season_short) |>
		dplyr::distinct()

		# And now we can make a table one where the statistical unit is a unique
		# person since we now only have one record per person.
	tbl_demographics <-
		tbl_demographics_data |>
		gtsummary::tbl_summary(
			by = study,
			include = -subject_id,
			label = list(
				aafe ~ "Age at First Enrollment",
				race ~ "Race/Ethnicity",
				sex ~ "Sex Assigned at Birth",
				birth_year ~ "Year of Birth",
				meas ~ "Contributed HAI assays",
				py ~ "Contributed person-years"
			),
			statistic = list(
				gtsummary::all_continuous() ~ "{median} ({min} - {max})"
			),
			digits = list(
				gtsummary::all_continuous() ~ 0,
				gtsummary::all_categorical() ~ 0
			)
		) |>
		gtsummary::add_overall(last = TRUE) |>
		gtsummary::as_flex_table()

	save_file_default(
		tbl_demographics,
		file_path
	)

	return(file_path)
}

# Strain information tables ####
## Table of strains used in the heterologous panel each year
create_annual_panel_table <- function(model_data, file_path) {
	strain_panel_dat <-
		model_data |>
		dplyr::transmute(
			`Season` = season,
			`Subtype` = hgp::replace_strain_names(
				strain_name,
				from = "short",
				to = "type-subtype"
			),
			`Strain` = strain_name
		) |>
		dplyr::arrange(`Season`, `Subtype`) |>
		dplyr::distinct() |>
		#tidyr::complete(Season, tidyr::nesting(Subtype, `Short name`)) |>
		dplyr::mutate(present = "X") |>
		tidyr::pivot_wider(
			names_from = Season,
			values_from = present,
			values_fill = list(present = "")
		) |>
		dplyr::rename_with(
			.cols = tidyr::starts_with("2"),
			.fn = \(x) gsub(" - 20(.{2})", "/\\1", x)
		) |>
		dplyr::arrange(Subtype, Strain)

	col_matrix <- ifelse(
		strain_panel_dat[, -c(1, 2)] == "X",
		"gray",
		"white"
	)

	strain_panel_table <-
		strain_panel_dat |>
		flextable::flextable() |>
		flextable::bg(
			j = -c(1, 2),
			bg = col_matrix
		) |>
		flextable::merge_v(1) |>
		flextable::valign(j = 1, valign = "top") |>
		flextable::vline_left() |>
		flextable::vline_right() |>
		flextable::border_inner_v() |>
		flextable::border_inner_h()

	save_file_default(
		strain_panel_table,
		file_path
	)

	return(file_path)

}

## Table of strain names ####
# Make a table of strain names with the associated short names to include
# in the supplement.
create_strain_name_translation_table <- function(model_data, file_path) {
	strain_names_tbl <-
		model_data |>
		dplyr::mutate(
			`Subtype` = hgp::replace_strain_names(
				strain_name,
				from = "short",
				to = "type-subtype"
			),
			`Strain name` = hgp::replace_strain_names(
				strain_name,
				from = "short",
				to = "full"
			),
			`Short name` = strain_name,
			.keep = "none"
		) |>
		dplyr::arrange(`Strain name`) |>
		dplyr::distinct() |>
		flextable::flextable() |>
		flextable::merge_v(j = 1) |>
		flextable::valign(j = 1, valign = "top") |>
		flextable::fontsize(size = 10, part = "all") |>
		flextable::padding(padding.top = 0.5, padding.bottom = 0.5)

	save_file_default(
		strain_names_tbl,
		file_path
	)

	return(file_path)
}

## Vaccine Strains for Each Year ####
create_fluzone_vaccine_table <- function(model_data, file_path) {
	season_vaccines_tbl <-
		model_data |>
		dplyr::select(type_subtype, season_short, vaccine_name) |>
		dplyr::distinct() |>
		tidyr::complete(type_subtype, season_short) |>
		tidyr::pivot_wider(
			names_from = c(type_subtype),
			values_from = vaccine_name
		) |>
		# Replace NAs with a dash
		dplyr::mutate(
			dplyr::across(dplyr::everything(), as.character),
			dplyr::across(dplyr::everything(), \(x) tidyr::replace_na(x, "—"))
		) |>
		dplyr::rename(
			Season = season_short
		) |>
		flextable::flextable() |>
		# flextable::hline(
		# 	i = 1,
		# 	j = c(2:5, 6:9),
		# 	part = "header"
		# ) |>
		flextable::vline(
			j = c(1),
			part = "all"
		) |>
		flextable::fix_border_issues()

	save_file_default(
		season_vaccines_tbl,
		file_path
	)

	return(file_path)
}

# Correlation between antigenic distances figure ###############################
# Correlation plot of only vaccine distances using normalized dists
clean_normalized_correlation_data <- function(model_data, metric_set) {
	if(missing(metric_set)) {
		metric_set <- c("Temporal", "Cartographic", "p-Epitope", "Grantham")
	}

	dist_cor_dat_norm <-
		model_data |>
		dplyr::filter(
			metric %in% metric_set
		) |>
		dplyr::distinct(d_norm, metric, strain_type, strain_name, vaccine_name) |>
		tidyr::pivot_wider(names_from = metric, values_from = d_norm)

	return(dist_cor_dat_norm)
}

make_normalized_correlation_pairplot <- function(
		normalized_correlation_data, file_name
	) {
	cor_plot_raw_norm <- normalized_correlation_data |>
		dplyr::select(
			strain_type, strain_name, vaccine_name,
			Cartographic, Grantham, `p-Epitope`, Temporal
		) |>
		GGally::ggpairs(
		mapping = ggplot2::aes(color = strain_type),
		columns = 4:7,
		upper = list(continuous = GGally::wrap(
			GGally::ggally_cor,
			stars = FALSE, method = "spearman", digits = 2,
			size = 6, title = "Overall"
		)),
		lower = list(continuous = GGally::wrap(
			GGally::ggally_points,
			size = 3, shape = 16, alpha = 0.5
		)),
		diag = list(continuous = GGally::wrap(
			GGally::ggally_barDiag,
			alpha = 0.5, bins = 20, position = "identity"
		))
	)

	# The autocalculated CI for the correlation will produce a warning because
	# there are ties in the data, we don't use the CI but there is no way to
	# turn it off so mute the warning instead.
	suppressWarnings(
		cor_plot_pub_norm <-
			cor_plot_raw_norm +
			ggplot2::scale_color_viridis_d(
				begin = 0.1, end = 0.8, aesthetics = c("color", "fill")
			 ) +
			hgp::theme_ms() +
			ggplot2::theme(
				strip.text = ggplot2::element_text(
					margin = ggplot2::margin(0, 0, 2, 2),
					size = 18
				)
			) +
			ggplot2::labs(
				x = "Distance metric specified by column",
				y = "Distance metric specified by row"
			)
	)

	ggplot2::ggsave(
		file_name,
		plot = cor_plot_pub_norm,
		width = 6.5 * 1.5,
		height = 6.5 * 1.5
	)

	return(file_name)
}

## Calculate the correlations and CIs using bayesian bootstrap #####
calculate_bayesboot_correlation <- function(x, y) {
	input_mat <- cbind(x, y)
	out <- bayesboot::bayesboot(
		input_mat,
		\(d) cor(x = d[, 1], y = d[, 2], method = "spearman")
	)

	return(out)
}


calculate_correlations_with_cis <- function(
		normalized_correlation_data, metrics
	) {
	require(bayesboot, quietly = TRUE)

	# First get all the column names that represent metrics
	if (missing(metrics)) {
		metrics <-
			normalized_correlation_data |>
			dplyr::select(-strain_type, -strain_name, -vaccine_name) |>
			colnames()
	}

	# Now get all of the unique pairwise combinations of those metrics
	metric_combs <- combn(metrics, 2)
	clean_mc <- metric_combs |>
		t() |>
		`colnames<-`(c("Metric1", "Metric2")) |>
		tibble::as_tibble()

	# First we nest the data by subtype, and use the cbind trick to also estimate
	# the overall correlations across subtype
	overall_data <- normalized_correlation_data |>
		dplyr::mutate(strain_type = "Overall")

	nested_data <-
		normalized_correlation_data |>
		dplyr::bind_rows(overall_data) |>
		tidyr::nest(cor_dat = -strain_type)

	# Outer loop: iterate over each of the nested data frames
	cors_list <- purrr::map(
		nested_data$cor_dat,
		\(dat) {
			# Inner map: iterate over each of the pairwise metric combinations
			cor_res <- purrr::map(
				1:ncol(metric_combs),
				# Actual function: get bootstrapped correlation between two metrics
				\(i) {
					these_metrics <- metric_combs[, i]
					calculate_bayesboot_correlation(
						dat[[ these_metrics[[1]] ]],
						dat[[ these_metrics[[2]] ]]
					)
				}
			)
			out <- tibble::add_column(clean_mc, "res" = cor_res)
		}
	)

	cors_df <- cors_list |>
		rlang::set_names(nested_data$strain_type) |>
		dplyr::bind_rows(.id = "strain_type") |>
		dplyr::mutate(
			res_ci = purrr::map(res, \(d) ggdist::mean_hdci(d$V1))
		) |>
		dplyr::select(-res) |>
		tidyr::unnest(res_ci)

	return(cors_df)
}

create_cor_ci_table <- function(cors_df, metrics = NULL) {
	if(!is.null(metrics)) {
		cors_df <-
			cors_df |>
			dplyr::filter(
				Metric1 %in% metrics,
				Metric2 %in% metrics
			)
	}

	cors_formatted <-
		cors_df |>
		dplyr::transmute(
			strain_type = clean_subtype(strain_type),
			Metric1, Metric2,
			dplyr::across(
				c(y, ymin, ymax),
				\(x) formatC(x, digits = 2, format = "f", flag = "0- ")
			)
		) |>
		dplyr::mutate(
			ci = paste0(y, " (", ymin, ",", ymax, ")"),
			.keep = "unused"
		)

	cors_wide <-
		cors_formatted |>
		tidyr::pivot_wider(
			names_from = Metric2,
			values_from = ci
		) |>
		dplyr::rename(" " = Metric1, "Subtype" = strain_type)

	cors_tbl <-
		cors_wide |>
		flextable::flextable()

	k <- flextable::ncol_keys(cors_tbl)

	out <-
		cors_tbl |>
		flextable::merge_v(j = 1) |>
		flextable::valign(j = 1, valign = "top") |>
		flextable::align(j = 3:k, align = "right") |>
		flextable::fix_border_issues() |>
		flextable::autofit()

	return(out)
}

## Make a table of the correlations

# ---- ICC for metrics ----
create_icc_input_data <- function(model_data, metric_set) {
	if(missing(metric_set)) {
		metric_set <- c("Temporal", "Cartographic", "p-Epitope", "Grantham")
	}

	dist_cor_dat <-
		model_data |>
		dplyr::filter(
			metric %in% metric_set
		) |>
		dplyr::distinct(d, metric, strain_type, strain_name, vaccine_name) |>
		tidyr::pivot_wider(names_from = metric, values_from = d)

	icc_data <-
		dist_cor_dat |>
		tidyr::pivot_longer(
			-c(strain_type, strain_name, vaccine_name),
			names_to = "method",
			values_to = "distance"
		)

	return(icc_data)
}

nest_icc_input_data <- function(icc_input_data) {
	icc_data_nested <- tidyr::nest(icc_input_data, data = -strain_type)
	return(icc_data_nested)
}

fit_icc_model <- function(icc_data_subset) {
	icc_model <-
		brms::brm(
			distance ~ 0 + method + (1 | strain_name) + (1 | vaccine_name),
			data = icc_data_subset,
			family = brms::brmsfamily("gaussian"),
			prior = c(
				brms::prior(student_t(3, 0, 5), class = "b"),
				brms::prior(student_t(3, 0, 1), class = "sd"), # for random effect sd_subject
				brms::prior(student_t(3, 0, 1), class = "sigma") # for residual sd
			),
			iter = 2000,
			warmup = 1000,
			chains = 12,
			cores = 12,
			seed = 12590,
			backend = "cmdstanr",
			control = list(
				adapt_delta = 0.99,
				max_treedepth = 10
			)
		)

	return(icc_model)
}

summarize_icc_from_model <- function(fitted_icc_model) {
	# Extract the posterior samples of the relevant variance components
	post <- brms::as_draws_df(fitted_icc_model) |> as.data.frame()
	sd_strain_name <- post[,"sd_strain_name__Intercept"] |> as.numeric()
	sd_vaccine_name <- post[,"sd_vaccine_name__Intercept"] |> as.numeric()
	sigma_res   <- post[,"sigma"] |> as.numeric()

	# Calculate the posterior samples of the ICC using the canonical formula
	icc_samples <-
		(sd_strain_name^2 + sd_vaccine_name^2) /
		(sd_strain_name^2 + sd_vaccine_name^2 + sigma_res^2)

	# Compute the ICC based on the posterior samples of variances
	icc_est_var <- ggdist::mean_hdci(icc_samples)

	# Compute the ICC based on the posterior predictive distribution
	# We could do this with performance::variance_decomposition but probably best
	# to make sure we do the CI's the same way
	# Get the posterior prediction samples for each one, note that
	# rows = 12k are samples
	# cols = 480 are values of dataframe to predict on
	fe_preds <- brms::posterior_predict(
		fitted_icc_model,
		re_formula = NA,
		summary = FALSE
	)
	me_preds <- brms::posterior_predict(
		fitted_icc_model,
		re_formula = NULL,
		summary = FALSE
	)

	icc_est_ppd <-
		ggdist::mean_hdci(1 - (apply(fe_preds, 1, var) / apply(me_preds, 1, var)))

	# Cleanup and export
	out <- dplyr::bind_rows(
		"var" = icc_est_var,
		"ppd" = icc_est_ppd,
		.id = "method"
	)
	return(out)
}

calculate_icc_for_subtype <- function(subtype_data) {
	icc_model <- fit_icc_model(subtype_data)
	icc_summary <- summarize_icc_from_model(icc_model)

	return(icc_summary)
}

combine_icc_summaries <- function(subtypes_list, summaries_list) {
	out <- summaries_list |>
		rlang::set_names(subtypes_list) |>
		dplyr::bind_rows(.id = "strain_type")

	return(out)
}

create_icc_summary_table <- function(icc_summary, file_path, column) {
	icc_table_alt <-
		icc_summary |>
		dplyr::mutate(
			dplyr::across(
				c(y, ymin, ymax),
				\(x) formatC(x, digits = 2, format = "f")
			)
		) |>
		dplyr::transmute(
			"Strain Type" = strain_type,
			method,
			res = paste0(y, " (", ymin, ", ", ymax, ")")
		) |>
		tidyr::pivot_wider(
			names_from = method,
			values_from = res
		) |>
		dplyr::rename(
			"ICC" = var,
			"PPD Ratio" = ppd
		) |>
		dplyr::select(`Strain Type`, {{column}}) |>
		flextable::flextable() |>
		flextable::fontsize(size = 12, part = "all") |>
		flextable::autofit()

	save_file_default(
		icc_table_alt,
		file_path
	)

	return(file_path)
}

# ---- Measures of Spread ----
normalize_icc_data <- function(icc_data) {
	norm_icc_dat <-
		icc_data |>
		dplyr::group_by(strain_type, method) |>
		dplyr::mutate(distance = minmax(distance))

	return(norm_icc_dat)
}

create_metrics_pc_plot <- function(normalized_icc_data, file_path) {
	# Paralell coordinates plot
	pc_plot <-
		normalized_icc_data |>
		ggplot2::ggplot() +
		ggplot2::aes(
			x = method,
			y = distance
		) +
		ggplot2::geom_line(
			ggplot2::aes(group = interaction(strain_name, vaccine_name)),
			alpha = 0.5
		) +
		ggplot2::geom_point(
			alpha = 0.5,
			stroke = 1
		) +
		ggplot2::facet_wrap(ggplot2::vars(strain_type), nrow = 1) +
		ggplot2::labs(
			x = NULL,
			y = "Normalized antigenic distance"
		) +
		hgp::theme_ms() +
		ggplot2::theme(
			axis.text.x = ggplot2::element_text(
				angle = 45, vjust = 1, hjust = 1, size = 14
			)
		)

	ggplot2::ggsave(
		file_path,
		plot = pc_plot,
		width = 13,
		height = 8
	)

	return(pc_plot)
}

gap_sd <- function(x, weights = NULL) {
	# Sort x (and weights, if provided)
	ord <- order(x)
	s <- x[ord]

	# Compute gaps between successive values
	gaps <- diff(s)

	# If no weights are given, return the regular SD of the gaps
	if(is.null(weights)) {
		if(length(s) < 2) return(NA)
		return(sd(gaps))
	}

	# Otherwise, sort weights correspondingly
	w <- weights[ord]
	# Define a weight for each gap as the average of the weights of its endpoints
	gap_w <- (w[-length(w)] + w[-1]) / 2

	# Compute the weighted mean of the gaps
	wmean <- sum(gap_w * gaps) / sum(gap_w)
	# Compute the weighted variance of the gaps
	wvar  <- sum(gap_w * (gaps - wmean)^2) / sum(gap_w)
	sqrt(wvar)
}

bayes_bootstrap_wrapper_gap_sd <- function(x, ...) {
	bb <- bayesboot::bayesboot(x, gap_sd, ...)
	ci <- ggdist::mean_hdci(bb$V1)
	return(ci)
}

calculate_gap_sd_with_bootstrap <- function(normalized_icc_data) {
	gap_sd_boot <-
		normalized_icc_data |>
		dplyr::summarize(
			bayes_bootstrap_wrapper_gap_sd(
				distance,
				use.weights = TRUE,
				R2 = 1e5
			),
			.groups = "drop"
		)

	return(gap_sd_boot)
}

create_gap_sd_plot <- function(gap_sd_bootstrap_results, file_path) {
	disp_plot <-
		ggplot2::ggplot(gap_sd_bootstrap_results) +
		ggplot2::aes(
			x = method,
			y = y,
			ymin = ymin,
			ymax = ymax
		) +
		ggplot2::geom_errorbar(
			width = 0.1,
			linewidth = 1
		) +
		ggplot2::geom_point(
			shape = 95,
			color = "red",
			size = 10,
			stroke = 4
		) +
		ggplot2::facet_wrap(ggplot2::vars(strain_type), nrow = 1) +
		ggplot2::labs(
			x = NULL,
			y = "Gap standard deviation"
		) +
		hgp::theme_ms() +
		ggplot2::theme(
			axis.text.x = ggplot2::element_text(
				angle = 45, vjust = 1, hjust = 1, size = 14
			)
		)

	ggplot2::ggsave(
		file_path,
		plot = disp_plot,
		width = 13,
		height = 6
	)

	return(disp_plot)
}

create_combined_metrics_plot <- function(pc_plot, disp_plot, file_path) {
	library(patchwork)

	p1 <- pc_plot + ggplot2::theme(strip.text = ggplot2::element_text(size = 24))
	p2 <- disp_plot + ggplot2::theme(strip.text = ggplot2::element_blank())

	comb_metrics_plot <-
		(p1) +
		(p2) +
		patchwork::plot_layout(ncol = 1, axes = "collect") +
		patchwork::plot_annotation(
			tag_levels = "A",
			tag_suffix = ")"
		) &
		ggplot2::theme(
			axis.title.y = ggplot2::element_text(size = 24),
			axis.text.x = ggplot2::element_text(size = 18),
			axis.text.y = ggplot2::element_text(size = 20),
			plot.tag = ggplot2::element_text(size = 24)
		)

	ggplot2::ggsave(
		file_path,
		plot = comb_metrics_plot,
		width = 13,
		height = 12
	)
	return(file_path)
}

# Pre and post titer plots ####
# This code makes nice plots showing a summary of the pre and post titer to
# each strain.

clean_data_for_titer_plots <- function(model_data) {
	dat_hai_counts <-
		model_data |>
		dplyr::mutate(
			strain_type = ifelse(
				strain_type_is_b_pre,
				"B-Pre",
				as.character(strain_type)
			) |>
				factor(
					levels = c("H1N1", "H3N2", "B-Pre", "B-Vic", "B-Yam"),
					labels = c("A(H1N1)", "A(H3N2)", "B/Pre", "B/Victoria",
										 "B/Yamagata")
				)
		) |>
		dplyr::select(log_pretiter, log_posttiter, strain_name, strain_type) |>
		tidyr::pivot_longer(
			cols = c(log_pretiter, log_posttiter),
			names_to = "time",
			values_to = "titer"
		) |>
		dplyr::group_by(strain_type, strain_name, time) |>
		dplyr::summarize(
			ggplot2::median_hilow(titer, conf.int = 0.5),
			.groups = "drop"
		) |>
		dplyr::mutate(dplyr::across(
			c(y, ymin, ymax),
			hgp::hai_to_natural_scale
		))

	return(dat_hai_counts)
}

create_raw_data_plot <- function(titer_plot_data, timepoint, file_path) {
	# Validate the timepoint argument
	if (timepoint == "pre") {
		titer_var <- "log_pretiter"
		timepoint_label <- "Pre-vaccination"
	} else if (timepoint == "post") {
		titer_var <- "log_posttiter"
		timepoint_label <- "Post-vaccination"
	} else {
		stop(paste0(
			"Argument `timepoint` should be 'pre' or 'post', not: ",
			timepoint, "."
		))
	}

	pre_A <-
		titer_plot_data |>
		dplyr::filter(
			time == titer_var,
			startsWith(as.character(strain_type), "A")
		) |>
		ggplot2::ggplot() +
		ggplot2::aes(x = strain_name, y = y, ymin = ymin, ymax = ymax) +
		ggplot2::geom_pointrange() +
		ggplot2::scale_y_continuous(
			trans = 'log2',
			breaks = hgp::hai_to_natural_scale(seq(0, 12, 1))
		) +
		ggplot2::facet_wrap(~strain_type, scales = "free_x") +
		ggplot2::labs(
			x = NULL,
			y = paste0(timepoint_label, " HAI titer")
		) +
		hgp::theme_ms()

	pre_B <-
		titer_plot_data |>
		dplyr::filter(
			time == titer_var,
			startsWith(as.character(strain_type), "B")
		) |>
		ggplot2::ggplot() +
		ggplot2::aes(x = strain_name, y = y, ymin = ymin, ymax = ymax) +
		ggplot2::geom_pointrange() +
		ggplot2::scale_y_continuous(
			trans = 'log2',
			breaks = hgp::hai_to_natural_scale(seq(0, 12, 1))
		) +
		ggplot2::facet_wrap(~strain_type, scales = "free_x") +
		ggplot2::labs(
			x = NULL,
			y = paste0(timepoint_label, " HAI titer")
		) +
		hgp::theme_ms()


	pre_plot <-
		pre_A / pre_B +
		patchwork::plot_layout(guides = "collect", axes = "collect") &
		ggplot2::theme(
			legend.position = "bottom",
			legend.key.width = ggplot2::unit(0.5, "in"),
			axis.text.x = ggplot2::element_text(
				size = 16, hjust = 1, vjust = 1, angle = 45
			),
			axis.text.y = ggplot2::element_text(size = 18),
			strip.text = ggplot2::element_text(size = 28),
			axis.title = ggplot2::element_text(size = 28),
			legend.title = ggplot2::element_text(size = 24),
			legend.text = ggplot2::element_text(size = 18)
		)

	dir.create(dirname(file_path), showWarnings = FALSE, recursive = TRUE)
	ggplot2::ggsave(
		plot = pre_plot,
		filename = file_path,
		width = 13,
		height = 13
	)

	return(file_path)
}
