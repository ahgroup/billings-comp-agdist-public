###
# Code for supplement and misc stuff
# Zane
# 2025-01-27
# Targetized 2025-04-17
###

# Function to create the DAG figure
create_dag <- function(file_path) {
	suppressPackageStartupMessages({
		library(dagitty)
		library(ggdag)
	})

	# Code for the DAG figure in the supplement.
	dag <- dagitty::dagitty(
		"dag {
	{sv sa U} -> ag -> y
	{U C sv sa} -> y
	{U C} -> p -> y
	ag [exposure]
	y [outcome]
	}")
	dagitty::coordinates(dag) <- list(
		x = c(
			y = 3, ag = 1, sv = 1.5, sa = 2.5, U = 1.5,	p = 2, C = 2.5
		),
		y = c(
			y = 0, ag = 0, sv = -0.25, sa = -0.25, U = 0.25, p = 0.25, C = 0.25
		)
	)

	tidy_dag <- ggdag::tidy_dagitty(dag)
	tidy_dag$data <-
		tidy_dag$data |>
		dplyr::mutate(lab = dplyr::case_match(
			name,
			"y" ~ "Outcome",
			"ag" ~ "Antigenic Distance",
			"sv" ~ "Vaccine strain",
			"sa" ~ "Assay strain",
			"U" ~ "U",
			"p" ~ "Pre-vaccination titer",
			"r" ~ "Race/ethnicity",
			"s" ~ "Sex assigned at birth",
			"b" ~ "Birth year",
			"a" ~ "Age"
		))

	my_ggdag <-
		tidy_dag |>
		ggdag::node_status() |>
		dplyr::mutate(
			observed = factor(
				ifelse(name %in% c("U", "i"), 0, 1),
				levels = c(1, 0),
				labels = c("Observed", "Unobserved")
			),
			status = factor(
				ifelse(is.na(status), "bleh", status),
				levels = c("1", "2", "bleh"),
				labels = c("Exposure", "Outcome", "Covariate")
			)
		) |>
		ggdag::ggdag() +
		ggdag::geom_dag_point(ggplot2::aes(color = status, shape = observed)) +
		ggdag::geom_dag_text() +
		ggplot2::scale_color_manual(values = c("#56B4E9", "#E69F00", "black")) +
		ggplot2::scale_shape_manual(values = c(16, 15)) +
		ggdag::theme_dag() +
		ggplot2::theme(
			plot.background = ggplot2::element_rect(fill = "white", color = "white")
		)

	ggplot2::ggsave(
		file_path,
		plot = my_ggdag,
		width = 9,
		height = 5,
		dpi = 300
	)

	return(file_path)
}

# Correlations between different sequence-based measures
analyze_sequence_distance_correlation <- function(
		model_data
	) {
	metrics <- c(
		"Hamming", "p-Epitope", "p-All-Epitope", "Grantham", "FLU Substitution"
	)
	metric_combs <- combn(metrics, 2)
	clean_mc <- metric_combs |>
		t() |>
		`colnames<-`(c("Metric1", "Metric2")) |>
		tibble::as_tibble()

	# Normalize the antigenic distances (had to do it here to avoid outdating
	# all model targets that depend on model_data)
	normalized_correlation_data <-
		model_data |>
		dplyr::group_by(strain_type, metric) |>
		dplyr::mutate(d_norm = minmax(d)) |>
		dplyr::ungroup() |>
		# Drop the dose column since everyone is SD
		dplyr::select(-dose) |>
		dplyr::mutate(metric = forcats::fct_drop(metric)) |>
		clean_normalized_correlation_data(metric_set = metrics) |>
		tidyr::pivot_longer(
			-c(strain_type, strain_name, vaccine_name),
			names_to = "metric",
			values_to = "d"
		)

	dist_cor_dat <-
		normalized_correlation_data |>
		dplyr::mutate(
			type_subtype =
				hgp::replace_strain_names(vaccine_name, "short", "type-subtype"),
			StrainComb = purrr::map_chr(
				1:nrow(normalized_correlation_data),
				\(i) c(
					normalized_correlation_data$vaccine_name[[i]],
					normalized_correlation_data$strain_name[[i]]
				) |>
					as.character() |>
					sort() |>
					paste0(collapse = "|")
			)
		) |>
		dplyr::distinct(type_subtype, metric, d, StrainComb, .keep_all = TRUE) |>
		dplyr::select(-StrainComb) |>
		tidyr::pivot_wider(names_from = metric, values_from = d) |>
		dplyr::mutate(type_group = type_subtype)

	# NRowbinding trick to add an Overall category
	dist_cor_dat_all <-
		dplyr::bind_rows(
			dist_cor_dat,
			dplyr::mutate(dist_cor_dat, type_group = "Overall")
		)

	# Copied from calculate_correlations_with_cis() so I don't have to adjust
	# that function to take a non-nested argument
	nested_data <-
		dist_cor_dat_all |>
		dplyr::rename(
			Strain1 = vaccine_name,
			Strain2 = strain_name
		) |>
		dplyr::select(
			type_group, type_subtype, Strain1, Strain2,
			dplyr::all_of(metrics)
		) |>
		tidyr::nest(cor_dat = -c(type_group))

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
		rlang::set_names(nested_data$type_group) |>
		dplyr::bind_rows(.id = "strain_type") |>
		dplyr::mutate(
			strain_type = clean_subtype(strain_type),
			res_ci = purrr::map(res, \(d) ggdist::mean_hdci(d$V1))
		) |>
		dplyr::select(-res) |>
		tidyr::unnest(res_ci)

	return(cors_df)
}

# Generate software bibliography -- this should always be right before
# the manuscript and supplement
generate_software_bibliography <- function(file_path) {
	require(softbib, quietly = TRUE)

	unlink(file_path)
	softbib::softbib(
		output = file_path
	)

	return(file_path)
}

# Function to make the accession number table formatted nicely
format_accession_number_table <- function(
		accession_number_table_data, file_path
	) {
	out <-
		accession_number_table_data |>
		flextable::flextable() |>
		flextable::fontsize(size = 7.5, part = "all") |>
		flextable::padding(padding = 0.01, part = "all")

	save_file_default(out, file_path)

	return(file_path)
}
