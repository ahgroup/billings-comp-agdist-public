###
# Utility functions
# Zane
# 2025-01-23
###

# Clean up the strain type to use Ben's preferred notation
clean_subtype <- function(x) {
	y <- tolower(x)
	out <- dplyr::case_when(
		grepl("h1n1", y) ~ "A(H1N1)",
		grepl("h3n2", y) ~ "A(H3N2)",
		grepl("vic", y) ~ "B/Victoria",
		grepl("yam", y) ~ "B/Yamagata",
		grepl("pre", y) ~ "B/Pre",
		grepl("overall", y) ~ "Overall",
		TRUE ~ NA_character_
	)

	out <- factor(
		out,
		levels = c("A(H1N1)", "A(H3N2)", "B/Victoria", "B/Yamagata", "B/Pre",
							 "Overall")
	)

	if (any(is.na(out))) {
		stop(paste0(
			"Unsupported pattern(s) in x: ",
			paste(x[is.na(out)], collapse = "; ")
		))
	}

	return(out)
}

minmax <- function(x, ...) {
	return((x - min(x, ...)) / (max(x, ...) - min(x, ...)))
}

# ensure_directory_exists <- function(file_path) {
# 	dir.create(dirname(file_path), recursive = TRUE, showWarnings = FALSE)
#
# 	invisible(file_path)
# }

# format_td <- function(td, digits = 1L) {
# 	value <- round(as.numeric(td), digits)
# 	unit <- attr(td, "units")
#
# 	return(paste(value, unit))
# }

# https://www.tidyverse.org/blog/2020/03/dplyr-1-0-0-summarise/
# TODO move to zlib
quibble <- function(x, q = c(0, 0.25, 0.5, 0.75, 1)) {
	tibble::tibble(q = q, x = stats::quantile(x, q))
}

# Generic function for saving plots


# Function to set the gtsummary theme
set_gtsummary_theme <- function() {
	suppressMessages({
		gtsummary::theme_gtsummary_journal(journal = "jama")
		gtsummary::theme_gtsummary_compact(font_size = 10)
		gtsummary::theme_gtsummary_language(language = "en", big.mark = "")
	})

	invisible(TRUE)
}

# Function to clean up the metrics variable of a dataset for analysis
clean_metrics <- function(.data, keep_set = "models") {
	# Set up metric levels and labels
	metric_levels <- c(
		"cartographic",
		"cophenetic_mltree",
		"dl",
		"dominant_pepitope",
		"flu_sub",
		"forward_temporal",
		"grantham",
		"hamming",
		"p_all_epitope",
		"absolute_temporal"
	)

	metric_labels <- c(
		"Cartographic",
		"Tree",
		"Damerau-Levenshtein",
		"p-Epitope",
		"FLU Substitution",
		"Temporal (asymmetric)",
		"Grantham",
		"Hamming",
		"p-All-Epitope",
		"Temporal"
	)

	# Validate inputs
	if (!is.data.frame(.data)) {
		stop(paste0(
			"`.data` should be a data frame, but instead it is an object of class: ",
			class(.data), "."
		))
	}

	if (!("metric" %in% colnames(.data))) {
		stop("There is no column in `.data` named 'metric'.")
	}

	if (!(keep_set %in% c("models", "extended", "all"))) {
		stop(paste0(
			"`keep_set` should be one of the following: 'models', 'extended', or ",
			"'all'. Instead it is: ", keep_set, "."
		))
	} else {
		vars_out <- switch(
			keep_set,
			"models" = c(
				"absolute_temporal", "cartographic", "cophenetic_mltree",
				"dominant_pepitope", "grantham", "p_all_epitope"
			),
			"extended" = c(
				"absolute_temporal", "cartographic", "cophenetic_mltree",
				"dominant_pepitope", "flu_sub", "grantham",
				"hamming", "p_all_epitope"
			),
			"all" = metric_levels
		)
	}

	# Filter out the ones we don't ever want and then make the ones we want to
	# keep look nice and drop levels that we already removed.
	out <- .data |>
		dplyr::filter(metric %in% vars_out) |>
		dplyr::mutate(
			metric = factor(metric, levels = metric_levels, labels = metric_labels),
			metric = forcats::fct_drop(metric)
		)

	return(out)
}

# Convenience function for saving one output as a rds file
# for use with tar_file.
save_file_default <- function(file_to_save, file_path) {
	file_path |>
		dirname() |>
		dir.create(showWarnings = FALSE, recursive = TRUE)
	readr::write_rds(file_to_save, file_path, compress = "xz")
	return(file_path)
}

# Function for making the flextables fit correctly in a word document
# In quarto 1.6.40 there are still a lot of issues with tables especially from
# gtsummary that you have to manually fix, but this is still the best thing
# to try first.
fit_flextable_to_page <- function(tbl, page_width = 6.5){

	ft_out <- tbl |> flextable::autofit()

	ft_width <- (dim(ft_out)$widths * page_width) /
		(flextable::flextable_dim(ft_out)$widths)

	ft_out <- flextable::width(ft_out, width = ft_width)

	return(ft_out)
}
