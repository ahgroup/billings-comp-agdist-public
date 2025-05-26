###
# Generating the required metadata files for submission to the CIVR-HRP repo
# Zane
# 2025-04-23
# The required files for submission at this time are:
# Protocol
# Study Design
# Study Personnel
# Study Inclusion
# Study Arms
# Study Planned Visits
# Interventions
# Human Subjects Demographics
# Computational Model
# Additionally, we must submit all SOPs and raw files included in the
# manuscript.
###

# Setup and common stuff ####
write_metadata_file <- function(metadata_sheet, header, path) {
	# Make sure the file path exists
	dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)

	# Open the file for writing
	con <- file(path)
	open(con, "w")

	# Write the header and column names in csv format
	write(header, con)
	write(colnames(metadata_sheet) |> paste0(collapse = ","), con, append = TRUE)

	# Write each row of the dataset
	for (i in 1:nrow(metadata_sheet)) {
		cur_row <- metadata_sheet[i, ] |> as.character()
		cur_txt <- paste0("\"", cur_row, "\"") |> paste0(collapse = ",")
		write(cur_txt, con, append = TRUE)
	}

	# Close the file to avoid corruption
	close(con)
	invisible(path)
}

# Global variables for things that get repeated a lot
set_common_metadata_variables <- function() {
	out <- list(
		STUDY_CODE  = "Study-320_antigenic_dist",
		PROTOCOL_ID = "Handelgroup_BayesModel"
	)

	return(out)
}

# Metadata sheets ####
## Sheet 1: Protocol ####
generate_metadata_sheet_protocol <- function(
		mdv, file_path
	) {
	protocol_data <- tibble::tibble(
		"Study_Code" = mdv$STUDY_CODE,
		"Protocol_ID" = mdv$PROTOCOL_ID,
		"Protocol_Name" = paste(
			"Bayesian Analysis and Modeling Protocol for Handelgroup Computational",
			"Projects"
		),
		"Protocol_Description" = paste(
			"Protocol for building Bayesian models,using R and Stan, via rstan,",
			"cmdstanr, or brms."
		),
		"Protocol_Type" = "Statistical Analysis Documentation",
		"Associated_Protocols" = "Not Applicable"
	)

	write_metadata_file(
		protocol_data,
		"#SDMCCDataTemplate:Protocol_v1.0",
		file_path
	)

	return(file_path)
}

## Sheet 2: Study Design ####
generate_metadata_sheet_study_design <- function(
		reporting_data, mdv, file_path
	) {
	study_design_data <- tibble::tibble(
		Study_Code = mdv$STUDY_CODE,
		Hypothesis = paste(
			"Sequence-based antigenic distance metrics generate similar predictions",
			"of vaccine breadth to antigenic distance metrics from cartography when",
			"used for modeling, even though the different distance metrics have low",
			"reliability."
		),
		Objectives = paste(
			"Compute cartographic distances using the HAI titer matrix from the",
			"UGAFluVac data, along with multiple sequence based distances.",
			"Calculate the reliability across measurements, and then generate",
			"computational models to predict vaccine-induced immunogenicity at",
			"multiple levels of antigenic distance, for each distance metric."
		),
		Endpoints = paste(
			"Pre and post-vaccination HAI titer measurements to Fluzone SD and HD",
			"vaccine to homologous strains and a panel of historical strains."
		),
		Target_Enrollment = "Not Provided",
		Minimum_Age = min(reporting_data$Min_Age),
		Maximum_Age = max(reporting_data$Max_Age),
		Age_Unit = "Years",
		Actual_Start_Date = "2024-11-04", # Date of first git commit
		Research_Focus = "Vaccine Response",
		Embargo_Date = "Not Applicable",
		Protocol_ID = mdv$PROTOCOL_ID,
		Link_Name = "Not Applicable",
		Link_URL = "Not Applicable",
		National_Clinical_Trial_Identifier = "Not Applicable",
		PMID = "Not Applicable"
	)

	write_metadata_file(
		study_design_data,
		"#SDMCCDataTemplate:Study_Design_v1.0",
		file_path
	)

	return(file_path)
}

## Sheet 3: Study Personnel ####
generate_metadata_sheet_personnel <- function(mdv, file_path) {
	personnel_data <- tibble::tibble(
		"Study_Code" = mdv$STUDY_CODE,
		"Email" = c("wesley.billings@uga.edu", "ahandel@uga.edu"),
		"Title_in_Study" = c("Primary researcher", "Study supervisor"),
		"Role_in_Study" = c("Sub-Investigator", "Principal-Investigator"),
		"Contributing_Institution_Code" = "UGAN_CIVR-HRP"
	)

	write_metadata_file(
		personnel_data,
		"#SDMCCDataTemplate:Study_Personnel_v1.0",
		file_path
	)

	return(file_path)
}

## Sheet 4: Study Inclusion ####
generate_metadata_sheet_study_inclusion <- function(mdv, file_path) {
	inclusion_data <- tibble::tibble(
		"Study_Code" = mdv$STUDY_CODE,
		"Criterion_ID" = c(
			"SD", "Season", "AllVisits"
		),
		"Criterion_Category" = c(
			"Exclusion", "Inclusion", "Exclusion"
		),
		"Criterion_Description" = c(
			"Only standard dose vaccine recipients were recruited.",
			"Records from the 2013/14 to 2017/2018 seasons were used.",
			paste(
				"Only patients with data from both a pre- and post-vaccination visit",
				"were used."
			)
		)
	)

	write_metadata_file(
		inclusion_data,
		"#SDMCCDataTemplate:Study_Inclusion_v1.0",
		file_path
	)

	return(file_path)
}
## Sheet 5: Study Arms ####
format_study_arms_data <- function(
		mdv
) {
	# Make a grid of study sites and years
	pa_fl_yrs <- tidyr::expand_grid(study = c("PA", "FL"), season = 2013:2016)
	uga_yrs <- tidyr::expand_grid(study = "UGA", season = 2016:2017)
	study_yrs <- dplyr::bind_rows(pa_fl_yrs, uga_yrs)

	# Add the two doses to the grid
	study_yrs_dose <- tidyr::expand_grid(study_yrs, dose = c("SD"))

	# Format the data correctly and add the other columns needed
	study_yrs_format <-
		study_yrs_dose |>
		dplyr::mutate(
			Cohort_ID = paste(paste0(study, season), dose, sep = "_"),
			long_place = dplyr::case_match(
				study,
				"PA" ~ "Pittsburgh, PA",
				"FL" ~ "Port St. Lucie, FL",
				"UGA" ~ "Athens, GA"
			),
			long_dose = dplyr::case_when(
				dose == "SD" & season <= 2014 ~ "Standard Dose TIV",
				dose == "SD" & season >= 2015 ~ "Standard Dose QIV",
				dose == "HD" & season <= 2019 ~ "High Dose TIV",
				dose == "HD" & season >= 2020 ~ "High Dose QIV"
			),
			long_season = paste0(season, "-", season + 1),
			Cohort_Name = paste(study, season, "Cohort", long_dose),
			Cohort_Description = paste(
				long_place, "Cohort, individuals who received", long_dose,
				"FluZone vaccine during the", long_season, "influenza season."
			),
			Intervention_ID = gsub(
				"Standard Dose ", "SD_", long_dose
			),
			Intervention_ID = gsub(
				"High Dose ", "HD_", Intervention_ID
			),
			Intervention_Name = paste(long_dose, "FluZone vaccine"),
			Product_Code = paste0(
				"Fluzone ", ifelse(dose == "HD", "High-Dose ", ""),
				ifelse(dose == "HD" & season >= 2020, "Quadrivalent ", ""),
				ifelse(dose == "SD" & season >= 2015, "Quadrivalent ", ""),
				long_season
			)
		)

	return(study_yrs_format)
}

generate_metadata_sheet_study_arms <- function(
		study_arms_data_input, mdv, file_path
	) {
	# Select the correct data for the study arms template
	study_arms_data <-
		study_arms_data_input |>
		dplyr::transmute(
			"Study_Code" = mdv$STUDY_CODE,
			Cohort_ID, Cohort_Name, Cohort_Description,
			Cohort_Type = "Observational"
		)

	# Save the template
	write_metadata_file(
		study_arms_data,
		"#SDMCCDataTemplate:Study_Arms_v1.0",
		file_path
	)

	return(file_path)
}

## Sheet 6: Study Planned Visits ####
generate_metadata_sheet_planned_visits <- function(
		mdv, file_path
	) {
	planned_visits_data <- tibble::tibble(
		"Study_Code" = mdv$STUDY_CODE,
		"Planned_Visit_ID" = c("Pre-Vac-Visit", "Post-Vac-Visit"),
		"Planned_Visit_Name" = c("Pre-vaccination visit", "Post-vaccination visit"),
		"Planned_Visit_Order_Number" = c(1, 2),
		"Interval" = "Not Applicable",
		"Planned_Visit_Min_Start_Day" = c(0, 21),
		"Planned_Visit_Max_Start_Day" = c(0, 28),
		"Planned_Visit_Start_Rule" = c(
			paste(
				"Individual presents to center where vaccine is available and has",
				"no contraindictions for administration of FluZone vaccine."
			),
			paste(
				"Previously vaccinated individual presents to study site at appointment",
				"time with no contraindications for sample collection."
			)
		),
		"Planned_Visit_End_Rule" = c(
			paste(
				"Subject completes intake form with demographic questionnaire,",
				"completes all pre-vaccination samples (e.g. serum), and vaccine is",
				"administered."
			),
			paste(
				"Subject completes all required samples for day 28 time point."
			)
		)
	)

	write_metadata_file(
		planned_visits_data,
		"#SDMCCDataTemplate:Study_Planned_Visits_v1.1",
		file_path
	)

	return(file_path)
}

## Sheet 7: Interventions ####
generate_metadata_sheet_interventions <- function(
		study_arms_data, mdv, file_path
	) {
	interventions_data <-
		study_arms_data |>
		dplyr::transmute(
			Study_Code = mdv$STUDY_CODE,
			Cohort_ID,
			Planned_Visit_ID = "Post-Vac-Visit",
			Intervention_ID,
			Intervention_Name,
			Exposure_Process_Reported = "Vaccination",
			Product_Code,
			MAB_Code = "Not Applicable",
			Adjuvant_Code = "Not Applicable",
			Challenge_Strain_Name = "Not Applicable",
			Passage_History = "Not Applicable",
			Dose = ifelse(
				grepl("High-Dose", Product_Code, fixed = TRUE),
				15, 60
			),
			Dose_Units = "ug",
			Delivery_Volume = "Not Provided",
			Formulation = "Not Provided",
			Route_Of_Admin_Reporter = "Intramuscular"
		)

	write_metadata_file(
		interventions_data,
		"#SDMCCDataTemplate:Interventions_v1.2",
		file_path
	)

	return(file_path)
}

## Sheet 8: Human Subjects Demographics ####
generate_metadata_sheet_human_subjects_demographics <- function(
		reporting_data, mdv, file_path
	) {
	hsd_data <- reporting_data |>
		dplyr::mutate(Study_Code = mdv$STUDY_CODE)

	# And save to file
	write_metadata_file(
		hsd_data,
		"#SDMCCDataTemplate:Human_Subjects_Demographics_v1.0",
		file_path
	)

	return(file_path)
}

## Sheet 9: Computational Model ####
generate_metadata_sheet_computational_model <- function(
		mdv, file_path
	) {
	pkgs <-
		renv::dependencies() |>
		dplyr::pull(Package) |>
		unique() |>
		sort()

	# Silently attach all of them to the session.
	suppressPackageStartupMessages(
		lapply(
			pkgs,
			\(p) requireNamespace(p, quietly = TRUE)
		) |>
			invisible()
	)

	# Create and format a string with the package requirements and versions.
	si <- sessionInfo()
	nv <- sapply(si$loadedOnly, \(x) paste0(x$Package, " (", x$Version, ")"))
	packages_versions <- paste0(sort(nv), collapse = "; ")

	# Create the data for this template
	computational_model_data <- tibble::tibble(
		Study_Code = mdv$STUDY_CODE,
		agdist = c(
			c("Cartographic", "Grantham", "pEpitope", "Temporal"),
			c("Cartographic", "Grantham", "pEpitope", "Temporal")
		),
		model = c(
			rep("LMM", times = 4),
			rep("GAMM", times = 4)
		),
		Model_ID = paste0(model, "_", agdist),
		Model_Name = paste(
			"Post-vac titer vs.", agdist, "antigenic distance", model
		),
		Model_Description = paste(
			"Uses a Bayesian mixed-effects", model, "to estimate the effect of",
			agdist, "antigenic distance on post-vaccination HAI titer."
		),
		Model_Usage = paste(
			"The model is used to estimate post-vaccination titer from antigenic",
			"distance data."
		),
		Model_Type = "OTH-Statistical Model",
		Protocol_ID = "Handelgroup_BayesModel",
		Operating_System_Description = paste0(
			"Developed on Windows 10 with ", version$version.string, "."
		),
		Programs_and_Packages = packages_versions,
		Code_File_Names = "Not Applicable",
		Code_Repository_Name = "Not Applicable",
		Code_Repository_Accession = "Not Applicable",
		CIVICs_Template = "N",
		Primary_Data_Study_Code = "Not Applicable",
		Primary_Data_File_Name = "Not Applicable",
		Data_Repository_Name = "Not Applicable",
		Data_Repository_Accession = "Not Applicable",
		Derived_Data_File_Name = "Not Applicable",
		Example_Data_File_Name = "Not Applicable",
		README_File_Name = "README.md",
		Additional_File_Names = "Not Applicable"
	) |>
		dplyr::select(-agdist, -model)

	# And save to file
	write_metadata_file(
		computational_model_data,
		"#SDMCCDataTemplate:Computational_Model_v1.0",
		file_path
	)

	return(file_path)
}

# Getting figures and tables together ####
move_figure_files_to_metadata_location <- function(
		prefix, base_dir, file_list
	) {
	# Make sure both arguments are character vectors and prefix is length 1
	if (!is.character(prefix) || !(length(prefix) == 1)) {
		stop("`prefix` should be a length one character vector.")
	}
	if (!is.character(base_dir) || !(length(base_dir) == 1)) {
		stop("`base_dir` should be a length one character vector.")
	}

	if (!is.character("file_list")) {
		stop("`file_list` should be a character vector.")
	}

	# Make sure all the files exist
	file_check <- file.exists(file_list)
	if (!all(file_check)) {
		plural <- (sum(!file_check) != 1)
		stop(paste0(
			"The following file", ifelse(plural, "s", ""), " in `file_list` ",
			ifelse(plural, "do", "does"), " not exist:\n",
			paste("•", file_list[!file_check], collapse = "\n")
		))
	}

	# Create the new file names
	nums <-
		length(file_list) |>
		seq_len()

	max_len <- max(nchar(nums))

	file_names <- paste0(
		"Figure_", prefix,
		stringr::str_pad(nums, width = max_len, side = "left", pad = "0"),
		".", tools::file_ext(file_list)
	)

	file_paths <- here::here(base_dir, file_names)

	# Move the old file names to new locations
	purrr::walk(
		dirname(file_paths),
		\(f) dir.create(f, showWarnings = FALSE, recursive = TRUE)
	)

	success <- file.copy(file_list, file_paths, overwrite = TRUE)

	if (!all(success)) {
		stop("Not every file copy was a success!")
	}

	return(file_paths)
}

move_table_data_to_metadata_location <- function(
		prefix, base_dir, file_list
	) {
	# Make sure both arguments are character vectors and prefix is length 1
	if (!is.character(prefix) || !(length(prefix) == 1)) {
		stop("`prefix` should be a length one character vector.")
	}
	if (!is.character(base_dir) || !(length(base_dir) == 1)) {
		stop("`base_dir` should be a length one character vector.")
	}

	if (!is.character("file_list")) {
		stop("`file_list` should be a character vector.")
	}

	# Make sure all the files exist
	file_check <- file.exists(file_list)
	if (!all(file_check)) {
		plural <- (sum(!file_check) != 1)
		stop(paste0(
			"The following file", ifelse(plural, "s", ""), " in `file_list` ",
			ifelse(plural, "do", "does"), " not exist:\n",
			paste("•", file_list[!file_check], collapse = "\n")
		))
	}

	# Create the new file names
	n <- length(file_list)
	nums <- seq_len(n)

	max_len <- max(nchar(nums))

	file_names <- paste0(
		"Table_", prefix,
		stringr::str_pad(nums, width = max_len, side = "left", pad = "0"),
		".txt"
	)

	file_paths <- here::here(base_dir, file_names)

	# Extra step for tables -- load and rip out the table data.
	# This assumes every table is a flextable and saved as a rds file
	for (i in 1:n) {
		this_file <- file_list[[i]]
		this_tbl <- readr::read_rds(this_file)

		# Use flextable's built in rtf formatting
		tmp <- tempfile()
		flextable::save_as_rtf(this_tbl, path = tmp)

		# Use striprtf to read the raw text as character vector
		txt <- striprtf::read_rtf(tmp, row_start = "| ") |>
			stringr::str_squish()

		# Add the second row we need to make it a markdown table
		# first find the number of pipe bars
		n_cols <- stringr::str_count(txt[[1]], "\\|") - 1
		# Now make the line with the dashes
		line2 <- paste0(
			"| ",
			paste0(rep(c("--- | "), times = n_cols), collapse = "")
		)
		# And put it all together
		txt_md <- c(
			txt[1],
			stringr::str_trim(line2),
			txt[2:length(txt)]
		)

		# Now do the annoying bit of opening a connection and writing the txt
		con <- file(file_paths[[i]])
		open(con, "w")
		writeLines(txt_md, con = con)
		close(con)
	}

	return(file_paths)
}

# Bundle everything into a zip folder for submission ####
create_zipped_metadata_submission <- function(
	# First all of the metadata sheets
	protocol_sheet_file,
	study_design_sheet_file,
	study_personnel_sheet_file,
	study_inclusion_sheet_file,
	study_arms_sheet_file,
	study_planned_visits_sheet_file,
	interventions_sheet_file,
	human_subjects_demographics_sheet_file,
	computational_model_sheet_file,
	# Then the figures and tables
	manuscript_figure_files,
	supplement_figure_files,
	manuscript_table_files,
	supplement_table_files,
	# And the other files
	other_files_to_include,
	# Control info
	zip_path
) {
	# Require the zip package -- we use zip because it doesn't require an
	# underlying system utility like utils::zip() does so in theory this is a
	# more reproducible solution.
	if (!requireNamespace("zip", quietly = TRUE)) {
		stop(paste0(
			"Package 'zip' is required.",
			"Please install it."
		))
	}

	# Create the master list of paths to zip
	all_paths <- c(
		protocol_sheet_file,
		study_design_sheet_file,
		study_personnel_sheet_file,
		study_inclusion_sheet_file,
		study_arms_sheet_file,
		study_planned_visits_sheet_file,
		interventions_sheet_file,
		human_subjects_demographics_sheet_file,
		computational_model_sheet_file,
		manuscript_figure_files,
		supplement_figure_files,
		manuscript_table_files,
		supplement_table_files,
		other_files_to_include
	)

	# Ensure output directory exists
	dir.create(dirname(zip_path), recursive = TRUE, showWarnings = FALSE)

	# Create the zip archive
	zip::zipr(zip_path, all_paths)
	invisible(zip_path)
}
