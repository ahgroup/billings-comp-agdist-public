###
# Data Processing
# Zane Billings
# 2024-11-11
# Targetized 2025-04-02
# Starting with the cleaned data from ahgroup/UGAFluVac repo, we need to do
# some data processing to get the data ready for our analyses.
# Also processes the multiple sources of genetic/antigenic distances.
# And sets up the CIVR-HRP reporting data with the specific options that are
# required to be reported there.
###

# Data cleaning ================================================================
# First filter out the records we don't want
prep_cohort_data <- function(raw_data) {
	# First step is basic filtering with our study's inclusion/exclusion criteria
	dat_filtered <-
		raw_data |>
		dplyr::filter(
			# Only look at SD individuals here
			dose == "SD",
			# Remove individuals with infinite titerincrease (1 person who has a
			# pretiter of 0, probably a data entry error)
			is.finite(titerincrease),
			# We really only want to use seasons up to 2017 - 2018
			# cause after that the heterologous panel got really small so the
			# cartographic distances are unreliable for some strains.
			#
			season <= "2017 - 2018"
		) |>
		dplyr::mutate(
			# Drop the missing factor levels that remain after filtering
			dplyr::across(
				dplyr::where(is.factor),
				forcats::fct_drop
			),
			# Create an identifying variable for the B-pre records, since we'll need
			# this later.
			strain_type_is_b_pre = (strain_type == "B-Pre")
		)

	# Next we'll deal with condensing the vaccine strain information into a single
	# column. We do this to make sure we're only comparing H1N1 assay strains
	# with H1N1 vaccine strains, and likewise for all the other subtypes.
	# However, for the presplit type B lineages, we need to decide whether those
	# should be compared to the B-yamagata or the B-victoria vaccine strain in
	# order to get the antigenic distance. For now we'll use a trick to get
	# both in the dataset so we can do a sensitivity analysis later.
	# First we'll select only the B-pre records.
	dat_filtered_b_pre <-
		dat_filtered |>
		dplyr::filter(strain_type_is_b_pre)

	# Now we create a new dataset which has one B-yam record and one B-vic record
	# for each of those records.
	dat_filtered_b_sub <-
		dplyr::bind_rows(
			dplyr::mutate(dat_filtered_b_pre, strain_type = "B-Vic"),
			dplyr::mutate(dat_filtered_b_pre, strain_type = "B-Yam")
		)

	# Now we remove the original B_pre records from the dataset and add these new
	# duplicated records.
	dat_filtered_b_dup <-
		dat_filtered |>
		dplyr::filter(!strain_type_is_b_pre) |>
		dplyr::bind_rows(dat_filtered_b_sub)

	# Now we need to clean up the vaccine strain columns to match.
	dat_prepped <-
		dat_filtered_b_dup |>
		# First we need to "condense" the vaccine name variables. We just need one
		# vaccine name variable that has the name of the H1 vaccine strain for H1
		# assays, and the name of the H3 vaccine strain for H3 assays.
		dplyr::mutate(
			vaccine_name = dplyr::case_match(
				strain_type,
				"H1N1" ~ h1n1_vaccine_strain,
				"H3N2" ~ h3n2_vaccine_strain,
				"B-Vic" ~ bvic_vaccine_strain,
				"B-Yam" ~ byam_vaccine_strain
			),
			# Replace the "None" values with NA
			vaccine_name = vaccine_name |>
				dplyr::na_if("None") |>
				forcats::fct_drop(),
			.after = dose
		) |>
		# Now drop the *_vaccine_fullname columns
		dplyr::select(-dplyr::ends_with("vaccine_strain")) |>
		# Now remove the observations where there isn't a vaccine strain
		tidyr::drop_na(vaccine_name)

	return(dat_prepped)
}

clean_cohort_data <- function(prepped_cohort_data) {
	# Now that we have the filtering and manipulation steps done, we'll select
	# only the variables we need for our analysis.
	dat_clean <-
		prepped_cohort_data |>
		dplyr::select(
			subject_id = id,
			season, year, study, age, birth_year, sex, dose, race = race_collapsed,
			vaccine_name, strain_type, strain_name,
			pretiter, posttiter, log_pretiter, log_posttiter, fold_change,
			titer_increase = titerincrease, seroprotection, seroconversion,
			strain_type_is_b_pre
		) |>
		# TODO fix this in UGAFluVac-data
		# there are currently two people who have multiple different race_collapsed
		# variables so they need to be standardized, we'll do that manually.
		dplyr::mutate(
			race = dplyr::if_else(
				subject_id %in% c("0045", "0072"),
				"Hispanic or Latino",
				race
			) |>
				forcats::fct_infreq() |>
				forcats::fct_na_value_to_level("Unknown")
		)

	return(dat_clean)
}

join_antigenic_distance_to_cohort_data <- function(
		cleaned_cohort_data, antigenic_distance_raw_data
	) {
	# Nest the antigenic distance data to make joining easier -- this allows us
	# to avoid making sure a many-to-many merge is working correctly.
	nested_agdist_data <-
		antigenic_distance_raw_data |>
		tidyr::nest(antigenic_distances = c(metric, d)) |>
		dplyr::select(-type_subtype)

	# Now join the antigenic distances
	dat_joined <-
		cleaned_cohort_data |>
		# First replace the vaccine and strain names with the short form
		dplyr::mutate(
			# Fix the weird Shangdong typo
			strain_name = replace(
				as.character(strain_name),
				strain_name == "H3N2-Shangdong-1993",
				"H3N2-Shandong-1993"
			),
			# Now do the strain name replacement
			dplyr::across(c(vaccine_name, strain_name), hgp::replace_strain_names),
			type_subtype = hgp::replace_strain_names(
				vaccine_name,
				from = "short", to = "type-subtype"
			),
			strain_type = factor(strain_type)
		) |>
		# Now join the antigenic distances
		dplyr::left_join(
			nested_agdist_data,
			by = c("vaccine_name" = "Strain1", "strain_name" = "Strain2"),
			relationship = "many-to-one"
		) |>
		# Now unnest the antigenic distance data
		tidyr::unnest(antigenic_distances) |>
		# Clean up the metric names
		clean_metrics(keep_set = "extended")

	return(dat_joined)
}

# Now we need to make a few additional changes before we can pass the
# data to our models.
create_model_data <- function(joined_data) {
	dat_model <-
		joined_data |>
		dplyr::mutate(
			# We need to make versions of the time and birth year variables that are
			# closer to being scale-free than the current versions -- models that
			# have those numbers in the thousands have worse conditioning problems that
			# can lead to numerical issues in an otherwise fine model. So we'll scale
			# the year variable by subtracting 2013, the first year, and we'll scale
			# the birth_year variable with minmax transformation.
			year_c = year - 2013,
			birth_year_c = minmax(birth_year),
			# We'll also minmax scale the age.
			age_c = minmax(age),
			# Finally we'll create indicator variables for the categorical covariates
			# that will go in the model
			# Indicator variable for sex, 0 is Male, 1 is female.
			sex_i = as.numeric(sex == "Female"),
			# Indicator variable race/ethnicity
			race_i = as.numeric(race != "White"),
			# Create a version of type_subtype that combines all the B strains
			type_group = factor(dplyr::if_else(
				startsWith(as.character(type_subtype), "B"),
				"B",
				as.character(type_subtype)
			)),
			# Short version of the season for plotting
			season_short = gsub(" - 20([0-9]{2})", "/\\1", as.character(season)) |>
				factor(ordered = TRUE)
		) |>
		# We also want to minmax scale the antigenic distances, but we want to
		# do it separately by subtype and method.
		dplyr::group_by(strain_type, metric) |>
		dplyr::mutate(d_norm = minmax(d)) |>
		dplyr::ungroup() |>
		# Drop the dose column since everyone is SD
		dplyr::select(-dose) |>
		# Keep only the distance metrics for the models
		dplyr::filter(
			metric %in% c(
				"Temporal",
				"Cartographic",
				"Tree",
				#"Damerau-Levenshtein",
				"p-Epitope",
				#"FLU Substitution",
				#"Temporal (asymmetric)",
				"Grantham",
				#"Hamming",
				"p-All-Epitope"
			)
		) |>
		dplyr::mutate(metric = forcats::fct_drop(metric))

	return(dat_model)
}

# Save the final model data to file
write_model_data_to_file <- function(model_data, filename_base) {
	filenames <- paste0(filename_base, c(".Rds", ".csv"))
	filenames[[1]] |>
		dirname() |>
		dir.create(showWarnings = FALSE, recursive = TRUE)

	readr::write_rds(
		model_data,
		file = filenames[[1]],
		compress = "xz"
	)
	readr::write_csv(
		model_data,
		file = filenames[[2]]
	)

	return(filenames)
}

# Create the CIVICs reporting data that we need for the data submission
create_civics_reporting_data <- function(prepped_cohort_data) {
	dat_used <- prepped_cohort_data |>
		dplyr::distinct(
			id, season, study, age, birth_year, gender, sex, bmi, dose, ethnicity,
			race_civics_standard, year, subject_id
		)

	# Get the minimum and maximum age in the study
	dat_ages <- dat_used |>
		dplyr::group_by(id, study) |>
		dplyr::mutate(
			Min_Age = min(age),
			Max_Age = max(age)
		) |>
		dplyr::ungroup()

	dat_reporting <-
		dat_used |>
		dplyr::transmute(
			Study_Code = "TBD",
			Subject_ID = subject_id,
			Cohort_ID = paste(study, year, dose, sep = "_"),
			Sex_Assigned_at_Birth = ifelse(
				is.na(sex), "Unknown", as.character(sex)
			),
			Gender = gender,
			Min_Age = dat_ages$Min_Age,
			Max_Age = dat_ages$Max_Age,
			Subject_Age_Unit = "Years",
			Birth_Year = birth_year,
			Subject_Age_Event = "Age at enrollment",
			Subject_Phenotype = "Not Collected",
			Subject_Location = ifelse(
				study == "UGA", "GA", study
			),
			Ethnicity = ifelse(
				is.na(ethnicity),
				"Unknown",
				ethnicity
			),
			Race = race_civics_standard,
			Subject_Description = "Not Provided"
		) |>
		dplyr::mutate(dplyr::across(dplyr::where(is.character), as.factor))

	return(dat_reporting)
}

# END OF FILE ==================================================================
