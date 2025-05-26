# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Add a header or whatever here

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes)
library(crew)
library(crew.cluster)

suppressPackageStartupMessages({
	library(renv)
	library(readr)
	library(markdown)
	library(rmarkdown)
	library(rstudioapi)
	library(visNetwork)
	library(quarto)
	library(brms)
	library(rstan)
	library(cmdstanr)
})

hpc <- !is.na(Sys.getenv("SLURM_JOB_ID", unset = NA))

# Set up crew controllers
controller_hpc_small <- crew.cluster::crew_controller_slurm(
	name = "hpc_small",
	workers = 4,
	seconds_idle = 120,  # time until workers are shut down after idle
	options_cluster = crew.cluster::crew_options_slurm(
		script_lines = c(
			"#SBATCH --ntasks=1",
			"module load R/4.4.1-foss-2022b"
			#add additional lines to the SLURM job script as necessary here
		),
		log_output = "logs/crew_small_log_%A.out",
		log_error = "logs/crew_small_log_%A.err",
		memory_gigabytes_required = 32,
		cpus_per_task = 2, #total 20gb RAM
		time_minutes = 1200, # wall time for each worker
		partition = "batch",
		verbose = TRUE
	)
)

controller_hpc_ll <- crew.cluster::crew_controller_slurm(
	name = "hpc_ll",
	workers = 12,
	seconds_idle = 120,  # time until workers are shut down after idle
	options_cluster = crew.cluster::crew_options_slurm(
		script_lines = c(
			"#SBATCH --ntasks=1",
			"module load R/4.4.1-foss-2022b"
			#add additional lines to the SLURM job script as necessary here
		),
		log_output = "logs/crew_ll_log_%A.out",
		log_error = "logs/crew_ll_log_%A.err",
		memory_gigabytes_required = 128,
		cpus_per_task = 1,
		time_minutes = 120,
		partition = "batch",
		verbose = TRUE
	)
)

controller_hpc_vd <- crew.cluster::crew_controller_slurm(
	name = "hpc_vd",
	workers = 12,
	seconds_idle = 120,
	options_cluster = crew.cluster::crew_options_slurm(
		script_lines = c(
			"#SBATCH --ntasks=1",
			"module load R/4.4.1-foss-2022b"
			#add additional lines to the SLURM job script as necessary here
		),
		log_output = "logs/crew_vd_log_%A.out",
		log_error = "logs/crew_vd_log_%A.err",
		memory_gigabytes_required = 64,
		cpus_per_task = 1,
		time_minutes = 120,
		partition = "batch",
		verbose = TRUE
	)
)

controller_brms <- crew.cluster::crew_controller_slurm(
	name = "brms",
	workers = 12,
	seconds_idle = 300,  # time until workers are shut down after idle
	options_cluster = crew.cluster::crew_options_slurm(
		script_lines = c(
			"#SBATCH --ntasks=1",
			"module load R/4.4.1-foss-2022b"
			#add additional lines to the SLURM job script as necessary here
		),
		log_output = "logs/crew_large_log_%A.out",
		log_error = "logs/crew_large_log_%A.err",
		memory_gigabytes_required = 32,
		cpus_per_task = 36,
		time_minutes = 7200, # wall time for each worker
		partition = "batch",
		verbose = TRUE
	)
)

n_local_cores <- ifelse(
	isTRUE(hpc),
	NA,
	parallelly::availableCores(omit = 2, constraints = "connections")
)
controller_local <- crew::crew_controller_local(
	name = "local",
	workers = n_local_cores,
	options_local = crew::crew_options_local(log_directory = "logs")
	# options_metrics = crew_options_metrics(
	# 	path = "logs",
	# 	seconds_interval = 1
	# )
)

# Set target options:
tar_option_set(
	format = "rds",
	# Set up parallel options via crew
	controller = crew::crew_controller_group(
		controller_hpc_ll,
		controller_hpc_vd,
		controller_hpc_small,
		controller_brms,
		controller_local
	),
	resources = targets::tar_resources(
		# If we detect that we're on slurm, by default use the small hpc controller
		# otherwise use the local controller
		# The large hpc controller must be manually specified
		crew = tar_resources_crew(controller = ifelse(hpc, "hpc_small", "local"))
	)
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source()

# Targets to run
list(
	# This target tells us whether to do a test run (very bad model fits with
	# only 20 or so samples, for debugging), or an actual run
	# FALSE = actual run; TRUE = test run
	targets::tar_target(
		name = is_test_run,
		command = FALSE
	),
	# UGAFluVac Cohort data set
	# This data is the output of the cleaning process from the repo
	# https://github.com/ahgroup/UGAFluVac-data/, which includes the raw excel
	# data files. We pulled the cleaned data from that repo on 2024-01-29.
	tarchetypes::tar_file_read(
		name = UGAFluVac_raw_data,
		command = here::here("data", "raw", "UGAFluVac-cohort-data.Rds"),
		read = readr::read_rds(!!.x)
	),
	# Antigenic distance data set
	# Load the clean antigenic distance data which is exported from the
	# ahgroup/influenza-antigenic-distance-data repo.
	tarchetypes::tar_file_read(
		name = antigenic_distance_raw_data,
		command = here::here("data", "raw", "antigenic-distance-data.Rds"),
		read = readr::read_rds(!!.x)
	),
	# Targets for data cleaning and processing
	targets::tar_target(
		name = prepped_cohort_data,
		command = prep_cohort_data(UGAFluVac_raw_data)
	),
	targets::tar_target(
		name = cleaned_cohort_data,
		command = clean_cohort_data(prepped_cohort_data)
	),
	targets::tar_target(
		name = joined_data,
		command = join_antigenic_distance_to_cohort_data(
			cleaned_cohort_data,
			antigenic_distance_raw_data
		)
	),
	targets::tar_target(
		name = model_data,
		command = create_model_data(joined_data)
	),
	# Save datasets to file
	tarchetypes::tar_file(
		name = model_data_files,
		command = write_model_data_to_file(
			model_data,
			here::here("data", "processed", "model-data")
		)
	),
	targets::tar_target(
		name = civics_reporting_data,
		command = create_civics_reporting_data(prepped_cohort_data)
	),
	tarchetypes::tar_file(
		name = civics_reporting_data_files,
		command = write_model_data_to_file(
			civics_reporting_data,
			here::here("data", "processed", "reporting-data")
		)
	),
	# Descriptive analysis targets
	tar_target(
		name = distinct_measurements,
		command = get_distinct_measurements(model_data)
	),
	tar_target(
		name = distinct_personyears,
		command = get_distinct_personyears(model_data)
	),
	tar_target(
		name = distinct_new_participants,
		command = get_distinct_new_participants(model_data)
	),
	tar_file(
		name = counts_by_study_table,
		command = create_count_by_study_table(
			distinct_measurements, distinct_personyears, distinct_new_participants,
			here::here("results", "tables", "tbl-counts.Rds")
		)
	),
	tar_file(
		name = data_counts,
		command = write_data_counts(
			distinct_measurements, distinct_personyears, distinct_new_participants,
			here::here("results", "data-counts.Rds")
		)
	),
	tar_file(
		name = measurements_per_person_year,
		command = calculate_measurements_per_person_year(
			distinct_measurements,
			here::here("results", "avg-measurements-per-py.Rds")
		)
	),
	tar_file(
		name = demographics_table,
		command = make_demographics_table(
			model_data,
			here::here("results", "tables", "tbl-demographics.Rds")
		)
	),
	# Strain information tables
	tar_file(
		name = strain_panel_table,
		command = create_annual_panel_table(
			model_data,
			here::here("results", "tables", "tbl-strain-panel.Rds")
		)
	),
	tar_file(
		name = strain_names_table,
		command = create_strain_name_translation_table(
			model_data,
			here::here("results", "tables", "tbl-strain-names.Rds")
		)
	),
	tar_file(
		name = annual_vaccine_table,
		command = create_fluzone_vaccine_table(
			model_data,
			here::here("results", "tables", "tbl-vaccine-strains.Rds")
		)
	),
	# Correlation between distances
	tar_target(
		name = normalized_correlation_data,
		command = clean_normalized_correlation_data(model_data)
	),
	tarchetypes::tar_file(
		name = vaccine_normalized_correlation_pairplot,
		command = make_normalized_correlation_pairplot(
			normalized_correlation_data,
			here::here("results", "figures", "dist-corr-norm.png")
		)
	),
	# Correlation table
	tar_target(
		name = correlation_cis,
		command = calculate_correlations_with_cis(normalized_correlation_data)
	),
	tar_target(
		name = correlation_cis_table,
		command = create_cor_ci_table(correlation_cis)
	),
	tar_file(
		name = correlation_cis_table_file,
		command = save_file_default(
			correlation_cis_table,
			here::here("results", "tables", "norm-corr-ci.Rds")
		)
	),
	# ICC analysis
	tar_target(
		name = icc_data,
		command = create_icc_input_data(model_data)
	),
	tar_target(
		name = icc_data_nested,
		command = nest_icc_input_data(icc_data)
	),
	# We have to use dynamic branching here because we want to iterate over a
	# target computed previously, which won't usually be defined at
	# the time the pipeline is created
	tar_target(
		name = icc_results_per_strain_type,
		command = calculate_icc_for_subtype(icc_data_nested$data),
		pattern = map(icc_data_nested),
		iteration = "list"
	),
	tar_target(
		name = icc_summary,
		command = combine_icc_summaries(
			icc_data_nested$strain_type,
			icc_results_per_strain_type
		)
	),
	tar_file(
		name = icc_summary_table,
		command = create_icc_summary_table(
			icc_summary,
			here::here("results", "tables", "icc-tab.Rds"),
			column = ICC
		)
	),
	tar_file(
		name = ppd_summary_table,
		command = create_icc_summary_table(
			icc_summary,
			here::here("results", "tables", "icc-tab-alt.Rds"),
			column = `PPD Ratio`
		)
	),
	# Targets for the spread/evenness analysis of metrics
	tar_target(
		name = normalized_icc_data,
		command = normalize_icc_data(icc_data)
	),
	tar_target(
		name = metrics_pc_plot,
		command = create_metrics_pc_plot(
			normalized_icc_data,
			here::here("results", "figures", "metrics-pc-plot.png")
		)
	),
	tar_target(
		name = gap_sd_bootstraps,
		command = calculate_gap_sd_with_bootstrap(normalized_icc_data)
	),
	tar_target(
		name = gap_sd_plot,
		command = create_gap_sd_plot(
			gap_sd_bootstraps,
			here::here("results", "figures", "metrics-ds-plot.png")
		)
	),
	tar_file(
		name = combined_metrics_plot,
		command = create_combined_metrics_plot(
			metrics_pc_plot,
			gap_sd_plot,
			here::here("results", "figures", "metrics-comb-plot.png")
		)
	),
	# Summary plots of pre and post titers
	tar_target(
		name = titer_plot_data,
		command = clean_data_for_titer_plots(model_data)
	),
	tar_file(
		name = pretiter_summary_plot,
		command = create_raw_data_plot(
			titer_plot_data,
			"pre",
			here::here("results", "figures", "pretiter-plot.png")
		)
	),
	tar_file(
		name = posttiter_summary_plot,
		command = create_raw_data_plot(
			titer_plot_data,
			"post",
			here::here("results", "figures", "posttiter-plot.png")
		)
	),
	# Model fitting prep
	tar_target(
		name = nested_model_data,
		command = nest_model_data_for_fitting(model_data)
	),
	tar_file_read(
		name = model_fitting_seeds,
		command = here::here("data", "raw", "model-fitting-seeds.txt"),
		read = scan(!!.x, quiet = TRUE)
	),
	tar_target(
		name = model_metadata,
		command = define_model_metadata(nested_model_data, model_fitting_seeds)
	),
	tar_file(
		name = model_metadata_file,
		command = save_file_default(
			model_metadata,
			here::here("results", "output", "model-metadata.Rds")
		)
	),
	# Modeling Fitting
	## Sample from the posterior (this one takes the longest by far)
	tar_target(
		name = brms_posterior_samples,
		command = fit_brms_model(
			data = model_metadata$model_data,
			model = model_metadata$model_forms,
			sample_priors = "no",
			# Add another random number from random dot org to get new random seeds
			seed = model_metadata$fitting_seed,
			test = is_test_run
		),
		pattern = map(model_metadata),
		iteration = "list",
		resources = tar_resources(
			crew = tar_resources_crew(controller = ifelse(hpc, "brms", "local"))
		)
	),
	# Now we save the fitted models to files
	tar_target(
		name = brms_posterior_fit_file_names,
		command = generate_brms_model_filenames(model_metadata, "posterior")
	),
	tar_target(
		name = brms_posterior_fit_files,
		command = save_file_default(
			brms_posterior_samples,
			brms_posterior_fit_file_names
		),
		pattern = map(brms_posterior_samples, brms_posterior_fit_file_names)
	),
	## Sample from the priors
	tar_target(
		name = brms_prior_samples,
		command = fit_brms_model(
			data = model_metadata$model_data,
			model = model_metadata$model_forms,
			sample_priors = "only",
			seed = model_metadata$fitting_seed + 877424L,
			test = is_test_run
		),
		pattern = map(model_metadata),
		iteration = "list",
		resources = tar_resources(
			crew = tar_resources_crew(controller = ifelse(hpc, "brms", "local"))
		)
	),
	# And save prior samples to file
	tar_target(
		name = brms_prior_fit_file_names,
		command = generate_brms_model_filenames(model_metadata, "prior")
	),
	tar_target(
		name = brms_prior_fit_files,
		command = save_file_default(
			brms_prior_samples,
			brms_prior_fit_file_names
		),
		pattern = map(brms_prior_samples, brms_prior_fit_file_names)
	),
	# Model diagnostics
	## First we need to get the cleaned diagnostics for each fit
	tar_target(
		name = posterior_fit_diagnostics,
		command = get_model_diagnostics(brms_posterior_samples),
		pattern = map(brms_posterior_samples),
		iteration = "list"
	),
	tar_target(
		name = prior_fit_diagnostics,
		command = get_model_diagnostics(brms_prior_samples),
		pattern = map(brms_prior_samples),
		iteration = "list"
	),
	## Now combine them into a nice table
	tar_target(
		name = posterior_fit_diagnostics_table,
		command = make_diagnostics_table(
			posterior_fit_diagnostics,
			model_metadata
		)
	),
	tar_target(
		name = prior_fit_diagnostics_table,
		command = make_diagnostics_table(
			prior_fit_diagnostics,
			model_metadata
		)
	),
	## And save them to files
	tar_file(
		name = posterior_fit_diagnostics_table_file,
		command = save_file_default(
			posterior_fit_diagnostics_table,
			here::here("results", "tables", "posterior-diagnostics.Rds")
		)
	),
	tar_file(
		name = prior_fit_diagnostics_table_file,
		command = save_file_default(
			prior_fit_diagnostics_table,
			here::here("results", "tables", "prior-diagnostics.Rds")
		)
	),
	# ELPD and WAIC calculation
	## First get the WAIC cause it's easy
	tar_target(
		name = model_waic,
		command = brms::waic(brms_posterior_samples),
		pattern = map(brms_posterior_samples),
		iteration = "list",
		resources = tar_resources(
			crew = tar_resources_crew(controller = ifelse(hpc, "hpc_ll", "local"))
		)
	),
	## Now find the problematic points across all models
	tar_target(
		name = ll_nan_values,
		command = find_ll_nan_values(model_waic),
		resources = tar_resources(
			crew = tar_resources_crew(controller = ifelse(hpc, "hpc_ll", "local"))
		)
	),
	## Correct the WAIC for NaN ll values
	tar_target(
		name = model_waic_corrected,
		command = correct_waic(model_waic, ll_nan_values),
		pattern = map(model_waic),
		iteration = "list",
		resources = tar_resources(
			crew = tar_resources_crew(controller = ifelse(hpc, "hpc_ll", "local"))
		)
	),
	## Next get the ELPD
	tar_target(
		name = model_log_likelihood,
		command = brms::log_lik(brms_posterior_samples),
		pattern = map(brms_posterior_samples),
		iteration = "list",
		resources = tar_resources(
			crew = tar_resources_crew(controller = ifelse(hpc, "hpc_ll", "local"))
		)
	),
	tar_target(
		name = model_elpd,
		command = calculate_elpd_from_ll(model_log_likelihood, ll_nan_values),
		pattern = map(model_log_likelihood),
		iteration = "list",
		resources = tar_resources(
			crew = tar_resources_crew(controller = ifelse(hpc, "hpc_ll", "local"))
		)
	),
	## Now make the comparisons
	tar_target(
		name = model_waic_comparisons,
		command = make_ic_comparisons(
			model_waic_corrected, model_metadata,
			ic = "waic"
		),
		resources = tar_resources(
			crew = tar_resources_crew(controller = ifelse(hpc, "hpc_ll", "local"))
		)
	),
	tar_target(
		name = model_elpd_comparisons,
		command = make_ic_comparisons(model_elpd, model_metadata, ic = "loo"),
		resources = tar_resources(
			crew = tar_resources_crew(controller = ifelse(hpc, "hpc_ll", "local"))
		)
	),
	## Now make nicely formatted tables for the elpd
	tar_target(
		name = formatted_elpd_comparisons,
		command = purrr::map(model_elpd_comparisons, format_elpd)
	),
	tar_target(
		name = elpd_comparison_table,
		command = make_elpd_comparison_table(formatted_elpd_comparisons)
	),
	tar_file(
		name = elpd_comparison_table_file,
		command = save_file_default(
			elpd_comparison_table,
			here::here("results", "tables", "loo-tbl-diff.Rds")
		)
	),
	## Now making the IC tables with diagnostics and all models
	tar_target(
		name = supp_elpd_data,
		command = clean_supp_elpd_data(formatted_elpd_comparisons)
	),
	tar_target(
		name = all_model_elpd_table,
		command = make_all_model_elpd_table(supp_elpd_data)
	),
	tar_target(
		name = elpd_diagnostic_table,
		command = make_elpd_diagnostic_table(supp_elpd_data)
	),
	tar_file(
		name = all_model_elpd_table_file,
		command = save_file_default(
			all_model_elpd_table,
			here::here("results", "tables", "model-elpd-table.Rds")
		)
	),
	tar_file(
		name = elpd_diagnostic_table_file,
		command = save_file_default(
			elpd_diagnostic_table,
			here::here("results", "tables", "elpd-diagnostic-table.Rds")
		)
	),
	# Model predictions
	## Get the prediction grid for each model
	tar_target(
		name = prediction_grids,
		command = generate_mem_grid(brms_posterior_samples),
		pattern = map(brms_posterior_samples),
		iteration = "list"
	),
	## Get the epred summaries all in one step
	tar_target(
		name = epred_draws,
		command = get_predictions(
			brms_posterior_samples,
			prediction_grids
		),
		pattern = map(brms_posterior_samples, prediction_grids),
		iteration = "list"
	),
	## Make the main predictions plot
	tar_file(
		name = model_prediction_plot,
		command = create_model_prediction_plot(
			epred_draws,
			model_metadata,
			here::here("results", "figures", "model-preds-gam-dose.png")
		)
	),
	## Calculate the contrast predictions and make the plot for both
	## model types
	tar_target(
		name = contrast_preds_lmm,
		command = calculate_contrast_preds(
			model_metadata,
			"lmm",
			brms_posterior_samples,
			prediction_grids
		)
	),
	tar_target(
		name = contrast_preds_gamm,
		command = calculate_contrast_preds(
			model_metadata,
			"gamm",
			brms_posterior_samples,
			prediction_grids
		)
	),
	tar_file(
		name = contrast_preds_lmm_plot,
		command = make_contrast_preds_plot(
			contrast_preds_lmm,
			here::here("results", "figures", "contrast-preds-plot-lmm.png")
		)
	),
	tar_file(
		name = contrast_preds_gamm_plot,
		command = make_contrast_preds_plot(
			contrast_preds_gamm,
			here::here("results", "figures", "contrast-preds-plot-gamm.png")
		)
	),
	## Plotting slopes and intercepts for LMM
	tar_file(
		name = lmm_parameters_plot,
		command = make_lmm_parameters_plot(
			model_metadata,
			brms_posterior_samples,
			here::here("results", "figures", "lmm-parm-plot.png")
		)
	),
	## Vaccine-specific predictions sensitivity analysis
	tar_file(
		name = vaccine_specific_predictions,
		command = create_vaccine_specific_prediction_plot(
			model_metadata,
			brms_posterior_samples,
			here::here("results", "figures", "vaccine-specific-preds-plot.png")
		)
	),
	## Fixed effects summaries
	tar_target(
		name = fixed_effect_summaries,
		command = summarize_fixed_effects(
			model_metadata,
			brms_posterior_samples
		)
	),
	tar_target(
		name = fixed_effect_summaries_table,
		command = create_fixed_effects_summary_table(fixed_effect_summaries)
	),
	tar_file(
		name = fixed_effect_summary_table_file,
		command = save_file_default(
			fixed_effect_summaries_table,
			here::here("results", "tables", "fixed-effects-coefs.Rds")
		)
	),
	## Variance decomposition
	tar_target(
		name = variance_decompositions_list,
		command = decompose_model_variance(
			brms_posterior_samples,
			model_metadata$model_forms
		),
		resources = tar_resources(
			crew = tar_resources_crew(controller = ifelse(hpc, "hpc_vd", "local"))
		),
		pattern = map(brms_posterior_samples, model_metadata),
		iteration = "list"
	),
	tar_target(
		name = variance_decompositions,
		command = bind_variance_decompositions(
			model_metadata,
			variance_decompositions_list
		)
	),
	tar_target(
		name = variance_decomposition_table,
		command = create_variance_ratio_table(variance_decompositions)
	),
	tar_file(
		name = variance_decomposition_table_file,
		command = save_file_default(
			variance_decomposition_table,
			here::here("results", "tables", "variance-ratio-table.Rds")
		)
	),
	# Supplementary analyses
	tar_file(
		name = dag_file,
		command = create_dag(here::here("results", "figures", "dag-model.png"))
	),
	tar_target(
		name = sequence_distance_correlations,
		command = analyze_sequence_distance_correlation(joined_data)
	),
	tar_target(
		name = sequence_correlation_table,
		command = create_cor_ci_table(
			sequence_distance_correlations,
			c("Hamming", "p-Epitope", "p-All-Epitope")
		)
	),
	tar_target(
		name = substitution_correlation_table,
		command = create_cor_ci_table(
			sequence_distance_correlations,
			c("Hamming", "Grantham", "FLU Substitution")
		)
	),
	tar_file(
		name = sequence_correlation_table_file,
		command = save_file_default(
			sequence_correlation_table,
			here::here("results", "tables", "sequence-dist-correlations.Rds")
		)
	),
	tar_file(
		name = substitution_correlation_table_file,
		command = save_file_default(
			substitution_correlation_table,
			here::here("results", "tables", "substitution-dist-correlations.Rds")
		)
	),
	# CIVICs metadata submission generate
	## Globals and setup
	tar_target(
		name = civics_metadata_globals,
		command = set_common_metadata_variables()
	),
	## Writing metadata sheets
	## This entire section could be rewritten with static branching
	## but i'm too lazy
	tar_file(
		name = civics_metadata_sheet_protocol,
		command = generate_metadata_sheet_protocol(
			civics_metadata_globals,
			here::here("assets", "CIVICs-metadata", "Protocol.csv")
		)
	),
	tar_file(
		name = civics_metadata_sheet_study_design,
		generate_metadata_sheet_study_design(
			civics_reporting_data,
			civics_metadata_globals,
			here::here("assets", "CIVICs-metadata", "Study_Design.csv")
		)
	),
	tar_file(
		name = civics_metadata_sheet_study_personnel,
		command = generate_metadata_sheet_personnel(
			civics_metadata_globals,
			here::here("assets", "CIVICs-metadata", "Study_Personnel.csv")
		)
	),
	tar_file(
		name = civics_metadata_sheet_study_inclusion,
		command = generate_metadata_sheet_study_inclusion(
			civics_metadata_globals,
			here::here("assets", "CIVICs-metadata", "Study_Inclusion.csv")
		)
	),
	tar_target(
		name = metadata_study_arms_dataset,
		command = format_study_arms_data(civics_metadata_globals)
	),
	tar_file(
		name = civics_metadata_sheet_study_arms,
		command = generate_metadata_sheet_study_arms(
			metadata_study_arms_dataset,
			civics_metadata_globals,
			here::here("assets", "CIVICs-metadata", "Study_Arms.csv")
		)
	),
	tar_file(
		name = civics_metadata_sheet_planned_visits,
		command = generate_metadata_sheet_planned_visits(
			civics_metadata_globals,
			here::here("assets", "CIVICs-metadata", "Study_Planned_Visits.csv")
		)
	),
	tar_file(
		name = civics_metadata_sheet_interventions,
		command = generate_metadata_sheet_interventions(
			metadata_study_arms_dataset,
			civics_metadata_globals,
			here::here("assets", "CIVICs-metadata", "Interventions.csv")
		)
	),
	tar_file(
		name = civics_metadata_sheet_human_subjects_demographics,
		command = generate_metadata_sheet_human_subjects_demographics(
			civics_reporting_data,
			civics_metadata_globals,
			here::here("assets", "CIVICs-metadata", "Human_Subjects_Demographics.csv")
		)
	),
	tar_file(
		name = civics_metadata_sheet_computational_model,
		command = generate_metadata_sheet_computational_model(
			civics_metadata_globals,
			here::here("assets", "CIVICs-metadata", "Computational_Model.csv")
		)
	),
	## Copying figures with approved names
	tar_file(
		name = civics_metadata_manuscript_figure_files,
		command = move_figure_files_to_metadata_location(
			prefix = "",
			base_dir = here::here("assets", "CIVICs-metadata"),
			c(
				vaccine_normalized_correlation_pairplot,
				model_prediction_plot,
				lmm_parameters_plot,
				vaccine_specific_predictions
			)
		)
	),
	tar_file(
		name = civics_metadata_supplement_figure_files,
		command = move_figure_files_to_metadata_location(
			prefix = "S",
			base_dir = here::here("assets", "CIVICs-metadata"),
			c(
				dag_file,
				pretiter_summary_plot,
				posttiter_summary_plot,
				combined_metrics_plot,
				contrast_preds_lmm_plot,
				contrast_preds_gamm_plot
			)
		)
	),
	## Copying/writing tables with approved format and names
	tar_file(
		name = civics_metadata_manuscript_table_files,
		command = move_table_data_to_metadata_location(
			prefix = "",
			base_dir = here::here("assets", "CIVICs-metadata"),
			c(
			counts_by_study_table,
			icc_summary_table,
			elpd_comparison_table_file,
			fixed_effect_summary_table_file,
			variance_decomposition_table_file
			)
		)
	),
	tar_file(
		name = civics_metadata_supplement_table_files,
		command = move_table_data_to_metadata_location(
			prefix = "S",
			base_dir = here::here("assets", "CIVICs-metadata"),
			c(
				annual_vaccine_table,
				strain_panel_table,
				strain_names_table,
				demographics_table,
				ppd_summary_table,
				correlation_cis_table_file,
				posterior_fit_diagnostics_table_file,
				prior_fit_diagnostics_table_file,
				elpd_diagnostic_table_file,
				substitution_correlation_table_file,
				sequence_correlation_table_file
			)
		)
	),
	## Zipping everything together
	tar_file(
		name = zipped_metadata_submission,
		command = create_zipped_metadata_submission(
			# Metadata sheets
			civics_metadata_sheet_protocol,
			civics_metadata_sheet_study_design,
			civics_metadata_sheet_study_personnel,
			civics_metadata_sheet_study_inclusion,
			civics_metadata_sheet_study_arms,
			civics_metadata_sheet_planned_visits,
			civics_metadata_sheet_interventions,
			civics_metadata_sheet_human_subjects_demographics,
			civics_metadata_sheet_computational_model,
			# Figures
			civics_metadata_manuscript_figure_files,
			civics_metadata_supplement_figure_files,
			# Tables
			civics_metadata_manuscript_table_files,
			civics_metadata_supplement_table_files,
			# Misc
			c(
				here::here("assets", "Protocol-Bayesian-Modeling.pdf")
			),
			# Location to save
			here::here("assets", "CIVICs-metadata-submission.zip")
		)
	),
	# Manuscript and supplement creation
	## First generate the software bibliography
	tar_file(
		name = software_bibliography,
		command = generate_software_bibliography(
			here::here("assets", "package-refs.bib")
		)
	),
	# Now render the manuscript
	tarchetypes::tar_quarto(
		name = manuscript,
		path = here::here("products", "manuscript", "Manuscript.qmd"),
		extra_files = c(
			here::here("assets", "refs.bib"),
			here::here("assets", "package-refs.bib"),
			here::here("assets", "vancouver.csl"),
			here::here("assets", "word-template.docx")
		),
		quiet = FALSE
	),
	## And the supplement
	tarchetypes::tar_quarto(
		name = supplement,
		path = here::here("products", "manuscript", "Supplement.qmd"),
		extra_files = c(
			here::here("assets", "refs.bib"),
			here::here("assets", "package-refs.bib"),
			here::here("assets", "vancouver.csl"),
			here::here("assets", "word-template.docx")
		),
		quiet = FALSE
	),
	## And the README
	tarchetypes::tar_quarto(
		name = README,
		path = here::here("README.qmd"),
		quiet = FALSE
	)
)
