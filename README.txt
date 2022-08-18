Scripts use the "here" package so users do not need to set paths. To use the here package, users must add the contents of this repository to an R project.

folder descriptions are as follows:

figures/: contains figures, for both the main text and SI, in .tiff and .png formats

model_fits/: contains a printoff of submodel diagnostics (submodel_diagnostics_printoff.docx), PSIS-LOOCV results (loos_D[1--3].rds), PSIS-LOOCV comparison tables (loo_table_D[1--3]].txt). 
	This folder normally contains model fits, but they are too large for GitHub. To obtain submodel fits, run scripts/write_submodels_D1.Rmd, scripts/write_submodels_D2.Rmd, and then scripts/fit_submodels.Rmd.
	To obtain full model fits, first get the submodel fits (see above), then run scripts/fit_full_model_D[1--3].Rmd.
	
models/: contains Stan models, both submodels and full models

scripts/: contains R and R markdown file for writing models, fitting them, creating figures, etc., see scripts/README.txt for individual script descriptions.

