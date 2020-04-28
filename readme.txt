The analyses of Roslin et al. manuscript "Large-scale variation in the strength and direction of phenological change"
can be replicated as follows:

1) Download the data from "Chronicles of nature calendar, a long-term and large-scale multitaxon database on phenology"Ovaskainen et al. 2020 (Nature Scientific Data 7, 47)

2) Place these data into the folder "chronicle-of-nature-calendar.v1.0.5" under your working directory. Create also folders "data" and "models" under the working directory.

3) Run the R-script S0_prepare_data.R. This script selects the data following the description of the methods section in Roslin et al. manuscript, and places the selected data to the folder "data"

4) Run the R-script S1_fit_model.R. This script fits the HMSC-models, and places the fitted models to the "models" folder.

5) Run the R-script S2_parameter.estimates.R. This scripts produces the results presented in the manuscript. 