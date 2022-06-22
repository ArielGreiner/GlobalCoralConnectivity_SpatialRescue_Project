Download the Beyer et al 2018 scores from: https://espace.library.uq.edu.au/view/UQ:0928a6a.
To calculate the scores used in this study from those, use the code in 'rescoring_50reefsdata.R'.
I then took this data (which included a score for all of the Beyer et al 2018 grid cells) and calculated the average score for every Wood et al. 2014 grid cell, using ArcMap and the process described in the 'Melding 50 reefs scores into WoodConnMat' file. 
The revised Beyer et al 2018 scores for each Wood et al. 2014 grid cell are included in the "scaled_scores_latlong_reproj_SJ.csv" dataframe (alongside other things generated in the MakingMasterDataframe_AGComputerVersion.R file found in the ConvertingWoodetal2014Data folder), which can be read directly into R.
