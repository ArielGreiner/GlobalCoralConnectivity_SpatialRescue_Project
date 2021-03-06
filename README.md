# GlobalCoralConnectivity_SpatialRescue_Project
Code used to generate results for Greiner et al. 2022 Global Ecology and Biogeography. All code was generated by me (Ariel Greiner) unless stated explicitly in the file itself. All outputs were generated between 2017-2021. Please contact me if you have any issues running the code. R was used to generate and manipulate the majority of the data used in this project. The revised Beyer et al 2018 scores included in DATAFRAME were generated using ArcMap (v. 10.4.1, ESRI 2016) from **DATAFRAME**, as ArcMap enabled us to calculate the average score for each of Wood et al. 2014's grid cells. The results from the random replicates were generated in a server farm (compute canada, https://www.computecanada.ca/, server farm accessible to Canadian researchers) and BASH code was used to run the R code in those instances.

HOW TO USE:
- The data used to generate the results in this paper may be found: The original worldwide reef scores determined in Beyer et al. (2018) are available at https://espace.library.uq.edu.au/view/UQ:0928a6a and https://github.com/WCS-Marine/local-reef-pressures/tree/main/data-raw/50-reefs (from Andrello et al. 2021). The original worldwide coral reef connectivity matrix generated for Wood et al. (2014) and the grid cell coordinates are available at the University of Bristol data repository, data.bris, at https://doi.org/10.5523/bris.2s0fn0bc89omq2kj2rol7iolwt.
-  If wish to start from the original data (see above), start from the code found in the 'OriginalDataConversion' folder, which will allow you to convert the data from its original form to the form used in the remainder of the code provided.
-  The data and plots for the Random Scenario were generated on compute canada because the data were too large to host on my personal computer and are also too large to be uploaded here, but the code for generating them on a similar server farm is included in the 'ComputeCanadaCode' folder.
- The code used to generate the results from the study can be found in the 'GeneratingFigures_Code' which is split up into code used to generate the main text results and the supplementary text results.

Please feel free to use and borrow from this code base according to the licence below.

Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
