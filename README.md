# A. bruennichi Data Analysis for Manuscript
This repository contains the R Project, scripts and output of analyses done in R pertaining to the publication, "Rapid ecological and evolutionary divergence during a poleward range expansion" by Sheffer et al. Due to the size of the files, the raw bioclim tif files are not included in this repository, but you can contact the corresponding author for that data. That data is only required for the 00_climate_data_wrangling.R script. For the main purposes of the analyses, you would not need to run the climate data wrangling script; the output is available and usable for all downstream analyses (used mainly in script 01_import_data.R and 01_maternalphenotype_latitude.R)

Scripts should be run in numerical order as named, with the exception of the script "myplot.R" which is a function for plotting that is called in the 00_climate_data_wrangling.R script. 

It also contains a python code snippet used for checking whether SNPs fall inside or outside of exons in the genome (inexon.py), and a csv file (Abruennichi_FASTAfilenames.csv) that can be used to match FASTA files with relevant metadata. 
