WNV_Mechanistic_Model

This repository is designed to be used with Kain and Bolker 2019: Predicting West Nile virus transmission in North American bird communities using phylogenetic mixed effects models and eBird citizen science data. <br />

This repository provides a copy of the code that is also available in the online supplementary material of Kain and Bolker 2019. All of the data needed to run these models can be obtained in the online supplementary material of Kain and Bolker 2019 at: URL HERE upon publication, or upon request (morganpkain@gmail.com) <br />

This code can be run using two methods:

	1) Clean eBird .zip file and run R code "manually" in two separate steps
		A) Run sh ebird_data_clean.sh in the command prompt once all eBird .zip files have been placed in the folder titled "ebird_zip_fresh"
		B) Open "top_level_script.R"
			i) Read the instructions at the top of this script
			ii) Adjust parameters and options as desired
			iii) Run the script as desired
	2) Clean eBird .zip file and run R code in a single step
		A) Run sh ebird_bash_run.sh 
			i) Will require some additional setup on the part of the user (see Step 4 in "ebird_bash_run.sh")
			ii) Cant be run if there are missing (see "saved_matching" below) 
			
	This repository contains a number of folders, many of which are empty. The folders are one of three types:
	1) Folders that need to have data placed in them prior to running the code. Data available as the supplementary material to Kain and Bolker 2019 at: URL HERE is organized in the appropriate folders
		A) trees -- Contains phylogenetic tree data
		B) ebird_zip_fresh -- Contains the .zip ebird file
		C) data -- Contains all other data (bird responses, county data etc.)
	2) Folders that contain model components
		A) stan -- Contains stan model definitions
        3) Folders that start empty but get filled as part of the automated workflow when output is saved to disk
        		A) ebird_data_for_R -- Will contain all of the eBird data once it is extracted
        		B) ebird_zip_dump -- Will contain .zip files after they get extracted
		C) ebird_pieces -- One of two intermediate folders to momentarily house extracted eBird pieces
		D) ebird_unzip -- Two of two intermediate folders to momentarily house extracted eBird pieces
		E) saved_fits -- Houses intermediate fitted model results to expedite code in future runs
		F) saved_matching -- Houses matched scientific names and bird body sizes
		G) saved_output -- Houses all other fits including final "product"

Details for the purpose of each R script are in "top_level_script.R"
	
