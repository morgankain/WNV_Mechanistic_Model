WNV_Mechanistic_Model

This repository is designed to be used with Kain and Bolker 2019: Predicting West Nile virus transmission in North American bird communities using phylogenetic mixed effects models and eBird citizen science data. (https://parasitesandvectors.biomedcentral.com/articles/10.1186/s13071-019-3656-8) <br />

This repository provides a copy of the code that is also available in the online supplementary material of Kain and Bolker 2019. All of the raw data (apart from eBird data) needed to run these models is present in this repository. eBird data can be obtained here (https://science.ebird.org/en/use-ebird-data) <br />

This code can be run using two methods: <br />

1) Clean eBird .zip file and run R code "manually" in two separate steps <br />
		A) Run sh ebird_data_clean.sh in the command prompt once all eBird .zip files have been placed in the folder titled "ebird_zip_fresh" <br />
		B) Open "top_level_script.R" <br />
			i) Read the instructions at the top of this script <br />
			ii) Adjust parameters and options as desired <br />
			iii) Run the script as desired <br />
			
2) Clean eBird .zip file and run R code in a single step <br />
		A) Run sh ebird_bash_run.sh  <br />
			i) Will require some additional setup on the part of the user (see Step 4 in "ebird_bash_run.sh") <br />
			ii) Cant be run if there are missing scientific names (see "saved_matching" below) <br />
	
This repository contains a number of folders, many of which are empty. The folders are one of three types: <br />		

(1) Folders that need to have data placed in them prior to running the code. Data available as the supplementary material to Kain and Bolker 2019 at: URL HERE is organized in the appropriate folders <br />
		A) trees -- Contains phylogenetic tree data <br />
		B) ebird_zip_fresh -- Contains the .zip ebird file <br />
		C) data -- Contains all other data (bird responses, county data etc.) <br />

(2) Folders that contain model components <br />
		A) stan -- Contains stan model definitions <br />

(3) Folders that start empty but get filled as part of the automated workflow when output is saved to disk <br />
        		A) ebird_data_for_R -- Will contain all of the eBird data once it is extracted <br />
        		B) ebird_zip_dump -- Will contain .zip files after they get extracted <br />
		C) ebird_pieces -- One of two intermediate folders to momentarily house extracted eBird pieces <br />
		D) ebird_unzip -- Two of two intermediate folders to momentarily house extracted eBird pieces <br />
		E) saved_fits -- Contains intermediate fitted model results to expedite code in future runs <br />
		F) saved_matching -- Contains matched scientific names and bird body sizes <br />
		G) saved_output -- Contains all other fits including final "product" <br />

Details for the purpose of each R script are in "top_level_script.R"
	
