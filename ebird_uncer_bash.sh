### Run over all sources of uncertainty, restarting R after each run to clear RAM

for i in `seq 1 24`

do

	R CMD BATCH top_level_script.R

done




