#### Step 1: Unzip all downloaded ebird data

## List of all of the files in the ebird zip folder
FILES=ebird_zip_fresh/*

## Initiate counter
c=1

## Unzip and subset each downloaded ebird file
for f in $FILES

do
    echo "Processing $f file"

## Unzip only the data file. It always has the prefix ebd (ebird data)
    7z x "$f" -r ebd* -oebird_unzip

## move the zip, enter and edit the data file in the resulting folder
    mv $f ebird_zip_dump/


#### Step 2: Clean/strip each data file

    echo "Processing file $c"

## rename the data file (there should only ever be one file in this folder at a time, so this should be ok)
    mv ebird_unzip/*.txt ebird_unzip/ebird_data.tsv

## not the best way to integrate perl splitting script, but ok...
    cd ebird_unzip
    make

## remove the large, unfiltered ebird data file
    rm ebird_data.tsv

## back up a folder
    cd ..

## Split the intermediate file. Store the desired columns as select.csv
    split -l 500000 ebird_unzip/ebird_stripped.tsv ebird_pieces/ebird_file_

## clean up, so naming of the sole file in the folder isn't screwed up the next time through the loop
    rm ebird_unzip/ebird_stripped.tsv
    rm ebird_unzip/ebird_data.tsv

#### Step 3: Rename and move each 500,000 line segment

## new counter for each file piece
cc=1
FILES=ebird_pieces/*

## Unzip and subset each downloaded ebird file
for ff in $FILES

do

    echo "Processing $ff file"

## rename number ebird_data_X
    mv $ff ebird_data_for_R/"ebird_data_$cc.txt"

## add to inner loop counter
    cc=$((cc+1))

done

#### Step 4: Run the R scripts to calculate community competence on these data

## Run R, loading the csv files contained within ebird_data_for_R/"ebird_data_$cc.txt"
#    Rscript top_level_script.R

## remove the csv pieces when done
#    rm ebird_data_for_R/*

## Rename the output from R to correspond to the counter
#    mv model_out/model_out.rds model_out/"model_out_$c.rds"

## Add to outer loop coutner
    c=$((c+1))

done





