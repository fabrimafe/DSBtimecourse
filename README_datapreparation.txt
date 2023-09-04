cp ~/workspace/daniela/csv/Psy1_Types.csv ~/workspace/daniela/csv/Psy1_june2023.csv
cat ~/workspace/daniela/csv/Psy1_june2023.csv | sed 's/highres_//g' | sed 's/control/0h_/' | sed 's/$/,0/' > temp
#then edit temp to add +0.5 on time and place controls as 0 (first)
mv temp ~/workspace/daniela/csv/Psy1_june2023.csv

#error because no MH_del. Warn Daniela that problem with the formatting. Also ask about time 0
Rscript ./csv2timecourse.R -i ~/workspace/daniela/csv/Psy1_june2023.csv -o Psy1.june2023.txt

cat Psy1.june2023.txt | awk '{if (NR==1 || ( NR >=6 && NR<=11) ){print}}' > Psy1.june2023_control.txt  
cat Psy1.june2023.txt | awk '{if (NR<6 || NR>11){print}}' > ~/workspace/daniela/input_datasets/timecourse_RNP_Psy1.june2023.txt
cat Psy1.june2023.txt | awk '{if (NR==1 || ( NR>1 && $1<25)){print}}' > ~/workspace/daniela/input_datasets/timecourse_RNP_Psy1.june2023.24h.txt

IGENE=Psy1.june2023
./calculate_error_matrix.R -i Psy1.june2023_control.txt -o ~/workspace/daniela/error_matrices/error_matrix4_${IGENE}_ -f 1


#create inputs for bootstraps
module load R/4.1.0
mkdir -p ~/workspace/daniela/input_datasets/stationarybootstraps
mkdir -p ~/workspace/daniela/input_datasets/MEbootstraps
mkdir -p ~/workspace/daniela/input_datasets/stratifiedbootstraps
#mkdir -p ~/workspace/daniela/input_datasets/stratifiedbootstraps2
myfiles=$( ls ~/workspace/daniela/input_datasets/timecourse_*Psy1x5.txt )
myfiles=$( ls ~/workspace/daniela/input_datasets/timecourse_RNP_CRTISO.cleanedandnov2021*.txt )
myfiles=$( ls ~/workspace/daniela/input_datasets/timecourse_*allb.txt )
myfiles=$( ls ~/workspace/daniela/input_datasets/timecourse_RNP_CRTISO.49and50bp_all.txt )
myfiles=$( ls ~/workspace/daniela/input_datasets/timecourse_RNP_*no72h.txt )
myfiles=$( ls /home/labs/alevy/fabrizio/workspace/daniela/input_datasets/timecourse_RNP_CRTISO.nov2021.49and50_mydelay0.25.txt )
myfiles=$( ls ~/workspace/daniela/input_datasets/timecourse_RNP_*9timepoint24h.txt )
myfiles=$( ls ~/workspace/daniela/input_datasets/timecourse_RNP_Psy1_real0.txt )
myfiles=$( ls ~/workspace/daniela/input_datasets/timecourse_RNP_*june2023*.txt )
myfiles=$( ls ~/workspace/daniela/input_datasets/timecourse_RNP_Psy1.24h.txt )
for i in $myfiles;do
echo $i
xname=$(basename $i .txt )
newpath=~/workspace/daniela/input_datasets/stratifiedbootstraps/${xname}
rm -r $newpath;mkdir -p $newpath
./timeseriesbootstraps.R -i $i -o ${newpath} -n 100 -m 2
done

