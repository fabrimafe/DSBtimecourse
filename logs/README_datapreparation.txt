cp ~/workspace/daniela/csv/Psy1_Types.csv ~/workspace/daniela/csv/Psy1_june2023.csv
cat ~/workspace/daniela/csv/Psy1_june2023.csv | sed 's/highres_//g' | sed 's/control/0h_/' | sed 's/$/,0/' > temp
#then edit temp to add +0.5 on time and place controls as 0 (first)
mv temp ~/workspace/daniela/csv/Psy1_june2023.csv

#error because no MH_del. Warn Daniela that problem with the formatting. Also ask about time 0
Rscript ./csv2timecourse.R -i ~/workspace/daniela/csv/Psy1_june2023.csv -o ~/workspace/daniela/input_dataset/Psy1.june2023.txt
cat Psy1.june2023.txt | awk '{if (NR==1 || ( NR >=6 && NR<=11) ){print}}' > Psy1.june2023_control.txt  
cat Psy1.june2023.txt | awk '{if (NR<6 || NR>11){print}}' > ~/workspace/daniela/input_datasets/timecourse_RNP_Psy1.june2023.txt
cat Psy1.june2023.txt | awk '{if (NR==1 || ( NR>1 && $1<25)){print}}' > ~/workspace/daniela/input_datasets/timecourse_RNP_Psy1.june2023.24h.txt
IGENE=Psy1.june2023
./calculate_error_matrix.R -i Psy1.june2023_control.txt -o ~/workspace/daniela/error_matrices/error_matrix4_${IGENE}_ -f 1

Rscript ./csv2timecourse.R -i ~/workspace/daniela/csv/Psy1.20230903_Types_MH_df.csv -o ~/workspace/daniela/input_datasets/timecourse_RNP_Psy1.20230903.txt
Rscript ./csv2timecourse.R -i ~/workspace/daniela/csv/CRTISO.20230903_RNP_Types_MH_df.csv -o ~/workspace/daniela/input_datasets/timecourse_RNP_CRTISO.20230903.txt
Rscript ./csv2timecourse.R -i ~/workspace/daniela/csv/CRTISO.49and50.20230903_Types_MH_df.csv -o ~/workspace/daniela/input_datasets/timecourse_RNP_CRTISO.49and50.20230903.txt
Rscript ./csv2timecourse.R -i ~/workspace/daniela/csv/PhyB2.9timepoint.20230903_RNP_Types_MH_df.csv -o ~/workspace/daniela/input_datasets/timecourse_RNP_PhyB2.20230903.txt
Rscript ./csv2timecourse.R -i ~/workspace/daniela/csv/Psy1.20231005_RNP_Types_MH_df.csv -o ~/workspace/daniela/input_datasets/timecourse_RNP_Psy1.20231005.txt

#REFRAME DATASET BATCH I AS RESUBMISSION
MYGENE=PhyB2.2.nov2021.24h;cp ~/workspace/daniela/input_datasets/timecourse_RNP_${MYGENE}.txt ~/workspace/daniela/input_datasets/timecourse_RNP_PhyB2.20231005.I.txt
Rscript ./csv2timecourse.R -i ~/workspace/daniela/csv/Psy1.20231005.I_RNP_Types_MH_df_timecourse.csv -o ~/workspace/daniela/input_datasets/timecourse_RNP_Psy1.20231005.I.txt
Rscript ./csv2timecourse.R -i /home/labs/alevy/fabrizio/workspace/daniela/csv/CRTISO.nov2021_RNP_Types_MH_df.csv -o ~/workspace/daniela/input_datasets/timecourse_RNP_CRTISO.20231005.I72h.txt
Rscript ./csv2timecourse.R -i /home/labs/alevy/fabrizio/workspace/daniela/csv/CRTISO.nov2021.49and50_RNP_Types_MH_df.csv -o ~/workspace/daniela/input_datasets/timecourse_RNP_CRTISO.49and50.20231005.I72h.txt

for MYGENE in CRTISO.49and50.20231005.I72h CRTISO.20231005.I72h;do
cat ~/workspace/daniela/input_datasets/timecourse_RNP_${MYGENE}.txt | awk -v OFS='\t' '{if ($1=="time"){print} else {$1=$1+0.5;print}}' > temp; mv temp ~/workspace/daniela/input_datasets/timecourse_RNP_${MYGENE}.txt  done
for MYGENE in Psy1.20231005.I72h PhyB2.20231005.I72h;do
Rscript ./csv2timecourse.R -i ~/workspace/daniela/csv/${MYGENE}.csv -o ~/workspace/daniela/input_datasets/timecourse_RNP_${MYGENE}.txt
done
for MYGENE in CRTISO.20231005.I72h CRTISO.49and50.20231005.I72h;do
Rscript ./csv2timecourse.R -i ~/workspace/daniela/csv/${MYGENE}_control_Types_MH_df.csv -o ~/workspace/daniela/input_datasets/control_${MYGENE}.txt
Rscript ./calculate_error_matrix.R -i ~/workspace/daniela/input_datasets/control_${MYGENE}.txt -o ~/workspace/daniela/error_matrices/error_matrix4_${MYGENE}_ -f 1
done


for MYGENE in Psy1.20231005.I72h PhyB2.20231005.I72h;do #Psy1.20231005.I;do #PhyB2.20231005.I Psy1.20231005.I;do # Psy1.20231005;do #PhyB2.2.nov2021.24h;do #CRTISO.49and50.20230903 Psy1.20230903 CRTISO.20230903 PhyB2.20230903;do   
cat ~/workspace/daniela/input_datasets/timecourse_RNP_${MYGENE}.txt | awk '{if (NR==1 || $1==0){print}}' > ~/workspace/daniela/input_datasets/${MYGENE}_control.txt
./calculate_error_matrix.R -i ~/workspace/daniela/input_datasets/${MYGENE}_control.txt -o ~/workspace/daniela/error_matrices/error_matrix4_${MYGENE}_ -f 1
cat ~/workspace/daniela/input_datasets/timecourse_RNP_${MYGENE}.txt | awk '{if ($1!=0){print}}' > temp.txt; mv temp.txt  ~/workspace/daniela/input_datasets/timecourse_RNP_${MYGENE}.txt 
done




#PhyB2 72h
cat ~/workspace/daniela/input_datasets/timecourse_RNP_PhyB2.2.nov2021.txt | awk -v OFS='\t' '{if ($1!="time"){$1=$1+0.5};print}' > ~/workspace/daniela/input_datasets/timecourse_RNP_PhyB2.20231005.I72h.txt
#cp ~/workspace/daniela/csv/control_PhyB2.2.nov2021_Types_MH_df.csv ~/workspace/daniela/csv/control_PhyB2.20231005.I72h.csv
#Psy1 72h
cat ~/workspace/daniela/input_datasets/timecourse_RNP_Psy1.txt | awk -v OFS='\t' '{if ($1!="time"){$1=$1+0.5};print}' > ~/workspace/daniela/input_datasets/timecourse_RNP_Psy1.20231005.I72h.txt
cp ~/workspace/daniela/csv/control_Psy1_Types_MH_df.csv ~/workspace/daniela/csv/control_Psy1.20231005.I72h.csv

for MYGENE in Psy1.20231005.I72h PhyB2.20231005.I72h;do 
./calculate_error_matrix.R -i ~/workspace/daniela/csv/control_${MYGENE}.csv -o ~/workspace/daniela/error_matrices/error_matrix4_${MYGENE}_ -f 1
done

cp ~/workspace/daniela/csv/PhyB2.9timepoint.20230909_RNP_Types_MH_df.csv ~/workspace/daniela/csv/PhyB2.20230909_RNP_Types_MH_df.csv
for MYGENE in Psy1 CRTISO CRTISO.49and50 PhyB2;do
Rscript ./csv2timecourse.R -i ~/workspace/daniela/csv/${MYGENE}.20230909_RNP_Types_MH_df.csv -o ~/workspace/daniela/input_datasets/timecourse_RNP_${MYGENE}.20230909.txt
cp ~/workspace/daniela/error_matrices/error_matrix4_${MYGENE}.20230903_errorsfromunbroken.tsv ~/workspace/daniela/error_matrices/error_matrix4_${MYGENE}.20230909_errorsfromunbroken.tsv
done


#PARSE TIME COURSES FOR NO DELAY
for MYGENE in Psy1.20231005 CRTISO.20230903 CRTISO.49and50.20230903 PhyB2.20230903;do
MYGENEOUT=$( echo $MYGENE | sed 's/20231005/20231026/' | sed 's/20230903/20231026/')
cat ~/workspace/daniela/input_datasets/timecourse_RNP_${MYGENE}.txt | awk -v OFS='\t' '{if ($1=="time"){print} else {$1=$1-0.5;print}}' > ~/workspace/daniela/input_datasets/timecourse_RNP_${MYGENEOUT}.txt
done



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
myfiles=$( ls ~/workspace/daniela/input_datasets/timecourse_RNP_*20230903.txt )
myfiles=$( ls ~/workspace/daniela/input_datasets/timecourse_RNP_PhyB2.2.nov2021.24h.txt )
myfiles=$( ls ~/workspace/daniela/input_datasets/timecourse_RNP_*20230909.txt )
myfiles=$( ls ~/workspace/daniela/input_datasets/timecourse_RNP_*20231005.txt )
myfiles=$( ls ~/workspace/daniela/input_datasets/timecourse_RNP_PhyB2.20231005.I.txt )
myfiles=$( ls ~/workspace/daniela/input_datasets/timecourse_RNP_Psy1.20231005.I.txt )
myfiles=$( ls ~/workspace/daniela/input_datasets/timecourse_RNP_Psy1.20231005.I72h.txt ~/workspace/daniela/input_datasets/timecourse_RNP_PhyB2.20231005.I72h.txt )
myfiles=$( ls  ~/workspace/daniela/input_datasets/timecourse_RNP_CRTISO.20231005.I72h.txt ~/workspace/daniela/input_datasets/timecourse_RNP_CRTISO.49and50.20231005.I72h.txt )
for i in $myfiles;do
echo $i
xname=$(basename $i .txt )
newpath=~/workspace/daniela/input_datasets/stratifiedbootstraps/${xname}
rm -r $newpath;mkdir -p $newpath
./timeseriesbootstraps.R -i $i -o ${newpath} -n 100 -m 2
done

