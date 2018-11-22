FinalReport="FinalReport_no_header.txt"
SNPTable="SNPTable_no_header.txt"
FinalReportHeader="FinalReportHeader.txt"

echo "Step 4. Generating tped file"
echo "Step 4.1 Generating genotype1.txt 'cut -f2- ${FinalReport}'"
cut -f2- ${FinalReport} >genotype1.txt
echo "Step 4.2 Replacing '-' with '0' and saving in genotype2.txt"
sed 's/-/0/g' genotype1.txt |sed 's/./& /g' >genotype2.txt
echo "Step 4.3 Replacing tab spaces with no spaces and saving into genotype3.txt"
sed 's/\t//g' genotype2.txt >genotype3.txt
echo "Combining Step3's mapfile with genotype3.txt and saving as genodata.tped"
paste -d' ' Step3_mapfile genotype3.txt >genodata.tped
echo "Created genodata.tped file"
echo " "
#echo "Remove intermediary (genotype[1-3].txt) files?"
#read cont
#if [[ $cont == "y" ]] ; then
#  rm genotype[1-3].txt
#else
#  echo "Did not remove intermediary files"
#fi
echo "Done; may want to remove intermediary genotype[1-3].txt files"
