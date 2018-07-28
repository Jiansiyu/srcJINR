#!/bin/bash
RUNLIST="runList.txt"

while read run_num
do
	FILE="/nobackup1/segarrae/data/dubna/raw/mpd_run_trigCode_"${run_num}".data.gz"
	echo "Working on file ${FILE}..."
	if [ -f $FILE ]
	then
		echo "Unpacking File..."
		gunzip $FILE
		FILE="/nobackup1/segarrae/data/dubna/raw/mpd_run_trigCode_"${run_num}".data"
		echo "Decoding file..."
		root -l -b -q "BmnDataToRoot.C(\"$FILE\")"
		rm bmn_run${run_num}_raw.root
		mv bmn_run${run_num}_digi.root /nobackup1/segarrae/data/dubna/root/H2-Target/
		echo "Re-packing file..."
		gzip --fast $FILE
	elif [ -f "/nobackup1/segarrae/data/dubna/raw/mpd_run_trigCode_"${run_num}".data" ]
	then
		echo "File already unpacked..."	
		FILE="/nobackup1/segarrae/data/dubna/raw/mpd_run_trigCode_"${run_num}".data"
		echo "Decoding file..."
		root -l -b -q "BmnDataToRoot.C(\"$FILE\")"
		rm bmn_run${run_num}_raw.root
		mv bmn_run${run_num}_digi.root /nobackup1/segarrae/data/dubna/root/H2-Target/
		echo "Re-packing file..."
		gzip --fast $FILE
	else
		echo "No file exists..."
	fi
	echo "Done with file ${FILE}"
	echo ""
done < $RUNLIST

