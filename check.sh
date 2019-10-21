#!/bin/bash

path[0]="/data/smechbal/Fluka/NonUniB/V4"
path[1]="/home/smechbal/ANALYSIS_SOFT"
path[2]="/data/smechbal/Data/FlightData/18A1_SplitBPD"
path[3]="/home/smechbal/Documents/AESOPLITE/Analysis/Macros"

reco[0]="KFone"
reco[1]="KFtwo"
reco[2]="RKfit"

r=1

for ((i=12;i<=22;i++))
do
     echo " $i " >> ${path[2]}/log.txt
     if [ "$i" -lt 10 ]
     then
        echo " $i less than 10 " >> ${path[2]}/log.txt
	if [ -e ${path[2]}/18A1_00$i.BPD.EVENT_${reco[$r]}.root ] && [ -e ${path[2]}/18A1_00$i.BPD.txt ]
	then
    		echo " 18A1_00$i.BPD.EVENT_${reco[$r]}.root reconstructed "  >  ${path[2]}/Checked_18A1_00$i.BPD.txt
                echo " $i BPD.root exist " >> ${path[2]}/log.txt
                rm ${path[2]}/18A1_00$i.BPD.root
       elif [ -e ${path[2]}/18A1_00$i.BPD.root ] 
               then
               echo " $i BPD.root exist " >> ${path[2]}/log.txt
               break
       else
    		echo " 18A1_00$i.BPD.EVENT_${reco[$r]}.root, NOT reconstructed" >  ${path[2]}/Checked_18A1_00$i.BPD.txt
	        echo " $i BPD.EVENT.root does not exist " >> ${path[2]}/log.txt
         	cd ${path[1]}
                echo "In ${PWD}"
                echo "${path[2]}/18A1_00$i.BPD 1 2 " > ${path[1]}/Datafilepaths.txt
                source ./setup
                cd prod
                echo "In ${PWD}" >> ${path[2]}/Checked_18A1_00$i.BPD.txt
                echo "nohup ./MainRawBPDEvent ../Datafilepaths.txt 6 0 1 ${reco[$r]} & " >>  ${path[2]}/Checked_18A1_00$i.BPD.txt            
                nohup ${path[1]}/prod/MainRawBPDEvent ${path[1]}/Datafilepaths.txt 6 0 1 ${reco[$r]} >> ${path[2]}/Checked_18A1_00$i.BPD.txt
                break
	fi
    else
        if [ -e ${path[2]}/18A1_0$i.BPD.EVENT_${reco[$r]}.root ] && [ -e ${path[2]}/18A1_0$i.BPD.txt ]
        then
                echo " 18A1_0$i.BPD.EVENT_${reco[$r]}.root reconstructed "  >  ${path[2]}/Checked_18A1_0$i.BPD.txt
                rm ${path[2]}/18A1_0$i.BPD.root
       elif [ -e ${path[2]}/18A1_0$i.BPD.root ]
       then
               echo " $i BPD.root exist " >> ${path[2]}/log.txt
               break
       else
                echo " 18A1_0$i.BPD.EVENT_${reco[$r]}.root, NOT reconstructed" >  ${path[2]}/Checked_18A1_0$i.BPD.txt
                echo " $i BPD.EVENT.root does not exist " >> ${path[2]}/log.txt    
                cd ${path[1]}
                echo "In ${PWD}"
                echo "${path[2]}/18A1_0$i.BPD 1 2 " > ${path[1]}/Datafilepaths.txt
                source ./setup
                cd prod
                echo "In ${PWD}" >> ${path[2]}/Checked_18A1_0$i.BPD.txt
                echo "nohup ./MainRawBPDEvent ../Datafilepaths.txt 6 0 1 ${reco[$r]} & " >>  ${path[2]}/Checked_18A1_0$i.BPD.txt
                nohup ${path[1]}/prod/MainRawBPDEvent ${path[1]}/Datafilepaths.txt 6 0 1 ${reco[$r]} >> ${path[2]}/Checked_18A1_0$i.BPD.txt
                break
        fi
   fi
done	
