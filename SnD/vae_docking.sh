#!/bin/bash


nligand=$1
nligand=$((nligand-1))
nligand=`seq 0 ${nligand}`

# Record output file directory
cdir=$2
# path to receptor file
receptor=$3
# goto working directory
#cd $2/$3

# run docking
echo "Docking in process..."
for i in ${nligand}
do
	ligand=design_${i}.sd
	if  [[  -f ${cdir}/${ligand} ]]
	then
		rbdock -i ${cdir}/${ligand} -o ${cdir}/pose_org_${i} -r ${receptor} -p dock.prm -T 1  -n 100 > ${cdir}/poses_org_${i}.out & 
		### docking without solvation terms, for effeciency use dock_solv.prm for solvation term
	fi
done
wait

# write scores
echo "Writing scores..."
rm -f ${cdir}/scores.txt
for i in ${nligand}
do
	ligand=design_${i}.sd
	
	if  [[  -f ${cdir}/${ligand} ]]
	then
		sdsort -n -f'SCORE.INTER' ${cdir}/pose_org_${i}.sd > ${cdir}/pose_org_sorted_${i}.sd
		score=`grep -A 1 '<SCORE.INTER>' ${cdir}/pose_org_sorted_${i}.sd | grep -v SCORE | grep -v '\-\-' | head -1 `
		echo "${i} ${score}" | tee -a ${cdir}/scores.txt
	fi
done
