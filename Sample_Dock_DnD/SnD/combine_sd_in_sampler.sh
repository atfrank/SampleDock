## Combine all .sd files in VAEdesigns_${date}/ by folders
## takes date as argument. Date Format Exp: Oct19
## Output file name is "all_desings.sd"
directory=$1
p=`ls -dv ${directory}/design_*`
for f in ${p}
do
    alldesigns=`ls -v ${f}/design*.sd`
    obabel -i sd ${alldesigns} -O ${f}/all_designs.sd
done
echo 'Last Cycle'
ls -dv ${directory}/design_* | tail -1