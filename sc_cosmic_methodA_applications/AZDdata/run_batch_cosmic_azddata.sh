module load gcc/6.2.0
module load R
taskfile=/n/groups/price/kushal/Causal/data/batch_joblist_AD.txt
cell=/n/groups/price/kushal/Causal/codes
IFS="
"

for step in `cat $taskfile`;
do
    code=`echo $step | awk '{print $1}'`
    context1=`echo $step | awk '{print $2}'`
    context2=`echo $step | awk '{print $3}'`
    celltype=`echo $step | awk '{print $4}'`
    method=`echo $step | awk '{print $5}'`
    output_cell=`echo $step | awk '{print $6}'`
    output_name=`echo $step | awk '{print $7}'`
    echo $context1 $context2 $celltype $method $output_cell $output_name
    if [ ! -f $output_cell/Causal_VIDE_${output_name}_${celltype}_${context1}_${context2}_${method}.twosided.txt ]
    then
	cmd="Rscript $code  $context1 $context2 $celltype $method $output_cell $output_name"
	sbatch --time=100:00 --mem=120000 --output=causal.out --error=causal.err -p short -c 1 --wrap="$cmd"
    fi
done


