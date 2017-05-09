prefix="exp_geneEC_adapt"
logdir="out/"
outdir="out/info/"

if [ ! -d $outdir ]; then
  mkdir -p $outdir
else
  echo "Cleaning everything in ${outdir}"
  rm -r ${outdir}* 2> /dev/null
fi

for exptype in 17_kl 18_kendall 19_rankdiff; do
  echo "Experiment type $exptype"
  expid=0
  for d in 1 2 3 4 5; do
    echo "Dimension no. $d"

    if [ $d -eq 5 ]; then
      dimstr=20D
      expid=0
    else
      dimstr=10D
    fi

    for fun in `seq 1 24`; do
      for inst in `seq 1 15`; do
        expid=`expr $expid + 1`
        suffix=${d}_${fun}_${inst}
        expfile=`find ${logdir}${prefix}_${exptype}_${dimstr}__$expid.o*`
        # echo "Processing file $expfile for dim $d, fun $fun, inst $inst"
        grep --color=never "Model.generationUpdate(): We are using the model" $expfile | wc -l > $outdir${exptype}_gmodel.$suffix
        grep --color=never "model trained" $expfile | wc -l > $outdir${exptype}_gorig.$suffix
        sed -n -e 's/^.*err = \([0-9.]*\).*$/\1/p' $expfile > $outdir${exptype}_err.$suffix

        if [ ! -s $outdir${exptype}_err.$suffix ]; then
          echo "0" > $outdir${exptype}_err.$suffix
        fi
      done
    done
  done
done
