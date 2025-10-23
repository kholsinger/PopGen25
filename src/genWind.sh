mapfile -t chr < chrs.txt
mapfile -t lengths < lengths.txt
wnd=$2
for k in `seq 1 ${#chr[@]}`; do
        let j=$k-1
        c=${chr[${j}]}
        START=0
        END=${lengths[${j}]}/$wnd+1
        MAX=${lengths[${j}]}
        for ((i=$START;i<$END;i++)); do
                        let start=$i*$wnd+1
                        let end=$start+$wnd
                        end=$((end<MAX ? end : MAX))
                        echo "##${c}_${i}_wnd${wnd}"
                        echo "bcftools view -r ${c}:${start}-${end} ${1} > ${c}_${i}_wnd${wnd}.vcf"
                        echo "plink2 --vcf ${c}_${i}_wnd${wnd}.vcf --snps-only --allow-extra-chr --export phylip-phased --out ${c}_${i}_wnd${wnd}"
                        echo "rm ${c}_${i}_wnd${wnd}.vcf"
                        echo "iqtree3 -redo -pre ${c}_${i}_wnd${wnd} -nt AUTO -ntmax 4 -bb 1000 -s ${c}_${i}_wnd${wnd}.phy"
                done
        done