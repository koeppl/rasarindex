vars=$1
seed_length=20

for i in `seq 1 1 5`;
do 
    echo $vars ${seed_length} $(r-index/build/ri-time indexes/chr19.${vars}.fa.ri reads/chr19.1.fa.seeds_${seed_length}.pizza 2>&1 | grep "Took" | cut -d' ' -f 2 | cut -d's' -f 1)
done 
