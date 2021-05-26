if [ ! -d output ]
then
    mkdir output
fi

for p in 0.0 0.1 0.25 0.5
do
    SEED=$RANDOM
    ./target/release/overlapping_generations -N 500 -n 20000 -P $p -S $SEED -x 1e-3 \
        -t output/treefile_$p --nthreads 12 --nreps 500
done
