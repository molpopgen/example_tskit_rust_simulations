for p in 0.0 0.1 0.25 0.5
do
    SEED=$RANDOM
    ./target/release/overlapping_generations -N 500 -n 10000 -P $p -S $SEED -x 1e-3 -t treefile --nthreads 4 --nreps 100
done
