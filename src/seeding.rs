use rand::rngs::StdRng;
use rand::Rng;
use rand::SeedableRng;

/// Generate a vector of `nseeds` unique [`u64`] values.
/// These values can be used as unique seeds for different replicates
/// of a simulation.
pub fn make_unique_seeds(initial_seed: u64, nseeds: usize) -> Vec<u64> {
    let rseed = rand_distr::Uniform::new(0_u64, u64::MAX);
    let mut rng = StdRng::seed_from_u64(initial_seed);
    let mut rv = vec![];
    let mut used_seeds = std::collections::HashSet::new();

    for _ in 0..nseeds {
        let mut repseed = rng.sample(rseed);
        while used_seeds.contains(&repseed) {
            repseed = rng.sample(rseed);
        }
        used_seeds.insert(repseed);
        rv.push(repseed);
    }
    rv
}
