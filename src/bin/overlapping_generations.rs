use clap::{value_t, App, Arg};
use rand::rngs::StdRng;
use rand::Rng;
use rand::SeedableRng;
use std::thread;
use std::sync::Arc;
use tskit_rust_example_programs::diploid::*;
use tskit_rust_example_programs::seeding as seeding;

struct ProgramOptions {
    params: SimParams,
    treefile: String,
    seed: u64,
    nthreads: i32,
    nreps: i32,
}

impl Default for ProgramOptions {
    fn default() -> Self {
        Self {
            params: SimParams::default(),
            treefile: String::from("treefile"),
            seed: 0,
            nthreads: 1,
            nreps: 1,
        }
    }
}

struct RunParams {
    params: SimParams,
    seeds: Vec<u64>,
    first_rep_id: usize,
    prefix: String,
}

#[derive(Debug, Clone)]
struct BadParameter {
    msg: String,
}

impl std::fmt::Display for BadParameter {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.msg)
    }
}

impl ProgramOptions {
    fn new() -> Self {
        let mut options = Self::default();

        let matches = App::new("overlapping_generations")
            .arg(
                Arg::with_name("popsize")
                    .short("N")
                    .long("popsize")
                    .help("Diploid population size. Default = 1,000.")
                    .takes_value(true),
            )
            .arg(
                Arg::with_name("nsteps")
                    .short("n")
                    .long("nsteps")
                    .help("Number of birth steps to simulate. For non-overlapping generations, this is the number of generations to simulate. Default = 1,000.")
                    .takes_value(true),
            )
            .arg(
                Arg::with_name("xovers")
                    .short("x")
                    .long("xovers")
                    .help("Mean number of crossovers per meiosis. The number of crossovers is Poisson-distributed with this value. Default = 0.0.")
                    .takes_value(true),
            )
            .arg(
                Arg::with_name("genome_length")
                    .short("L")
                    .long("genome_length")
                    .help("Genome length (continuous units).  Default = 1e6.")
                    .takes_value(true),
            )
            .arg(
                Arg::with_name("simplification_interval")
                    .short("s")
                    .long("simplify")
                    .help("Number of birth steps between simplifications. Default = 100.")
                    .takes_value(true),
            )
            .arg(
                Arg::with_name("treefile")
                    .short("t")
                    .long("treefile")
                    .help("Pref of output file. The format is a tskit \"trees\" file. Default = \"treefile\".")
                    .takes_value(true),
            )
            .arg(
                Arg::with_name("seed")
                    .short("S")
                    .long("seed")
                    .help("Random number seed. Default = 0.")
                    .takes_value(true),
            )
            .arg(
                Arg::with_name("psurvival")
                    .short("P")
                    .long("psurvival")
                    .help("Survival probability. A value of 0.0 is the Wright-Fisher model of non-overlapping generations.  Values must b 0.0 <= p < 1.0.  Default = 0.0.")
                    .takes_value(true),
            )
            .arg(Arg::with_name("nthreads").short("T").long("nthreads").help("Number of threads to use. Default = 1").takes_value(true))
            .arg(Arg::with_name("nreps").short("r").long("nreps").help("Number replicates to run. Default = 1").takes_value(true))
            .get_matches();

        options.params.popsize =
            value_t!(matches.value_of("popsize"), u32).unwrap_or(options.params.popsize);
        options.params.nsteps =
            value_t!(matches.value_of("nsteps"), u32).unwrap_or(options.params.nsteps);
        options.params.xovers =
            value_t!(matches.value_of("xovers"), f64).unwrap_or(options.params.xovers);
        options.params.genome_length = value_t!(matches.value_of("genome_length"), f64)
            .unwrap_or(options.params.genome_length);
        options.params.simplification_interval =
            value_t!(matches.value_of("simplification_interval"), u32)
                .unwrap_or(options.params.simplification_interval);
        options.params.psurvival =
            value_t!(matches.value_of("psurvival"), f64).unwrap_or(options.params.psurvival);
        options.seed = value_t!(matches.value_of("seed"), u64).unwrap_or(options.seed);
        options.treefile =
            value_t!(matches.value_of("treefile"), String).unwrap_or(options.treefile);
        options.nthreads = value_t!(matches.value_of("nthreads"), i32).unwrap_or(options.nthreads);
        options.nreps = value_t!(matches.value_of("nreps"), i32).unwrap_or(options.nreps);

        options.validate().unwrap();
        options
    }

    // NOTE: This function is incomplete.
    fn validate(&self) -> Result<(), BadParameter> {
        match self.params.psurvival.partial_cmp(&0.0) {
            Some(std::cmp::Ordering::Less) => {
                return Err(BadParameter {
                    msg: String::from("psurvival must be 0 <= p < 1.0"),
                });
            }
            Some(_) => (),
            None => (),
        }

        match self.params.psurvival.partial_cmp(&1.0) {
            Some(std::cmp::Ordering::Less) => (),
            Some(_) => {
                return Err(BadParameter {
                    msg: String::from("psurvival must be 0 <= p < 1.0"),
                });
            }
            None => (),
        }

        Ok(())
    }
}

fn overlapping_generations(params: SimParams, seed: u64) -> tskit::TableCollection {
    let mut tables = match tskit::TableCollection::new(params.genome_length) {
        Ok(x) => x,
        Err(e) => panic!("{}", e),
    };

    let mut rng = StdRng::seed_from_u64(seed);

    let mut alive: Vec<Diploid> = vec![];
    for _ in 0..params.popsize {
        let node0 = match tables.add_node(0, params.nsteps as f64, tskit::TSK_NULL, tskit::TSK_NULL)
        {
            Ok(x) => x,
            Err(e) => panic!("{}", e),
        };
        let node1 = match tables.add_node(0, params.nsteps as f64, tskit::TSK_NULL, tskit::TSK_NULL)
        {
            Ok(x) => x,
            Err(e) => panic!("{}", e),
        };
        alive.push(Diploid { node0, node1 });
    }

    let mut parents: Vec<Parents> = vec![];

    for step in (0..params.nsteps).rev() {
        parents.clear();
        death_and_parents(&alive, &params, &mut parents, &mut rng);
        births(&parents, &params, step, &mut tables, &mut alive, &mut rng);

        if step % params.simplification_interval == 0 {
            simplify(&mut alive, &mut tables);
        }
    }

    tables.build_index().unwrap();

    tables
}

fn run(params: SimParams, seed: u64, nreps: i32, prefix: String) {
    let mut rng = StdRng::seed_from_u64(seed);
    let rseed = rand_distr::Uniform::new(0_u64, u64::MAX);
    for _ in 0..nreps {
        let repseed = rng.sample(rseed);
        let mut tables = overlapping_generations(params, repseed);
        let mut outfile = prefix.to_string();
        outfile.push_str(&"_".to_string());
        outfile.push_str(&params.psurvival.to_string());
        outfile.push_str(&"_".to_string());
        outfile.push_str(&repseed.to_string());
        outfile.push_str(&".trees".to_string());
        tables.build_index().unwrap();
        tables
            .dump(&outfile, tskit::TableOutputOptions::empty())
            .unwrap();
    }
}

fn main() {
    let options = ProgramOptions::new();

    let seeds = seeding::make_unique_seeds(options.seed, options.nreps);

    let reps_per_thread = seeds.len() / options.nthreads as usize;
    let mut repid = 0_usize;

    for _ in 0..options.nthreads-1 {
        if repid + reps_per_thread < seeds.len() {
            println!("this thread gets seeds {} to {}", repid, repid+reps_per_thread);
        } 
        repid += reps_per_thread;
    }
    println!("last thread gets {} {}", repid, seeds.len());


    let mut rng = StdRng::seed_from_u64(options.seed);
    let rseed = rand_distr::Uniform::new(0_u64, u64::MAX);
    let mut handles = vec![];
    for _ in 0..options.nthreads {
        // rust won't let us borrow values
        // to pass to threads, so we make local copies
        // that we then move into the lambda function send to
        // the spawned thread.
        let threadseed = rng.sample(rseed);
        let pclone = options.params.clone();
        let nr = options.nreps;
        let fname = options.treefile.to_string();
        let h = thread::spawn(move || run(pclone, threadseed, nr, fname));
        handles.push(h);
    }

    // If you don't join a thread, it never runs.
    // Unlike C++, not joining a thread is not a runtime error.
    for h in handles {
        h.join().unwrap();
    }
}
