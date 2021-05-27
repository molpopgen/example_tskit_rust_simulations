use clap::{value_t, App, Arg};
use rand::rngs::StdRng;
use rand::Rng;
use rand::SeedableRng;
use rand_distr::Uniform;
use std::sync::Arc;
use std::thread;
use tskit_rust_example_programs::diploid::*;
use tskit_rust_example_programs::seeding;

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

// Replace nodes at positions
// 2i and 2i + 1 with node1 and node2,
// respectively
struct Replacement {
    index: usize,
    node1: tskit::tsk_id_t,
    node2: tskit::tsk_id_t,
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

    let mut alive = vec![];

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
        alive.push(node0);
        alive.push(node1);
    }

    let mut replacements = vec![];

    // Used to pick the parents for a Replacement
    let picker = Uniform::new(0, params.popsize as usize);

    for step in (0..params.nsteps).rev() {
        replacements.clear();

        // Generate deaths, record replacement nodes
        for index in 0..params.popsize as usize {
            let x: f64 = rng.gen();
            match x.partial_cmp(&params.psurvival) {
                Some(std::cmp::Ordering::Greater) => {
                    // Generate two offspring nodes
                    let node1 = tables
                        .add_node(0, step as f64, tskit::TSK_NULL, tskit::TSK_NULL)
                        .unwrap();
                    let node2 = tables
                        .add_node(0, step as f64, tskit::TSK_NULL, tskit::TSK_NULL)
                        .unwrap();
                    // Record that individual i will be replaced
                    // by the two new nodes
                    replacements.push(Replacement {
                        index,
                        node1,
                        node2,
                    });
                }
                Some(_) => (),
                None => panic!("bad floating point comparison"),
            }
        }

        // For each replacement, pick parents and add edges
        for rep in &replacements {
            for offspring_node_ in &[rep.node1, rep.node2] {
                let parent_index = rng.sample(picker);
                let mut node1 = alive[2 * parent_index];
                let mut node2 = alive[2 * parent_index + 1];

                // FIXME: use crossover code in lib
                // Pick which gamete to pass on
                let x: f64 = rng.gen();
                match x.partial_cmp(&0.5) {
                    Some(std::cmp::Ordering::Less) => {
                        std::mem::swap(&mut node1, &mut node2);
                    }
                    Some(_) => (),
                    None => panic!("Unexpected None"),
                }
                // record the edge
                tables
                    .add_edge(0., tables.sequence_length(), node1, *offspring_node_)
                    .unwrap();
            }
        }

        // Finally, replace the parent nodes with the new births
        for rep in &replacements {
            alive[2 * rep.index] = rep.node1;
            alive[2 * rep.index + 1] = rep.node1;
        }

        if step % params.simplification_interval == 0 {
            match tables.full_sort(tskit::TableSortOptions::default()) {
                Ok(_) => (),
                Err(e) => panic!("{}", e),
            }
            match tables.simplify(&alive, tskit::SimplificationOptions::empty(), true) {
                Ok(x) => match x {
                    Some(idmap) => {
                        for a in alive.iter_mut() {
                            *a = idmap[*a as usize];
                        }
                    }
                    None => panic!("expected an id map!"),
                },
                Err(e) => panic!("{}", e),
            }
        }
    }

    tables.build_index().unwrap();

    tables
}

fn finalise_tables_and_output(
    params: SimParams,
    seed: u64,
    repid: usize,
    tables: tskit::TableCollection,
    outfile_prefix: &str,
) {
    let mut tables = tables; // this is a simple move
    use tskit::provenance::Provenance;
    let provenance = format!(
        "{{\"seed\": {}, \"N\": {}, \"psurvival\": {}, \"nsteps\": {}, \"recrate\": {}}}",
        seed, params.popsize, params.psurvival, params.nsteps, params.xovers,
    );
    tables.add_provenance(&provenance).unwrap();
    let mut outfile = outfile_prefix.to_string();
    outfile.push_str(&"_".to_string());
    outfile.push_str(&repid.to_string());
    outfile.push_str(&".trees".to_string());
    tables.build_index().unwrap();
    tables
        .dump(&outfile, tskit::TableOutputOptions::empty())
        .unwrap();
}

fn run_from_seeds(params: SimParams, seeds: &[u64], first_rep_id: usize, outfile_prefix: &str) {
    for (idx, seed) in seeds.iter().enumerate() {
        let tables = overlapping_generations(params, *seed);
        finalise_tables_and_output(params, *seed, first_rep_id + idx, tables, outfile_prefix);
    }
}

fn run_in_thread(run_params_arc: Arc<RunParams>) {
    let run_params = &*run_params_arc;

    run_from_seeds(
        run_params.params,
        &run_params.seeds,
        run_params.first_rep_id,
        &run_params.prefix,
    );
}

// This function handles using many threads to run the simulations.
// rust requires "safe" sharing of data, so no "just sent a pointer
// and promise" stuff will work here.
// We use Arc (https://doc.rust-lang.org/std/sync/struct.Arc.html)
// to manage the data being sent.
// This function creates the Arc, and run_in_thread (see above) will consume it.
fn run_threaded(options: ProgramOptions, seeds: Vec<u64>) {
    let mut handles = vec![];
    let reps_per_thread = seeds.len() / options.nthreads as usize;
    let mut repid = 0_usize;
    for _ in 0..options.nthreads - 1 {
        assert!(repid + reps_per_thread < seeds.len());
        let run_params = Arc::new(RunParams {
            params: options.params,
            seeds: seeds[repid..repid + reps_per_thread].to_vec(),
            first_rep_id: repid,
            prefix: options.treefile.to_string(),
        });
        let h = thread::spawn(|| run_in_thread(run_params));
        handles.push(h);
        repid += reps_per_thread;
    }
    let run_params = Arc::new(RunParams {
        params: options.params,
        seeds: seeds[repid..seeds.len()].to_vec(),
        first_rep_id: repid,
        prefix: options.treefile,
    });
    let h = thread::spawn(|| run_in_thread(run_params));
    handles.push(h);

    // If you don't join a thread, it never runs.
    // Unlike C++, not joining a thread is not a runtime error.
    for h in handles {
        h.join().unwrap();
    }
}

fn main() {
    let options = ProgramOptions::new();

    if options.nreps > 1 {
        // Then the input seed is an initial seed used to generate
        // unique seeds for each replicate
        let seeds = seeding::make_unique_seeds(options.seed, options.nreps);

        if options.nthreads > 1 {
            run_threaded(options, seeds);
        } else {
            run_from_seeds(options.params, &seeds, 0, &options.treefile);
        }
    } else {
        // The input seed is the seed for the replicate.
        assert_eq!(options.nreps, 1);
        let tables = overlapping_generations(options.params, options.seed);
        finalise_tables_and_output(options.params, options.seed, 0, tables, &options.treefile);
    }
}
