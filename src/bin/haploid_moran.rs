use clap::{value_t, App, Arg};
use rand::rngs::StdRng;
use rand::Rng;
use rand::SeedableRng;
use rand_distr::Uniform;
use std::sync::Arc;
use std::thread;
use tskit::TableAccess;
use tskit_rust_example_programs::seeding;

#[derive(Clone)]
struct ProgramOptions {
    treefile: String,
    popsize: i32,
    nsteps: i32,
    seed: u64,
    nthreads: i32,
    nreps: i32,
}

impl Default for ProgramOptions {
    fn default() -> Self {
        Self {
            treefile: String::from("treefile_moran"),
            popsize: 1000,
            nsteps: 10000,
            seed: 0,
            nthreads: 1,
            nreps: 1,
        }
    }
}

impl ProgramOptions {
    fn new() -> Self {
        let mut options = Self::default();

        let matches = App::new("haploid_moran")
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
                Arg::with_name("treefile")
                    .short("t")
                    .long("treefile")
                    .help("Pref of output file. The format is a tskit \"trees\" file. Default = \"treefile_moran\".")
                    .takes_value(true),
            )
            .arg(
                Arg::with_name("seed")
                    .short("S")
                    .long("seed")
                    .help("Random number seed. Default = 0.")
                    .takes_value(true),
            )
            .arg(Arg::with_name("nthreads").short("T").long("nthreads").help("Number of threads to use. Default = 1").takes_value(true))
            .arg(Arg::with_name("nreps").short("r").long("nreps").help("Number replicates to run. Default = 1").takes_value(true))
            .get_matches();

        options.popsize = value_t!(matches.value_of("popsize"), i32).unwrap_or(options.popsize);
        options.nsteps = value_t!(matches.value_of("nsteps"), i32).unwrap_or(options.nsteps);
        options.seed = value_t!(matches.value_of("seed"), u64).unwrap_or(options.seed);
        options.treefile =
            value_t!(matches.value_of("treefile"), String).unwrap_or(options.treefile);
        options.nthreads = value_t!(matches.value_of("nthreads"), i32).unwrap_or(options.nthreads);
        options.nreps = value_t!(matches.value_of("nreps"), i32).unwrap_or(options.nreps);

        options
    }
}

fn moran(popsize: i32, nsteps: i32, seed: u64) -> tskit::TableCollection {
    let mut rng = StdRng::seed_from_u64(seed);

    let mut tables = tskit::TableCollection::new(1.0).unwrap();

    let mut alive = vec![];
    for _ in 0..popsize {
        let id = tables
            .add_node(0, nsteps as f64, tskit::TSK_NULL, tskit::TSK_NULL)
            .unwrap();
        alive.push(id);
    }

    let picker = Uniform::new(0, popsize);

    for step in (0..nsteps).rev() {
        let dead = rng.sample(picker);
        let replace = rng.sample(picker);

        if dead != replace {
            // birth (and coalescence)
            let new_birth = tables
                .add_node(0, step as f64, tskit::TSK_NULL, tskit::TSK_NULL)
                .unwrap();
            assert!(
                tables.nodes().time(new_birth).unwrap()
                    < tables.nodes().time(alive[replace as usize]).unwrap()
            );
            tables
                .add_edge(0., 1., alive[replace as usize], new_birth)
                .unwrap();
            alive[dead as usize] = new_birth;
        }

        if step % 100 == 0 {
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

    tables
}

fn finalise_tables_and_output(
    options: ProgramOptions,
    repid: usize,
    seed: u64,
    tables: tskit::TableCollection,
) {
    let mut tables = tables; // this is a simple move
    use tskit::provenance::Provenance;
    let provenance = format!(
        "{{\"seed\": {}, \"N\": {},  \"nsteps\": {} }}",
        seed, options.popsize, options.nsteps,
    );
    tables.add_provenance(&provenance).unwrap();
    let mut outfile = options.treefile.to_string();
    outfile.push_str(&"_".to_string());
    outfile.push_str(&repid.to_string());
    outfile.push_str(&".trees".to_string());
    tables.build_index().unwrap();
    tables
        .dump(&outfile, tskit::TableOutputOptions::empty())
        .unwrap();
}

fn run_from_seeds(params: ProgramOptions, seeds: &[u64], first_rep_id: usize) {
    for (idx, seed) in seeds.iter().enumerate() {
        let tables = moran(params.popsize, params.nsteps, *seed);
        finalise_tables_and_output(params.clone(), first_rep_id + idx, *seed, tables);
    }
}

fn run_in_thread(run_params_arc: Arc<(ProgramOptions, usize, Vec<u64>)>) {
    let run_params = &*run_params_arc;

    run_from_seeds(run_params.0.clone(), &run_params.2, run_params.1);
}

fn run_threaded(options: ProgramOptions, seeds: Vec<u64>) {
    let mut handles = vec![];
    let reps_per_thread = seeds.len() / options.nthreads as usize;
    let mut repid = 0_usize;
    for _ in 0..options.nthreads - 1 {
        assert!(repid + reps_per_thread < seeds.len());
        let run_params = Arc::new((
            options.clone(),
            repid,
            seeds[repid..repid + reps_per_thread].to_vec(),
        ));
        let h = thread::spawn(|| run_in_thread(run_params));
        handles.push(h);
        repid += reps_per_thread;
    }
    let run_params = Arc::new((options, repid, seeds[repid..seeds.len()].to_vec()));
    let h = thread::spawn(|| run_in_thread(run_params));
    handles.push(h);

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
            run_from_seeds(options, &seeds, 0);
        }
    } else {
        // The input seed is the seed for the replicate.
        assert_eq!(options.nreps, 1);
        let tables = moran(options.popsize, options.nsteps, options.seed);
        let seed = options.seed;
        finalise_tables_and_output(options, 0, seed, tables);
    }
}
