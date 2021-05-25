use clap::{value_t, App, Arg};
use example_tskit_rust_simulations::diploid::*;
use rand::rngs::StdRng;
use rand::SeedableRng;

struct ProgramOptions {
    params: SimParams,
    treefile: String,
    seed: u64,
}

impl Default for ProgramOptions {
    fn default() -> Self {
        Self {
            params: SimParams::default(),
            treefile: String::from("treefile.trees"),
            seed: 0,
        }
    }
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
                    .help("Name of output file. The format is a tskit \"trees\" file. Default = \"treefile.trees\".")
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

fn main() {
    let options = ProgramOptions::new();

    let tables = overlapping_generations(options.params, options.seed);

    tables
        .dump(&options.treefile, tskit::TableOutputOptions::empty())
        .unwrap();
}
