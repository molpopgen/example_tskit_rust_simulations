use rand::rngs::StdRng;
use rand::Rng;
use rand_distr::{Exp, Uniform};

#[derive(Copy, Clone)]
pub struct SimParams {
    pub popsize: u32,
    pub nsteps: u32,
    pub xovers: f64,
    pub psurvival: f64,
    pub genome_length: f64,
    pub simplification_interval: u32,
}

impl Default for SimParams {
    fn default() -> Self {
        Self {
            popsize: 1000,
            nsteps: 1000,
            xovers: 0.,
            psurvival: 0.0,
            genome_length: 1e6,
            simplification_interval: 100,
        }
    }
}

#[derive(Copy, Clone)]
pub struct Diploid {
    pub node0: tskit::tsk_id_t,
    pub node1: tskit::tsk_id_t,
}

pub struct Parents {
    pub index: usize,
    pub parent0: Diploid,
    pub parent1: Diploid,
}

pub fn death_and_parents(
    alive: &[Diploid],
    params: &SimParams,
    parents: &mut Vec<Parents>,
    rng: &mut StdRng,
) {
    let random_parents = Uniform::new(0_usize, params.popsize as usize);
    for index in 0..alive.len() {
        let x: f64 = rng.gen();
        match x.partial_cmp(&params.psurvival) {
            Some(std::cmp::Ordering::Greater) => {
                let parent0 = alive[rng.sample(random_parents)];
                let parent1 = alive[rng.sample(random_parents)];
                parents.push(Parents {
                    index,
                    parent0,
                    parent1,
                });
            }
            Some(_) => (),
            None => (),
        }
    }
}

fn mendel(pnodes: &mut (tskit::tsk_id_t, tskit::tsk_id_t), rng: &mut StdRng) {
    let x: f64 = rng.gen();
    match x.partial_cmp(&0.5) {
        Some(std::cmp::Ordering::Less) => {
            std::mem::swap(&mut pnodes.0, &mut pnodes.1);
        }
        Some(_) => (),
        None => panic!("Unexpected None"),
    }
}

pub fn crossover_and_record_edges_details(
    parent: Diploid,
    offspring_node: tskit::tsk_id_t,
    params: &SimParams,
    tables: &mut tskit::TableCollection,
    rng: &mut StdRng,
) {
    let mut pnodes = (parent.node0, parent.node1);
    mendel(&mut pnodes, rng);

    if params.xovers == 0.0 {
        match tables.add_edge(0., tables.sequence_length(), pnodes.0, offspring_node) {
            Ok(_) => (),
            Err(e) => panic!("{}", e),
        }
    } else {
        let exp = match Exp::new(params.xovers / tables.sequence_length()) {
            Ok(e) => e,
            Err(e) => panic!("{}", e),
        };
        let mut current_pos = 0.0;
        loop {
            let next_length = rng.sample(exp);
            match (current_pos + next_length).partial_cmp(&tables.sequence_length()) {
                Some(std::cmp::Ordering::Less) => {
                    match tables.add_edge(
                        current_pos,
                        current_pos + next_length,
                        pnodes.0,
                        offspring_node,
                    ) {
                        Ok(_) => (),
                        Err(e) => panic!("{}", e),
                    }
                    std::mem::swap(&mut pnodes.0, &mut pnodes.1);
                    current_pos += next_length;
                }
                Some(_) => {
                    match tables.add_edge(
                        current_pos,
                        tables.sequence_length(),
                        pnodes.0,
                        offspring_node,
                    ) {
                        Ok(_) => (),
                        Err(e) => panic!("{}", e),
                    }
                    break;
                }
                None => panic!("Unexpected None"),
            }
        }
    }
}

pub fn crossover_and_record_edges(
    parents: &Parents,
    offspring_nodes: (tskit::tsk_id_t, tskit::tsk_id_t),
    params: &SimParams,
    tables: &mut tskit::TableCollection,
    rng: &mut StdRng,
) {
    crossover_and_record_edges_details(parents.parent0, offspring_nodes.0, params, tables, rng);
    crossover_and_record_edges_details(parents.parent1, offspring_nodes.1, params, tables, rng);
}

pub fn births(
    parents: &[Parents],
    params: &SimParams,
    birth_time: u32,
    tables: &mut tskit::TableCollection,
    alive: &mut [Diploid],
    rng: &mut StdRng,
) {
    for p in parents {
        // Register the two nodes for our offspring
        let node0 = match tables.add_node(
            0,                 // flags
            birth_time as f64, // time
            tskit::TSK_NULL,   // population
            // individual
            tskit::TSK_NULL,
        ) {
            Ok(x) => x,
            Err(e) => panic!("{}", e),
        };
        let node1 = match tables.add_node(0, birth_time as f64, tskit::TSK_NULL, tskit::TSK_NULL) {
            Ok(x) => x,
            Err(e) => panic!("{}", e),
        };

        // Replace a dead individual
        // with our newborn.
        alive[p.index] = Diploid { node0, node1 };

        crossover_and_record_edges(p, (node0, node1), params, tables, rng);
    }
}

pub fn simplify(alive: &mut [Diploid], tables: &mut tskit::TableCollection) {
    let mut samples = vec![];
    for a in alive.iter() {
        assert!(a.node0 != a.node1);
        samples.push(a.node0);
        samples.push(a.node1);
    }

    match tables.full_sort(tskit::TableSortOptions::default()) {
        Ok(_) => (),
        Err(e) => panic!("{}", e),
    }

    match tables.simplify(&samples, tskit::SimplificationOptions::empty(), true) {
        Ok(x) => match x {
            Some(idmap) => {
                for a in alive.iter_mut() {
                    a.node0 = idmap[a.node0 as usize];
                    assert!(a.node0 != tskit::TSK_NULL);
                    a.node1 = idmap[a.node1 as usize];
                    assert!(a.node1 != tskit::TSK_NULL);
                }
            }
            None => panic!("Unexpected None"),
        },
        Err(e) => panic!("{}", e),
    };
}
