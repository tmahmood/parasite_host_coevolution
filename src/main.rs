#[macro_use]
extern crate serde;
extern crate serde_derive;
extern crate serde_ini;
extern crate rand;
extern crate core;

use std::collections::HashMap;
use std::fmt::{Display, Formatter};
use std::fs;
use futures::executor::block_on;
use ndarray::{Array, Array1, ArrayView, Axis, Ix, Ix1};
use ndarray_rand::RandomExt;
use rand::distributions::Uniform;
use rand::{Rng, thread_rng};
use rand::rngs::ThreadRng;
use rand::seq::SliceRandom;
use crate::hosts::{Host, HostTypes};
use crate::simulation::{new_simulation, Simulation};
use crate::simulation_pref::SimulationPref;

// It's important to shuffle lists before the random processes. In the past when I've hired for
// these sorts of programs, the biggest problem has been that individuals of one type were more
// likely to get chosen for reproduction or death because lists weren't shuffled, and that causes
// the random processes to be biased.

pub mod a2d;
pub mod simulation_pref;
pub mod parasites;
pub mod hosts;
pub mod simulation;

pub mod world {
    use crate::a2d::A2D;
    use crate::Host;

    pub struct World {
        hosts: A2D<Host>,
    }
}

pub fn generate_individual(f: usize, len: usize) -> Array1<usize> {
    Array::random(len, Uniform::new(0, f))
    // let mut rng = rand::thread_rng();
    // let range = ;
    // Array::from((0..len).map(|_| rng.sample(&range)).collect::<usize>())
}

#[derive(Debug, Clone)]
pub struct ParasiteSpeciesIndex {
    species_index: usize,
    parasite_index: usize,
    match_count: usize,
}

impl ParasiteSpeciesIndex {
    fn species(&self) -> usize {
        self.species_index
    }
    fn parasite(&self) -> usize {
        self.parasite_index
    }
    fn match_count(&self) -> usize {
        self.match_count
    }
}

impl Display for ParasiteSpeciesIndex {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "S{: <3} P{: <3} M{: <3} ", self.species_index, self.parasite_index, self.match_count)
    }
}

/**
In the instructions, the thing that distinguishes between the two types is the part that begins
here: "If 1) the current generation is after the Lth generation and if 2) less than an M fraction of
host individuals (a total of reservation and wild host individuals) have been killed, then AA% of
reservation hosts get exposed to one individual each from I additional parasite species...." So, if
L is some large number like 10000 or M=0, then there shouldn't be a difference between reservation
and wild. They should be winning at roughly the same rates. Also P goes against reservation type in
reproduction but that's a minor difference.
 **/
pub async fn expose_hosts(simulation: &mut Simulation) {
    let mut killed_hosts_list = vec![];
    let mut unmatched_parasites = vec![];
    let mut faced_parasites = HashMap::new();
    let mut no_of_reservation_host_died = 0;
    let mut no_of_wild_host_died = 0;
    // If a host individual has a match score with at least X parasite individuals that is lower than N,
    // then that host individual is considered killed.
    let hosts = simulation.hosts();
    for i in 0..simulation.pref().a() + simulation.pref().b() {
        let host = &hosts[i];
        let mut match_score_bellow_threshold = 0;
        let selected_parasites = parasites::chose_parasites(simulation).await;
        let all_parasites = simulation.parasites();
        for parasite_k in selected_parasites.clone() {
            let s = all_parasites.index_axis(Axis(0), parasite_k.species());
            let parasite = s.index_axis(Axis(0), parasite_k.parasite());
            let (bellow_threshold, match_count) = find_match_score(&host, &parasite, parasite_k.species(), simulation.pref());
            if bellow_threshold {
                match_score_bellow_threshold += 1;
            } else {
                unmatched_parasites.push(ParasiteSpeciesIndex { match_count, ..parasite_k });
            }
        }
        if match_score_bellow_threshold >= simulation.pref().x() {
            match host.host_type() {
                HostTypes::Reservation => no_of_reservation_host_died += 1,
                HostTypes::Wild => no_of_wild_host_died += 1
            }
            killed_hosts_list.push(i)
        } else {
            faced_parasites.insert(i, selected_parasites);
        }
    }
    // clear the killed hosts
    // simulation.remove_hosts(&killed_hosts_list);
    // secondary exposer
    let secondary_allowed = ((simulation.pref().a() + simulation.pref().b()) as f32 * simulation.pref().m());
    if simulation.current_generation() > simulation.pref().l() && killed_hosts_list.len() < secondary_allowed as usize {
        println!("Secondary exposer used");
    }

    let no_of_reservation_host_alive = simulation.pref().a() - no_of_reservation_host_died;
    let no_of_wild_host_alive = simulation.pref().b() - no_of_wild_host_died;

    let chance = get_chance(
        HostTypes::Reservation,
        no_of_wild_host_died as f32,
        no_of_reservation_host_alive as f32,
        no_of_reservation_host_died as f32,
        no_of_wild_host_alive as f32,
        no_of_wild_host_alive as f32,
        simulation.pref().o(),
        simulation.pref().p(),
    );
    println!("{} hosts got killed, {} left {}", killed_hosts_list.len(), simulation.hosts().len(), chance);
}

pub fn get_chance(host_type: HostTypes, v: f32, u: f32, w: f32, y: f32, t: f32, o: f32, p: f32) -> f32 {
    match host_type {
        HostTypes::Reservation => (1. + o * y * w / u + (1. - o) * y * w / (t + u) + (1. - o) * y * v / (t + u) - p) as f32,  // for reservation host individuals
        HostTypes::Wild => (1. + o * y * v / t + (1. - o) * y * v / (t + u) + (1. - o) * y * w / (t + u)) as f32, // for wild host individuals
    }
}


async fn main_async() {
    // TODO Need better way to feed the input file
    let s: SimulationPref = serde_ini::from_str(&fs::read_to_string("params.conf").unwrap()).unwrap();
    let mut simulation = new_simulation(s).await;
    println!("{}", simulation);
}

pub fn find_match_score(host: &Host, parasite: &ArrayView<usize, Ix1>, index: usize, pref: &SimulationPref) -> (bool, usize) {
    let mut match_count = 0;
    // for ii in index..pref.g() {
    //     if host.number_set()[ii] == parasite[ii - index] {
    //         match_count += 1
    //     }
    // }
    (match_count < pref.n(), match_count)
}


fn main() {
    block_on(main_async());
}


