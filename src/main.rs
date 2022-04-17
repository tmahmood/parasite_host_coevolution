#[macro_use]
extern crate serde;
extern crate serde_derive;
extern crate serde_ini;
extern crate rand;
extern crate core;

use std::fs;
use futures::executor::block_on;
use rand::distributions::Uniform;
use rand::{Rng, thread_rng};
use rand::rngs::ThreadRng;
use rand::seq::SliceRandom;
use crate::hosts::Host;
use crate::parasites::Parasite;
use crate::simulation::{new_simulation, Simulation};
use crate::simulation_pref::SimulationPref;

// It's important to shuffle lists before the random processes. In the past when I've hired for
// these sorts of programs, the biggest problem has been that individuals of one type were more
// likely to get chosen for reproduction or death because lists weren't shuffled, and that causes
// the random processes to be biased.

pub mod simulation_pref;
pub mod parasites;
pub mod hosts;
pub mod simulation;

pub fn generate_individual(f: i32, len: i32) -> Vec<i32> {
    let mut rng = rand::thread_rng();
    let range = Uniform::new(0, f);
    (0..len).map(|_| rng.sample(&range)).collect()
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
async fn main_async() {
    // TODO Need better way to feed the input file
    let s: SimulationPref = serde_ini::from_str(
        &fs::read_to_string("params.conf").unwrap()
    ).unwrap();
    let mut simulation = new_simulation(s).await;
    println!("{}", simulation);
    expose_hosts(&mut simulation).await;
}

pub fn find_match_score(host: &Host, parasite: &Parasite, index: usize, pref: &SimulationPref) -> bool {
    let number_set = host.number_set();
    let parasite_numer_set = parasite.number_set();
    let mut match_count = 0;
    for ii in index..pref.g() as usize {
        if number_set[ii] == parasite_numer_set[ii - index] { match_count += 1 }
    }
    match_count < 
}

pub async fn expose_hosts(simulation: &mut Simulation) {
    let mut killed_list = vec![];
    let hosts = simulation.hosts();
    // If a host individual has a match score with at least X parasite individuals that is lower than N,
    // then that host individual is considered killed.
    for (i, host) in hosts.iter().enumerate() {
        let mut how_many_parasites_matched = 0;
        let parasites = parasites::chose_parasites(simulation).await;
        for (index, parasite) in parasites {
            find_match_score(host, &parasite, index as usize, simulation.pref());
        }
        if how_many_parasites_matched >= simulation.pref().x() {
            killed_list.push(how_many_parasites_matched)
        }
    }
    println!("{} of {} hosts got killed", killed_list.len(), simulation.hosts().len());
}

fn main() {
    block_on(main_async());
}


