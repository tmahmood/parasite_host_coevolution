use std::fmt::{Display, Formatter};
use futures::join;
use rand::distributions::Uniform;
use rand::prelude::SliceRandom;
use rand::Rng;
use crate::hosts::{create_random_hosts, Host};
use crate::parasites::{create_random_parasites, Parasite};
use crate::{generate_individual, SimulationPref};

pub struct Simulation {
    pref: SimulationPref,
    no_of_simulation_run: i32,
    no_of_generations: i32,
    hosts: Vec<Host>,
    parasites: Vec<Vec<Parasite>>,
}

pub async fn new_simulation(pref: SimulationPref) -> Simulation {
    let hosts = create_random_hosts(&pref);
    let parasites = create_random_parasites(&pref);
    let result = join!(hosts, parasites);
    Simulation {
        pref: pref.clone(),
        no_of_simulation_run: pref.gg(),
        no_of_generations: pref.ff(),
        hosts: result.0,
        parasites: result.1,
    }
}

impl Simulation {
    pub fn pref(&self) -> &SimulationPref {
        &self.pref
    }
    pub fn no_of_generations(&self) -> i32 {
        self.no_of_generations
    }
    pub fn hosts(&self) -> &Vec<Host> {
        &self.hosts
    }

    pub fn parasites(&self) -> &Vec<Vec<Parasite>> {
        &self.parasites
    }
}


impl Display for Simulation {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let mut s = String::new();
        let hosts = self.hosts();
        for host in hosts {
            s.push_str(&host.to_string());
            s.push_str("\n");
        }
        s.push_str("\n");
        let parasite_species = self.parasites();
        for (k, species) in parasite_species.iter().enumerate() {
            s.push_str(&format!("+ {}", k));
            for (j, parasite) in species.iter().enumerate() {
                if j % 10 == 0 { s.push_str("\n") }
                s.push_str(&format!("| {: >2}", parasite))
            }
            s.push_str("\n");
        }
        write!(f, "{}", s)
    }
}
