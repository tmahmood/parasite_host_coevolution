use std::fmt::{Display, Formatter};
use futures::join;
use ndarray::{Array, Ix, Ix1, Ix2, Ix3};
use ndarray_rand::RandomExt;
use rand::distributions::Uniform;
use rand::prelude::SliceRandom;
use rand::Rng;
use crate::hosts::{create_random_hosts, Host};
use crate::{generate_individual, SimulationPref};
use crate::a2d::A2D;

pub struct Simulation {
    pref: SimulationPref,
    no_of_simulation_run: usize,
    no_of_generations: usize,
    hosts: Array<Host, Ix1>,
    parasites: Array<usize, Ix3>,
    current_generation: usize,
}

pub async fn new_simulation(pref: SimulationPref) -> Simulation {
    let total = pref.a() + pref.b();
    // Array::random((total), Uniform::new(0, pref.f()));
    // create random hosts
    let hosts = create_random_hosts(&pref);
    // create parasites for each species
    let parasites = async {
        Array::random(
            (pref.d(), pref.e(), pref.g()),
            Uniform::new(0, pref.f()),
        )
    };
    let result = join!(hosts, parasites);
    Simulation {
        pref: pref.clone(),
        no_of_simulation_run: pref.gg(),
        no_of_generations: pref.ff(),
        hosts: result.0,
        parasites: result.1,
        current_generation: 0,
    }
}

impl Simulation {
    pub fn pref(&self) -> &SimulationPref {
        &self.pref
    }
    pub fn no_of_generations(&self) -> usize {
        self.no_of_generations
    }
    pub fn hosts(&self) -> &Array<Host, Ix1> {
        &self.hosts
    }

    pub fn kill_host(&mut self, index: usize) {
        self.hosts[index].set_alive(false);
    }

    pub fn hosts_mut(&mut self) -> &mut Array<Host, Ix1> {
        &mut self.hosts
    }

    pub fn parasites(&self) -> &Array<usize, Ix3> {
        &self.parasites
    }

    pub fn no_of_simulation_run(&self) -> usize {
        self.no_of_simulation_run
    }

    pub fn current_generation(&self) -> usize {
        self.current_generation
    }

    pub fn inc_generation(&mut self) {
        self.current_generation += 1;
    }
}


impl Display for Simulation {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let mut s = String::new();
        let hosts = self.hosts();
        s.push_str(&hosts.to_string());
        s.push_str("\n");
        s.push_str(&format!("{}", self.parasites()));
        write!(f, "{}", s)
    }
}
