use std::collections::HashMap;
use std::fmt::{Formatter};
use ndarray::{Array, Axis, Ix, Ix1, Ix2, Ix3};
use ndarray_rand::RandomExt;
use rand::distributions::Uniform;
use crate::hosts::{create_random_hosts, Host};
use crate::{HostTypes, SimulationPref};

pub struct Simulation {
    pref: SimulationPref,
    simulation_state: SimulationState,
    generations: Vec<SimulationState>,
}

impl Simulation {
    pub(crate) fn update_host_match_score(&mut self, host_match_score: HashMap<usize, usize>) {
        self.simulation_state.host_match_score = host_match_score;
    }
}

impl Simulation {
    pub(crate) fn species_match_score(&self) -> &HashMap<usize, usize> {
        &self.simulation_state.species_match_score
    }

    pub(crate) fn update_species_match_score(&mut self, species_match_score: HashMap<usize, usize>) {
        self.simulation_state.species_match_score = species_match_score;
    }
}

#[derive(Clone, Debug)]
pub struct SimulationState {
    match_scores: HashMap<(usize, usize), usize>,
    species_match_score: HashMap<usize, usize>,
    current_generation: usize,
    hosts: Array<Host, Ix1>,
    parasites: Array<usize, Ix3>,
    host_match_score: HashMap<usize, usize>,
}

impl Default for SimulationState {
    fn default() -> Self {
        SimulationState {
            match_scores: Default::default(),
            species_match_score: Default::default(),
            current_generation: 0,
            hosts: Default::default(),
            parasites: Default::default(),
            host_match_score: Default::default()
        }
    }
}

impl SimulationState {
    pub fn update_match_store(&mut self, k: (usize, usize), v: usize) -> Option<usize> {
        if self.match_scores.contains_key(&k) {
            println!("{}", self.match_scores.get(&k).unwrap());
        }
        self.match_scores.insert(k, v)
    }

    pub fn match_scores(&self) -> &HashMap<(usize, usize), usize> {
        &self.match_scores
    }
}

pub fn new_simulation(pref: SimulationPref) -> Simulation {
    let hosts = create_random_hosts(&pref);
    let parasites: Array<usize, Ix3> = Array::random(
        (pref.d(), pref.e(), pref.g()),
        Uniform::new(0, pref.f()),
    );
    Simulation {
        pref: pref.clone(),
        simulation_state: SimulationState {
            hosts,
            parasites,
            ..SimulationState::default()
        },
        generations: vec![],
    }
}

impl Simulation {
    /**
    update parasite with parent parasite's number set
     */
    pub(crate) fn update_parasites(&mut self, species: usize, parasite: usize, parent: usize) {
        let base = self.simulation_state.parasites.index_axis(Axis(0), species);
        let parent_view = base.index_axis(Axis(0), parent).to_owned();
        for v in 0..parent_view.len() {
            self.simulation_state.parasites[[species, parasite, v]] = parent_view[v];
        }
    }

    pub(crate) fn update_parasites_exposed_to(&mut self, parasites_exposed_to: HashMap<(usize, usize), usize>) {
        self.simulation_state.match_scores = parasites_exposed_to;
    }

    pub(crate) fn parasites_exposed_to(&mut self) -> HashMap<(usize, usize), usize> {
        self.simulation_state.match_scores.clone()
    }

    pub fn pref(&self) -> SimulationPref {
        self.pref.clone()
    }

    pub fn hosts(&self) -> &Array<Host, Ix1> {
        &self.simulation_state.hosts
    }

    pub fn update_dead_host(&mut self, index: usize, parent_host_index: usize) -> (usize, usize, usize) {
        let (host_type, number_set) = {
            let parent_host = &self.simulation_state.hosts[parent_host_index];
            let host_type = parent_host.host_type();
            (host_type, parent_host.number_set().clone())
        };
        self.simulation_state.hosts[index]
            .set_host_type(host_type)
            .set_number_set(number_set)
            .set_alive(true);
        self.count_dead_hosts()
    }

    pub fn kill_host(&mut self, index: usize) {
        let host = &mut self.simulation_state.hosts[index];
        host.set_alive(false);
    }

    pub fn count_alive_hosts(&self) -> (usize, usize, usize) {
        self.count_alive_hosts_from_generation(self.simulation_state())
    }

    pub fn count_alive_hosts_from_generation(&self, simulation_state: &SimulationState) -> (usize, usize, usize) {
        let mut count_hosts_alive_reservation = 0;
        let mut count_hosts_alive_wild = 0;
        for host in simulation_state.hosts.iter() {
            if host.alive() {
                if host.host_type() == HostTypes::Reservation {
                    count_hosts_alive_reservation += 1;
                } else {
                    count_hosts_alive_wild += 1;
                }
            }
        }
        (count_hosts_alive_reservation + count_hosts_alive_wild, count_hosts_alive_reservation, count_hosts_alive_wild)
    }

    pub fn count_dead_hosts(&mut self) -> (usize, usize, usize) {
        let mut count_hosts_dead_reservation = 0;
        let mut count_hosts_dead_wild = 0;
        for host in self.simulation_state.hosts.iter() {
            if !host.alive() {
                if host.host_type() == HostTypes::Reservation {
                    count_hosts_dead_reservation += 1;
                } else {
                    count_hosts_dead_wild += 1;
                }
            }
        }
        (count_hosts_dead_wild + count_hosts_dead_reservation, count_hosts_dead_reservation, count_hosts_dead_wild)
    }

    pub fn hosts_mut(&mut self) -> &mut Array<Host, Ix1> {
        &mut self.simulation_state.hosts
    }

    pub fn parasites(&self) -> &Array<usize, Ix3> {
        &self.simulation_state.parasites
    }

    pub fn parasites_mut(&mut self) -> &mut Array<usize, Ix3> {
        &mut self.simulation_state.parasites
    }

    pub fn current_generation(&self) -> usize {
        self.simulation_state.current_generation
    }

    pub fn next_generation(&mut self) {
        self.generations.push(self.simulation_state.clone());
        self.simulation_state = SimulationState {
            hosts: self.simulation_state.hosts.to_owned(),
            parasites: self.simulation_state.parasites.to_owned(),
            ..Default::default()
        };
        self.simulation_state.current_generation += 1;
    }
    pub fn simulation_state(&self) -> &SimulationState {
        &self.simulation_state
    }

    pub fn set_simulation_state(&mut self, simulation_state: SimulationState) {
        self.simulation_state = simulation_state;
    }

    pub fn generate_report(&self) {
        let mut wild_loner = 0;
        let mut resv_loner = 0;
        let mut more_wild = 0;
        let mut more_reserv = 0;
        let mut tied = 0;
        for generation in &self.generations {
            let (_, r, w) = self.count_alive_hosts_from_generation(generation);
            if r > w { more_reserv += 1 }
            if r < w { more_wild += 1 }
            if r == w { tied += 1 }
            if r == 0 { wild_loner += 1 }
            if w == 0 { resv_loner += 1 }
        }

        println!("- {} runs ended with wild individuals the lone type remaining", wild_loner);
        println!("- {} runs ended with reservation individuals the lone type remaining", resv_loner);
        println!("- {} runs ended with wild individuals a higher quantity than reservation individuals", more_wild);
        println!("- {} runs ended with reservation individuals a higher quantity than wild individuals", more_reserv);
        println!("- {} runs ended with the quantities of the two types tied", tied);
    }
}