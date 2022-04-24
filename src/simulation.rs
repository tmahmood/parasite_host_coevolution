use std::collections::HashMap;
use std::fmt::{Display, Formatter};
use futures::join;
use ndarray::{Array, Axis, Ix, Ix1, Ix2, Ix3};
use ndarray_rand::RandomExt;
use rand::distributions::Uniform;
use crate::hosts::{create_random_hosts, Host};
use crate::{HostTypes, SimulationPref};

pub struct Simulation {
    pref: SimulationPref,
    no_of_simulation_run: usize,
    no_of_generations: usize,
    hosts: Array<Host, Ix1>,
    parasites: Array<usize, Ix3>,
    simulation_state: SimulationState,
}


pub struct SimulationState {
    current_generation: usize,
    count_dead_reservation_hosts: usize,
    count_dead_wild_hosts: usize,
    dead_reservation_hosts: Vec<usize>,
    dead_wild_hosts: Vec<usize>,
    match_scores: HashMap<(usize, usize), usize>,
    match_scores_by_species: HashMap<usize, Vec<(usize, usize)>>
}

impl SimulationState {

    pub fn default() -> Self {
        SimulationState {
            current_generation: 0,
            count_dead_reservation_hosts: 0,
            count_dead_wild_hosts: 0,
            dead_reservation_hosts: Default::default(),
            dead_wild_hosts: Default::default(),
            match_scores: Default::default(),
            match_scores_by_species: Default::default(),
        }
    }

    pub fn update_match_store(&mut self, k: (usize, usize), v: usize) -> Option<usize> {
        if self.match_scores.contains_key(&k) {
            println!("{}", self.match_scores.get(&k).unwrap());
        }
        self.match_scores.insert(k, v)
    }
    pub fn current_generation(&self) -> usize {
        self.current_generation
    }
    pub fn count_dead_reservation_hosts(&self) -> usize {
        self.count_dead_reservation_hosts
    }
    pub fn count_dead_wild_hosts(&self) -> usize {
        self.count_dead_wild_hosts
    }
    pub fn dead_reservation_hosts(&self) -> &Vec<usize> {
        &self.dead_reservation_hosts
    }
    pub fn dead_wild_hosts(&self) -> &Vec<usize> {
        &self.dead_wild_hosts
    }
    pub fn match_scores(&self) -> &HashMap<(usize, usize), usize> {
        &self.match_scores
    }
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
        simulation_state: SimulationState::default(),
    }
}

impl Simulation {

    /**
    update parasite with parent parasite's number set
    */
    pub(crate) fn update_parasites(&mut self, species: usize, parasite: usize, parent: usize) {
        let base = self.parasites.index_axis(Axis(0), species);
        let parent_view = base.index_axis( Axis(0), parent).to_owned();
        for v in 0..parent_view.len() {
            self.parasites[[species, parasite, v]] = parent_view[v];
        }
    }

    pub(crate) fn update_parasites_exposed_to(&mut self, parasites_exposed_to: HashMap<(usize, usize), usize>) {
        self.simulation_state.match_scores = parasites_exposed_to;
    }

    pub(crate) fn parasites_exposed_to(&mut self) -> HashMap<(usize, usize), usize>{
        self.simulation_state.match_scores.clone()
    }

    pub fn count_dead_hosts(&self) -> (usize, usize, usize) {
        (
            self.simulation_state.count_dead_reservation_hosts,
            self.simulation_state.count_dead_wild_hosts,
            self.simulation_state.count_dead_wild_hosts + self.simulation_state.count_dead_reservation_hosts
        )
    }

    pub fn pref(&self) -> SimulationPref {
        self.pref.clone()
    }

    pub fn no_of_generations(&self) -> usize {
        self.no_of_generations
    }

    pub fn hosts(&self) -> &Array<Host, Ix1> {
        &self.hosts
    }

    pub fn update_dead_host(&mut self, index: usize, parent_host_index: usize) -> (usize, usize, usize) {
        let (host_type, number_set) = {
            let parent_host = &self.hosts[parent_host_index];
            let host_type = parent_host.host_type();
            match host_type {
                HostTypes::Reservation => self.simulation_state.count_dead_reservation_hosts -= 1,
                HostTypes::Wild => self.simulation_state.count_dead_wild_hosts -= 1,
            }
            (host_type, parent_host.number_set().clone())
        };
        self.hosts[index]
            .set_host_type(host_type)
            .set_number_set(number_set)
            .set_alive(true);
        self.count_dead_hosts()
    }

    pub fn kill_host(&mut self, index: usize) {
        let host = &mut self.hosts[index];
        host.set_alive(false);
        match host.host_type() {
            HostTypes::Reservation =>  {
                self.simulation_state.count_dead_reservation_hosts += 1;
                self.simulation_state.dead_reservation_hosts.push(index);
            },
            HostTypes::Wild => {
                self.simulation_state.count_dead_wild_hosts += 1;
                self.simulation_state.dead_wild_hosts.push(index);
            }
        }
    }

    pub fn hosts_mut(&mut self) -> &mut Array<Host, Ix1> {
        &mut self.hosts
    }

    pub fn parasites(&self) -> &Array<usize, Ix3> {
        &self.parasites
    }


    pub fn parasites_mut(&mut self) -> &mut Array<usize, Ix3> {
        &mut self.parasites
    }
    pub fn no_of_simulation_run(&self) -> usize {
        self.no_of_simulation_run
    }

    pub fn current_generation(&self) -> usize {
        self.simulation_state.current_generation
    }

    pub fn next_generation(&mut self) {
        self.simulation_state = SimulationState::default();
        self.simulation_state.current_generation += 1;
    }
    pub fn simulation_state(&self) -> &SimulationState {
        &self.simulation_state
    }

    pub fn set_simulation_state(&mut self, simulation_state: SimulationState) {
        self.simulation_state = simulation_state;
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
