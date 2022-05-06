use std::collections::HashMap;
use std::fmt::{Display, Formatter};
use std::fs::{create_dir_all, File, OpenOptions};
use std::io::Write;
use std::num::ParseIntError;

use ndarray::{Array, Array3, Axis, Ix, Ix1, Ix2, Ix3};
use ndarray_rand::RandomExt;
use rand::distributions::Uniform;
use rayon::prelude::IntoParallelRefIterator;

use crate::{HostTypes, SimulationPref};
use crate::hosts::{create_random_hosts, Host};

#[derive(Copy, Clone)]
pub enum ProgramVersions {
    One,
    Two,
    Three,
    Four,
}

pub type SpeciesParasite = (usize, usize);

impl Display for ProgramVersions {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            ProgramVersions::One => "Version 1",
            ProgramVersions::Two => "Version 2",
            ProgramVersions::Three => "Version 3",
            ProgramVersions::Four => "Version 4",
        };
        write!(f, "{}", s)
    }
}

impl From<String> for ProgramVersions {
    fn from(s: String) -> Self {
        let i = s.parse::<i32>();
        let v = match i {
            Ok(v) => v,
            Err(_) => 1
        };
        if v == 4 { ProgramVersions::Four } else if v == 3 { ProgramVersions::Three } else if v == 2 { ProgramVersions::Two } else { ProgramVersions::One }
    }
}

pub struct Simulation {
    pref: SimulationPref,
    simulation_state: SimulationState,
    generations: Vec<SimulationState>,
    program_version: ProgramVersions,
    gg: usize,
}

impl Simulation {
    pub(crate) fn parasites_possible(&self) -> Vec<Vec<usize>> {
        self.simulation_state.parasites_possible.clone()
    }

    pub(crate) fn set_parasites_possible(&mut self, parasites_possible: Vec<Vec<usize>>) {
        self.simulation_state.parasites_possible = parasites_possible;
    }
}

impl Simulation {
}


#[derive(Clone, Debug)]
pub struct SimulationState {
    /** parasite scores, key is (Species, Parasite) */
    match_scores: HashMap<SpeciesParasite, usize>,
    /** */
    species_match_score: HashMap<usize, usize>,
    current_generation: usize,
    hosts: Array<Host, Ix1>,
    parasites: Array<usize, Ix3>,
    host_match_score: HashMap<usize, usize>,
    species_left: HashMap<usize, Vec<usize>>,
    parasites_possible: Vec<Vec<usize>>
}

impl Default for SimulationState {
    fn default() -> Self {
        SimulationState {
            match_scores: Default::default(),
            species_match_score: Default::default(),
            current_generation: 0,
            hosts: Default::default(),
            parasites: Default::default(),
            host_match_score: Default::default(),
            species_left: Default::default(),
            parasites_possible: vec![]
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

    /** parasite scores, key is (Species, Parasite) => Match score */
    pub fn match_scores(&self) -> &HashMap<(usize, usize), usize> {
        &self.match_scores
    }

    pub fn host_match_score(&self) -> &HashMap<usize, usize> {
        &self.host_match_score
    }

}

pub fn new_simulation(pref: SimulationPref, program_version: ProgramVersions, gg: usize) -> Simulation {
    let hosts = create_random_hosts(&pref);
    let parasites: Array<usize, Ix3> = Array::random(
        (pref.d(), pref.e(), pref.g()),
        Uniform::new(0, pref.f()),
    );
    // create log folder
    if gg < 3 {
        let g_folder = format!("report/sim_{}", gg);
        create_dir_all(&g_folder);
        let mut f = File::create(format!("{}/hosts", g_folder)).expect("Unable to create file");
        f.write_all(&format!("{:#?}", hosts).to_string().as_bytes()).expect("Unable to write data");
        let mut f = File::create(format!("{}/parasites", g_folder)).expect("Unable to create file");
        f.write_all(print_parasites(&parasites).as_bytes()).expect("Unable to write data");
        for ff in 0..3 {
            create_dir_all(format!("{}/{}", g_folder, ff));
        }
    }
    Simulation {
        pref: pref.clone(),
        simulation_state: SimulationState {
            hosts,
            parasites,
            ..SimulationState::default()
        },
        generations: vec![],
        program_version,
        gg,
    }
}

impl Simulation {
    pub fn set_species_left(&mut self, species_left: HashMap<usize, Vec<usize>>) {
        self.simulation_state.species_left = species_left;
    }

    pub(crate) fn species_left(&self) -> HashMap<usize, Vec<usize>> {
        self.simulation_state.species_left.clone()
    }

    pub(crate) fn species_left_mut(&mut self) -> &mut HashMap<usize, Vec<usize>> {
        &mut self.simulation_state.species_left
    }

    pub(crate) fn set_host_match_score(&mut self, host_match_score: HashMap<usize, usize>) {
        self.simulation_state.host_match_score = host_match_score;
    }

    pub(crate) fn update_host_match_score(&mut self, host_index: usize, inc: usize) {
        *self.simulation_state.host_match_score.entry(host_index).or_insert(0) += inc;
    }

    pub(crate) fn species_match_score(&self) -> &HashMap<usize, usize> {
        &self.simulation_state.species_match_score
    }
    /** Updates species total match score (species index) => match score) */
    pub(crate) fn set_species_match_score(&mut self, species_match_score: HashMap<usize, usize>) {
        self.simulation_state.species_match_score = species_match_score;
    }

    pub fn update_species_match_score(&mut self, species_index: usize, match_score: usize) {
        *self.simulation_state.species_match_score.entry(species_index).or_insert(0) += match_score;
    }

    pub fn pv(&self, file_name: &str, content: &str, append: bool) {
        if self.gg() < 3 && self.current_generation() < 3 {
            let f_folder = format!("report/sim_{}/{}", self.gg(), self.current_generation());
            let file_path = format!("{}/{}.txt", f_folder, file_name);
            let mut f = if append {
                OpenOptions::new()
                    .append(true)
                    .create(true) // Optionally create the file if it doesn't already exist
                    .open(file_path)
                    .expect("Unable to open file")
            } else {
                File::create(file_path).expect("Unable to create file")
            };
            f.write_all(content.as_bytes()).expect("unable to write data");
        }
    }
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

    /** set parasites scores, key is (Species, Parasite) => Match score */
    pub(crate) fn set_parasites_exposed_to(&mut self, parasites_exposed_to: HashMap<SpeciesParasite, usize>) {
        self.simulation_state.match_scores = parasites_exposed_to;
    }

    /** update parasites scores, key is (Species, Parasite) => Match score */
    pub(crate) fn update_parasites_exposed_to(&mut self, key: SpeciesParasite, parasite: usize) {
        self.simulation_state.match_scores.insert(key, parasite);
    }

    /** get a clone of parasites scores, key is (Species, Parasite) => Match score */
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
        simulation_state.hosts.for_each(|host| {
            if host.alive() {
                if host.host_type() == HostTypes::Reservation {
                    count_hosts_alive_reservation += 1;
                } else {
                    count_hosts_alive_wild += 1;
                }
            }
        });
        (count_hosts_alive_reservation + count_hosts_alive_wild, count_hosts_alive_reservation, count_hosts_alive_wild)
    }

    pub fn count_dead_hosts(&mut self) -> (usize, usize, usize) {
        let mut count_hosts_dead_reservation = 0;
        let mut count_hosts_dead_wild = 0;
        self.simulation_state.hosts.for_each(|host| {
            if !host.alive() {
                if host.host_type() == HostTypes::Reservation {
                    count_hosts_dead_reservation += 1;
                } else {
                    count_hosts_dead_wild += 1;
                }
            }
        });
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
        let generation = self.simulation_state.current_generation + 1;
        self.simulation_state = SimulationState {
            hosts: self.simulation_state.hosts.to_owned(),
            parasites: self.simulation_state.parasites.to_owned(),
            current_generation: generation,
            ..Default::default()
        };
    }

    pub fn simulation_state(&self) -> &SimulationState {
        &self.simulation_state
    }

    pub fn set_simulation_state(&mut self, simulation_state: SimulationState) {
        self.simulation_state = simulation_state;
    }

    pub fn generate_report(&self) -> HostsCount {
        let (_, r, w) = self.count_alive_hosts_from_generation(self.generations.last().unwrap());
        HostsCount {
            wild_host: w,
            reservation_host: r,
        }
    }
    pub fn program_version(&self) -> ProgramVersions {
        self.program_version.clone()
    }

    pub fn gg(&self) -> usize {
        self.gg
    }
}

#[derive(Debug, Clone)]
pub struct HostsCount {
    pub wild_host: usize,
    pub reservation_host: usize,
}

impl Display for HostsCount {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "w: {} r: {}", self.wild_host, self.reservation_host)
    }
}

#[derive(Debug)]
pub struct GGRunReport {
    pub hosts_count: HostsCount,
    pub wild_loner: usize,
    pub reservation_loner: usize,
    pub high_wild: usize,
    pub high_reservation: usize,
    pub tied: usize,
}

impl Display for GGRunReport {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let mut s = String::new();
        s.push_str(&format!("- {} runs ended with wild individuals the lone type remaining", self.wild_loner));
        s.push_str(&format!("- {} runs ended with reservation individuals the lone type remaining", self.reservation_loner));
        s.push_str(&format!("- {} runs ended with wild individuals a higher quantity than reservation individuals", self.high_wild));
        s.push_str(&format!("- {} runs ended with reservation individuals a higher quantity than wild individuals", self.high_reservation));
        s.push_str(&format!("- {} runs ended with the quantities of the two types tied", self.tied));
        write!(f, "{}", s)
    }
}

pub fn print_parasites(all_parasites: &Array3<usize>) -> String {
    let mut _s = String::new();
    let mut i = 0;
    let l1 = all_parasites.len_of(Axis(0));
    let l2 = all_parasites.len_of(Axis(1));
    for ii in 0..l1 {
        _s.push_str(&format!("Species: {}\n", ii));
        for jj in 0..l2 {
            let v = all_parasites.index_axis(Axis(0), ii);
            _s.push_str(&format!("({}, {}) {}\n", ii, jj, v.index_axis(Axis(0), jj)));
        }
        _s.push_str("\n");
    }
    format!("PARASITES:\n{}", _s)
}