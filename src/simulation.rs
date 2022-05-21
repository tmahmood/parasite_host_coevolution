use std::collections::HashMap;
use std::fmt::{Display, Formatter};
use std::fs::{create_dir_all, File, OpenOptions};
use std::io::Write;
use std::ops::Div;

use ndarray::{Array, Array1, Array3, Axis, Ix1, Ix3};

use crate::{generate_individual, HostTypes, SimulationPref};
use crate::hosts::{create_random_hosts, Host, print_hosts};

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
    log_files: HashMap<String, String>,
}

impl Simulation {
    pub(crate) fn additional_exposure(&self) -> bool {
        self.simulation_state.additional_exposure
    }

    pub(crate) fn has_additional_exposure(&mut self) {
        self.simulation_state.additional_exposure = true;
    }

    pub(crate) fn set_additional_exposure(&mut self, s: bool) {
        self.simulation_state.additional_exposure = s;
    }
}

impl Simulation {
    pub(crate) fn parasites_possible(&self) -> Vec<Vec<usize>> {
        self.simulation_state.parasites_possible.clone()
    }

    pub(crate) fn set_parasites_possible(&mut self, parasites_possible: Vec<Vec<usize>>) {
        self.simulation_state.parasites_possible = parasites_possible;
    }
}

impl Simulation {}


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
    parasites_possible: Vec<Vec<usize>>,
    additional_exposure: bool,
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
            parasites_possible: vec![],
            additional_exposure: false
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

pub fn create_random_parasites(pref: &SimulationPref) -> Array3<usize> {
    let mut v = vec![];
    for _ in 0..pref.d() {
        v.push(generate_individual(pref.f(), pref.g()));
    }
    Array::from_shape_fn([pref.d(), pref.e(), pref.g()], |(i, _, k)| {
        v[i][k]
    })
}

pub fn new_simulation(pref: SimulationPref, program_version: ProgramVersions, gg: usize) -> Simulation {
    let hosts = create_random_hosts(&pref);
    let parasites = create_random_parasites(&pref);
    // create log folder
    if gg < 3 {
        let g_folder = format!("report/sim_{}", gg);
        create_dir_all(&g_folder).expect("Failed to create directory");
        let mut f = File::create(format!("{}/hosts", g_folder)).expect("Unable to create file");
        let s = print_hosts(&hosts);
        f.write_all(s.as_bytes()).expect("Unable to write data");
        // f.write_all(&format!("{:#?}", hosts).to_string().as_bytes()).expect("Unable to write data");
        let mut f = File::create(format!("{}/parasites", g_folder)).expect("Unable to create file");
        f.write_all(print_parasites(&parasites).as_bytes()).expect("Unable to write data");
        for ff in 0..3 {
            create_dir_all(format!("{}/{}", g_folder, ff)).expect("Failed to create directory");
        }
    }
    let initial_simulation_state = SimulationState {
        hosts,
        parasites,
        ..SimulationState::default()
    };
    Simulation {
        pref: pref.clone(),
        simulation_state: initial_simulation_state.clone(),
        generations: vec![initial_simulation_state],
        program_version,
        gg,
        log_files: HashMap::new(),
    }
}

impl Simulation {
    pub fn generations(&self) -> &Vec<SimulationState> {
        &self.generations
    }

    pub fn set_species_left(&mut self, species_left: HashMap<usize, Vec<usize>>) {
        self.simulation_state.species_left = species_left;
    }

    pub(crate) fn species_left(&self) -> HashMap<usize, Vec<usize>> {
        self.simulation_state.species_left.clone()
    }

    pub(crate) fn update_host_match_score(&mut self, host_index: usize, inc: usize) {
        *self.simulation_state.host_match_score.entry(host_index).or_insert(0) += inc;
    }

    pub(crate) fn species_match_score(&self) -> &HashMap<usize, usize> {
        &self.simulation_state.species_match_score
    }

    pub fn update_species_match_score(&mut self, species_index: usize, match_score: usize) {
        *self.simulation_state.species_match_score.entry(species_index).or_insert(0) += match_score;
    }

    pub fn pp(&self, file_name: &str, content: &str, txt: bool) {
        let file_path = self.get_file_path(file_name, txt, false);
        let mut f = OpenOptions::new()
            .append(true)
            .create(true) // Optionally create the file if it doesn't already exist
            .open(file_path)
            .expect("Unable to open file");
        f.write_all(content.as_bytes()).expect("unable to write data");
    }

    pub fn write_all(&self) {
        self.log_files.iter().for_each(|(file_path, content)| {
            let mut f = OpenOptions::new()
                .append(true)
                .create(true) // Optionally create the file if it doesn't already exist
                .open(file_path)
                .expect("Unable to open file");
            f.write_all(content.as_bytes()).expect("unable to write data");
        });
    }

    pub fn pv(&mut self, file_name: &str, content: &str, txt: bool) {
        if self.gg() < 3 && self.current_generation() < 3 {
            let file_path = self.get_file_path(file_name, txt, true);
            self.log_files.entry(file_path.clone()).or_insert(String::new()).push_str(content);
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


    /** update parasites scores, key is (Species, Parasite) => Match score */
    pub(crate) fn update_parasites_exposed_to(&mut self, key: SpeciesParasite, parasite: usize) {
        self.simulation_state.match_scores.insert(key, parasite);
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

    pub fn next_generation(&mut self) -> HostsCount {
        self.generations.push(self.simulation_state.clone());
        let generation = self.simulation_state.current_generation + 1;
        self.simulation_state = SimulationState {
            hosts: self.simulation_state.hosts.to_owned(),
            parasites: self.simulation_state.parasites.to_owned(),
            current_generation: generation,
            ..Default::default()
        };
        self.get_hosts_count()
    }

    pub fn simulation_state(&self) -> &SimulationState {
        &self.simulation_state
    }

    pub fn set_simulation_state(&mut self, simulation_state: SimulationState) {
        self.simulation_state = simulation_state;
    }


    pub fn get_hosts_count(&self) -> HostsCount {
        let (_, r, w) = self.count_alive_hosts_from_generation(
            self.generations.last().unwrap()
        );
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

    fn get_file_path(&self, file_name: &str, txt: bool, inc_gen: bool) -> String {
        let f_folder = if inc_gen {
            format!("report/sim_{}/{}", self.gg(), self.current_generation())
        } else {
            format!("report/sim_{}", self.gg())
        };
        let extension = if txt {
            "txt"
        } else {
            "csv"
        };
        format!("{}/{}.{}", f_folder, file_name, extension)
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

pub struct ReportHostType {
    pub hosts_count: Array1<usize>,
    total_host: usize,
    mean: f32,
    standard_deviation: f32,
    standard_error: f32,
    confidence_interval_high_point: f32,
    confidence_interval_low_point: f32,
}

impl ReportHostType {
    pub fn new(hosts_count: Array1<usize>) -> Self {
        ReportHostType {
            total_host: hosts_count.sum(),
            hosts_count,
            mean: 0.0,
            standard_deviation: 0.0,
            standard_error: 0.0,
            confidence_interval_high_point: 0.0,
            confidence_interval_low_point: 0.0,
        }
    }

    pub fn calculate(&mut self, pref: SimulationPref) {
        self.mean = self.total_host as f32 / self.hosts_count.len() as f32;
        // Sqrt[ ( (100-100.4)^2 + (80-100.4)^2 + (92-100.4)^2 + (110-100.4)^2 + (120-100.4)^2 ) /4
        //   ] = 15.5177
        self.standard_deviation = self.hosts_count.iter().map(|v| {
            (*v as f32 - self.mean).powf(2.)
        }).sum::<f32>()
            .div(pref.gg() as f32 - 1.)
            .sqrt();
        self.standard_error = self.standard_deviation / (pref.gg() as f32).sqrt();
        self.confidence_interval_high_point = self.mean + pref.z() * self.standard_error;
        self.confidence_interval_low_point = self.mean - pref.z() * self.standard_error;
    }
    pub fn total_host(&self) -> usize {
        self.total_host
    }
    pub fn mean(&self) -> f32 {
        self.mean
    }
    pub fn standard_deviation(&self) -> f32 {
        self.standard_deviation
    }
    pub fn standard_error(&self) -> f32 {
        self.standard_error
    }
    pub fn confidence_interval_high_point(&self) -> f32 {
        self.confidence_interval_high_point
    }
    pub fn confidence_interval_low_point(&self) -> f32 {
        self.confidence_interval_low_point
    }
}

impl Display for ReportHostType {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let mut _s = String::new();
        _s.push_str(&format!("standard deviation:\n        {:10.3}\n", self.standard_deviation));
        _s.push_str(&format!("Confidence Interval:\n {: >24} {: >11}{: >12}\n", "Means", "High Point", "Low Point"));
        _s.push_str(&format!("             {:12.3}{:12.3}{:12.3}\n", self.mean, self.confidence_interval_high_point, self.confidence_interval_low_point));
        write!(f, "{}", _s)
    }
}

#[derive(Debug)]
pub struct GGRunReport {
    pub hosts: (Array1<usize>, Array1<usize>),
    simulations: usize,
    wild_loner: usize,
    reservation_loner: usize,
    high_wild: usize,
    high_reservation: usize,
    tied: usize,
}

impl GGRunReport {
    pub fn new(hosts: (Array1<usize>, Array1<usize>), ff: usize) -> Self {
        GGRunReport {
            hosts,
            simulations: ff,
            wild_loner: 0,
            reservation_loner: 0,
            high_wild: 0,
            high_reservation: 0,
            tied: 0,
        }
    }

    pub fn calculations(&mut self) {
        for i in 0..self.simulations {
            let r = self.hosts.0[i];
            let w = self.hosts.1[i];
            if r > w { self.high_reservation += 1 }
            if r < w { self.high_wild += 1 }
            if r == w { self.tied += 1 }
            if r == 0 { self.wild_loner += 1 }
            if w == 0 { self.reservation_loner += 1 }
        }
    }
}

impl Display for GGRunReport {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let mut _s = String::new();
        _s.push_str(&format!("- {} simulation runs ended with wild individuals the lone type remaining\n", self.wild_loner));
        _s.push_str(&format!("- {} simulation runs ended with reservation individuals the lone type remaining\n", self.reservation_loner));
        _s.push_str(&format!("- {} simulation runs ended with wild individuals a higher quantity than reservation individuals\n", self.high_wild));
        _s.push_str(&format!("- {} simulation runs ended with reservation individuals a higher quantity than wild individuals\n", self.high_reservation));
        _s.push_str(&format!("- {} simulation runs ended with the quantities of the two types tied\n", self.tied));
        write!(f, "{}", _s)
    }
}

pub fn print_parasites(all_parasites: &Array3<usize>) -> String {
    let mut _s = String::new();
    let l1 = all_parasites.len_of(Axis(0));
    let l2 = all_parasites.len_of(Axis(1));
    for ii in 0..l1 {
        _s.push_str(&format!("Species: {}\n", ii));
        for jj in 0..l2 {
            let v = all_parasites.index_axis(Axis(0), ii);
            _s.push_str(&format!("({: >3}, {: >3}) {}\n", ii, jj, v.index_axis(Axis(0), jj)));
        }
        _s.push_str("\n");
    }
    format!("PARASITES:\n{}", _s)
}