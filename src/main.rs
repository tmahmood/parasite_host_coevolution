#![feature(iter_collect_into)]
#![feature(map_first_last)]
#![feature(drain_filter)]
extern crate core;
extern crate rand;
extern crate serde;
extern crate serde_derive;
extern crate serde_ini;

use std::{fs, time};
use std::cmp::min;
use std::collections::HashMap;
use std::env::args;
use std::fmt::{Display, Formatter};
use std::fs::{File, remove_dir_all};
use std::io::Write;
use std::path::Path;

use indicatif::{MultiProgress, ParallelProgressIterator, ProgressBar, ProgressStyle};
use log;
use log::{info, LevelFilter, SetLoggerError};
use log4rs::{Config, Handle};
use log4rs::append::console::ConsoleAppender;
use log4rs::append::file::FileAppender;
use log4rs::config::{Appender, Root};
use log4rs::encode::pattern::PatternEncoder;
use ndarray::{Array, Array1, Array2, Array3, ArrayBase, Axis, Ix2, OwnedRepr};
use ndarray_rand::RandomExt;
use rand::{Rng, thread_rng};
use rand::distributions::{Uniform, WeightedIndex};
use rand::prelude::*;
use rayon::prelude::*;

use crate::hosts::{Host, HostTypes, print_hosts};
use crate::mutations::{mutate_hosts, mutate_parasites};
use crate::simulation::{GGRunReport, HostsCount, new_simulation, print_parasites, ProgramVersions, ReportHostType, Simulation};
use crate::simulation_pref::SimulationPref;

pub mod simulation_pref;
pub mod hosts;
pub mod simulation;
pub mod mutations;

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
}

impl Display for ParasiteSpeciesIndex {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "S{: <3} P{: <3} M{: <3} ", self.species_index, self.parasite_index, self.match_count)
    }
}

fn main() {
    let multi_progress_bar = MultiProgress::new();
    let sty = ProgressStyle::with_template(
        "[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7} {msg}",
    ).unwrap().progress_chars("##-");
    if Path::new("report").exists() {
        remove_dir_all("report").expect("Failed to remove report directory");
    };
    let program = ProgramVersions::from(args().nth(2).unwrap_or(format!("1")));
    config_logger(1).unwrap();
    let param_file = args().nth(1).unwrap_or(format!("conf/params.conf"));
    let pref: SimulationPref = serde_ini::from_str(&fs::read_to_string(param_file).unwrap()).unwrap();
    let pb_all_simulation = multi_progress_bar.add(ProgressBar::new(pref.gg() as u64));
    pb_all_simulation.set_style(sty.clone());
    let pb_generations = multi_progress_bar.insert_after(&pb_all_simulation, ProgressBar::new(pref.ff() as u64));
    pb_generations.set_style(sty.clone());
    let now = time::Instant::now();
    println!("Running version {}, build 0.1.36_wrong_sort_when_doing_parasite_replacement", program);
    let program_clone = program.clone();
    let pref_clone = pref.clone();
    let mut wild_hosts: Vec<Vec<usize>> = vec![];
    let mut reservation_hosts: Vec<Vec<usize>> = vec![];
    if option_env!("DEBUG").is_some() {
        (0..pref.gg()).into_iter().for_each(|gg| {
            let k = run_generation_step(program_clone, gg.clone(), pref_clone);
            wild_hosts.push(k.0);
            reservation_hosts.push(k.1);
        });
    } else {
        (0..pref.gg()).into_par_iter().progress_with(pb_all_simulation).map(move |gg| {
            run_generation_step(program_clone, gg.clone(), pref_clone)
        }).unzip_into_vecs(
            &mut wild_hosts,
            &mut reservation_hosts,
        );
    }

    let flat_wild_hosts: Vec<usize> = wild_hosts.iter().flatten().cloned().collect();
    let flat_reservation_hosts: Vec<usize> = reservation_hosts.iter().flatten().cloned().collect();
    let reservation_hosts_ar = Array2::from_shape_vec((pref.gg(), pref.ff()), flat_reservation_hosts).unwrap();
    let wild_hosts_ar = Array2::from_shape_vec((pref.gg(), pref.ff()), flat_wild_hosts).unwrap();
    //
    for k in 0..wild_hosts.len() {
        info!("\n{}\nr: {: >3?}\nw: {: >3?}\n", k,
                 reservation_hosts_ar.index_axis(Axis(0), k).as_slice().unwrap(),
                 wild_hosts_ar.index_axis(Axis(0), k).as_slice().unwrap());
    }
    //
    let mut f = File::create("report/confidence.csv").expect("Unable to create file");
    f.write_all("Generation, Standard Deviation (R), Means (R), High Point (R), Low Point (R), \
    Standard Deviation (W), Means (W), High Point (W), Low Point (W)\n".as_bytes())
        .expect("Failed to write to file");
    //
    let l1 = generate_excel(&reservation_hosts_ar, &pref);
    let l2 = generate_excel(&wild_hosts_ar, &pref);
    for ii in 0..pref.ff() {
        f.write_all(format!("{}, {}, {}\n", ii, l1[ii], l2[ii]).as_bytes()).expect("Failed to write row");
    }
    let last_row_r = reservation_hosts_ar.index_axis(Axis(1), pref.ff() - 1).to_owned();
    let last_row_w = wild_hosts_ar.index_axis(Axis(1), pref.ff() - 1).to_owned();
    let r = calculate_result((last_row_r, last_row_w), pref.gg());
    println!("{}", r);
    info!("\n{}", r);
    println!("took {} secs", now.elapsed().as_secs_f32())
}

fn generate_excel(hosts: &ArrayBase<OwnedRepr<usize>, Ix2>, pref: &SimulationPref) -> Vec<String> {
    hosts.columns().into_iter().map(|v| {
        let mut r = ReportHostType::new(v.to_owned());
        r.calculate(pref.clone());
        format!("{},{},{},{}",
                r.standard_deviation(),
                r.mean(),
                r.confidence_interval_high_point(),
                r.confidence_interval_low_point()
        )
    }).collect()
}

fn calculate_result(result: (Array1<usize>, Array1<usize>), ff: usize) -> GGRunReport {
    let mut report_gen = GGRunReport::new(result, ff);
    report_gen.calculations();
    report_gen
}

fn fill_host(simulation: &mut Simulation) -> HostsCount {
    let mut k = simulation.next_generation();
    if k.wild_host == 0 && k.reservation_host == 0 {
        k.reservation_host = 0;
        k.wild_host = 0;
    } else if k.wild_host == 0 {
        k.reservation_host = simulation.pref().a() + simulation.pref().b();
    } else {
        k.wild_host = simulation.pref().a() + simulation.pref().b();
    }
    return k;
}

fn run_generation_step(program: ProgramVersions, gg: usize, pref: SimulationPref) -> (Vec<usize>, Vec<usize>) {
    let mut simulation = new_simulation(pref.clone(), program, gg);
    let mut lines = vec![];
    let mut simulation_ended = false;
    (0..pref.ff()).into_iter().map(|_| {
        if simulation_ended {
            return fill_host(&mut simulation);
        }
        expose_all_hosts_to_parasites(&mut simulation);
        additional_exposure(&mut simulation);
        if !should_continue(&mut simulation) {
            simulation_ended = true;
            return fill_host(&mut simulation);
        }
        birth_hosts(&mut simulation);
        parasite_truncation_and_birth(&mut simulation);
        mutation(&mut simulation);
        parasite_replacement(&mut simulation);
        let _s = print_parasites(&simulation.parasites());
        simulation.pv("parasites_at_end", &_s, true);
        let _s = print_hosts(simulation.hosts());
        simulation.pv("hosts_at_end", &_s, true);
        simulation.next_generation()
    }).collect_into(&mut lines);
    simulation.write_all();
    let mut wild = vec![];
    let mut reservation = vec![];
    lines.clone()
        .into_par_iter()
        .map(|v| (v.wild_host, v.reservation_host))
        .unzip_into_vecs(&mut wild, &mut reservation);
    (wild, reservation)
}


pub fn generate_individual(f: usize, len: usize) -> Array1<usize> {
    Array::random(len, Uniform::new(0, f))
}


pub fn expose_all_hosts_to_parasites(simulation: &mut Simulation) {
    let mut rng = thread_rng();
    let file_name = "host_exposed_to";
    simulation.pv(file_name, "Initial Exposure\n", true);
    let all_parasites = simulation.parasites().clone();
    let hosts = simulation.hosts().clone();
    // we are calculating all the random species for each hosts beforehand, that way we don't
    // have to calculate it in the loop and most likely it will be optimized, we remove the species
    // that is used, that way we don't have to check for used parasite individual
    let mut species_possible = HashMap::new();
    for i in 0..hosts.len() {
        let mut k = (0..simulation.pref().d()).collect::<Vec<usize>>();
        k.shuffle(&mut rng);
        species_possible.insert(i, k);
    }
    // same with parasites
    let mut parasites_possible = vec![];
    for _ in 0..simulation.pref().d() {
        let mut k = (0..simulation.pref().e()).collect::<Vec<usize>>();
        k.shuffle(&mut rng);
        parasites_possible.push(k);
    }
    for host_index in 0..simulation.pref().a() + simulation.pref().b() {
        let host = &hosts[host_index];
        simulation.pv(file_name, &format!("{: <3} {:12}{}\n(species, parasite): match count -> code, \n", host_index, &host.host_type().to_string(), &host.number_set()), true);
        let mut match_score_bellow_threshold = 0;
        let species = species_possible.get_mut(&host_index).unwrap();
        for _ in 0..simulation.pref().h() {
            let species_index = species.pop().unwrap();
            let parasite_index = parasites_possible[species_index].pop().unwrap();
            let mut p_idx = ParasiteSpeciesIndex {
                species_index,
                parasite_index,
                match_count: 0,
            };
            let match_score = find_match_score(&host, &all_parasites, &mut p_idx, simulation);
            simulation.update_parasites_exposed_to((p_idx.species(), p_idx.parasite()), match_score);
            simulation.update_species_match_score(p_idx.species(), match_score);
            simulation.update_host_match_score(host_index, 1);
            if match_score < simulation.pref().n() {
                match_score_bellow_threshold += 1;
            }
            simulation.update_host_match_score_bellow_dd(host_index, if match_score < simulation.pref().dd() { 1 } else { 0 });
            simulation.update_host_match_score_bellow_j(host_index, if match_score < simulation.pref().j() { 1 } else { 0 });
            //
            let d = parasite_row(&all_parasites, p_idx.species(), p_idx.parasite_index);
            let p_grid = print_matching_number_sets(
                d,
                p_idx.species(),
            );
            simulation.pv(file_name, &format!("{: >3},{: >3} : {: >2} -> {: >3}\n", p_idx.species(), p_idx.parasite(), match_score, p_grid), true);
            //
        }
        simulation.pv(file_name, &format!("-------------------------\n"), true);
        if match_score_bellow_threshold >= simulation.pref().x() {
            simulation.kill_host(host_index);
            simulation.pv("host_dying_initial_exposure", &format!("{:3} {}\n", host_index, &simulation.hosts()[host_index].to_string()), true);
        }
    }
    let (t, r, w) = simulation.count_dead_hosts();
    simulation.pv("host_dying_initial_exposure", &format!("{} R({}) W({})", t, r, w), true);
    simulation.set_species_left(species_possible);
    simulation.set_parasites_possible(parasites_possible);
}

pub fn additional_exposure(simulation: &mut Simulation) {
    let file_name = "host_additional_exposure";
    let (total_dead_hosts, _, _) = simulation.count_dead_hosts();
    simulation.pv(file_name, "Additional Exposure\n", true);
    // secondary exposure
    if !additional_exposure_selected(simulation, total_dead_hosts) {
        return;
    }
    simulation.has_additional_exposure();
    simulation.pv(file_name, &format!("current generation: {}\n", simulation.current_generation()), true);
    let all_parasites = simulation.parasites().clone();
    let (_, reservation_hosts_alive, _) = simulation.count_alive_hosts();
    let mut rng = thread_rng();
    let mut hosts_to_try = 0;
    let mut species_par_host = simulation.species_left();
    let mut parasites_possible = simulation.parasites_possible();
    let mut alive_hosts: Vec<usize> = simulation.hosts_alive().drain_filter(|v| {
        simulation.host_type(*v) == HostTypes::Reservation
    }).collect();
    alive_hosts.shuffle(&mut rng);
    let mut no_of_additional_host = (reservation_hosts_alive as f32 * simulation.pref().aa()).ceil() as usize;
    no_of_additional_host = min(alive_hosts.len(), no_of_additional_host);
    simulation.pv(file_name, &format!("Additional Exposure candidate {}\n", no_of_additional_host), true);
    while hosts_to_try < no_of_additional_host {
        let host_index: usize = alive_hosts.pop().unwrap();
        let host = simulation.hosts()[host_index].clone();
        simulation.pv(file_name, &format!("{: <3} {:12}{}\n(species, parasite): match count -> code, \n", host_index, &host.host_type().to_string(), &host.number_set()), true);
        hosts_to_try += 1;
        // expose to parasite
        for _ in 0..simulation.pref().i() {
            let species_index = species_par_host.get_mut(&host_index).unwrap().pop().unwrap();
            let parasite_index = parasites_possible[species_index].pop().unwrap();
            let mut p_idx = ParasiteSpeciesIndex {
                species_index,
                parasite_index,
                match_count: 0,
            };
            let match_score = find_match_score(&host, &all_parasites, &mut p_idx, &simulation);
            // update simulation state
            simulation.update_parasites_exposed_to((p_idx.species(), p_idx.parasite()), match_score);
            simulation.update_species_match_score(p_idx.species(), match_score);
            simulation.update_host_match_score(host_index, 1);
            simulation.update_host_match_score_bellow_dd(host_index, if match_score < simulation.pref().dd() { 1 } else { 0 });
            simulation.update_host_match_score_bellow_j(host_index, if match_score < simulation.pref().j() { 1 } else { 0 });
            let d = parasite_row(&all_parasites, p_idx.species(), p_idx.parasite_index);
            let p_grid = print_matching_number_sets(d, p_idx.species());
            simulation.pv(file_name, &format!("{: >3},{: >3} : {: >2} -> {: >3}\n", p_idx.species(), p_idx.parasite(), match_score, p_grid), true);
            //
        }
        simulation.pv(file_name, &format!("-------------------------\n"), true);
        if *simulation.ss().host_match_scores_bellow_dd().get(&host_index).unwrap() >= simulation.pref().cc() {
            simulation.kill_host(host_index);
            simulation.pv("host_dying_additional_exposure", &format!("{:3} {}\n", host_index, &host.to_string()), true);
        } else {
            simulation.ss_mut().add_hosts_tried(host_index);
        }
    }
    let (t, r, w) = simulation.count_dead_hosts();
    simulation.pv("host_dying_additional_exposure", &format!("{} R({}) W({})", t, r, w), true);
}

/**
If
  1. the current generation is after the Lth generation and if
  2. less than an M fraction of host individuals (a total of reservation and wild host individuals)
     have been killed
 */
fn additional_exposure_selected(simulation: &mut Simulation, total_dead_hosts: usize) -> bool {
    // is the current generation is after the Lth generation
    let k = simulation.current_generation() as i32 > simulation.pref().l();
    // less than an M fraction of host individuals (a total of reservation and wild host
    // individuals) have been killed
    let m_fraction = (simulation.pref().a() + simulation.pref().b()) as f32 * simulation.pref().m();
    let l = total_dead_hosts < m_fraction as usize;
    simulation.pv("additional_exposure_cond_check",
                  &format!("gen({}), l({}), dead_hosts({}) m fraction({})\n",
                           simulation.current_generation(),
                           simulation.pref().l(),
                           total_dead_hosts,
                           m_fraction),
                  true);
    if k && l {
        simulation.pv("additional_exposure_cond_check",
                      "additional exposure selected", true);
        return true;
    }
    simulation.pv("additional_exposure_cond_check",
                  "additional exposure not selected", true);
    return false;
}

pub fn birth_hosts(simulation: &mut Simulation) {
    let file_name = "hosts_birth";
    let (dist, choices) = match simulation.program_version() {
        ProgramVersions::One | ProgramVersions::Three => {
            birth_generation_version_1(simulation)
        }
        ProgramVersions::Two | ProgramVersions::Four => {
            calculate_qi(simulation);
            birth_generation_version_2(simulation)
        }
    };
    let mut rng = thread_rng();
    simulation.pv(file_name,
                  &format!("{}\n{: >4} {: <width$} {}\n",
                           simulation.program_version().to_string(),
                           "index", "host", "target",
                           width = simulation.pref().c() + 64),
                  true);
    simulation.hosts().clone().into_iter().enumerate()
        .filter(|(_, host)| !host.alive())
        .map(|v| v.0)
        .collect::<Vec<usize>>().into_iter()
        .for_each(|host_index| {
            let parent_index = choices[dist.sample(&mut rng)];
            let parent_host = simulation.hosts()[parent_index].clone();
            simulation.pv(file_name, &format!(" {: >4} {: >width$} {:4}\n",
                                              parent_index,
                                              parent_host, host_index,
                                              width = simulation.pref().c() + 13), true);
            simulation.update_dead_host(host_index, parent_index);
        });
}

pub fn birth_generation_version_1(simulation: &mut Simulation) -> (WeightedIndex<f32>, Vec<usize>) {
    let (_, no_of_dead_reservation_host, no_of_dead_wild_host) = simulation.count_dead_hosts();
    let (_, no_of_reservation_host_alive, no_of_wild_host_alive) = simulation.count_alive_hosts();
    let (chance_reservation_p_0, _w) = get_chances_v1(
        no_of_dead_wild_host as f32,
        no_of_reservation_host_alive as f32,
        no_of_dead_reservation_host as f32,
        simulation.pref().y(),
        no_of_wild_host_alive as f32,
        simulation.pref().o(),
        0.,
    );
    let (chance_reservation_p, chance_wild) = get_chances_v1(
        no_of_dead_wild_host as f32,
        no_of_reservation_host_alive as f32,
        no_of_dead_reservation_host as f32,
        simulation.pref().y(),
        no_of_wild_host_alive as f32,
        simulation.pref().o(),
        simulation.pref().p(),
    );
    // select index of the hosts that are alive
    simulation.pv("hosts_birth", "Birth V1\n", true);
    let choices: Vec<usize> = simulation.hosts_alive();
    // chances of each hosts that are alive for birth
    simulation.pv("hosts_birth", &format!("host_index,     host_type,    p,    chance\n"), true);
    let pv = |hi, ht, p, c| -> String {
        format!("{:10},{: >14},{:5},{:10}\n", hi, ht, p, c)
    };

    let chances: Vec<f32> = choices.iter().map(|host_index| {
        match simulation.hosts()[*host_index].host_type() {
            HostTypes::Reservation => {
                return if simulation.ss().hosts_tried().contains(host_index) {
                    simulation.pv("hosts_birth", &pv(host_index, "Reservation", simulation.pref().p(), chance_reservation_p), true);
                    chance_reservation_p
                } else {
                    simulation.pv("hosts_birth", &pv(host_index, "Reservation", 0., chance_reservation_p_0), true);
                    chance_reservation_p_0
                };
            }
            HostTypes::Wild => {
                simulation.pv("hosts_birth", &pv(host_index, "Wild", 0., chance_wild), true);
                chance_wild
            }
        }
    }).collect();
    (WeightedIndex::new(chances).unwrap(), choices)
}

pub fn get_chances_v1(v: f32, u: f32, w: f32, y: f32, t: f32, o: f32, p: f32) -> (f32, f32) {
    (
        (1. + o * y * w / u + (1. - o) * y * w / (t + u) + (1. - o) * y * v / (t + u) - p),  // for reservation host individuals
        (1. + o * y * v / t + (1. - o) * y * v / (t + u) + (1. - o) * y * w / (t + u)), // for wild host individuals
    )
}

pub fn birth_generation_version_2(simulation: &mut Simulation) -> (WeightedIndex<f32>, Vec<usize>) {
    let (_, no_of_dead_reservation_host, no_of_dead_wild_host) = simulation.count_dead_hosts();
    let hh = simulation.pref().hh();
    let qr = simulation.ss().qr();
    let qw = simulation.ss().qw();
    let mut qi_values = vec![];
    // select index of the hosts that are alive
    let choices: Vec<usize> = simulation.hosts_alive();
    // chances of each hosts that are alive for birth
    simulation.pv("hosts_birth", "Birth V2\n", true);
    simulation.pv("hosts_birth",
                  &format!("host_index,      host_type,    p,    qi,     chance\n"), true);
    let pv = |hi, ht, p, c, qi| -> String {
        format!("{: >10}, {: >14}, {: >4}, {: >5}, {: >10}\n", hi, ht, p, qi, c)
    };
    let chances: Vec<f32> = choices.iter().map(|host_index| {
        let score_bellow_j = *simulation.ss().host_match_scores_bellow_j().get(&host_index).unwrap();
        let qi = (hh - score_bellow_j) as f32;
        qi_values.push(qi);
        let p = if simulation.ss().hosts_tried().contains(&host_index) {
            simulation.pref().p()
        } else {
            0.
        };
        let vc = get_chances_v2(
            simulation.host_type(*host_index),
            qi, qr, qw,
            simulation.pref().r() as f32 * no_of_dead_reservation_host as f32,
            simulation.pref().s() as f32 * no_of_dead_wild_host as f32,
            simulation.pref().o(),
            simulation.pref().y(),
            p,
        );
        let vc = if vc < 0. || vc.is_nan() || vc.is_infinite() { 0. } else { vc };
        simulation.pv("hosts_birth", &pv(host_index,
                                         simulation.host_type(*host_index).to_string(), p, vc, qi), true);
        vc
    }).collect();
    (WeightedIndex::new(chances).unwrap(), choices)
}

fn calculate_qi(simulation: &mut Simulation) {
    let mut qr = 0.;
    let mut qw = 0.;
    let hosts = simulation.hosts();

    simulation.ss().host_match_scores_bellow_j().iter()
        .filter(|(index, _)| hosts[**index].alive())
        .for_each(|(index, score)| {
            match hosts[*index].host_type() {
                HostTypes::Reservation => qr += simulation.pref().hh() as f32 - *score as f32,
                HostTypes::Wild => qw += simulation.pref().hh() as f32 - *score as f32
            }
        });
    simulation.ss_mut().set_qr(qr);
    simulation.ss_mut().set_qw(qw);
}

/**
- Qi=Q_subscript_i=Qi for each surviving host individual is the variable HH minus the number
      of match scores it has that is below a match threshold,
      J (J is input at the beginning of the program)
- Qr=Q_subscript_r=The sum of Qi for every reservation individual remaining alive in the
      population.
- Qw=Q_subscript_w=The sum of Qi for every wild individual remaining alive in the population.
- Lr=L_subscript_r=R\*Kr
- Lw=L_subscript_w=S\*Kw
- R=reserve_constant
- S=wild_constant
- Kr=K_subscript_r=number of reservation individuals killed this generation
- Kw=K_subscript_w=number of wild individuals killed this generation
 */
pub fn get_chances_v2(host_type: HostTypes, qi: f32, qr: f32, qw: f32, lr: f32, lw: f32, o: f32, y: f32, p: f32) -> f32 {
    match host_type {
        HostTypes::Reservation => qi + (qi / qr) * o * y * lr + qi / (qr + qw) * (1. - o) * y * lr + qi / (qr + qw) * (1. - o) * y * lw - p,
        HostTypes::Wild => qi + qi / (qr + qw) * (1. - o) * y * lr + (qi / qw) * o * y * lw + qi / (qr + qw) * (1. - o) * y * lw
    }
}

fn parasite_truncation_and_birth(simulation: &mut Simulation) {
    let mut _s = String::new();

    let match_scores = simulation.ss().match_scores().clone();
    // match score -> total found
    let mut frequency = HashMap::<usize, usize>::new();
    // addition of all the frequency of scores lower than current
    let mut cumulative_frequency = Vec::with_capacity(simulation.pref().g() + 1);
    // match_score -> Vec<(species_index, parasite_index)>
    let mut individuals_with_score = HashMap::<usize, Vec<(usize, usize)>>::new();
    // prepare data store for frequency, cumulative freq.
    for k in 0..simulation.pref().g() + 1 {
        cumulative_frequency.push(0);
        frequency.insert(k, 0);
        individuals_with_score.insert(k, Vec::new());
    }
    // calculate frequencies and cumulative frequency of each score
    match_scores.iter().for_each(|((species, parasites), match_score)| {
        *frequency.entry(*match_score).or_insert(0) += 1;
        individuals_with_score.get_mut(&match_score).unwrap().push((species.clone(), parasites.clone()));
        cumulative_frequency[*match_score] += 1;
    });

    //
    let mut rng = thread_rng();
    //
    _s.push_str(&format!("{:10} {:10}\n", "killed_individual", "parent_individual"));
    // now calculate the percentile of each match scores
    for i in 1..cumulative_frequency.len() {
        cumulative_frequency[i] += cumulative_frequency[i - 1];
        // percentile rank parasite based on highest match scores
        let percentile = calculate_percentile(
            cumulative_frequency[i],
            *frequency.get(&i).unwrap(),
            match_scores.len(),
        );
        //
        if percentile >= simulation.pref().bb() {
            // get all the parasites that were used
            let parasites = individuals_with_score.get(&i).unwrap();
            let mut already_tried = vec![];
            for (s, p) in parasites {
                // get random existing parasite from the same species (s), excluding this parasite (i)
                let parent_parasite_index = loop {
                    let i1 = rng.gen_range(0..simulation.pref().e());
                    if i1 != *p && !already_tried.contains(&(s, i1)) {
                        already_tried.push((s, i1));
                        break i1;
                    }
                };
                _s.push_str(&format!("({:3}, {:3})         ({:3}, {:3})\n", s, p, s, parent_parasite_index));
                simulation.update_parasites(*s, *p, parent_parasite_index);
            }
        }
    }
    simulation.pv(
        "parasite_birth",
        &format!("PARASITE_BIRTH\n{}", _s),
        true,
    );
}

pub fn print_matching_number_sets(n2: Array1<usize>, species: usize) -> String {
    let mut s = String::new();
    for _ in 0..species * n2.len() {
        s.push_str("   ");
    }
    s.push_str(&format!("{}", n2));
    s
}

pub fn find_match_score(host: &Host, all_parasites: &Array3<usize>, p_idx: &mut ParasiteSpeciesIndex, simulation: &Simulation) -> usize {
    let mut match_count = 0;
    let number_set = host.number_set();
    let start_at = p_idx.species() * simulation.pref().g();

    for ii in start_at..start_at + simulation.pref().g() {
        match simulation.program_version() {
            ProgramVersions::One | ProgramVersions::Two => {
                if number_set[ii] == all_parasites[[p_idx.species(), p_idx.parasite_index, ii - start_at]] {
                    match_count += 1
                }
            }
            ProgramVersions::Three | ProgramVersions::Four => {
                if number_set[ii] != all_parasites[[p_idx.species(), p_idx.parasite_index, ii - start_at]] {
                    match_count += 1
                }
            }
        }
    }
    p_idx.match_count = match_count;
    match_count
}

fn parasite_row(all_parasites: &Array3<usize>, species: usize, parasite: usize) -> Array1<usize> {
    //
    let d = all_parasites.index_axis(Axis(0), species).to_owned();
    d.index_axis(Axis(0), parasite).to_owned()
}


fn mutation(simulation: &mut Simulation) {
    //
    mutate_hosts(simulation);
    mutate_parasites(simulation);
}


/**
Replacement of parasite species ->

For each parasite species during a generation, a total match score is calculated which is the sum
of all match scores for individuals in that parasite species during that generation, including those
match scores resulting from the extra exposure above. The Q species with the highest total match
scores are eliminated (that is, all the E individuals from each of the Q species with the highest
total match scores are eliminated).

In their place come Q new species. For each new species, a set of G numbers is randomly generated,
and these numbers can be either of F possible values. For each individual in the new species,
that individual’s set of G numbers is the same as this set that was just generated for that species.
(This sameness of the set of G numbers for each individual in a parasite species is only
temporary. Differences in the values between each individual in a parasite species will probably
be caused in the Mutation step of each generation.)
 */
fn parasite_replacement(simulation: &mut Simulation) {
    let mut replaced = 0;
    let to_be_replaced = simulation.pref().q();
    // find the max value
    let mut s = String::new();
    let species_match_score = simulation.species_match_score().clone();
    let mut max_keys: Vec<(&usize, &usize)> = species_match_score.iter().map(|(k, v)| (v, k)).collect();
    max_keys.sort();
    while replaced < to_be_replaced {
        let ky = max_keys.pop();
        if ky.is_none() { break }
        let (_, species_index) = ky.unwrap();
        let f = simulation.pref().f();
        let g = simulation.pref().g();
        let species = simulation.parasites_mut().index_axis_mut(Axis(0), *species_index);
        let i = loop {
            let i = generate_individual(f, g);
            if species.index_axis(Axis(0), 1) != i {
                break i;
            }
        };
        let mut species = simulation.parasites_mut()
            .index_axis_mut(Axis(0), *species_index);
        s.push_str(&format!("REPLACE:\n{}\n{:#?}\n", *species_index, species));
        for mut row in species.rows_mut() {
            row.assign(&i);
        }
        replaced += 1;
        s.push_str(&format!("REPLACED_WITH:\n{}\n{:#?}\n", i, species));
    }
    simulation.pv("replaced_parasites", &s, true)
}


fn should_continue(simulation: &mut Simulation) -> bool {
    let (_, r, w) = simulation.count_alive_hosts();
    !(r == 0 || w == 0)
}

#[inline]
/**
Calculates Percentile rank using the following formula

PR = (CF - (0.5 * F)) * 100 / N

where
**CF**—the cumulative frequency—is the count of all scores less than or equal to the score of interest,
**F** is the frequency for the score of interest, and
**N** is the number of scores in the distribution

https://en.wikipedia.org/wiki/Percentile_rank
 */
fn calculate_percentile(cf: usize, f: usize, total: usize) -> f32 {
    // cf, f, total
    ((cf as f32 - f as f32 / 2.) / total as f32) * 100.
}

#[test]
fn test_calculate_percentile() {
    assert_eq!(
        calculate_percentile(10, 2, 10), 90.
    );
    assert_eq!(
        calculate_percentile(8, 3, 10), 65.
    );
    assert_eq!(
        calculate_percentile(5, 5, 10), 25.
    );
}

fn config_logger(parallel_run: i32) -> Result<Handle, SetLoggerError> {
    let logfile = FileAppender::builder()
        .encoder(Box::new(PatternEncoder::new("{l} - {m}\n")))
        .build("report/output.log").unwrap();
    let console = ConsoleAppender::builder()
        .encoder(Box::new(PatternEncoder::new("{l} - {m}\n")))
        .build();
    let mut root_builder = Root::builder();
    let mut config = Config::builder()
        .appender(Appender::builder().build("logfile", Box::new(logfile)));
    if parallel_run == 0 {
        config = config.appender(Appender::builder().build("console", Box::new(console)));
        root_builder = root_builder.appender("console");
    }
    let c = config.build(root_builder.appender("logfile").build(LevelFilter::Info));
    log4rs::init_config(c.unwrap())
}

