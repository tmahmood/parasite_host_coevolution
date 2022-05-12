#![feature(iter_collect_into)]
extern crate core;
extern crate rand;
extern crate serde;
extern crate serde_derive;
extern crate serde_ini;

use std::{fs, time};
use std::collections::HashMap;
use std::env::args;
use std::fmt::{Display, Formatter};
use std::fs::{File, remove_dir_all};
use std::io::Write;
use std::ops::Div;

use indicatif::{MultiProgress, ParallelProgressIterator, ProgressBar, ProgressIterator, ProgressStyle};
use log;
use log::{LevelFilter, SetLoggerError};
use log4rs::{Config, Handle};
use log4rs::append::console::ConsoleAppender;
use log4rs::append::file::FileAppender;
use log4rs::append::rolling_file::RollingFileAppender;
use log4rs::config::{Appender, Root};
use log4rs::encode::pattern::PatternEncoder;
use ndarray::{arr2, Array, Array1, Array2, Array3, ArrayBase, ArrayView1, Axis, Ix, Ix1, Ix2, Ix3, OwnedRepr};
use ndarray_rand::RandomExt;
use rand::{Rng, thread_rng};
use rand::distributions::{Uniform, WeightedIndex};
use rand::prelude::*;
use rayon::prelude::*;

use crate::hosts::{create_random_hosts, Host, HostTypes};
use crate::simulation::{create_random_parasites, GGRunReport, HostsCount, new_simulation, print_parasites, ProgramVersions, ReportHostType, Simulation};
use crate::simulation_pref::SimulationPref;

pub mod simulation_pref;
pub mod parasites;
pub mod hosts;
pub mod simulation;

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

fn main() {
    let multi_progress_bar = MultiProgress::new();
    let sty = ProgressStyle::with_template(
        "[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7} {msg}",
    )
        .unwrap()
        .progress_chars("##-");


    remove_dir_all("report").expect("Failed to remove report directory");
    let param_file = args().nth(1).unwrap_or(format!("params.conf"));
    let program = ProgramVersions::from(args().nth(2).unwrap_or(format!("1")));
    config_logger(1).unwrap();
    let pref: SimulationPref = serde_ini::from_str(&fs::read_to_string(param_file).unwrap()).unwrap();
    //
    let pb_all_simulation = multi_progress_bar.add(ProgressBar::new(pref.gg() as u64));
    pb_all_simulation.set_style(sty.clone());
    //
    let pb_generations = multi_progress_bar.insert_after(&pb_all_simulation, ProgressBar::new(pref.ff() as u64));
    pb_generations.set_style(sty.clone());
    //

    type UsizeVec = Vec<(Vec<usize>, Vec<usize>)>;
    let now = time::Instant::now();

    println!("Running parallel version, will log to file. {}", program);
    let program_clone = program.clone();
    let pref_clone = pref.clone();
    let mut wild_hosts = vec![];
    let mut reservation_hosts = vec![];
    (0..pref.gg()).into_par_iter().progress_with(pb_all_simulation).map(move |gg| {
        run_generation_step(program_clone, gg.clone(), pref_clone)
    }).unzip_into_vecs(
        &mut wild_hosts,
        &mut reservation_hosts,
    );

    let flat_wild_hosts: Vec<usize> = wild_hosts.iter().flatten().cloned().collect();
    let flat_reservation_hosts: Vec<usize> = reservation_hosts.iter().flatten().cloned().collect();
    let reservation_hosts_ar = Array2::from_shape_vec((pref.gg(), pref.ff()), flat_reservation_hosts).unwrap();
    let wild_hosts_ar = Array2::from_shape_vec((pref.gg(), pref.ff()), flat_wild_hosts).unwrap();

    generate_excel(&reservation_hosts_ar, "reservation_host_confidence", &pref);
    generate_excel(&wild_hosts_ar, "wild_host_confidence", &pref);
    let last_row_r = reservation_hosts_ar.index_axis(Axis(1), pref.ff() - 1).to_owned();
    let last_row_w = wild_hosts_ar.index_axis(Axis(1), pref.ff() - 1).to_owned();
    let r = calculate_result((last_row_r, last_row_w), pref.gg());
    println!("{}", r);
    println!("took {} secs", now.elapsed().as_secs())
}

fn generate_excel(hosts: &ArrayBase<OwnedRepr<usize>, Ix2>, x: &str, pref: &SimulationPref) {
    let mut f = File::create(format!("report/{}.csv", x)).expect("Unable to create file");
    f.write_all(format!("Standard Deviation, Means, High Point, Low Point\n").as_bytes()).expect("Failed to write to file");
    hosts.columns().into_iter().for_each(|v| {
        let mut r = ReportHostType::new(v.to_owned());
        r.calculate(pref.clone());
        f.write_all(
            format!("{},{},{},{}\n",
                    r.standard_deviation(),
                    r.mean(),
                    r.confidence_interval_high_point(),
                    r.confidence_interval_low_point()
            ).as_bytes()).expect("Failed to write to file");
    });
}

fn calculate_result(result: (Array1<usize>, Array1<usize>), ff: usize) -> GGRunReport {
    let mut report_gen = GGRunReport::new(result, ff);
    report_gen.calculations();
    report_gen
}

fn run_generation_step(program: ProgramVersions, gg: usize, pref: SimulationPref) -> (Vec<usize>, Vec<usize>) {
    let mut simulation = new_simulation(pref.clone(), program, gg);
    let mut lines = vec![];
    (0..pref.ff()).into_iter().map(|_| {
        expose_all_hosts_to_parasites(&mut simulation);
        additional_exposure(&mut simulation);
        if !should_continue(&mut simulation) {
            return simulation.next_generation();
        }
        birth_hosts(&mut simulation);
        parasite_truncation_and_birth(&mut simulation);
        mutation(&mut simulation);
        parasite_replacement(&mut simulation);
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
    // we are calculating all the random species for each hosts before hand, that way we don't
    // have to calculate it in the loop and most likely it will be optimized, we remove the species
    // that is used, that way we don't have to check for used parasite individual
    let _species_possible = (0..simulation.pref().d()).collect::<Vec<usize>>();
    let mut species_possible = HashMap::new();
    for i in 0..hosts.len() {
        let mut k = _species_possible.clone();
        k.shuffle(&mut rng);
        species_possible.insert(i, k);
    }
    // same with parasites
    let _parasites_possible = (0..simulation.pref().e()).collect::<Vec<usize>>();
    let mut parasites_possible = vec![];
    for _ in 0..simulation.pref().d() {
        let mut k = _parasites_possible.clone();
        k.shuffle(&mut rng);
        parasites_possible.push(k);
    }
    for host_index in 0..simulation.pref().a() + simulation.pref().b() {
        let host = &hosts[host_index];
        simulation.pv(file_name, &format!("{: <3} {:12}{}\n(species, parasite): match count -> code, \n", host_index, &host.host_type().to_string(), &host.number_set()), true);
        let mut match_score_bellow_threshold = 0;
        let mut species = species_possible.get_mut(&host_index).unwrap();
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
            //
            let d = parasite_row(&all_parasites, p_idx.species(), p_idx.parasite_index);
            let p_grid = print_matching_number_sets(
                host.number_set().clone(),
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
    let (_, alive_reservation_hosts, _) = simulation.count_alive_hosts();
    // secondary exposure
    let secondary_allowed = (simulation.pref().a() + simulation.pref().b()) as f32 * simulation.pref().m();
    if simulation.current_generation() as i32 <= simulation.pref().l() {
        return;
    }
    if total_dead_hosts > secondary_allowed as usize {
        return;
    }
    simulation.pv(file_name, &format!("{}\n", simulation.current_generation()), true);
    let all_parasites = simulation.parasites().clone();
    let no_of_additional_host = (alive_reservation_hosts as f32 * simulation.pref().aa()).ceil() as usize;
    simulation.pv(file_name, &format!("Additional Exposure candidate {}\n", no_of_additional_host), true);
    let mut rng = thread_rng();
    let mut hosts_to_try = 0;
    let mut species_par_host = simulation.species_left();
    let mut parasites_possible = simulation.parasites_possible();
    let mut tried = HashMap::new();
    while hosts_to_try < no_of_additional_host {
        let host_index: usize = rng.gen_range(0..simulation.pref().a() + simulation.pref().b());
        if tried.contains_key(&host_index) { continue; }
        tried.insert(host_index, 0);
        let host = simulation.hosts()[host_index].clone();
        // only reservation hosts
        if host.host_type() == HostTypes::Wild { continue; }
        // skip dead hosts
        if !host.alive() { continue; }
        simulation.pv(file_name, &format!("{: <3} {:12}{}\n(species, parasite): match count -> code, \n", host_index, &host.host_type().to_string(), &host.number_set()), true);
        hosts_to_try += 1;
        // expose to parasite
        let mut match_score_bellow_threshold = 0;
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
            if match_score < simulation.pref().dd() {
                match_score_bellow_threshold += 1;
            }
            //
            let d = parasite_row(&all_parasites, p_idx.species(), p_idx.parasite_index);
            let p_grid = print_matching_number_sets(
                host.number_set().clone(),
                d,
                p_idx.species(),
            );
            simulation.pv(file_name, &format!("{: >3},{: >3} : {: >2} -> {: >3}\n", p_idx.species(), p_idx.parasite(), match_score, p_grid), true);
            //
        }
        simulation.pv(file_name, &format!("-------------------------\n"), true);
        if match_score_bellow_threshold >= simulation.pref().cc() {
            simulation.kill_host(host_index);
            simulation.pv("host_dying_additional_exposure", &format!("{:3} {}\n", host_index, &host.to_string()), true);
        }
    }
    let (t, r, w) = simulation.count_dead_hosts();
    simulation.pv("host_dying_additional_exposure", &format!("{} R({}) W({})", t, r, w), true);
}

pub fn birth_hosts(simulation: &mut Simulation) {
    let file_name = "hosts_birth";
    let (dist, choices) = match simulation.program_version() {
        ProgramVersions::One | ProgramVersions::Three => birth_generation_version_1(simulation),
        ProgramVersions::Two | ProgramVersions::Four => birth_generation_version_2(simulation),
    };
    let mut rng = thread_rng();
    simulation.pv(file_name,
                  &format!("{}\n{: >4} {: <width$} {}\n",
                           simulation.program_version().to_string(),
                           "index", "host", "target",
                           width = simulation.pref().c() + 57),
                  true);
    //simulation.pv(file_name, , true);
    loop {
        // pick up parent host
        let random_host_index = choices[dist.sample(&mut rng)];
        let parent_index = match simulation.program_version() {
            ProgramVersions::One | ProgramVersions::Three => {
                random_host_selection_v1(&simulation, random_host_index)
            }
            ProgramVersions::Two | ProgramVersions::Four => {
                choices[dist.sample(&mut rng)]
            }
        };

        let parent_host = simulation.hosts()[parent_index].clone();
        let mut index = None;
        for (ii, host) in simulation.hosts().iter().enumerate() {
            if !host.alive() {
                index = Some(ii);
                break;
            }
        }
        if index == None {
            panic!("I couldn't find any dead [{}] hosts!", parent_host.host_type());
        }
        let host_index = index.unwrap();
        simulation.pv(file_name, &format!(" {: >4} {: >width$} {:4}\n",
                                          parent_index,
                                          parent_host, host_index,
                                          width = simulation.pref().c() + 13), true);
        // now we get the host
        let (total_dead_hosts, _, _) = simulation.update_dead_host(host_index, parent_index);
        if total_dead_hosts == 0 { break; }
    }
}

fn random_host_selection_v1(simulation: &Simulation, random_host_index: usize) -> usize {
    let mut rng = thread_rng();
    loop {
        let p = rng.gen_range(0..simulation.hosts().len());
        let p_host = &simulation.hosts()[p];
        if p_host.host_type() == HostTypes::Reservation && random_host_index == 1 && p_host.alive() {
            break p;
        }
        if random_host_index == 0 && p_host.host_type() == HostTypes::Wild && p_host.alive() {
            break p;
        }
    }
}

pub fn birth_generation_version_1(simulation: &mut Simulation) -> (WeightedIndex<f32>, Vec<usize>) {
    let (_, no_of_dead_reservation_host, no_of_dead_wild_host) = simulation.count_dead_hosts();
    let (_, no_of_reservation_host_alive, no_of_wild_host_alive) = simulation.count_alive_hosts();
    let (chance_reservation, chance_wild) = get_chances_v1(
        no_of_dead_wild_host as f32,
        no_of_reservation_host_alive as f32,
        no_of_dead_reservation_host as f32,
        simulation.pref().y(),
        no_of_wild_host_alive as f32,
        simulation.pref().o(),
        simulation.pref().p(),
    );
    let choices = vec![0, 1];
    simulation.pv("hosts_birth", &format!("Birth V1: {}, {}, {:?}\n", chance_reservation, chance_wild, choices), true);
    (WeightedIndex::new(vec![chance_wild, chance_reservation]).unwrap(), choices)
}

pub fn birth_generation_version_2(simulation: &mut Simulation) -> (WeightedIndex<f32>, Vec<usize>) {
    let (_, no_of_dead_reservation_host, no_of_dead_wild_host) = simulation.count_dead_hosts();
    let host_match_score = simulation.simulation_state().host_match_score().clone();
    let qr: f32 = find_sum_of_match_score(&host_match_score, simulation.hosts(), HostTypes::Reservation);
    let qw: f32 = find_sum_of_match_score(&host_match_score, simulation.hosts(), HostTypes::Wild);
    let choices: Vec<usize> = simulation.hosts().iter().enumerate()
        .filter(|v| v.1.alive())
        .map(|v| v.0)
        .collect();
    let chances: Vec<f32> = choices.iter().map(|v| {
        let qi = *host_match_score.get(&v).unwrap();
        get_chances_v2(
            simulation.hosts()[*v].host_type(),
            qi as f32, qr, qw,
            simulation.pref().r() as f32,
            simulation.pref().s() as f32,
            no_of_dead_wild_host as f32,
            no_of_dead_reservation_host as f32,
            simulation.pref().o(),
            simulation.pref().y(),
            simulation.pref().p(),
        )
    }).collect();
    simulation.pv("hosts_birth", &format!("Birth V2: {:#?}, {:#?}\n", chances, choices), true);
    (WeightedIndex::new(chances).unwrap(), choices)
}


pub fn get_chances_v1(v: f32, u: f32, w: f32, y: f32, t: f32, o: f32, p: f32) -> (f32, f32) {
    (
        (1. + o * y * w / u + (1. - o) * y * w / (t + u) + (1. - o) * y * v / (t + u) - p),  // for reservation host individuals
        (1. + o * y * v / t + (1. - o) * y * v / (t + u) + (1. - o) * y * w / (t + u)), // for wild host individuals
    )
}

pub fn get_chances_v2(host_type: HostTypes, qi: f32, qr: f32, qw: f32, r: f32, s: f32, kr: f32, kw: f32, o: f32, y: f32, p: f32) -> f32 {
    let lr = r * kr;
    let lw = s * kw;
    match host_type {
        HostTypes::Reservation => qi + (qi / qr) * o * y * lr + qi / (qr + qw) * (1. - o) * y * lr + qi / (qr + qw) * (1. - o) * y * lw - p,
        HostTypes::Wild => qi + qi / (qr + qw) * (1. - o) * y * lr + (qi / qw) * o * y * lw + qi / (qr + qw) * (1. - o) * y * lw
    }
}

fn parasite_truncation_and_birth(simulation: &mut Simulation) {
    let mut _s = String::new();

    let match_scores = simulation.simulation_state().match_scores().clone();
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
    _s.push_str(&format!("{: >6} {: >12} {: >12}\n", "parasite", "killed", "new parent"));
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
            for (s, p) in parasites {
                // get random existing parasite from the same species (s), excluding this parasite (i)
                let parent_parasite_index = loop {
                    let i1 = rng.gen_range(0..simulation.pref().e());
                    if i1 != *p { break i1; }
                };
                // only for print
                let b = simulation.parasites().index_axis(Axis(0), *s);
                let d = b.index_axis(Axis(0), *p);
                //
                let b = simulation.parasites().index_axis(Axis(0), *s);
                let v = b.index_axis(Axis(0), parent_parasite_index);
                _s.push_str(&format!("{:?}   {: >2} {: >2}\n", (s, i), d, v));
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

pub fn print_matching_number_sets(n1: Array1<usize>, n2: Array1<usize>, species: usize) -> String {
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

fn find_sum_of_match_score(host_match_score: &HashMap<usize, usize>, hosts: &Array<Host, Ix1>, host_type: HostTypes) -> f32 {
    host_match_score.into_par_iter().fold(|| 0f32, |mut acc, v| {
        if hosts[*v.0].host_type() == host_type {
            acc += *v.1 as f32;
            acc
        } else {
            acc
        }
    }).sum::<f32>()
}

fn mutation(simulation: &mut Simulation) {
    let c = simulation.pref().c();
    let f = simulation.pref().f();
    let g = simulation.pref().g();
    let weights: Vec<f32> = (0..simulation.pref().f()).map(|_| simulation.pref().ee()).collect();
    let dist = WeightedIndex::new(&weights).unwrap();
    let choices = [0, 1];
    let mut rng = rand::thread_rng();
    let mut mutated_hosts = String::new();
    for (i, host) in simulation.hosts_mut().iter_mut().enumerate() {
        let mut m = host.number_set().clone();
        let mut changes = 0;
        for cc in 0..c {
            let k = choices[dist.sample(&mut rng)];
            if k == 1 {
                changes += 1;
                m[cc] = rng.gen_range(0..f);
            }
        }
        host.set_number_set(m);
        mutated_hosts.push_str(&format!("{:4} {}\n", i, host.number_set().clone()));
    }
    simulation.pv("mutated_hosts", &mutated_hosts, true);
    let weights: Vec<f32> = (0..simulation.pref().f()).map(|_| simulation.pref().k()).collect();
    let dist = WeightedIndex::new(&weights).unwrap();
    for mut parasite in simulation.parasites_mut().rows_mut() {
        for cc in 0..g {
            let k = choices[dist.sample(&mut rng)];
            if k == 1 {
                parasite[cc] = rng.gen_range(0..f);
            }
        }
    }
    simulation.pv("mutated_parasites",
                  &format!("{}",
                           print_parasites(simulation.parasites())),
                  true);
}

fn parasite_replacement(simulation: &mut Simulation) {
    let mut replaced = 0;
    let to_be_replaced = simulation.pref().q();
    // find the max value
    let species_match_score = *simulation.species_match_score().iter()
        .max_by(|a, b| a.1.cmp(&b.1))
        .map(|(k, v)| v)
        .unwrap();
    let mut max_keys: Vec<usize> = simulation.species_match_score().iter()
        .filter(|v| *v.1 == species_match_score)
        .map(|v| *v.0)
        .collect();
    let i = generate_individual(simulation.pref().f(), simulation.pref().g());
    while replaced < to_be_replaced && max_keys.len() > 0 {
        let ky = max_keys.pop().unwrap();
        let mut species = simulation.parasites_mut().index_axis_mut(Axis(0), ky);
        for (ii, iv) in i.iter().enumerate() {
            for mut row in species.rows_mut() {
                row[ii] = *iv;
            }
        }
        replaced += 1;
    }
}


fn should_continue(simulation: &mut Simulation) -> bool {
    let (t, r, w) = simulation.count_alive_hosts();
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

