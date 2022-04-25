extern crate core;
extern crate ndarray_parallel;
extern crate rand;
extern crate serde;
extern crate serde_derive;
extern crate serde_ini;

use std::fs;
use std::collections::HashMap;
use std::fmt::{Display, Formatter};
use std::fs::{create_dir_all, OpenOptions};
use std::io::Write;
use std::ops::Div;
use std::path::Path;

use futures::executor::block_on;
use ndarray::{Array, Array1, Array3, ArrayBase, Axis, Ix, Ix1, Ix3, OwnedRepr};
use ndarray_rand::RandomExt;
use rand::{Rng, thread_rng};
use rand::distributions::{Uniform, WeightedIndex};
use rand::prelude::*;
use rayon::prelude::*;

use parasites::random_parasite;

use crate::hosts::{Host, HostTypes};
use crate::simulation::{HostsCount, new_simulation, ProgramVersions, Simulation};
use crate::simulation_pref::SimulationPref;

// It's important to shuffle lists before the random processes. In the past when I've hired for
// these sorts of programs, the biggest problem has been that individuals of one type were more
// likely to get chosen for reproduction or death because lists weren't shuffled, and that causes
// the random processes to be biased.

pub mod simulation_pref;
pub mod parasites;
pub mod hosts;
pub mod simulation;


async fn main_async() {
    let pref: SimulationPref = serde_ini::from_str(&fs::read_to_string("params.conf").unwrap()).unwrap();
    let mut result: Vec<HostsCount> = vec![];
    (0..pref.gg()).into_par_iter().map(|gg| {
        println!("simulation running: {}", gg);
        let mut simulation = new_simulation(pref.clone());
        let mut lines = vec![];
        for ff in 0..pref.ff() {
            expose_all_hosts_to_parasites(&mut simulation);
            additional_exposure(&mut simulation);
            if !should_continue(&mut simulation) {
                simulation.next_generation();
                break;
            }
            birth_hosts(&mut simulation);
            parasite_truncation_and_birth(&mut simulation);
            mutation(&mut simulation);
            parasite_replacement(&mut simulation);
            simulation.next_generation();
            lines.push(simulation.generate_report());
            if lines.len() == 3 {
                let s = calculate_result(lines.clone(), pref.clone());
                write_to_file(gg, ff, &s);
            }
        }
        simulation.generate_report()
    }).collect_into_vec(&mut result);
    println!("{:#?}", result);
    let s = calculate_result(result, pref.clone());
    write_to_file(10000, 1000, &s);
}

fn write_to_file(gg: usize, ff: usize, content: &str) {
    let k = gg / 100;
    let p = format!("reports/{}", k);
    let path = Path::new(&p);
    if !std::path::Path::exists(path) {
        std::fs::create_dir_all(path).unwrap();
    }
    let mut f = OpenOptions::new()
        .append(true)
        .create(true) // Optionally create the file if it doesn't already exist
        .open(format!("{}/{}_{}.txt", path.to_str().unwrap(), gg, ff))
        .expect("Unable to open file");
    f.write_all(content.as_bytes()).expect("Unable to write data");
}

fn calculate_result(result: Vec<HostsCount>, pref: SimulationPref) -> String {
    let mut wild_loner = 0;
    let mut reservation_loner = 0;
    let mut high_wild = 0;
    let mut high_reservation = 0;
    let mut tied = 0;
    let mut total_reservation_host = 0;
    let mut total_wild_host = 0;
    for generation in result.iter() {
        let r = generation.reservation_host;
        let w = generation.wild_host;
        if r > w { high_reservation += 1 }
        if r < w { high_wild += 1 }
        if r == w { tied += 1 }
        if r == 0 { wild_loner += 1 }
        if w == 0 { reservation_loner += 1 }
        total_reservation_host += r;
        total_wild_host += w
    }

    let mut _s = String::new();
    _s.push_str(&format!("- {} runs ended with wild individuals the lone type remaining\n", wild_loner));
    _s.push_str(&format!("- {} runs ended with reservation individuals the lone type remaining\n", reservation_loner));
    _s.push_str(&format!("- {} runs ended with wild individuals a higher quantity than reservation individuals\n", high_wild));
    _s.push_str(&format!("- {} runs ended with reservation individuals a higher quantity than wild individuals\n", high_reservation));
    _s.push_str(&format!("- {} runs ended with the quantities of the two types tied\n", tied));
    let mean_wild = total_wild_host as f32 / result.len() as f32;
    let mean_reservation = total_reservation_host as f32 / result.len() as f32;
    _s.push_str(&format!("{},{}\n", mean_wild, mean_reservation));
    let standard_deviation_w: f32 = result.par_iter().map(|v| {
        (v.wild_host as f32 - mean_wild).powf(2.)
    }).sum::<f32>()
        .div(pref.gg() as f32 - 1.)
        .sqrt();
    let standard_deviation_r: f32 = result.par_iter().map(|v| {
        (v.reservation_host as f32 - mean_reservation).powf(2.)
    }).sum::<f32>()
        .div(pref.gg() as f32 - 1.)
        .sqrt();
    _s.push_str(&format!("standard deviation, {}, {}\n", standard_deviation_w, standard_deviation_r));

    let standard_error_w = standard_deviation_w / (pref.gg() as f32).sqrt();
    let standard_error_r = standard_deviation_r / (pref.gg() as f32).sqrt();

    let confidence_interval_w_high_point = mean_wild + pref.z() * standard_error_w;
    let confidence_interval_w_low_point = mean_wild - pref.z() * standard_error_w;

    let confidence_interval_r_high_point = mean_wild + pref.z() * standard_error_r;
    let confidence_interval_r_low_point = mean_wild - pref.z() * standard_error_r;

    _s.push_str(&format!("{},{},{},{}\n", confidence_interval_r_high_point, confidence_interval_r_low_point, confidence_interval_w_high_point, confidence_interval_w_low_point));
    _s
}


fn should_continue(simulation: &mut Simulation) -> bool {
    let (t, r, w) = simulation.count_alive_hosts();
    !(r == 0 || w == 0)
}

fn parasite_replacement(simulation: &mut Simulation) {
    let mut replaced = 0;
    let to_be_replaced = simulation.pref().q();
    // find the max value
    let species_match_score = *simulation.species_match_score().iter().max_by(|a, b| a.1.cmp(&b.1)).map(|(k, v)| v).unwrap();
    let mut max_keys: Vec<usize> = simulation.species_match_score().iter().filter(|v| *v.1 == species_match_score).map(|v| *v.0).collect();
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

fn mutation(simulation: &mut Simulation) {
    let c = simulation.pref().c();
    let f = simulation.pref().f();
    let g = simulation.pref().g();
    let weights: Vec<f32> = (0..simulation.pref().f()).map(|_| simulation.pref().ee()).collect();
    let dist = WeightedIndex::new(&weights).unwrap();
    let choices = [0, 1];
    let mut rng = rand::thread_rng();
    for host in simulation.hosts_mut() {
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
    }
    let weights: Vec<f32> = (0..simulation.pref().f()).map(|_| simulation.pref().k()).collect();
    let dist = WeightedIndex::new(&weights).unwrap();
    for mut parasite in simulation.parasites_mut().rows_mut() {
        let old = parasite.to_owned();
        for cc in 0..g {
            let k = choices[dist.sample(&mut rng)];
            if k == 1 {
                parasite[cc] = rng.gen_range(0..f);
            }
        }
    }
}

/**
All parasite individuals have a match score for a particular generation that is the number of
matching digits it has with the host it matched with in the earlier step. For each parasite species,
all parasite individuals with a match score that is higher than *BB%* of all parasite individuals in
the same parasite species become eliminated. To replace these individuals and restore the number of
parasite individuals in the species to the original quantity of *E*, the parasite individuals in that
species that have not been eliminated have an equal likelihood of producing each new individual in
the species, in a random process. A parasite is able to produce more than one birthed individual
if it is chosen more than once during this random process. A parasite individual has the same
*G*-digit series of numbers as its parent.
 */
fn parasite_truncation_and_birth(simulation: &mut Simulation) {
    let match_scores = simulation.simulation_state().match_scores().clone();
    let mut frequency = HashMap::<usize, usize>::new();
    let mut cumulative_frequency = Vec::with_capacity(simulation.pref().g() + 1);
    let mut individuals_with_score = HashMap::<usize, Vec<(usize, usize)>>::new();
    for k in 0..simulation.pref().g() + 1 {
        cumulative_frequency.push(0);
        frequency.insert(k, 0);
        individuals_with_score.insert(k, Vec::new());
    }
    // calculate frequencies of each score
    for ((s, p), v) in match_scores.iter() {
        *frequency.entry(*v).or_insert(0) += 1;
        individuals_with_score.get_mut(&v).unwrap().push((s.clone(), p.clone()));
        cumulative_frequency[*v] += 1;
    }
    let mut percentiles = HashMap::<usize, f32>::new();
    // calculate cumulative frequency of each match
    let mut rng = thread_rng();
    for i in 1..cumulative_frequency.len() {
        cumulative_frequency[i] += cumulative_frequency[i - 1];
        percentiles.insert(i, calculate_percentile(cumulative_frequency[i], *frequency.get(&i).unwrap(), match_scores.len()));
        if percentiles.get(&i).unwrap() >= &simulation.pref().bb() {
            let parasites = individuals_with_score.get(&i).unwrap();
            let mut _s = String::new();
            for (s, i) in parasites {
                // get random existing parasite from the same species (s), excluding this parasite (i)
                let parent_parasite_index = loop {
                    let i1 = rng.gen_range(0..simulation.pref().e());
                    if i1 != *i { break i1; }
                };
                let b = simulation.parasites().index_axis(Axis(0), *s);
                let v = b.index_axis(Axis(0), *i);
                _s.push_str(&format!("{: >2} {: >3}: {: >2} : ", s, i, v));
                let v = b.index_axis(Axis(0), parent_parasite_index);
                _s.push_str(&format!("{: >2} -> ", v));
                simulation.update_parasites(*s, *i, parent_parasite_index);
                let b = simulation.parasites().index_axis(Axis(0), *s);
                let v = b.index_axis(Axis(0), *i);
                _s.push_str(&format!("{: >2}\n", v));
            }
        }
    }
}

#[inline]
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

fn main() {
    block_on(main_async());
}


pub fn generate_individual(f: usize, len: usize) -> Array1<usize> {
    Array::random(len, Uniform::new(0, f))
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

pub fn expose_all_hosts_to_parasites(simulation: &mut Simulation) {
    let mut parasites_exposed_to = HashMap::new();
    let mut species_match_score = HashMap::<usize, usize>::new();
    let mut host_match_score = HashMap::<usize, usize>::new();
    let all_parasites = simulation.parasites().clone();
    let hosts = simulation.hosts().clone();
    for i in 0..simulation.pref().a() + simulation.pref().b() {
        let host = &hosts[i];
        let mut match_score_bellow_threshold = 0;
        for _ in 0..simulation.pref().h() {
            let mut p_idx = get_random_parasite(simulation, &parasites_exposed_to);
            let match_score = find_match_score(&host, &all_parasites, &mut p_idx, simulation.pref(), simulation.program_version());
            parasites_exposed_to.insert((p_idx.species(), p_idx.parasite()), match_score);
            *species_match_score.entry(p_idx.species()).or_insert(0) += match_score;
            if match_score < simulation.pref().n() {
                match_score_bellow_threshold += 1;
            }
            if match_score < simulation.pref().j() {
                *host_match_score.entry(i).or_insert(0) += 1;
            }
        }
        if match_score_bellow_threshold >= simulation.pref().x() {
            simulation.kill_host(i);
        }
    }
    simulation.update_parasites_exposed_to(parasites_exposed_to);
    simulation.update_species_match_score(species_match_score);
    simulation.update_host_match_score(host_match_score);
}

fn get_random_parasite(simulation: &mut Simulation, parasites_exposed_to: &HashMap<(usize, usize), usize>) -> ParasiteSpeciesIndex {
    loop {
        let p_idx = random_parasite(simulation.pref());
        let ky = (p_idx.species(), p_idx.parasite());
        if parasites_exposed_to.contains_key(&ky) { continue; }
        break p_idx;
    }
}

pub fn additional_exposure(simulation: &mut Simulation) {
    let mut species_match_score = simulation.species_match_score().clone();
    let mut host_match_score = simulation.simulation_state().host_match_score().clone();
    let (total_dead_hosts, _, _) = simulation.count_dead_hosts();
    let (_, alive_reservation_hosts, _) = simulation.count_alive_hosts();
    // secondary exposure
    let secondary_allowed = (simulation.pref().a() + simulation.pref().b()) as f32 * simulation.pref().m();
    if simulation.current_generation() >= simulation.pref().l() && total_dead_hosts < secondary_allowed as usize {
        return;
    }
    let all_parasites = simulation.parasites().clone();
    let no_of_additional_host = (alive_reservation_hosts as f32 * simulation.pref().aa()).ceil() as usize;
    let mut rng = thread_rng();
    let mut hosts_to_try = 0;
    let mut parasites_exposed_to = simulation.parasites_exposed_to();
    while hosts_to_try < no_of_additional_host {
        let index: usize = rng.gen_range(0..simulation.pref().a() + simulation.pref().b());
        // let host = &simulation.hosts()[index];
        // only reservation hosts
        if simulation.hosts()[index].host_type() == HostTypes::Wild { continue; }
        // skip dead hosts
        if !simulation.hosts()[index].alive() { continue; }
        hosts_to_try += 1;
        // expose to parasite
        let mut match_score_bellow_threshold = 0;
        for _ in 0..simulation.pref().i() {
            let mut p_idx = get_random_parasite(simulation, &parasites_exposed_to);
            let match_score = find_match_score(&simulation.hosts()[index], &all_parasites, &mut p_idx, simulation.pref(), simulation.program_version());
            parasites_exposed_to.insert((p_idx.species(), p_idx.parasite()), match_score);
            *species_match_score.entry(p_idx.species()).or_insert(0) += match_score;
            if match_score < simulation.pref().dd() {
                match_score_bellow_threshold += 1;
            }
            if match_score < simulation.pref().j() {
                *host_match_score.entry(index).or_insert(0) += 1;
            }
        }
        if match_score_bellow_threshold >= simulation.pref().cc() {
            simulation.kill_host(index);
        }
    }
    simulation.update_parasites_exposed_to(parasites_exposed_to);
    simulation.update_species_match_score(species_match_score);
    simulation.update_host_match_score(host_match_score);
}

pub fn birth_hosts(simulation: &mut Simulation) {
    let (dist, choices) = match simulation.program_version() {
        ProgramVersions::One => birth_generation_version_1(simulation),
        ProgramVersions::Two => birth_generation_version_2(simulation),
        ProgramVersions::Three => birth_generation_version_1(simulation),
        ProgramVersions::Four => birth_generation_version_1(simulation),
    };
    let mut rng = rand::thread_rng();
    loop {
        // pick up parent host
        let pick_reservation_host = choices[dist.sample(&mut rng)];
        let parent_index = loop {
            let p = rng.gen_range(0..simulation.hosts().len());
            let p_host = &simulation.hosts()[p];
            if p_host.host_type() == HostTypes::Reservation && pick_reservation_host == 1 && p_host.alive() {
                break p;
            }
            if pick_reservation_host == 0 && p_host.host_type() == HostTypes::Wild && p_host.alive() {
                break p;
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
        // now we get the host
        let (total_dead_hosts, _, _) = simulation.update_dead_host(host_index, parent_index);
        if total_dead_hosts == 0 { break; }
    }
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
    (WeightedIndex::new(vec![chance_wild, chance_reservation]).unwrap(), choices)
}

pub fn birth_generation_version_2(simulation: &mut Simulation) -> (WeightedIndex<f32>, Vec<usize>) {
    let (_, no_of_dead_reservation_host, no_of_dead_wild_host) = simulation.count_dead_hosts();
    let (_, no_of_reservation_host_alive, no_of_wild_host_alive) = simulation.count_alive_hosts();
    let host_match_score = simulation.simulation_state().host_match_score();
    let qr: f32 = find_sum_of_match_score(host_match_score, simulation.hosts(), HostTypes::Reservation);
    let qw: f32 = find_sum_of_match_score(host_match_score, simulation.hosts(), HostTypes::Wild);
    let choices: Vec<usize> = simulation.hosts().iter().enumerate()
        .filter(|v| v.1.alive())
        .map(|v| v.0)
        .collect();
    let chances: Vec<f32> = choices.iter().map(|v| {
        get_chances_v2(
            simulation.hosts()[*v].host_type(),
            *host_match_score.get(&v).unwrap() as f32,
            qr,
            qw,
            simulation.pref().r() as f32,
            simulation.pref().s() as f32,
            no_of_dead_wild_host as f32,
            no_of_dead_reservation_host as f32,
            simulation.pref().o(),
            simulation.pref().y(),
            simulation.pref().p(),
        )
    }).collect();
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

pub fn find_match_score(host: &Host, all_parasites: &Array3<usize>, p_idx: &mut ParasiteSpeciesIndex, pref: SimulationPref, version: ProgramVersions) -> usize {
    let mut match_count = 0;
    let number_set = host.number_set();
    let index = p_idx.species();
    for ii in index..index + pref.g() {
        match version {
            ProgramVersions::One | ProgramVersions::Two => {
                if number_set[ii] == all_parasites[[index, ii - index, 0]] {
                    match_count += 1
                }
            }
            ProgramVersions::Three | ProgramVersions::Four => {
                if number_set[ii] != all_parasites[[index, ii - index, 0]] {
                    match_count += 1
                }
            }
        }
    }
    p_idx.match_count = match_count;
    match_count
}