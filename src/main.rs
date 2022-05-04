extern crate core;
extern crate rand;
extern crate serde;
extern crate serde_derive;
extern crate serde_ini;

use std::{fs, time};
use std::collections::HashMap;
use std::env::args;
use std::fmt::{Display, format, Formatter};
use std::fs::{create_dir_all, File, OpenOptions, remove_dir_all, remove_file};
use std::io::Write;
use std::ops::Div;
use std::path::Path;

use futures::executor::block_on;
use indicatif::{MultiProgress, ParallelProgressIterator, ProgressBar, ProgressIterator, ProgressStyle};
use log;
use log::{info, LevelFilter, SetLoggerError};
use log4rs::{Config, Handle};
use log4rs::append::console::ConsoleAppender;
use log4rs::append::file::FileAppender;
use log4rs::append::rolling_file::RollingFileAppender;
use log4rs::config::{Appender, Root};
use log4rs::encode::pattern::PatternEncoder;
use ndarray::{Array, Array1, Array3, ArrayBase, ArrayView1, Axis, Ix, Ix1, Ix3, OwnedRepr};
use ndarray_rand::RandomExt;
use rand::{Rng, thread_rng};
use rand::distributions::{Uniform, WeightedIndex};
use rand::prelude::*;
use rayon::prelude::*;

use parasites::random_parasite;

use crate::hosts::{Host, HostTypes};
use crate::simulation::{HostsCount, new_simulation, print_parasites, ProgramVersions, Simulation};
use crate::simulation_pref::SimulationPref;
use crate::simulation_runner::SimulationRunner;

// It's important to shuffle lists before the random processes. In the past when I've hired for
// these sorts of programs, the biggest problem has been that individuals of one type were more
// likely to get chosen for reproduction or death because lists weren't shuffled, and that causes
// the random processes to be biased.

pub mod simulation_pref;
pub mod parasites;
pub mod hosts;
pub mod simulation;

pub mod simulation_runner {
    use crate::{Simulation, SimulationPref};

    pub struct SimulationRunner {
        pref: SimulationPref,
    }

    impl SimulationRunner {
        pub fn new(pref: SimulationPref) -> SimulationRunner {
            SimulationRunner {
                pref
            }
        }
    }
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
    let mean_wild = total_wild_host as f32 / result.len() as f32;
    let mean_reservation = total_reservation_host as f32 / result.len() as f32;
    // Sqrt[
    //   (
    //      (100-100.4)^2 +
    //      (80-100.4)^2 +
    //      (92-100.4)^2 +
    //      (110-100.4)^2 +
    //      (120-100.4)^2
    //   ) /4
    //   ] = 15.5177
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

    let standard_error_w = standard_deviation_w / (pref.gg() as f32).sqrt();
    let standard_error_r = standard_deviation_r / (pref.gg() as f32).sqrt();

    let confidence_interval_w_high_point = mean_wild + pref.z() * standard_error_w;
    let confidence_interval_w_low_point = mean_wild - pref.z() * standard_error_w;

    let confidence_interval_r_high_point = mean_reservation + pref.z() * standard_error_r;
    let confidence_interval_r_low_point = mean_reservation - pref.z() * standard_error_r;

    _s.push_str(&format!("- {} runs ended with wild individuals the lone type remaining\n", wild_loner));
    _s.push_str(&format!("- {} runs ended with reservation individuals the lone type remaining\n", reservation_loner));
    _s.push_str(&format!("- {} runs ended with wild individuals a higher quantity than reservation individuals\n", high_wild));
    _s.push_str(&format!("- {} runs ended with reservation individuals a higher quantity than wild individuals\n", high_reservation));
    _s.push_str(&format!("- {} runs ended with the quantities of the two types tied\n", tied));
    // _s.push_str(&format!("means:\n WILD: {: <10}\n RESERVATION: {: <10}\n", mean_wild, mean_reservation));
    _s.push_str(&format!("standard deviation:\n WILD        {:10.3}\n RESERVATION {:10.3}\n", standard_deviation_w, standard_deviation_r));
    _s.push_str(&format!("Confidence Interval:\n {: >24} {: >11}{: >12}\n", "Means", "High Point", "Low Point"));
    _s.push_str(&format!(" WILD        {:12.3}{:12.3}{:12.3}\n RESERVATION {:12.3}{:12.3}{:12.3}\n",
                         mean_wild, confidence_interval_w_high_point, confidence_interval_w_low_point,
                         mean_reservation, confidence_interval_r_high_point, confidence_interval_r_low_point));
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


fn main() {

    let multi_progress_bar = MultiProgress::new();
    let sty = ProgressStyle::with_template(
        "[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7} {msg}",
    )
        .unwrap()
        .progress_chars("##-");



    remove_dir_all("report");
    let param_file = args().nth(1).unwrap_or(format!("params.conf"));
    let program = ProgramVersions::from(args().nth(2).unwrap_or(format!("1")));
    let parallel_run = args().nth(3).unwrap_or(format!("0"));
    config_logger(
        if parallel_run == "0" {
            0
        } else {
            1
        }
    ).unwrap();
    let pref: SimulationPref = serde_ini::from_str(&fs::read_to_string(param_file).unwrap()).unwrap();

    //
    let pb_all_simulation = multi_progress_bar.add(ProgressBar::new(pref.gg() as u64));
    pb_all_simulation.set_style(sty.clone());
    //
    let pb_generations = multi_progress_bar.insert_after(&pb_all_simulation, ProgressBar::new(pref.ff() as u64));
    pb_generations.set_style(sty.clone());
    //

    let now = time::Instant::now();
    let mut result: Vec<HostsCount> = vec![];

    if parallel_run == "1" {
        println!("Running parallel version, will log to file. {}", program);
        let program_clone = program.clone();
        let pref_clone = pref.clone();
        (0..pref.gg()).into_par_iter().progress_with(pb_all_simulation).map(move |gg| {
            run_generation_step(program_clone, gg.clone(), pref_clone)
        }).collect_into_vec(&mut result);
    } else {
        info!("Running normal version. {}", program);
        let program_clone = program.clone();
        let pref_clone = pref.clone();
        result = (0..pref.gg()).into_iter().progress_with(pb_all_simulation).map(|gg| {
            run_generation_step(program_clone, gg.clone(), pref_clone)
        }).collect();
    }
    let r = calculate_result(result.clone(), pref.clone());
    info!("{}", r);
    println!("took {} secs", now.elapsed().as_secs())
}

fn run_generation_step(program: ProgramVersions, gg: usize, pref: SimulationPref) -> HostsCount {
    let bar = ProgressBar::new(pref.ff() as u64);
    let mut simulation = new_simulation(pref.clone(), program, gg);
    let mut lines = vec![];
    (0..pref.ff()).into_iter().progress().for_each(|_| {
        expose_all_hosts_to_parasites(&mut simulation);
        additional_exposure(&mut simulation);
        if !should_continue(&mut simulation) {
            simulation.next_generation();
            return;
        }
        birth_hosts(&mut simulation);
        parasite_truncation_and_birth(&mut simulation);
        mutation(&mut simulation);
        parasite_replacement(&mut simulation);
        simulation.next_generation();
        lines.push(simulation.generate_report());
    });
    simulation.generate_report()
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
    let mut rng = thread_rng();
    let file_name = "host_exposed_to";
    let mut parasites_exposed_to = HashMap::new();
    let mut species_match_score = HashMap::<usize, usize>::new();
    let mut host_match_score = HashMap::<usize, usize>::new();
    // let mut host_parasite_score = HashMap::new();
    simulation.pv(file_name, "Initial Exposure\n", true);
    let all_parasites = simulation.parasites().clone();
    let hosts = simulation.hosts().clone();

    let mut species_possible = (0..simulation.pref().d()).collect::<Vec<usize>>();
    let mut parasites_possible = (0..simulation.pref().e()).collect::<Vec<usize>>();
    for i in 0..simulation.pref().a() + simulation.pref().b() {
        let host = &hosts[i];
        simulation.pv(file_name, &format!("{: <3} {:12}{}\n(species, parasite): match count -> code, \n", i, &host.host_type().to_string(), &host.number_set()), true);
        let mut match_score_bellow_threshold = 0;
        let mut species = species_possible.clone();
        let mut parasites = parasites_possible.clone();
        species.shuffle(&mut rng);
        parasites.shuffle(&mut rng);
        for _ in 0..simulation.pref().h() {
            let mut p_idx = ParasiteSpeciesIndex {
                species_index: species.pop().unwrap(),
                parasite_index: parasites.pop().unwrap(),
                match_count: 0
            };
            let match_score = find_match_score(&host, &all_parasites, &mut p_idx, simulation);
            parasites_exposed_to.insert((p_idx.species(), p_idx.parasite()), match_score);
            *species_match_score.entry(p_idx.species()).or_insert(0) += match_score;
            *host_match_score.entry(i).or_insert(0) += 1;
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
        simulation.update_species_not_exposed_to(i, species);
        simulation.pv(file_name, &format!("-------------------------\n"), true);
        if match_score_bellow_threshold >= simulation.pref().x() {
            simulation.kill_host(i);
            simulation.pv("host_dying_initial_exposure", &format!("{}\n", &simulation.hosts()[i].to_string()), true);
        }
    }
    let (t, r, w) = simulation.count_dead_hosts();
    simulation.pv("host_dying_initial_exposure", &format!("{} R({}) W({})", t, r, w), true);
    simulation.update_parasites_exposed_to(parasites_exposed_to);
    simulation.update_species_match_score(species_match_score);
    simulation.update_host_match_score(host_match_score);
}

pub fn additional_exposure(simulation: &mut Simulation) {
    let file_name = "host_additional_exposure";
    let mut species_match_score = simulation.species_match_score().clone();
    let mut host_match_score = simulation.simulation_state().host_match_score().clone();
    let (total_dead_hosts, _, _) = simulation.count_dead_hosts();
    simulation.pv(file_name, "Additional Exposure\n", true);
    let (_, alive_reservation_hosts, _) = simulation.count_alive_hosts();
    // secondary exposure
    let secondary_allowed = (simulation.pref().a() + simulation.pref().b()) as f32 * simulation.pref().m();
    if simulation.current_generation() <= simulation.pref().l()  {
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
    let mut parasites_exposed_to = simulation.parasites_exposed_to();
    let mut species_par_host = simulation.species_left();
    let mut parasites = (0..simulation.pref().e()).collect::<Vec<usize>>();
    let mut tried = vec![];
    while hosts_to_try < no_of_additional_host {
        let index: usize = rng.gen_range(0..simulation.pref().a() + simulation.pref().b());
        if tried.contains(&index) { continue }
        tried.push(index);
        simulation.pv(file_name,
                      &format!("host: {: >3}{: >4}      {}\nMATCHED PARASITES\n",
                               index, &simulation.hosts()[index].host_type(),
                               &simulation.hosts()[index].number_set()), true);
        // only reservation hosts
        if simulation.hosts()[index].host_type() == HostTypes::Wild { continue; }
        // skip dead hosts
        if !simulation.hosts()[index].alive() { continue; }
        hosts_to_try += 1;
        // expose to parasite
        let mut match_score_bellow_threshold = 0;
        for _ in 0..simulation.pref().i() {
            let mut p_idx = ParasiteSpeciesIndex {
                species_index: species_par_host.get_mut(&index).unwrap().pop().unwrap(),
                parasite_index: parasites.pop().unwrap(),
                match_count: 0
            };
            // let mut p_idx = get_random_parasite(simulation, &parasites_exposed_to);
            let match_score = find_match_score(&simulation.hosts()[index], &all_parasites, &mut p_idx, &simulation);
            // update simulation state
            parasites_exposed_to.insert((p_idx.species(), p_idx.parasite()), match_score);
            *species_match_score.entry(p_idx.species()).or_insert(0) += match_score;
            *host_match_score.entry(index).or_insert(0) += 1;
            //
            if match_score < simulation.pref().dd() {
                match_score_bellow_threshold += 1;
            }
            //
            let d = parasite_row(&all_parasites, p_idx.species(), p_idx.parasite_index);
            let p_grid = print_matching_number_sets(
                simulation.hosts()[index].number_set().clone(),
                d,
                p_idx.species(),
            );
            simulation.pv(file_name, &format!("{: >3},{: >3} : {: >2} -> {: >3}\n", p_idx.species(), p_idx.parasite(), match_score, p_grid), true);
            //
        }
        if match_score_bellow_threshold >= simulation.pref().cc() {
            simulation.kill_host(index);
            simulation.pv("host_dying_additional_exposure", &format!("{}\n", &simulation.hosts()[index].to_string()), true);
        }
        simulation.update_species_not_exposed_to(index, species_par_host.get(&index).unwrap().clone());
    }
    let (t, r, w) = simulation.count_dead_hosts();
    simulation.pv("host_dying_additional_exposure", &format!("{} R({}) W({})", t, r, w), true);
    simulation.update_parasites_exposed_to(parasites_exposed_to);
    simulation.update_species_match_score(species_match_score);
    simulation.update_host_match_score(host_match_score);
}

pub fn birth_hosts(simulation: &mut Simulation) {
    let file_name = "hosts_birth";
    let (dist, choices) = match simulation.program_version() {
        ProgramVersions::One => birth_generation_version_1(simulation),
        ProgramVersions::Two => birth_generation_version_2(simulation),
        ProgramVersions::Three => birth_generation_version_1(simulation),
        ProgramVersions::Four => birth_generation_version_1(simulation),
    };
    let mut rng = thread_rng();
    loop {
        // pick up parent host
        let random_host_index = choices[dist.sample(&mut rng)];
        let parent_index = match simulation.program_version() {
            ProgramVersions::One => {
                random_host_selection_v1(&simulation, random_host_index)
            }
            ProgramVersions::Two => {
                rng.gen_range(0..simulation.hosts().len())
            }
            ProgramVersions::Three => {
                random_host_selection_v1(&simulation, random_host_index)
            }
            ProgramVersions::Four => {
                rng.gen_range(0..simulation.hosts().len())
            }
        };

        let parent_host = simulation.hosts()[parent_index].clone();
        simulation.pv(file_name, &format!("{:3} {}\n", parent_index, parent_host), true);
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
    simulation.pv("hosts_birth", &format!("Birth V1: {:#?}, {:#?}\n", chances, choices), true);
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
    for ((species, parasites), match_score) in match_scores.iter() {
        *frequency.entry(*match_score).or_insert(0) += 1;
        // group parasites by score
        individuals_with_score.get_mut(&match_score).unwrap().push((species.clone(), parasites.clone()));
        // update cumulative frequency
        cumulative_frequency[*match_score] += 1;
    }
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
    let index = p_idx.species() * simulation.pref().g();

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

fn get_random_parasite(simulation: &mut Simulation, parasites_exposed_to: &HashMap<(usize, usize), usize>) -> ParasiteSpeciesIndex {

    loop {
        let p_idx = random_parasite(simulation.pref());
        let ky = (p_idx.species(), p_idx.parasite());
        if parasites_exposed_to.contains_key(&ky) { continue; }
        break p_idx;
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