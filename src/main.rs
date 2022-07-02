#![feature(iter_collect_into)]
#![feature(map_first_last)]
#![feature(drain_filter)]
#![feature(hash_drain_filter)]
#![feature(is_some_with)]
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
use std::path::Path;

use indicatif::{MultiProgress, ParallelProgressIterator, ProgressBar, ProgressStyle};
use log;
use log::{debug, info, LevelFilter, SetLoggerError, warn};
use log4rs::{Config, Handle};
use log4rs::append::console::ConsoleAppender;
use log4rs::append::file::FileAppender;
use log4rs::config::{Appender, Root};
use log4rs::encode::pattern::PatternEncoder;
use ndarray::{Array, Array1, Array2, Array3, ArrayBase, Axis, OwnedRepr};
use ndarray_rand::RandomExt;
use rand::{Rng, thread_rng};
use rand::distributions::{Uniform, WeightedIndex};
use rand::prelude::*;
use rayon::prelude::*;

use crate::exposure::{additional_exposure, expose_all_hosts_to_parasites};
use crate::hosts::{Host, HostTypes, print_hosts};
use crate::mutations::{mutate_hosts, mutate_parasites};
use crate::qi_calculations::calculate_qi;
use crate::simulation::{DeathRule, GGRunReport, HostsCount, new_simulation, print_parasites, ProgramVersions, ReportHostType, Simulation};
use crate::simulation_pref::SimulationPref;

pub mod simulation_pref;
pub mod hosts;
pub mod simulation;
pub mod mutations;
pub mod exposure;
pub mod qi_calculations;
pub mod host_death_rules;

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
    if Path::new("report").exists() {
        remove_dir_all("report").expect("Failed to remove report directory");
    };
    config_logger(1).unwrap();
    let param_file = args().nth(1).unwrap_or(format!("conf/params.conf"));
    let _ss = &fs::read_to_string(param_file).unwrap();
    let _pref = serde_ini::from_str(_ss);
    if _pref.is_err() {
        println!("{:?}", _pref.err());
        panic!(r#"Invalid input, please check.
        Floats should be: 0.XX not .XX
        Comments should not be in the same line with input variables
        "#);
    }
    let pref: SimulationPref = _pref.unwrap();
    //
    let program = if let Some(x) = args().nth(2) {
        ProgramVersions::from(x)
    } else {
        ProgramVersions::from(pref.program_version())
    };
    let death_rule = if let Some(x) = args().nth(3) {
        DeathRule::from(x)
    } else {
        DeathRule::from(pref.death_rule())
    };
    //
    // progressbar setup
    let pb_all_simulation = setup_progressbar(pref.gg());
    // timer
    let now = time::Instant::now();
    // starting up
    println!("Running version {} with Death Step {}, build 0.1.42_new_kill_algo", program, death_rule);
    let program_clone = program.clone();
    let pref_clone = pref.clone();
    let mut wild_hosts: Vec<Vec<usize>> = vec![];
    let mut reservation_hosts: Vec<Vec<usize>> = vec![];
    if option_env!("DEBUG").is_some() {
        (0..pref.gg()).into_iter().for_each(|gg| {
            let k = run_generation_step(program_clone, gg.clone(), pref_clone, death_rule);
            wild_hosts.push(k.0);
            reservation_hosts.push(k.1);
        });
    } else {
        (0..pref.gg()).into_par_iter().progress_with(pb_all_simulation).map(move |gg| {
            run_generation_step(program_clone, gg.clone(), pref_clone, death_rule)
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

fn setup_progressbar(gg: usize) -> ProgressBar {
    let multi_progress_bar = MultiProgress::new();
    let sty = ProgressStyle::with_template(
        "[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7} {msg}",
    ).unwrap().progress_chars("##-");
    let pb_all_simulation = multi_progress_bar.add(ProgressBar::new(gg as u64));
    pb_all_simulation.set_style(sty.clone());
    pb_all_simulation
}

fn generate_excel(hosts: &ArrayBase<OwnedRepr<usize>, ndarray::Ix2>, pref: &SimulationPref) -> Vec<String> {
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

fn run_generation_step(program: ProgramVersions, gg: usize, pref: SimulationPref, death_rule: DeathRule) -> (Vec<usize>, Vec<usize>) {
    debug!("starting simulation");
    let mut simulation = new_simulation(pref.clone(), program, gg, death_rule);
    let mut lines = vec![];
    let mut simulation_ended = false;
    (0..pref.ff()).into_iter().map(|_| {
        if simulation_ended {
            return fill_host(&mut simulation);
        }
        expose_all_hosts_to_parasites(&mut simulation);
        host_death_rules::v2::initial(&mut simulation);
        additional_exposure(&mut simulation);
        host_death_rules::v2::additional(&mut simulation);
        if !should_continue(&mut simulation) {
            simulation_ended = true;
            return fill_host(&mut simulation);
        }
        if !birth_hosts(&mut simulation) {
            simulation_ended = true;
            return fill_host(&mut simulation);
        }
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

pub fn birth_hosts(simulation: &mut Simulation) -> bool {
    let file_name = "hosts_birth";
    let result = match simulation.program_version() {
        ProgramVersions::One | ProgramVersions::Three => {
            birth_generation_version_1(simulation)
        }
        ProgramVersions::Two | ProgramVersions::Four
        | ProgramVersions::Five | ProgramVersions::Six
        | ProgramVersions::Seven | ProgramVersions::Eight
        | ProgramVersions::Nine | ProgramVersions::Ten => {
            calculate_qi(simulation);
            let qq = simulation.ss().qi_host_individual();
            let avg = qq.into_iter().fold(0., |mut a, (_, v)| {
                a += v;
                a
            }) / qq.len() as f32;
            info!("average qi: {}", avg);
            birth_generation_version_2(simulation)
        }
    };
    if result.is_err() {
        match result {
            Err(0) => {
                let o = format!("All chances are Zero, Please check input parameters. Ending at Generation, {}", simulation.current_generation());
                simulation.pv(file_name, &o, true);
                warn!("{}", o)
            }
            Err(1) => {
                let o = format!("Tied due to all negative chances {}", simulation.current_generation());
                simulation.pv(file_name, &o, true);
                warn!("{}", o);
            }
            _ => {}
        }
        return false;
    }
    let (dist, choices) = result.unwrap();
    simulation.pv(file_name, &format!("{:#?}\n{:?}\n", dist, choices), true);
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
    true
}

pub fn birth_generation_version_1(simulation: &mut Simulation) -> Result<(WeightedIndex<f32>, Vec<usize>), usize> {
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
    Ok((WeightedIndex::new(chances).unwrap(), choices))
}

pub fn get_chances_v1(v: f32, u: f32, w: f32, y: f32, t: f32, o: f32, p: f32) -> (f32, f32) {
    (
        (1. + o * y * w / u + (1. - o) * y * w / (t + u) + (1. - o) * y * v / (t + u) - p),  // for reservation host individuals
        (1. + o * y * v / t + (1. - o) * y * v / (t + u) + (1. - o) * y * w / (t + u)), // for wild host individuals
    )
}

pub fn get_chances_v2(host_type: HostTypes, qi: f32, qr: f32, qw: f32, lr: f32, lw: f32, o: f32, y: f32, p: f32) -> f32 {
    match host_type {
        HostTypes::Reservation => qi + (qi / qr) * o * y * lr + qi / (qr + qw) * (1. - o) * y * lr + qi / (qr + qw) * (1. - o) * y * lw - p,
        HostTypes::Wild => qi + qi / (qr + qw) * (1. - o) * y * lr + (qi / qw) * o * y * lw + qi / (qr + qw) * (1. - o) * y * lw
    }
}

pub fn birth_generation_version_2(simulation: &mut Simulation) -> Result<(WeightedIndex<f32>, Vec<usize>), usize> {
    let (_, no_of_dead_reservation_host, no_of_dead_wild_host) = simulation.count_dead_hosts();
    let qr = simulation.ss().qr();
    let qw = simulation.ss().qw();
    let mut qi_values = vec![];
    // select index of the hosts that are alive
    let choices: Vec<usize> = simulation.hosts_alive();
    // chances of each hosts that are alive for birth
    simulation.pv("hosts_birth", "Birth V2\n", true);
    simulation.pv("hosts_birth",
                  &format!("host_index,      host_type,    p,             qi,     chance\n"), true);
    let pv = |hi, ht, p, c, qi| -> String {
        format!("{: >10}, {: >14}, {: >4}, {: >14}, {: >10}\n", hi, ht, p, qi, c)
    };
    let chances: Vec<f32> = choices.iter().map(|host_index| {
        let qi = simulation.ss().qi_host_individual().get(host_index).unwrap().clone();
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
    if chances.iter().sum::<f32>() == 0. {
        return Err(0);
    }
    // if all are less than zero
    if chances.iter().all(|x| *x < 0.) {
        return Err(1);
    }
    Ok((WeightedIndex::new(chances).unwrap(), choices))
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
        cumulative_frequency[*match_score] += 1;
        individuals_with_score.get_mut(&match_score).unwrap().push((species.clone(), parasites.clone()));
    });
    let mut rng = thread_rng();
    _s.push_str(&format!("{:width$} {:10}\n", "killed_individual", "parent_individual", width = 26 + simulation.pref().g()));
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
                let p1 = parasite_row(simulation.parasites(), *s, *p);
                let p2 = parasite_row(simulation.parasites(), *s, parent_parasite_index);
                _s.push_str(&format!("({:3}, {:3}) {}        ({:3}, {:3}) {}\n", s, p, p1, s, parent_parasite_index, p2));
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
    simulation.pv("replaced_parasites", &format!("{:?}\n", species_match_score), true);

    while replaced < to_be_replaced {
        let ky = max_keys.pop();
        if ky.is_none() { break; }
        let (score, species_index) = ky.unwrap();
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
        s.push_str(&format!("REPLACE:\n{}: Score({})\n{:#?}\n", *species_index, score, species));
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
