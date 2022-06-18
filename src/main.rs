#![feature(iter_collect_into)]
#![feature(map_first_last)]
#![feature(drain_filter)]
#![feature(is_some_with)]
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
use std::path::Path;

use indicatif::{MultiProgress, ParallelProgressIterator, ProgressBar, ProgressStyle};
use log;
use log::{debug, info, LevelFilter, SetLoggerError, warn};
use log4rs::{Config, Handle};
use log4rs::append::console::ConsoleAppender;
use log4rs::append::file::FileAppender;
use log4rs::config::{Appender, Root};
use log4rs::encode::pattern::PatternEncoder;
use ndarray::{Array, Array1, Array2, Array3, ArrayBase, Axis, Ix2, Ix3, OwnedRepr};
use ndarray_rand::RandomExt;
use rand::{Rng, thread_rng};
use rand::distributions::{Uniform, WeightedIndex};
use rand::prelude::*;
use rayon::prelude::*;

use crate::exposure::{additional_exposure, expose_all_hosts_to_parasites};
use crate::hosts::{Host, HostTypes, print_hosts};
use crate::mutations::{mutate_hosts, mutate_parasites};
use crate::simulation::{DeathRule, GGRunReport, HostsCount, new_simulation, print_parasites, ProgramVersions, ReportHostType, Simulation};
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
    if Path::new("report").exists() {
        remove_dir_all("report").expect("Failed to remove report directory");
    };
    config_logger(1).unwrap();
    let param_file = args().nth(1).unwrap_or(format!("conf/params.conf"));
    let pref: SimulationPref = serde_ini::from_str(&fs::read_to_string(param_file).unwrap()).unwrap();
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
    println!("Running version {} with Death Step {}, build 0.1.38_adding_another_4_versions", program, death_rule);
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

fn run_generation_step(program: ProgramVersions, gg: usize, pref: SimulationPref, death_rule: DeathRule) -> (Vec<usize>, Vec<usize>) {
    debug!("starting simulation");
    let mut simulation = new_simulation(pref.clone(), program, gg, death_rule);
    let mut lines = vec![];
    let mut simulation_ended = false;
    (0..pref.ff()).into_iter().map(|_| {
        if simulation_ended {
            debug!("all hosts, or one type of host killed");
            return fill_host(&mut simulation);
        }
        debug!("exposing hosts to parasite");
        expose_all_hosts_to_parasites(&mut simulation);
        additional_exposure(&mut simulation);
        if !should_continue(&mut simulation) {
            simulation_ended = true;
            return fill_host(&mut simulation);
        }
        debug!("host birth");
        if !birth_hosts(&mut simulation) {
            simulation_ended = true;
            return fill_host(&mut simulation);
        }
        debug!("parasite truncation and birth");
        parasite_truncation_and_birth(&mut simulation);
        debug!("parasite mutation");
        mutation(&mut simulation);
        debug!("parasite replacement");
        parasite_replacement(&mut simulation);
        let _s = print_parasites(&simulation.parasites());
        simulation.pv("parasites_at_end", &_s, true);
        let _s = print_hosts(simulation.hosts());
        simulation.pv("hosts_at_end", &_s, true);
        debug!("starting next generation");
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

pub mod exposure;

fn kill_host_condition_matched(match_score_bellow_threshold: usize, simulation: &mut Simulation, host_index: usize, is_additional: bool) -> bool {
    match simulation.death_rule() {
        DeathRule::Default => kill_host_condition_v1(match_score_bellow_threshold, simulation, is_additional, host_index),
        DeathRule::VersionOne => kill_host_condition_v2(match_score_bellow_threshold, simulation, is_additional, host_index)
    }
}

// Default kill condition
pub fn kill_host_condition_v1(match_score_bellow_threshold: usize, simulation: &mut Simulation, is_additional: bool, _: usize) -> bool {
    if is_additional {
        match_score_bellow_threshold >= simulation.pref().cc()
    } else {
        match_score_bellow_threshold >= simulation.pref().x()
    }
}

// Death to a reservation individual if the total number of unmatched digits is above SS (another
// variable) and death to a wild individual if the total number of unmatched digits is above
// TT (another variable).
pub fn kill_host_condition_v2(match_score_bellow_threshold: usize, simulation: &mut Simulation, _: bool, host_index: usize) -> bool {
    match simulation.host_type(host_index) {
        HostTypes::Reservation => match_score_bellow_threshold > simulation.pref().ss(),
        HostTypes::Wild => match_score_bellow_threshold > simulation.pref().tt()
    }
}

fn print_exposure_state(all_parasites: &ArrayBase<OwnedRepr<usize>, Ix3>, p_idx: &ParasiteSpeciesIndex) -> String {
    let d = parasite_row(&all_parasites, p_idx.species(), p_idx.parasite_index);
    print_matching_number_sets(d, p_idx.species())
}

fn update_exposure_states(simulation: &mut Simulation, p_idx: &ParasiteSpeciesIndex, host_index: usize, match_score: usize) {
    simulation.update_parasites_exposed_to((p_idx.species(), p_idx.parasite()), match_score);
    simulation.update_species_match_score(p_idx.species(), match_score);
    simulation.update_host_match_score(host_index, 1);
    simulation.update_host_match_score_bellow_dd(host_index, if match_score < simulation.pref().dd() { 1 } else { 0 });
    simulation.update_host_match_score_bellow_j(host_index, if match_score < simulation.pref().j() { 1 } else { 0 });
    simulation.update_host_match_score_bellow_oo(host_index, if match_score < simulation.pref().oo() { Some(match_score) } else { None });
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

struct QiParams<'a> {
    hh: usize,
    score: usize,
    simulation: &'a Simulation,
    host_index: usize,
}

fn calculate_qi(simulation: &mut Simulation) {
    let mut qr = 0.;
    let mut qw = 0.;
    let match_score_bellow_j = simulation.ss().host_match_scores_bellow_j().clone();
    let hosts = simulation.hosts().clone();

    let callback = match simulation.program_version() {
        ProgramVersions::Two | ProgramVersions::Four => calculate_qi_v1,
        ProgramVersions::Five | ProgramVersions::Six => calculate_qi_v2,
        ProgramVersions::Seven | ProgramVersions::Eight => calculate_qi_v3,
        ProgramVersions::Nine | ProgramVersions::Ten => calculate_qi_v4,
        _ => panic!("Should not come here!")
    };

    match_score_bellow_j.iter()
        .filter(|(index, _)| hosts[**index].alive())
        .for_each(|(index, score)| {
            let hh = simulation.pref().hh();
            let qi_params = QiParams {
                hh,
                score: *score,
                simulation,
                host_index: *index,
            };
            let mut qi = callback(qi_params);
            if qi < 0. { qi = 0. }
            simulation.set_host_individual_qi(*index, qi);
            match simulation.host_type(*index) {
                HostTypes::Reservation => qr += qi,
                HostTypes::Wild => qw += qi
            }
        });
    simulation.ss_mut().set_qr(qr);
    simulation.ss_mut().set_qw(qw);
}

fn calculate_qi_v1(qi_params: QiParams) -> f32 {
    qi_params.hh as f32 - qi_params.score as f32
}

// Suppose: An individual's match scores are 4, 6, and 9. And OO=7. HH=12. JJ=0.1. PP=10. And
// suppose there are C=36=length of set of numbers associated with a host.
//
// First way -> NN=(OO=7-4)+(OO=7-6)=4. Then Qi=(HH=12)*(1-JJ=.1*NN=4)=7.2.
// Second way -> 4 < 7, and 6 < 7, but 9 > 7. NN=2. Then Qi=(HH=12)*(1-JJ=.1*NN=2)=9.6.
// Third way -> NN=(C=36)-(4+6+9)-(PP=10)=7. Then Qi=(HH=12)*(1-JJ=.1*NN=7)=3.6.
//
// And the new death rule=Suppose SS=12 and TT=19
//
// if the individual is reservation, (C=36)- (4+6+9)=17 > SS=12. The individual dies.
// if the individual is wild, (C=36)-(4+6+9)=17 < TT=19. The individual lives.


/**
Same as version 2 except Qi ->Q_subscript_i=Qi for each surviving host individual is the
 variable HH times (1-JJ*NN), where JJ and NN are new variables different from J and N and
 JJ is a given input and NN is the total number of digits that each match score is under OO,
 another new variable.
 **/
fn calculate_qi_v2(qi_params: QiParams) -> f32 {
    let match_scores = qi_params.simulation.ss().match_scores_bellow_oo_for_host(qi_params.host_index);
    let oo = qi_params.simulation.pref().oo();
    let jj = qi_params.simulation.pref().jj();
    let nn = match_scores.iter()
        .filter(|v| **v < oo)
        .fold(0., |mut a, v| {
            a += oo as f32 - *v as f32;
            a
        });
    qi_params.hh as f32 * (1. - (jj as f32 * nn))
}

// 3) same as version 2 except Qi ->Q_subscript_i=Qi for each surviving host individual is the
// variable HH times (1-JJ*NN), where JJ is a given input and NN is the total number of match
// scores under OO.
fn calculate_qi_v3(qi_params: QiParams) -> f32 {
    let match_scores = qi_params.simulation.ss().match_scores_bellow_oo_for_host(qi_params.host_index);
    let nn = match_scores.iter()
        .filter(|v| **v < qi_params.simulation.pref().oo())
        .count();
    qi_params.hh as f32 * (1. - (qi_params.simulation.pref().jj() as f32 * nn as f32))
}

// A version with: Qi for each surviving host individual is the variable HH times (1-JJ*NN),
// where JJ is a given input and NN is the number of unmatched digits above PP (another variable).
fn calculate_qi_v4(qi_params: QiParams) -> f32 {
    let match_scores = qi_params.simulation.ss().match_scores_bellow_oo_for_host(qi_params.host_index);
    let pp = qi_params.simulation.pref().pp();
    let jj = qi_params.simulation.pref().jj();
    let c = qi_params.simulation.pref().c();

    let nn = c as f32 - match_scores.iter()
        .fold(0., |mut a, v| {
            a += *v as f32;
            a
        }) - pp as f32;
    qi_params.hh as f32 * (1. - (jj as f32 * nn))
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

pub fn find_match_score(host: &Host, all_parasites: &Array3<usize>, p_idx: &mut ParasiteSpeciesIndex, simulation: &Simulation) -> usize {
    let mut match_count = 0;
    let number_set = host.number_set();
    let start_at = p_idx.species() * simulation.pref().g();

    for ii in start_at..start_at + simulation.pref().g() {
        match simulation.program_version() {
            ProgramVersions::One | ProgramVersions::Two | ProgramVersions::Five | ProgramVersions::Seven | ProgramVersions::Nine => {
                if number_set[ii] == all_parasites[[p_idx.species(), p_idx.parasite_index, ii - start_at]] {
                    match_count += 1
                }
            }
            ProgramVersions::Three | ProgramVersions::Four | ProgramVersions::Six | ProgramVersions::Eight | ProgramVersions::Ten => {
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
