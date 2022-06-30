use std::cmp::min;
use std::collections::HashMap;

use ndarray::{Array3, ArrayBase, Ix3, OwnedRepr};
use rand::seq::SliceRandom;
use rand::thread_rng;

use crate::{Host, HostTypes, ParasiteSpeciesIndex, ProgramVersions, Simulation};

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
        simulation.pv(file_name,
                      &format!("{: <3} {:12}{}\n(species, parasite): match count -> code, \n",
                               host_index, &host.host_type().to_string(), &host.number_set()),
                      true);
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
            update_exposure_states(simulation, &p_idx, host_index, match_score);
            let p_grid = print_exposure_state(&all_parasites, &p_idx);
            simulation.pv(file_name, &format!("{: >3},{: >3} : {: >2} -> {: >3}\n", p_idx.species(), p_idx.parasite(), match_score, p_grid), true);
        }
        simulation.pv(file_name, &format!("-------------------------\n"), true);
    }
    simulation.set_species_left(species_possible);
    simulation.set_parasites_possible(parasites_possible);
}

pub fn additional_exposure(simulation: &mut Simulation) {
    let file_name = "host_additional_exposure";
    let (total_dead_hosts, _, _) = simulation.count_dead_hosts();
    simulation.pv(file_name,
                  &format!(
                      "Additional Exposure\ncurrent generation: {}\n",
                      simulation.current_generation()
                  ), true);
    // secondary exposure
    if !additional_exposure_selected(simulation, total_dead_hosts) {
        return;
    }
    simulation.has_additional_exposure();
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
    let no_of_additional_host = min(
        alive_hosts.len(),
        (reservation_hosts_alive as f32 * simulation.pref().aa()).ceil() as usize,
    );
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
            update_exposure_states(simulation, &p_idx, host_index, match_score);
            let p_grid = print_exposure_state(&all_parasites, &p_idx);
            simulation.pv(file_name, &format!("{: >3},{: >3} : {: >2} -> {: >3}\n", p_idx.species(), p_idx.parasite(), match_score, p_grid), true);
            //
        }
        // adding hosts tried to check with death rule
        simulation.ss_mut().add_hosts_tried(host_index);
        simulation.pv(file_name, &format!("-------------------------\n"), true);
    }
}


pub fn update_exposure_states(simulation: &mut Simulation, p_idx: &ParasiteSpeciesIndex, host_index: usize, match_score: usize) {
    simulation.update_parasites_exposed_to((p_idx.species(), p_idx.parasite()), match_score);
    simulation.update_species_match_score(p_idx.species(), match_score);
    simulation.update_host_match_score_all(host_index, match_score);
    simulation.update_host_match_score(host_index, 1);
    simulation.update_host_match_score_bellow_n(host_index, if match_score < simulation.pref().n() { 1 } else { 0 });
    simulation.update_host_match_score_bellow_dd(host_index, if match_score < simulation.pref().dd() { 1 } else { 0 });
    simulation.update_host_match_score_bellow_j(host_index, if match_score < simulation.pref().j() { 1 } else { 0 });
}

pub fn print_exposure_state(all_parasites: &ArrayBase<OwnedRepr<usize>, Ix3>, p_idx: &ParasiteSpeciesIndex) -> String {
    let d = crate::parasite_row(&all_parasites, p_idx.species(), p_idx.parasite_index);
    crate::print_matching_number_sets(d, p_idx.species())
}
