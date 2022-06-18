use std::cmp::min;
use std::collections::HashMap;

use log::debug;
use rand::seq::SliceRandom;
use rand::thread_rng;

use crate::{additional_exposure_selected, find_match_score, HostTypes, kill_host_condition_matched, ParasiteSpeciesIndex, print_exposure_state, Simulation, update_exposure_states};

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
            update_exposure_states(simulation, &p_idx, host_index, match_score);
            if match_score < simulation.pref().n() {
                match_score_bellow_threshold += 1;
            }
            //
            let p_grid = print_exposure_state(&all_parasites, &p_idx);
            simulation.pv(file_name, &format!("{: >3},{: >3} : {: >2} -> {: >3}\n", p_idx.species(), p_idx.parasite(), match_score, p_grid), true);
            //
        }
        simulation.pv(file_name, &format!("-------------------------\n"), true);
        if kill_host_condition_matched(match_score_bellow_threshold, simulation, host_index, false) {
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
        debug!("additional exposure not selected");
        return;
    }
    debug!("additional exposure selected");
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
    debug!("exposing additional hosts");
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
        debug!("done additional exposure");
        simulation.pv(file_name, &format!("-------------------------\n"), true);
        let v = *simulation.ss().host_match_scores_bellow_dd().get(&host_index).unwrap();
        if kill_host_condition_matched(v, simulation, host_index, true) {
            simulation.kill_host(host_index);
            simulation.pv("host_dying_additional_exposure", &format!("{:3} {}\n", host_index, &host.to_string()), true);
        } else {
            simulation.ss_mut().add_hosts_tried(host_index);
        }
    }
    let (t, r, w) = simulation.count_dead_hosts();
    simulation.pv("host_dying_additional_exposure", &format!("{} R({}) W({})", t, r, w), true);
}
