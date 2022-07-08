use crate::{DeathRule, HostsCount, HostTypes, ProgramVersions, Simulation, SimulationPref};
use crate::host_death_rules::v2::v1;
use crate::hosts::create_random_hosts;
use crate::simulation::{create_random_parasites, SimulationState};

fn create_test_simulation(conf: &str) -> Simulation {
    let pref: SimulationPref = serde_ini::from_str(conf).unwrap();
    let hosts = create_random_hosts(&pref);
    let parasites = create_random_parasites(&pref);
    let initial_simulation_state = SimulationState::new(
        hosts, parasites, HostsCount {
            reservation_host: pref.a(),
            wild_host: pref.b(),
        },
    );
    Simulation::new(
        pref, initial_simulation_state,
        ProgramVersions::from(pref.program_version()), 0,
        DeathRule::from(pref.death_rule()),
    )
}

fn set_host_match_scores(simulation: &mut Simulation, scores: &Vec<Vec<i32>>) {
    for host_index in 0..scores.len() {
        for k in scores[host_index].iter() {
            let match_score = *k as usize;
            simulation.update_host_match_score_all(host_index, match_score);
            simulation.update_host_match_score(host_index, 1);
            simulation.update_host_match_score_bellow_n(host_index, if match_score < simulation.pref().n() { 1 } else { 0 });
            simulation.update_host_match_score_bellow_dd(host_index, if match_score < simulation.pref().dd() { 1 } else { 0 });
            simulation.update_host_match_score_bellow_j(host_index, if match_score < simulation.pref().j() { 1 } else { 0 });
        }
    }
}

fn print_log(v: &Simulation) {
    for (f, s) in v.log_files().iter() {
        println!("{:14} ->\n{}", f, s);
    }
}

fn print_simulations_states_s(v1_s: &Simulation) {
    print_log(v1_s);
    println!("n:{:?}\ndd:{:?}\nj:{:?}\nall:{:?}",
             v1_s.ss().host_match_scores_bellow_n(),
             v1_s.ss().host_match_scores_bellow_dd(),
             v1_s.ss().host_match_scores_bellow_j(),
             v1_s.ss().host_match_scores_all(),
    );
    for i in 0..4 {
        println!("{}: {:14}", i, v1_s.hosts()[i].to_string());
    }
}


pub mod test_death_rules;
pub mod test_qi_calculations;