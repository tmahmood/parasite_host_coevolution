use std::collections::HashMap;

use ndarray::Axis;
use rand::distributions::WeightedIndex;
use rand::prelude::*;
use rand::Rng;

use crate::Simulation;

fn find_parasite_index(e: usize, g: usize, x: usize) -> (usize, usize, usize) {
    // row
    let pi = x / g;
    // which group of row, consisting of e rows per group
    let s = pi / e;
    //
    (s, pi % e, x % g)
}

pub fn mutate_parasites(simulation: &mut Simulation) {
    let mut rng = rand::thread_rng();
    let d = simulation.pref().d();
    let g = simulation.pref().g();
    let e = simulation.pref().e();
    let f = simulation.pref().f();
    let total_parasites = d * e;
    let total_digits = total_parasites * g;
    let weights: Vec<f32> = vec![simulation.pref().k(); total_digits];
    let dist = WeightedIndex::new(&weights).unwrap();
    let mut mutated_parasites = String::new();
    let mut already_changed = HashMap::new();
    for _ in 0..total_parasites {
        let index = dist.sample(&mut rng);
        let (spc, ind, col) = find_parasite_index(e, g, index);
        if already_changed.contains_key(&ind) {
            continue;
        }
        let k = simulation.parasites().index_axis(Axis(0), spc);
        mutated_parasites.push_str(&format!("({:3}, {:3}) {}\n", spc, ind, k.index_axis(Axis(0), ind)));
        // address k in hosts
        simulation.parasites_mut()[[spc, ind, col]] = loop {
            let ll = rng.gen_range(0..f);
            if ll != simulation.parasites_mut()[[spc, ind, col]] {
                break ll;
            }
        };
        already_changed.insert(ind, 1);
        let k = simulation.parasites().index_axis(Axis(0), spc);
        mutated_parasites.push_str(&format!(" {:3}  {:3}  {}\n\n", " ", " ", k.index_axis(Axis(0), ind)));
    }
    simulation.pv("mutated_parasites", &format!("Total Mutation: {}\n", already_changed.len()), true);
    simulation.pv("mutated_parasites", &format!("{}\n", mutated_parasites), true);
}

pub fn mutate_hosts(simulation: &mut Simulation) {
    let c = simulation.pref().c();
    let f = simulation.pref().f();
    let ee = simulation.pref().ee();
    let total_hosts = simulation.pref().a() + simulation.pref().b();
    let total_digits = total_hosts * c;
    // setting weights for all the possible digits,
    let weights: Vec<f32> = vec![ee; total_digits];
    let dist = WeightedIndex::new(&weights).unwrap();
    let mut rng = rand::thread_rng();
    let mut mutated_hosts = String::new();
    //
    let mut already_changed = HashMap::new();
    for _ in 0..total_hosts {
        let index = dist.sample(&mut rng);
        let host_index = (index as f32 / c as f32).floor() as usize;
        if already_changed.contains_key(&host_index) {
            continue;
        }
        mutated_hosts.push_str(&format!("{:3} {}\n", host_index, simulation.hosts()[host_index]));
        // address k in hosts
        let which_digit = (index as f32 % c as f32) as usize;
        let host_num_set = simulation.hosts_mut()[host_index].number_set_mut();
        host_num_set[which_digit] = loop {
            let ll = rng.gen_range(0..f);
            if ll != host_num_set[which_digit] {
                break ll;
            }
        };
        already_changed.insert(host_index, 1);
        mutated_hosts.push_str(&format!("{:3} {}\n\n", " ", simulation.hosts()[host_index]));
    }
    simulation.pv("mutated_hosts", &format!("total mutation: {}\n", already_changed.len()), true);
    simulation.pv("mutated_hosts", &mutated_hosts, true);
}
