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
    let mut rng = thread_rng();
    let d = simulation.pref().d();
    let g = simulation.pref().g();
    let e = simulation.pref().e();
    let f = simulation.pref().f();
    let total_parasites = d * e;
    let total_digits = total_parasites * g;
    let changed_digits = mutate_digits(simulation.pref().k(), total_digits);
    let mut mutated_parasites = String::new();
    for index in changed_digits.iter() {
        let (spc, ind, col) = find_parasite_index(e, g, index.clone());
        let k = simulation.parasites().index_axis(Axis(0), spc);
        mutated_parasites.push_str(&format!("({:3}, {:3}) {}\n", spc, ind, k.index_axis(Axis(0), ind)));
        simulation.parasites_mut()[[spc, ind, col]] = loop {
            let ll = rng.gen_range(0..f);
            if ll != simulation.parasites_mut()[[spc, ind, col]] {
                break ll;
            }
        };
        let k = simulation.parasites().index_axis(Axis(0), spc);
        mutated_parasites.push_str(&format!(" {:3}  {:3}  {}\n\n", " ", " ", k.index_axis(Axis(0), ind)));
    }
    let l = changed_digits.len();
    simulation.pv("mutated_parasites", &format!("Total Mutation: {} {}%\n", l, l as f64 * 100. / total_digits as f64), true);
    simulation.pv("mutated_parasites", &format!("{}\n", mutated_parasites), true);
}

pub fn mutate_hosts(simulation: &mut Simulation) {
    let mut rng = thread_rng();
    let c = simulation.pref().c();
    let f = simulation.pref().f();
    let ee = simulation.pref().ee();
    let total_hosts = simulation.pref().a() + simulation.pref().b();
    let total_digits = total_hosts * c;
    let changed_digits = mutate_digits(ee, total_digits);
    // setting weights for all the possible digits,
    let mut mutated_hosts = String::new();
    for index in changed_digits.iter() {
        let host_index = (*index as f32 / c as f32).floor() as usize;
        mutated_hosts.push_str(&format!("{:3} {}\n", host_index, simulation.hosts()[host_index]));
        // address k in hosts
        let which_digit = (*index as f32 % c as f32) as usize;
        let host_num_set = simulation.hosts_mut()[host_index].number_set_mut();
        host_num_set[which_digit] = loop {
            let ll = rng.gen_range(0..f);
            if ll != host_num_set[which_digit] {
                break ll;
            }
        };
        mutated_hosts.push_str(&format!("{:3} {}\n\n", " ", simulation.hosts()[host_index]));
    }
    let l = changed_digits.len();
    simulation.pv("mutated_hosts", &format!("total mutation: {} {}%\n", l, l as f64 * 100. / total_digits as f64), true);
    simulation.pv("mutated_hosts", &mutated_hosts, true);
}

fn mutate_digits(mutation_rate: f32, total_digits: usize) -> Vec<usize>{
    let mut rng = thread_rng();
    if mutation_rate >= 1. {
        panic!("Mutation rate should be bellow 1")
    }
    let other = 1. - mutation_rate;
    let weights: Vec<f32> = vec![other, mutation_rate];
    let dist = WeightedIndex::new(&weights).unwrap();
    let mut to_be_changed = vec![];
    for index in 0..total_digits {
        let change = dist.sample(&mut rng);
        if change > 0 {
            to_be_changed.push(index);
        }
    }
    to_be_changed
}
