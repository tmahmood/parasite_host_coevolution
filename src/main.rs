extern crate serde;
extern crate serde_derive;
extern crate serde_ini;
extern crate rand;
extern crate core;

use std::collections::HashMap;
use std::fmt::{Display, Formatter};
use std::fs;
use futures::executor::block_on;
use ndarray::{Array, Array1, ArrayBase, Axis, Ix, Ix1, Ix3, OwnedRepr};
use ndarray_rand::{RandomExt};
use rand::distributions::{Uniform, WeightedIndex};
use rand::prelude::*;
use rand::{Rng, thread_rng};
use parasites::random_parasite;
use crate::hosts::{Host, HostTypes};
use crate::simulation::{new_simulation, Simulation};
use crate::simulation_pref::SimulationPref;

// It's important to shuffle lists before the random processes. In the past when I've hired for
// these sorts of programs, the biggest problem has been that individuals of one type were more
// likely to get chosen for reproduction or death because lists weren't shuffled, and that causes
// the random processes to be biased.

pub mod a2d;
pub mod simulation_pref;
pub mod parasites;
pub mod hosts;
pub mod simulation;


async fn main_async() {
    // TODO Maybe Need better way to feed the input file
    let s: SimulationPref = serde_ini::from_str(&fs::read_to_string("params.conf").unwrap()).unwrap();
    let mut simulation = new_simulation(s.clone()).await;
    println!("#S1##############");
    expose_all_hosts_to_parasites(&mut simulation).await;
    println!("#S2##############");
    additional_exposure(&mut simulation);
    println!("#S3##############");
    birth_hosts(&mut simulation);
    println!("#S4##############");
    parasite_truncation_and_birth(&mut simulation);
    println!("#S5##############");
    mutation(&mut simulation);
    println!("#S6##############");
    parasite_replacement(&mut simulation);
}

fn parasite_replacement(simulation: &mut Simulation) {
    let mut replaced = 0;
    let to_be_replaced = simulation.pref().q();
    // find the max value
    let species_match_score = *simulation.species_match_score().iter().max_by(|a, b| a.1.cmp(&b.1)).map(|(k, v)| v).unwrap();
    let mut max_keys: Vec<usize> = simulation.species_match_score().iter().filter(|v| *v.1 == species_match_score).map(|v| *v.0).collect();
    let i = generate_individual(simulation.pref().f(), simulation.pref().g());
    while replaced < to_be_replaced && max_keys.len() > 0{
        let ky = max_keys.pop().unwrap();
        let mut species = simulation.parasites_mut().index_axis_mut(Axis(0), ky);
        for (ii, iv) in i.iter().enumerate() {
            for mut row in species.rows_mut() {
                row[ii] = *iv;
            }
        }
        replaced += 1;
    }
    println!("{:?} {:?}", simulation.species_match_score(), species_match_score);
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
        let old = host.number_set().clone();
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
        let new = host.number_set();
        println!("{} {} {}", old, new, changes);
    }

    let weights: Vec<f32> = (0..simulation.pref().f()).map(|_| simulation.pref().k()).collect();
    let dist = WeightedIndex::new(&weights).unwrap();
    for mut parasite in simulation.parasites_mut().rows_mut() {
        let old = parasite.to_owned();
        let mut changes = 0;
        for cc in 0..g {
            let k = choices[dist.sample(&mut rng)];
            if k == 1 {
                parasite[cc] = rng.gen_range(0..f);
                if old != parasite {
                    changes += 1;
                }
            }
        }
        println!("{} {} {}", old, parasite, changes);
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
        individuals_with_score.insert(k, Vec::new());
    }
    // calculate frequencies of each score
    for ((s, p), v) in match_scores.iter() {
        println!("{: >2} {: >3} {: >2}", s, p, v);
        if !frequency.contains_key(&v) {
            frequency.insert(*v, 0);
        }
        *frequency.get_mut(&v).unwrap() += 1usize;
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
            for (s, i) in parasites {
                // get random existing parasite from the same species (s), excluding this parasite (i)
                let parent_parasite_index = loop {
                    let i1 = rng.gen_range(0..simulation.pref().e());
                    if i1 != *i { break i1; }
                };
                let b = simulation.parasites().index_axis(Axis(0), *s);
                let v = b.index_axis(Axis(0), *i);
                print!("{: >2} {: >3}: {: >2} : ", s, i, v);
                let v = b.index_axis(Axis(0), parent_parasite_index);
                print!("{: >2} -> ", v);
                simulation.update_parasites(*s, *i, parent_parasite_index);
                let b = simulation.parasites().index_axis(Axis(0), *s);
                let v = b.index_axis(Axis(0), *i);
                println!("{: >2}", v);
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

pub async fn expose_all_hosts_to_parasites(simulation: &mut Simulation) {
    let mut parasites_exposed_to = HashMap::new();
    let mut species_match_score = HashMap::<usize, usize>::new();
    let all_parasites = simulation.parasites().clone();
    let hosts = simulation.hosts().clone();
    for i in 0..simulation.pref().a() + simulation.pref().b() {
        let host = &hosts[i];
        println!("            {}", host.number_set());
        let mut match_score_bellow_threshold = 0;
        for _ in 0..simulation.pref().h() {
            let mut p_idx = get_random_parasite(simulation, &parasites_exposed_to);
            let match_score = find_match_score(&host, &all_parasites, &mut p_idx, simulation.pref());
            parasites_exposed_to.insert((p_idx.species(), p_idx.parasite()), match_score);
            *species_match_score.entry(p_idx.species()).or_insert(0) += match_score;
            if match_score < simulation.pref().n() {
                match_score_bellow_threshold += 1;
            }
        }
        println!("----------------------");
        if match_score_bellow_threshold >= simulation.pref().x() {
            simulation.kill_host(i);
        }
    }
    println!("species match score: {:?}", species_match_score);
    simulation.update_parasites_exposed_to(parasites_exposed_to);
    simulation.update_species_match_score(species_match_score);
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
    let (dead_reservation_hosts, _, total_dead_hosts) = simulation.count_dead_hosts();
    let alive_reservation_hosts = simulation.pref().a() - dead_reservation_hosts;
    // secondary exposure
    let secondary_allowed = (simulation.pref().a() + simulation.pref().b()) as f32 * simulation.pref().m();
    if simulation.current_generation() >= simulation.pref().l() && total_dead_hosts < secondary_allowed as usize {
        println!("Skipping");
        return;
    }
    let all_parasites = simulation.parasites().clone();
    let no_of_additional_host = (alive_reservation_hosts as f32 * simulation.pref().aa()).ceil() as usize;
    println!("# host {}", no_of_additional_host);
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
        println!("            {}", simulation.hosts()[index].number_set());
        let mut match_score_bellow_threshold = 0;
        for _ in 0..simulation.pref().i() {
            // we should not expose the host to same parasite more than once in current generation
            // so we will pick random parasite till we get unused one
            let mut p_idx = get_random_parasite(simulation, &parasites_exposed_to);
            let match_score = find_match_score(&simulation.hosts()[index], &all_parasites, &mut p_idx, simulation.pref());
            parasites_exposed_to.insert((p_idx.species(), p_idx.parasite()), match_score);
            *species_match_score.entry(p_idx.species()).or_insert(0) += match_score;

            if match_score < simulation.pref().dd() {
                match_score_bellow_threshold += 1;
            }
        }
        println!("----------------------");
        if match_score_bellow_threshold >= simulation.pref().cc() {
            simulation.kill_host(index);
        }
    }
    simulation.update_parasites_exposed_to(parasites_exposed_to);
    simulation.update_species_match_score(species_match_score);
}

pub fn birth_hosts(simulation: &mut Simulation) {
    let (no_of_dead_reservation_host, no_of_dead_wild_host, _) = simulation.count_dead_hosts();
    let no_of_reservation_host_alive = simulation.pref().a() - no_of_dead_reservation_host;
    let no_of_wild_host_alive = simulation.pref().b() - no_of_dead_wild_host;

    let (chance_reservation, chance_wild) = get_chances(
        no_of_dead_wild_host as f32,
        no_of_reservation_host_alive as f32,
        no_of_dead_reservation_host as f32,
        simulation.pref().y(),
        no_of_wild_host_alive as f32,
        simulation.pref().o(),
        simulation.pref().p(),
    );
    println!("{} {}\n----------------------", chance_reservation, chance_wild);
    let mut rng = thread_rng();
    let mut new_born = 0;
    let choices: Vec<(usize, f32)> = simulation.hosts().iter().enumerate().filter(|v| v.1.alive()).map(|(ii, v)| {
        match v.host_type() {
            HostTypes::Reservation => (ii, chance_reservation),
            HostTypes::Wild => (ii, chance_wild)
        }
    }).collect();
    let dist = WeightedIndex::new(choices.iter().map(|item| item.1)).unwrap();
    let (mut count_dead_r, mut count_dead_w, mut total_dead_hosts) = simulation.count_dead_hosts();
    loop {
        // pick up parent host
        let parent_index: usize = choices[dist.sample(&mut rng)].0;
        let parent_host = simulation.hosts()[parent_index].clone();
        if !parent_host.alive() { continue; }
        if count_dead_r <= 0 && parent_host.host_type() == HostTypes::Reservation {
            continue;
        }
        if count_dead_w <= 0 && parent_host.host_type() == HostTypes::Wild {
            continue;
        }
        let mut index = None;
        for (ii, host) in simulation.hosts().iter().enumerate() {
            if !host.alive() && host.host_type() == parent_host.host_type() {
                index = Some(ii);
                break;
            }
        }
        if index == None {
            panic!("I couldn't find any dead [{}] hosts!", parent_host.host_type());
        }
        let host_index = index.unwrap();
        // now we get the host
        let k = simulation.update_dead_host(host_index, parent_index);
        count_dead_r = k.0;
        count_dead_w = k.1;
        total_dead_hosts = k.2;
        new_born += 1;
        println!("new_born: {: >3} {: >3}, total_dead: {: >3}, dead_r: {: >3}, dead_w: {: >3}", new_born, parent_host, total_dead_hosts, count_dead_r, count_dead_w);
        if total_dead_hosts == 0 { break; }
    }
}

pub fn get_chances(v: f32, u: f32, w: f32, y: f32, t: f32, o: f32, p: f32) -> (f32, f32) {
    (
        (1. + o * y * w / u + (1. - o) * y * w / (t + u) + (1. - o) * y * v / (t + u) - p) as f32,  // for reservation host individuals
        (1. + o * y * v / t + (1. - o) * y * v / (t + u) + (1. - o) * y * w / (t + u)) as f32, // for wild host individuals
    )
}

pub fn find_match_score(host: &Host, all_parasites: &ArrayBase<OwnedRepr<usize>, Ix3>, p_idx: &mut ParasiteSpeciesIndex, pref: SimulationPref) -> usize {
    let mut match_count = 0;
    let number_set = host.number_set();
    let mut p = String::new();
    let index = p_idx.species();
    print!(" ");
    for ii in index..index + pref.g() {
        p.push_str(&format!("{}, ", all_parasites[[index, ii - index, 0]]));
        if number_set[ii] == all_parasites[[index, ii - index, 0]] {
            match_count += 1
        }
    }
    print!("{} | {}->{} || ", match_count, index, index + pref.g());
    for _ in 0..index {
        print!("   ");
    }
    println!("{}", p.trim());
    p_idx.match_count = match_count;
    match_count
}