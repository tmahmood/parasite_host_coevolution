use rand::distributions::WeightedIndex;
use rand::prelude::*;
use rand::Rng;

use crate::Simulation;

pub fn mutate_parasites(simulation: &mut Simulation) {
    let mut rng = rand::thread_rng();
    let g = simulation.pref().g();
    let e = simulation.pref().e();
    let f = simulation.pref().f();
    let weights: Vec<f32> = (0..simulation.pref().f()).map(|_| simulation.pref().k()).collect();
    let dist = WeightedIndex::new(&weights).unwrap();
    let choices = [0, 1];
    let mut no_changes = 0;
    let bound = (1. / simulation.pref().k()).round() as i32;
    let mut cnt = 0;
    let mut mutated_parasites = String::new();
    simulation
        .parasites_mut()
        .rows_mut()
        .into_iter()
        .for_each(|mut parasite| {
            for cc in 0..g {
                let k = choices[dist.sample(&mut rng)];
                if k == 1 && no_changes > bound {
                    let spc = (cnt as f32 / e as f32).floor();
                    mutated_parasites.push_str(&format!("({:3}, {:3}) {} ", spc, cnt % e, parasite));
                    parasite[cc] = loop {
                        let ll = rng.gen_range(0..f);
                        if ll != parasite[cc] {
                            break ll;
                        }
                    };
                    mutated_parasites.push_str(&format!("{}\n", parasite));
                    no_changes = 0;
                    break;
                }
                no_changes += 1;
            }
            cnt += 1;
        });
    simulation.pv("mutated_parasites", &format!("{}\n", mutated_parasites), true);
}

pub fn mutate_hosts(simulation: &mut Simulation) {
    let c = simulation.pref().c();
    let f = simulation.pref().f();
    let weights: Vec<f32> = (0..simulation.pref().f()).map(|_| simulation.pref().ee()).collect();
    let dist = WeightedIndex::new(&weights).unwrap();
    let choices = [0, 1];
    let mut rng = rand::thread_rng();
    let mut mutated_hosts = String::new();
    let bound = (1. / simulation.pref().ee()).round() as i32;
    let mut no_changes = 0i32;
    //
    simulation
        .hosts_mut()
        .iter_mut()
        .enumerate()
        .for_each(|(i, host)| {
            //
            if (no_changes + c as i32) < 1000 {
                // jump to next host, as we are sure no mutation can occur
                no_changes += c as i32;
                return;
            }
            let mut changed = false;
            let m = host.number_set_mut();
            for cc in 0..c {
                let k = choices[dist.sample(&mut rng)];
                if k == 1 && no_changes > bound {
                    m[cc] = loop {
                        let ll = rng.gen_range(0..f);
                        if ll != m[cc] {
                            break ll;
                        }
                    };
                    changed = true;
                    no_changes = 0;
                    break;
                }
            }
            if changed {
                mutated_hosts.push_str(&format!("{:4} {}\n", i, host));
            }
        });
    //
    simulation.pv("mutated_hosts", &mutated_hosts, true);
}
