pub mod v1 {
    use crate::{DeathRule, HostTypes, Simulation};

    pub fn initial(simulation: &mut Simulation) {
        for host_index in 0..simulation.pref().a() + simulation.pref().b() {
            let v = *simulation.ss().host_match_scores_bellow_n().get(&host_index).unwrap();
            if match_condition(v, simulation, host_index, false) {
                simulation.kill_host(host_index);
                simulation.pv("host_dying_initial_exposure", &format!("{:3} {}\n", host_index, &simulation.hosts()[host_index].to_string()), true);
            }
        }
        let (t, r, w) = simulation.count_dead_hosts();
        simulation.pv("host_dying_initial_exposure", &format!(">> {} R({}) W({})\n", t, r, w), true);
    }

    pub fn additional(simulation: &mut Simulation) {
        let already_tried = simulation.ss().hosts_tried().clone();
        simulation.ss_mut().clear_hosts_tried();
        for host_index in already_tried {
            let v = *simulation.ss().host_match_scores_bellow_dd().get(&host_index).unwrap();
            if match_condition(v, simulation, host_index, true) {
                simulation.kill_host(host_index);
                simulation.pv("host_dying_additional_exposure", &format!("{:3} {}\n", host_index, simulation.hosts()[host_index].to_string()), true);
            } else {
                simulation.ss_mut().add_hosts_tried(host_index);
            }
        }
        let (t, r, w) = simulation.count_dead_hosts();
        simulation.pv("host_dying_additional_exposure", &format!(">> {} R({}) W({})\n", t, r, w), true);
    }

    fn match_condition(match_score_bellow_threshold: usize, simulation: &mut Simulation, host_index: usize, is_additional: bool) -> bool {
        match simulation.death_rule() {
            DeathRule::Default => v1(match_score_bellow_threshold, simulation, is_additional, host_index),
            DeathRule::VersionOne => v2(match_score_bellow_threshold, simulation, is_additional, host_index),
            DeathRule::VersionTwo => unimplemented!("Not available in this version")
        }
    }

    // Default kill condition
    pub fn v1(match_score_bellow_threshold: usize, simulation: &mut Simulation, is_additional: bool, host_index: usize) -> bool {
        let test = if is_additional {
            simulation.pref().cc()
        } else {
            simulation.pref().x()
        };
        simulation.pv("kill_condition_test", &format!("{: >3} -> {: >3} >= {: >3}\n", host_index, match_score_bellow_threshold, test), true);
        match_score_bellow_threshold >= test
    }

    /**
    Death to a reservation individual if the total number of unmatched digits is above SS (another
    variable) and death to a wild individual if the total number of unmatched digits is above
    TT (another variable).
     */
    pub fn v2(_: usize, simulation: &mut Simulation, _: bool, host_index: usize) -> bool {
        let g = simulation.pref().g();
        let match_score = simulation.ss().host_match_scores_all().get(&host_index).unwrap().iter().map(|v| g - v).sum::<usize>();
        let test = match simulation.host_type(host_index) {
            HostTypes::Reservation => simulation.pref().ss(),
            HostTypes::Wild => simulation.pref().tt()
        };
        simulation.pv("kill_condition_test", &format!("{: >3} -> {: >3} > {: >3}\n", host_index, match_score, test), true);
        match_score > test
    }
}


pub mod v2 {
    use std::collections::btree_map::BTreeMap;
    use std::collections::HashMap;

    use crate::{calculate_percentile, DeathRule, HostTypes, Simulation};

    pub fn initial(simulation: &mut Simulation) {
        let host_indices = (0..simulation.pref().a() + simulation.pref().b()).collect();
        run_death_rule(simulation, host_indices, false);
    }

    pub fn additional(simulation: &mut Simulation) {
        let host_indices = simulation.ss().hosts_tried().clone();
        simulation.ss_mut().clear_hosts_tried();
        run_death_rule(simulation, host_indices, true);
    }

    pub fn run_death_rule(simulation: &mut Simulation, host_indices: Vec<usize>, is_additional: bool) {
        let file_name = if is_additional {
            "host_dying_additional_exposure"
        } else {
            "host_dying_initial_exposure"
        };
        match simulation.death_rule() {
            DeathRule::Default => v1(simulation, is_additional, host_indices),
            DeathRule::VersionOne => v2(simulation, is_additional, host_indices),
            DeathRule::VersionTwo => v3(simulation, is_additional, host_indices),
        }
        let (t, r, w) = simulation.count_dead_hosts();
        simulation.pv(file_name, &format!(">> {} R({}) W({})\n", t, r, w), true);
    }

    // Default kill condition
    pub fn v1(simulation: &mut Simulation, is_additional: bool, host_indices: Vec<usize>) {
        for host_index in host_indices {
            let (should_kill, v, test, file_name) = if is_additional {
                let v = *simulation.ss().host_match_scores_bellow_dd().get(&host_index).unwrap();
                (v >= simulation.pref().cc(), v, simulation.pref().cc(), "host_dying_additional_exposure")
            } else {
                let v = *simulation.ss().host_match_scores_bellow_n().get(&host_index).unwrap();
                (v >= simulation.pref().x(), v, simulation.pref().x(), "host_dying_initial_exposure")
            };
            simulation.pv("kill_condition_test",
                          &format!("{} {: >3} -> {: >3} > {: >3}\n", file_name, host_index, v, test), true);
            if should_kill {
                simulation.kill_host(host_index);
                simulation.pv(file_name, &format!("{:3} {}\n", host_index, &simulation.hosts()[host_index].to_string()), true);
            }
        }
    }

    pub fn v2(simulation: &mut Simulation, is_additional: bool, host_indices: Vec<usize>) {
        let g = simulation.pref().g();
        let file_name = if is_additional {
            "host_dying_additional_exposure"
        } else {
            "host_dying_initial_exposure"
        };
        for host_index in host_indices {
            let match_score = simulation.ss().host_match_scores_all().get(&host_index).unwrap().iter().map(|v| g - v).sum::<usize>();
            let test = match simulation.host_type(host_index) {
                HostTypes::Reservation => simulation.pref().ss(),
                HostTypes::Wild => simulation.pref().tt()
            };
            simulation.pv("kill_condition_test", &format!("{} {: >3} -> {: >3} > {: >3}\n", file_name, host_index, match_score, test), true);
            if match_score > test {
                simulation.kill_host(host_index);
                simulation.pv(file_name, &format!("{:3} {}\n", host_index, &simulation.hosts()[host_index].to_string()), true);
            }
        }
    }

    /**
    I am proposing a third rule: calculate for each individual's match scores (including additiona
    exposure) G minus the match score and sum them for all the individual's match scores jus
    like in the second rule. But the reservation individual dies if the result is greater tha
    the result of ww% of all individuals and the wild individual dies if the result is greate
    than the result of vv% of individuals.
     */
    pub fn v3(simulation: &mut Simulation, is_additional: bool, host_indices: Vec<usize>) {
        if host_indices.len() == 0 { return; }
        let file_name = if is_additional {
            "host_dying_additional_exposure"
        } else {
            "host_dying_initial_exposure"
        };
        let vv = simulation.pref().vv();
        let ww = simulation.pref().ww();
        let mut _s = String::new();
        let g = simulation.pref().g();
        let d = simulation.pref().d();
        let match_scores: HashMap<usize, Vec<usize>> = simulation.ss().host_match_scores_all().clone()
            .drain_filter(|host_index, _| host_indices.contains(host_index))
            .collect();
        let mut frequency = BTreeMap::<usize, usize>::new();
        let mut cumulative_frequency = vec![0; g * d];
        let mut reservation_individuals_with_score = BTreeMap::<usize, Vec<usize>>::new();
        let mut wild_individuals_with_score = BTreeMap::<usize, Vec<usize>>::new();
        match_scores.iter().for_each(|(host_index, ind_match_scores)| {
            let match_score = ind_match_scores.iter().map(|v| g - v).sum::<usize>();
            *frequency.entry(match_score).or_insert(0) += 1;
            cumulative_frequency[match_score] += 1;
            if simulation.host_type(*host_index) == HostTypes::Wild {
                wild_individuals_with_score.entry(match_score).or_insert(vec![]).push(*host_index);
            } else {
                reservation_individuals_with_score.entry(match_score).or_insert(vec![]).push(*host_index);
            }
        });
        let mut freq = 0;
        let last = frequency.last_key_value().unwrap().0.clone();
        let mut kill_fn = |list: &mut BTreeMap<usize, Vec<usize>>, i| {
            if !list.contains_key(&i) {
                return;
            }
            // println!("{} {:?}", i, list[&i]);
            list.get(&i).unwrap().iter().for_each(|host_index| {
                simulation.kill_host(*host_index);
                simulation.pv(file_name, &format!("{:3} {}\n", host_index, &simulation.hosts()[*host_index].to_string()), true);
            });
        };
        // now calculate the percentile of each match scores
        for i in 0..last + 1 {
            freq += cumulative_frequency.get(i).unwrap_or(&0);
            cumulative_frequency[i] = freq;
            // percentile rank parasite based on highest match scores
            if !frequency.contains_key(&i) { continue; }
            let percentile = calculate_percentile(
                cumulative_frequency[i],
                *frequency.get(&i).unwrap(),
                match_scores.len(),
            );
            // no need to check if the percentile < vv or ww
            if percentile < vv && percentile < ww {
                continue;
            }
            // println!("{} {} {}", percentile, ww, vv);
            if percentile >= vv && !is_additional {
                kill_fn(&mut wild_individuals_with_score, i)
            }
            if percentile >= ww {
                kill_fn(&mut reservation_individuals_with_score, i)
            }
            // println!("---------------");
        }
    }
}

