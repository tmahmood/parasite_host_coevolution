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
        let (file_name, test) = if is_additional {
            ("host_dying_additional_exposure", simulation.pref().ww() / 100.)
        } else {
            ("host_dying_initial_exposure", simulation.pref().vv() / 100.)
        };
        let mut _s = String::new();
        let g = simulation.pref().g();
        let d = simulation.pref().d();
        let match_scores: HashMap<usize, Vec<usize>> = simulation.ss().host_match_scores_all().clone()
            .drain_filter(|host_index, _| host_indices.contains(host_index))
            .collect();
        let mut frequency = BTreeMap::<usize, usize>::new();
        let mut cumulative_frequency = vec![0; g * d];
        let mut individuals_with_score = BTreeMap::<usize, Vec<usize>>::new();
        match_scores.iter().for_each(|(host_index, ind_match_scores)| {
            let match_score = ind_match_scores.iter().map(|v| g - v).sum::<usize>();
            *frequency.entry(match_score).or_insert(0) += 1;
            cumulative_frequency[match_score] += 1;
            individuals_with_score.entry(match_score).or_insert(vec![]).push(*host_index);
        });
        let mut freq = 0;
        let last = frequency.last_key_value().unwrap().0.clone();
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
            if percentile < test {
                continue;
            }
            individuals_with_score.get(&i).unwrap().iter().for_each(|host_index| {
                simulation.kill_host(*host_index);
                simulation.pv(file_name, &format!("{:3} {}\n", host_index, &simulation.hosts()[*host_index].to_string()), true);
            });
        }
    }
}

#[cfg(test)]
mod testing_simulation {
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
            println!("{:14} -> {}", f, s);
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


    const TV1: &str = "death_rule=0|program_version=1|a=2|b=2|c=18|d=6|e=4|f=2|g=3|h=4|i=2|j=2|k=0.1|l=-1|m=0.5|o=0.04|p=0.05|q=1|r=6|s=6|n=2|x=1|y=0.71|z=1.96|aa=0.5|bb=50|cc=2|dd=1|ee=0.1|ff=1000|gg=115|hh=12|oo=3|jj=0.01|pp=10|ss=9|tt=12";

    #[test]
    fn test_v1_death_rule_v1_initial() {
        let s = TV1.replace("|", "\n");
        let mut s_v1 = create_test_simulation(&s);
        // initial
        // host individual that has X match scores less than N are killed
        // x = 1
        // n = 2
        let scores_initial = vec![
            ////////////////////
            vec![2, 4, 2, 1], // 1
            vec![1, 1, 1, 2], // 3
            vec![2, 4, 2, 3], // 0
            vec![2, 3, 3, 4], // 0
        ];
        set_host_match_scores(&mut s_v1, &scores_initial);
        crate::host_death_rules::v1::initial(&mut s_v1);
        print_simulations_states_s(&s_v1);
        assert_eq!(s_v1.is_host_alive(0), false);
        assert_eq!(s_v1.is_host_alive(1), false);
        assert_eq!(s_v1.is_host_alive(2), true);
        assert_eq!(s_v1.is_host_alive(3), true);
    }

    #[test]
    fn test_v2_death_rule_v1_initial() {
        let s = TV1.replace("|", "\n");
        let mut s_v2 = create_test_simulation(&s);
        let scores_initial = vec![
            ////////////////////
            vec![2, 4, 2, 1], // 1
            vec![1, 1, 1, 2], // 3
            vec![2, 4, 2, 3], // 0
            vec![2, 3, 3, 4], // 0
        ];
        set_host_match_scores(&mut s_v2, &scores_initial);
        crate::host_death_rules::v2::initial(&mut s_v2);
        // crate::host_death_rules::v2::v1(
        //     &mut s_v2, false,
        //     (0..scores_initial.len()).collect(),
        // );
        print_simulations_states_s(&s_v2);
        assert_eq!(s_v2.is_host_alive(0), false);
        assert_eq!(s_v2.is_host_alive(1), false);
        assert_eq!(s_v2.is_host_alive(2), true);
        assert_eq!(s_v2.is_host_alive(3), true);
    }

    #[test]
    fn test_v1_death_rule_v1_additional() {
        // cc = 2
        // dd = 1
        let scores_additional = vec![
            vec![],
            vec![0, 3, 4, 2, 2, 1], // 1
            vec![],
            vec![0, 3, 0, 4, 2, 2], // 2
        ];
        let s = TV1.replace("|", "\n");
        let mut v1_s = create_test_simulation(&s);
        v1_s.hosts_mut()[1].set_host_type(HostTypes::Reservation);
        v1_s.hosts_mut()[3].set_host_type(HostTypes::Reservation);
        set_host_match_scores(&mut v1_s, &scores_additional);
        v1_s.ss_mut().add_hosts_tried(1);
        v1_s.ss_mut().add_hosts_tried(3);
        crate::host_death_rules::v1::additional(&mut v1_s);
        print_simulations_states_s(&v1_s);
        assert_eq!(v1_s.is_host_alive(1), true);
        assert_eq!(v1_s.is_host_alive(3), false);
    }

    #[test]
    fn test_v2_death_rule_v1_additional() {
        // cc = 2
        // dd = 1
        let scores_additional = vec![
            vec![],
            vec![0, 3, 4, 2, 2, 1], // 1
            vec![],
            vec![0, 3, 0, 4, 2, 2], // 2
        ];
        let s = TV1.replace("|", "\n");
        let mut v1_s = create_test_simulation(&s);
        v1_s.hosts_mut()[1].set_host_type(HostTypes::Reservation);
        v1_s.hosts_mut()[3].set_host_type(HostTypes::Reservation);
        set_host_match_scores(&mut v1_s, &scores_additional);
        v1_s.ss_mut().add_hosts_tried(1);
        v1_s.ss_mut().add_hosts_tried(3);
        crate::host_death_rules::v2::additional(&mut v1_s);
        print_simulations_states_s(&v1_s);
        assert_eq!(v1_s.is_host_alive(1), true);
        assert_eq!(v1_s.is_host_alive(3), false);
    }

    const TV2: &str = "death_rule=1|program_version=1|a=2|b=2|c=30|d=4|e=4|f=2|g=10|h=3|i=1|j=2|k=0.1|l=-1|m=0.5|o=0.04|p=0.05|q=1|r=6|s=6|n=2|x=1|y=0.71|z=1.96|aa=0.5|bb=50|cc=2|dd=1|ee=0.1|ff=1000|gg=115|hh=12|oo=3|jj=0.01|pp=10|ss=9|tt=12";

    #[test]
    fn test_v1_death_rule_v2_initial() {
        let s = TV2.replace("|", "\n");
        let mut v1_s = create_test_simulation(&s);
        v1_s.hosts_mut()[0].set_host_type(HostTypes::Wild);
        let scores_initial = vec![
            vec![4, 6, 9],
            vec![0, 0, 0],
            vec![0, 0, 0],
            vec![0, 0, 0],
        ];
        set_host_match_scores(&mut v1_s, &scores_initial);
        crate::host_death_rules::v1::initial(&mut v1_s);
        print_log(&v1_s);
        println!("all:\n{:?}", v1_s.ss().host_match_scores_all());
        assert_eq!(v1_s.is_host_alive(0), true);
    }

    #[test]
    fn test_v1_death_rule_v2_additional() {
        let s = TV2.replace("|", "\n");
        let mut v1_s = create_test_simulation(&s);
        v1_s.hosts_mut()[0].set_host_type(HostTypes::Reservation);
        v1_s.hosts_mut()[1].set_host_type(HostTypes::Wild);
        let scores_additional = vec![
            vec![4, 6, 9, 3],
            vec![],
            vec![],
            vec![],
        ];
        set_host_match_scores(&mut v1_s, &scores_additional);
        v1_s.ss_mut().add_hosts_tried(0);
        crate::host_death_rules::v1::additional(&mut v1_s);
        print_log(&v1_s);
        println!("all:\n{:?}", v1_s.ss().host_match_scores_all());
        assert_eq!(v1_s.is_host_alive(0), false);
    }

    #[test]
    fn test_v2_death_rule_v2_initial() {
        let s = TV2.replace("|", "\n");
        let mut v1_s = create_test_simulation(&s);
        v1_s.hosts_mut()[0].set_host_type(HostTypes::Wild);
        let scores_initial = vec![
            vec![4, 6, 9],
            vec![0, 0, 0],
            vec![0, 0, 0],
            vec![0, 0, 0],
        ];
        set_host_match_scores(&mut v1_s, &scores_initial);
        crate::host_death_rules::v2::initial(&mut v1_s);
        print_log(&v1_s);
        println!("all:\n{:?}", v1_s.ss().host_match_scores_all());
        assert_eq!(v1_s.is_host_alive(0), true);
    }

    #[test]
    fn test_v2_death_rule_v2_additional() {
        let s = TV2.replace("|", "\n");
        let mut v1_s = create_test_simulation(&s);
        v1_s.hosts_mut()[0].set_host_type(HostTypes::Reservation);
        v1_s.hosts_mut()[1].set_host_type(HostTypes::Wild);
        let scores_additional = vec![
            vec![4, 6, 9, 3],
            vec![],
            vec![],
            vec![],
        ];
        set_host_match_scores(&mut v1_s, &scores_additional);
        v1_s.ss_mut().add_hosts_tried(0);
        crate::host_death_rules::v2::additional(&mut v1_s);
        print_log(&v1_s);
        assert_eq!(v1_s.is_host_alive(0), false);
    }

    const TV3: &str = "death_rule=2|program_version=1|ww=85|vv=90|a=6|b=6|c=40|d=4|e=12|f=2|g=10|h=3|i=1|j=2|k=0.1|l=-1|m=0.5|o=0.04|p=0.05|q=1|r=6|s=6|n=2|x=1|y=0.71|z=1.96|aa=0.5|bb=50|cc=2|dd=1|ee=0.1|ff=1000|gg=115|hh=12|oo=3|jj=0.01|pp=10|ss=9|tt=12";

    #[test]
    fn test_match_score_calculation() {
        let s = TV3.replace("|", "\n");
        let mut v1_s = create_test_simulation(&s);
        v1_s.hosts_mut()[0].set_host_type(HostTypes::Reservation);
        v1_s.hosts_mut()[9].set_host_type(HostTypes::Reservation);
        let scores = vec![
            vec![8, 4, 3], // 0
            vec![2, 3, 2], // 1
            vec![1, 3, 1], // 2
            vec![9, 4, 2], // 3
            vec![2, 9, 2], // 4
            vec![3, 5, 2], // 5
            vec![2, 4, 2], // 6
            vec![2, 2, 4], // 7
            vec![2, 3, 2], // 8
            vec![7, 3, 1], // 9
            vec![3, 1, 3], // 10
            vec![3, 2, 1], // 11
            vec![9, 5, 4], // 12
            vec![3, 8, 4], // 13
            vec![1, 1, 1], // 14
        ];
        set_host_match_scores(&mut v1_s, &scores);
        crate::host_death_rules::v2::initial(&mut v1_s);
        print_log(&v1_s);
        assert!(!v1_s.is_host_alive(2));
    }

    #[test]
    fn test_match_score_calculation_additional() {
        let s = TV3.replace("|", "\n");
        let mut v1_s = create_test_simulation(&s);
        v1_s.hosts_mut()[0].set_host_type(HostTypes::Reservation);
        v1_s.hosts_mut()[9].set_host_type(HostTypes::Reservation);
        let scores = vec![
            vec![8, 4, 3, 4], // 0
            vec![2, 3, 2], // 1
            vec![1, 3, 1], // 2
            vec![9, 4, 2], // 3
            vec![2, 9, 2], // 4
            vec![3, 5, 2], // 5
            vec![2, 4, 2], // 6
            vec![2, 2, 4], // 7
            vec![2, 3, 2], // 8
            vec![7, 3, 1, 5], // 9
            vec![3, 1, 3], // 10
            vec![3, 2, 1], // 11
            vec![9, 5, 4], // 12
            vec![3, 8, 4], // 13
            vec![1, 1, 1], // 14
        ];
        set_host_match_scores(&mut v1_s, &scores);
        v1_s.ss_mut().add_hosts_tried(0);
        v1_s.ss_mut().add_hosts_tried(9);
        crate::host_death_rules::v2::additional(&mut v1_s);
        print_log(&v1_s);
        assert!(v1_s.is_host_alive(0));
        assert!(!v1_s.is_host_alive(9));
        assert!(false)
    }
}

