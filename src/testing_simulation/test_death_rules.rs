use crate::testing_simulation::*;

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

const TV3: &str = "death_rule=2|program_version=1|ww=75|vv=49|a=7|b=7|c=40|d=4|e=14|f=2|g=10|h=3|i=1|j=2|k=0.1|l=-1|m=0.5|o=0.04|p=0.05|q=1|r=6|s=6|n=2|x=1|y=0.71|z=1.96|aa=0.5|bb=50|cc=2|dd=1|ee=0.1|ff=1000|gg=115|hh=12|oo=3|jj=0.01|pp=10|ss=9|tt=12";

fn setup_death_rule_v3(additional: bool) -> Simulation {
    let s = TV3.replace("|", "\n");
    let mut v1_s = create_test_simulation(&s);
    v1_s.hosts_mut()[0].set_host_type(HostTypes::Reservation);
    v1_s.hosts_mut()[1].set_host_type(HostTypes::Wild);
    v1_s.hosts_mut()[2].set_host_type(HostTypes::Wild);
    v1_s.hosts_mut()[3].set_host_type(HostTypes::Reservation);
    v1_s.hosts_mut()[4].set_host_type(HostTypes::Reservation);
    v1_s.hosts_mut()[5].set_host_type(HostTypes::Wild);
    v1_s.hosts_mut()[6].set_host_type(HostTypes::Wild);
    v1_s.hosts_mut()[7].set_host_type(HostTypes::Wild);
    v1_s.hosts_mut()[8].set_host_type(HostTypes::Reservation);
    v1_s.hosts_mut()[9].set_host_type(HostTypes::Reservation);
    v1_s.hosts_mut()[10].set_host_type(HostTypes::Reservation);
    v1_s.hosts_mut()[11].set_host_type(HostTypes::Wild);
    v1_s.hosts_mut()[12].set_host_type(HostTypes::Wild);
    v1_s.hosts_mut()[13].set_host_type(HostTypes::Reservation);
    let scores = vec![
        vec![8, 4, 3],    //  0 A R
        vec![2, 3, 2],    //  1 K W
        vec![1, 1, 2],    //  2 K W
        vec![9, 4, 2],    //  3 A R
        vec![2, 9, 2],    //  4 A R
        vec![3, 5, 2],    //  5 A W
        vec![2, 4, 2],    //  6 K W
        vec![2, 2, 4],    //  7 K W
        vec![1, 1, 3],    //  8 K R
        vec![1, 2, 3],    //  9 K R
        vec![3, 1, 3],    // 10 A R
        vec![3, 2, 1],    // 11 A W
        vec![1, 1, 1],    // 12 K W
        vec![1, 1, 1],    // 13 K R
    ];
    set_host_match_scores(&mut v1_s, &scores);
    if additional {
        let scores = vec![
            vec![4],
            vec![5],
        ];
        let mut i = 0;
        for host_index in vec![0, 9] {
            let match_score = scores[i][0];
            v1_s.update_host_match_score_all(host_index, match_score);
            v1_s.update_host_match_score(host_index, 1);
            v1_s.update_host_match_score_bellow_n(host_index, if match_score < v1_s.pref().n() { 1 } else { 0 });
            v1_s.update_host_match_score_bellow_dd(host_index, if match_score < v1_s.pref().dd() { 1 } else { 0 });
            v1_s.update_host_match_score_bellow_j(host_index, if match_score < v1_s.pref().j() { 1 } else { 0 });
            i += 1;
        }
        v1_s.ss_mut().add_hosts_tried(0);
        v1_s.ss_mut().add_hosts_tried(9);
    }
    for (host_index, scores) in v1_s.ss().host_match_scores_all() {
        println!("{:3} {:?}", host_index, scores);
    }
    v1_s
}

#[test]
fn test_match_score_calculation() {
    let mut v1_s = setup_death_rule_v3(false);
    crate::host_death_rules::v2::initial(&mut v1_s);
    print_log(&v1_s);
    assert!(v1_s.is_host_alive(0));
    assert!(!v1_s.is_host_alive(1));
    assert!(!v1_s.is_host_alive(2));
    assert!(v1_s.is_host_alive(3));
    assert!(v1_s.is_host_alive(4));
    assert!(v1_s.is_host_alive(5));
    assert!(v1_s.is_host_alive(6));
    assert!(v1_s.is_host_alive(7));
    assert!(!v1_s.is_host_alive(8));
    assert!(v1_s.is_host_alive(9));
    assert!(v1_s.is_host_alive(10));
    assert!(!v1_s.is_host_alive(11));
    assert!(!v1_s.is_host_alive(12));
    assert!(!v1_s.is_host_alive(13));
}

#[test]
fn test_match_score_calculation_additional() {
    let mut v1_s = setup_death_rule_v3(true);
    crate::host_death_rules::v2::additional(&mut v1_s);
    print_log(&v1_s);
    assert!(v1_s.is_host_alive(0));
    assert!(!v1_s.is_host_alive(9));
}
