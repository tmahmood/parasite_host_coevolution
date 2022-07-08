use crate::{HostTypes, ProgramVersions, Simulation};

pub struct QiParams<'a> {
    pub hh: usize,
    pub score: usize,
    pub simulation: &'a Simulation,
    pub host_index: usize,
    pub xx: f32,
}

pub fn calculate_qi(simulation: &mut Simulation) {
    let mut qr = 0.;
    let mut qw = 0.;
    let match_score_bellow_j = simulation.ss().host_match_scores_bellow_j().clone();
    let hosts = simulation.hosts().clone();

    let callback = match simulation.program_version() {
        ProgramVersions::Two | ProgramVersions::Four => calculate_qi_v1,
        ProgramVersions::Five | ProgramVersions::Six => calculate_qi_v2,
        ProgramVersions::Seven | ProgramVersions::Eight => calculate_qi_v3,
        ProgramVersions::Nine | ProgramVersions::Ten => calculate_qi_v4,
        _ => panic!("Should not come here!")
    };

    match_score_bellow_j.iter()
        .filter(|(index, _)| hosts[**index].alive())
        .for_each(|(index, score)| {
            let hh = simulation.pref().hh();
            let qi_params = QiParams {
                hh,
                score: *score,
                simulation,
                host_index: *index,
                xx: simulation.pref().xx(),
            };
            let mut qi = callback(qi_params);
            if qi < 0. { qi = 0. }
            simulation.set_host_individual_qi(*index, qi);
            match simulation.host_type(*index) {
                HostTypes::Reservation => qr += qi,
                HostTypes::Wild => qw += qi
            }
        });
    simulation.ss_mut().set_qr(qr);
    simulation.ss_mut().set_qw(qw);
}

pub fn calculate_qi_v1(qi_params: QiParams) -> f32 {
    qi_params.hh as f32 - qi_params.score as f32
}

// Suppose: An individual's match scores are 4, 6, and 9. And OO=7. HH=12. JJ=0.1. PP=10. And
// suppose there are C=36=length of set of numbers associated with a host.
//
// First way -> NN=(OO=7-4)+(OO=7-6)=4. Then Qi=(HH=12)*(1-JJ=.1*NN=4)=7.2.
// Second way -> 4 < 7, and 6 < 7, but 9 > 7. NN=2. Then Qi=(HH=12)*(1-JJ=.1*NN=2)=9.6.
// Third way -> NN=(C=36)-(4+6+9)-(PP=10)=7. Then Qi=(HH=12)*(1-JJ=.1*NN=7)=3.6.
//
// And the new death rule=Suppose SS=12 and TT=19
//
// if the individual is reservation, (C=36)- (4+6+9)=17 > SS=12. The individual dies.
// if the individual is wild, (C=36)-(4+6+9)=17 < TT=19. The individual lives.


/**
Same as version 2 except Qi ->Q_subscript_i=Qi for each surviving host individual is the
 variable HH times (1-JJ*NN), where JJ and NN are new variables different from J and N and
 JJ is a given input and NN is the total number of digits that each match score is under OO,
 another new variable.
 **/
pub fn calculate_qi_v2(qi_params: QiParams) -> f32 {
    let match_scores = qi_params.simulation.ss().host_match_scores_all().get(&qi_params.host_index).unwrap();
    let oo = qi_params.simulation.pref().oo();
    let jj = qi_params.simulation.pref().jj();
    let g = qi_params.simulation.pref().g();
    let mut nn: f32 = match_scores.iter()
        .filter(|v| g - **v > oo)
        .fold(0., |mut a, v| {
            a += g as f32 - *v as f32 - oo as f32;
            a
        });
    if nn < 0. { nn = 0. }
    nn = nn.powf(qi_params.xx);
    let jj_nn = (1. - (jj as f32 * nn));
    qi_params.hh as f32 * jj_nn
}

// 3) same as version 2 except Qi ->Q_subscript_i=Qi for each surviving host individual is the
// variable HH times (1-JJ*NN), where JJ is a given input and NN is the total number of match
// scores under OO.
pub fn calculate_qi_v3(qi_params: QiParams) -> f32 {
    let match_scores = qi_params.simulation.ss().host_match_scores_all().get(&qi_params.host_index).unwrap();
    let nn = (match_scores.iter()
        .filter(|v| **v < qi_params.simulation.pref().oo())
        .count() as f32).powf(qi_params.xx);
    qi_params.hh as f32 * (1. - (qi_params.simulation.pref().jj() as f32 * nn as f32))
}

// A version with: Qi for each surviving host individual is the variable HH times (1-JJ*NN),
// where JJ is a given input and NN is the number of unmatched digits above PP (another variable).
pub fn calculate_qi_v4(qi_params: QiParams) -> f32 {
    let match_scores = qi_params.simulation.ss().host_match_scores_all().get(&qi_params.host_index).unwrap();
    let pp = qi_params.simulation.pref().pp();
    let jj = qi_params.simulation.pref().jj();
    let g = qi_params.simulation.pref().g();
    let mut nn = match_scores.iter()
        .fold(0., |mut a, v| {
            a += g as f32 - *v as f32;
            a
        }) - pp as f32;
    if nn < 0. { nn = 0. }
    nn = nn.powf(qi_params.xx);
    qi_params.hh as f32 * (1. - (jj as f32 * nn))
}
