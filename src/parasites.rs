use std::fmt::Formatter;

use rand::distributions::Distribution;
use rand::distributions::Uniform;

use crate::{ParasiteSpeciesIndex, Simulation, SimulationPref};

// pub(crate) fn random_parasite(pref: SimulationPref) -> ParasiteSpeciesIndex {
//     let mut rng = rand::thread_rng();
//     // select random species
//     let range = Uniform::new(0, pref.d());
//     let species_index = range.sample(&mut rng);
//     // select random parasite
//     let range = Uniform::new(0, pref.e());
//     let parasite_index = range.sample(&mut rng);
//     ParasiteSpeciesIndex {
//         species_index,
//         parasite_index,
//         match_count: 0
//     }
// }