use std::fmt::{Display, Formatter};
use ndarray::Array3;
use rand::distributions::Uniform;
use rand::distributions::Distribution;
use rand::Rng;
use rand::prelude::SliceRandom;
use crate::{generate_individual, ParasiteSpeciesIndex, Simulation, SimulationPref};
use crate::a2d::A2D;


pub(crate) fn random_parasite(pref: SimulationPref) -> ParasiteSpeciesIndex {
    let mut rng = rand::thread_rng();
    // select random species
    let range = Uniform::new(0, pref.d());
    let species_index = range.sample(&mut rng);
    // select random parasite
    let range = Uniform::new(0, pref.e());
    let parasite_index = range.sample(&mut rng);
    ParasiteSpeciesIndex {
        species_index,
        parasite_index,
        match_count: 0
    }
}