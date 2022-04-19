use std::fmt::{Display, Formatter};
use rand::distributions::Uniform;
use rand::distributions::Distribution;
use rand::Rng;
use rand::prelude::SliceRandom;
use crate::{generate_individual, ParasiteSpeciesIndex, Simulation, SimulationPref};
use crate::a2d::A2D;


pub async fn chose_parasites(simulation: &Simulation) -> Vec<ParasiteSpeciesIndex> {
    let d = simulation.pref().d();
    let h = simulation.pref().h();
    // select random species
    let mut rng = rand::thread_rng();
    let range = Uniform::new(0, d);
    let species: Vec<usize> = (0..h).map(|_| rng.sample(&range)).collect();
    // select random parasite individuals
    let mut selected_parasites = vec![];
    for species_index in species {
        let mut range = Uniform::new(0, simulation.pref().e());
        let parasite_index = range.sample(&mut rng);
        selected_parasites.push(ParasiteSpeciesIndex {
            species_index,
            parasite_index,
            match_count: 0
        });
    }
    selected_parasites
}

