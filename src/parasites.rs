use std::fmt::{Display, Formatter};
use rand::distributions::Uniform;
use rand::Rng;
use rand::prelude::SliceRandom;
use crate::{generate_individual, Simulation, SimulationPref};

#[derive(Debug, Clone)]
pub struct Parasite {
    number_set: Vec<i32>,
    f: i32,
    c: i32,
}

impl Parasite {
    pub fn new(f: i32, c: i32) -> Self {
        Parasite {
            number_set: generate_individual(f, c),
            f,
            c,
        }
    }

    pub fn number_set(&self) -> &Vec<i32> {
        &self.number_set
    }

    pub fn f(&self) -> i32 {
        self.f
    }

    pub fn c(&self) -> i32 {
        self.c
    }

    pub fn set_number_set(&mut self, number_set: Vec<i32>) {
        self.number_set = number_set;
    }

    pub fn set_f(&mut self, f: i32) {
        self.f = f;
    }

    pub fn set_c(&mut self, c: i32) {
        self.c = c;
    }
}

impl Display for Parasite {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let mut s = String::new();
        for number in &self.number_set {
            s += &format!("{}", number);
        }
        write!(f, "{}", s)
    }
}

pub async fn create_random_parasites(pref: &SimulationPref) -> Vec<Vec<Parasite>> {
    let mut all_parasites = Vec::with_capacity(
        pref.d() as usize
    );
    for i in 0..pref.d() as usize {
        all_parasites.insert(i, Vec::with_capacity(pref.e() as usize));
        for _ in 0..pref.e() as usize {
            all_parasites[i].push(Parasite::new(pref.f(), pref.g()))
        }
    }
    all_parasites
}

pub async fn chose_parasites(simulation: &Simulation) -> Vec<(i32, Parasite)> {
    let mut rng = rand::thread_rng();
    let d = simulation.pref().d();
    let h = simulation.pref().h();
    let range = Uniform::new(0, d);
    let species: Vec<i32> = (0..h).map(|_| rng.sample(&range)).collect();
    let mut parasites = vec![];
    for species_index in species {
        let parasite = simulation.parasites()[species_index as usize].choose(&mut rng);
        parasites.push((species_index, parasite.unwrap().clone()));
    }
    parasites
}

