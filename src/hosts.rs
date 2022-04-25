use std::collections::HashMap;
use std::fmt::{Display, format, Formatter};
use std::ptr::write;

use ndarray::{Array, Array1, Ix, Ix1, Ix2};
use ndarray_rand::RandomExt;
use rand::{Rng, thread_rng};
use rand::distributions::{Distribution, Standard, Uniform};
use rand::prelude::SliceRandom;

use crate::{generate_individual, SimulationPref};

#[derive(Clone, Debug, PartialOrd, PartialEq, Copy)]
pub enum HostTypes {
    Reservation,
    Wild,
}

impl Display for HostTypes {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let r = match self {
            HostTypes::Reservation => "R",
            HostTypes::Wild => "W"
        };
        write!(f, "{}", r)
    }
}

impl From<usize> for HostTypes {
    fn from(num: usize) -> Self {
        if num == 0 { HostTypes::Reservation } else { HostTypes::Wild }
    }
}

#[derive(Debug, Clone)]
pub struct Host {
    host_type: HostTypes,
    number_set: Array1<usize>,
    alive: bool,
}

impl Default for Host {
    fn default() -> Self {
        Host {
            host_type: HostTypes::Reservation,
            number_set: Default::default(),
            alive: true
        }
    }
}

impl Display for Host {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let mut s = String::new();
        for number in &self.number_set {
            s += &format!("{}", number);
        }
        write!(f, "({}) {}:{}", self.alive, self.host_type, s)
    }
}

impl Host {
    pub fn new_reservation_host(number_set: Array1<usize>) -> Self {
        Host::new(HostTypes::Reservation, number_set)
    }

    pub fn new_wild_host(number_set: Array1<usize>) -> Self {
        Host::new(HostTypes::Wild, number_set)
    }

    fn new(host_type: HostTypes, number_set: Array1<usize>) -> Self {
        Host { host_type, number_set, alive: true }
    }

    pub fn host_type(&self) -> HostTypes {
        self.host_type.clone()
    }
    pub fn number_set(&self) -> &Array1<usize> {
        &self.number_set
    }

    pub fn alive(&self) -> bool {
        self.alive
    }

    pub fn set_host_type(&mut self, host_type: HostTypes) -> &mut Self{
        self.host_type = host_type;
        self

    }

    pub fn set_number_set(&mut self, number_set: Array1<usize>) -> &mut Self {
        self.number_set = number_set;
        self
    }

    pub fn set_alive(&mut self, alive: bool) -> &mut Self {
        self.alive = alive;
        self
    }
}

pub fn create_random_hosts(pref: &SimulationPref) -> Array1<Host> {
    let total = (pref.a() + pref.b()) as usize;
    let mut all_hosts = Vec::with_capacity(total);
    (0..pref.a()).for_each(|_| {
        all_hosts.push(Host::new_reservation_host(
            generate_individual(pref.f(), pref.c()))
        );
    });
    (0..pref.b()).for_each(|_|
        all_hosts.push(Host::new_wild_host(
            generate_individual(pref.f(), pref.c()))
        )
    );
    let mut rng = rand::thread_rng();
    all_hosts.shuffle(&mut rng);
    Array::from(all_hosts)
}

