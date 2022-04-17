use std::fmt::{Display, format, Formatter};
use std::ptr::write;
use rand::{Rng, thread_rng};
use rand::prelude::SliceRandom;
use crate::{generate_individual, SimulationPref};

#[derive(Clone, Debug)]
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

#[derive(Debug)]
pub struct Host {
    host_type: HostTypes,
    number_set: Vec<i32>,
    f: i32,
    c: i32,
}

impl Display for Host {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let mut s = String::new();
        for number in &self.number_set {
            s += &format!("{}", number);
        }
        write!(f, "{}:{}", self.host_type , s)
    }
}

impl Host {
    pub fn new_random(f: i32, c: i32) -> Self {
        let mut rng = thread_rng();
        let n: usize = rng.gen_range(0..1);
        Host::new(HostTypes::from(n), f, c)
    }

    pub fn new_reservation_host(f: i32, c: i32) -> Self {
        Host::new(HostTypes::Reservation, f, c)
    }

    pub fn new_wild_host(f: i32, c: i32) -> Self {
        Host::new(HostTypes::Wild, f, c)
    }

    fn new(host_type: HostTypes, f: i32, c: i32) -> Self {
        Host {
            host_type,
            number_set: generate_individual(f, c),
            f, c,
        }
    }

    pub fn host_type(&self) -> HostTypes {
        self.host_type.clone()
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

    pub fn set_host_type(&mut self, host_type: HostTypes) {
        self.host_type = host_type;
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

pub async fn create_random_hosts(pref: &SimulationPref) -> Vec<Host> {
    let total = (pref.a() + pref.b()) as usize;
    let mut all_hosts = Vec::with_capacity(total);
    (0..pref.a()).for_each(|_|
        all_hosts.push(Host::new_reservation_host(pref.f(), pref.c()))
    );
    (0..pref.b()).for_each(|_|
        all_hosts.push(Host::new_wild_host(pref.f(), pref.c()))
    );
    let mut rng = rand::thread_rng();
    all_hosts.shuffle(&mut rng);
    all_hosts
}

