use crate::serde_derive::Deserialize;
use serde_ini::{Deserializer, Parser};

#[derive(Deserialize, Debug, Clone)]
pub struct SimulationPref {
    // number of reservation hosts to start each simulation run
    a: i32,
    // number of wild hosts to start each simulation run
    b: i32,
    // length of the set of numbers associated with each host
    c: i32,
    // number of species of parasites
    d: i32,
    // number of individuals in each parasite species
    e: i32,
    // possible values of each number in the set of numbers associated with both parasites and hosts
    f: i32,
    // length of the set of numbers associated with each parasite
    g: i32,
    // number of parasite species a host is exposed to in each generation
    h: i32,
    // number of parasite species of additional exposure for some reservation individuals
    i: i32,
    // Match threshold
    j: i32,
    // mutation rate for parasites
    k: f32,
    // no additional exposure until after the Lth generation
    l: i32,
    // no additional exposure unless less than an M fraction of host (reservation and wild) individuals have been killed
    m: f32,
    // host killed if it has a match score with at least X parasite individuals that’s lower than N
    n: i32,
    // relatedness factor
    o: f32,
    // reservation factor
    p: f32,
    // the number of parasite species that are replaced each generation
    q: i32,
    // reserve_constant
    r: i32,
    // wild_constant
    s: i32,
    // If a host individual has a match score with at least X parasite individuals that is lower than N, then that host individual is considered killed.
    x: i32,
    // inefficiency factor
    y: f32,
    // T variables
    z: f32,
    // percentage of reservation hosts exposed to additional parasites
    aa: f32,
    // 0.5
    // all parasite individuals with a total match score that is higher than BB% of all parasite individuals in the same parasite species become eliminated
    bb: i32,
    // on the conditions that it is after the Lth generation and if less than an M fraction of host individuals (a total of reservation and wild host individuals) have been killed if a reservation host individual has a match score (including this additional exposure) with at least CC parasite individuals that is lower than DD, then that host individual is considered killed
    cc: i32,
    // on the conditions that it is after the Lth generation and if less than an M fraction of host individuals (a total of reservation and wild host individuals) have been killed if a reservation host individual has a match score (including this additional exposure) with at least CC parasite individuals that is lower than DD, then that host individual is considered killed
    dd: i32,
    // mutation rate for hosts
    ee: f32,
    // number of generations per simulation run
    ff: i32,
    // number of simulation runs
    gg: i32,
    // HH-number of match scores below the threshold
    hh: i32,
}

impl SimulationPref {
    pub fn a(&self) -> i32 {
        self.a
    }
    pub fn b(&self) -> i32 {
        self.b
    }
    pub fn c(&self) -> i32 {
        self.c
    }
    pub fn d(&self) -> i32 {
        self.d
    }
    pub fn e(&self) -> i32 {
        self.e
    }
    pub fn f(&self) -> i32 {
        self.f
    }
    pub fn g(&self) -> i32 {
        self.g
    }
    pub fn h(&self) -> i32 {
        self.h
    }
    pub fn i(&self) -> i32 {
        self.i
    }
    pub fn j(&self) -> i32 {
        self.j
    }
    pub fn k(&self) -> f32 {
        self.k
    }
    pub fn l(&self) -> i32 {
        self.l
    }
    pub fn m(&self) -> f32 {
        self.m
    }
    pub fn n(&self) -> i32 {
        self.n
    }
    pub fn o(&self) -> f32 {
        self.o
    }
    pub fn p(&self) -> f32 {
        self.p
    }
    pub fn q(&self) -> i32 {
        self.q
    }
    pub fn r(&self) -> i32 {
        self.r
    }
    pub fn s(&self) -> i32 {
        self.s
    }
    pub fn x(&self) -> i32 {
        self.x
    }
    pub fn y(&self) -> f32 {
        self.y
    }
    pub fn z(&self) -> f32 {
        self.z
    }
    pub fn aa(&self) -> f32 {
        self.aa
    }
    pub fn bb(&self) -> i32 {
        self.bb
    }
    pub fn cc(&self) -> i32 {
        self.cc
    }
    pub fn dd(&self) -> i32 {
        self.dd
    }
    pub fn ee(&self) -> f32 {
        self.ee
    }
    pub fn ff(&self) -> i32 {
        self.ff
    }
    pub fn gg(&self) -> i32 {
        self.gg
    }
    pub fn hh(&self) -> i32 {
        self.hh
    }
}
