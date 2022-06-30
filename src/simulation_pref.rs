use crate::serde_derive::Deserialize;

#[derive(Deserialize, Debug, Clone, Copy)]
pub struct SimulationPref {
    /** number of reservation hosts to start each simulation run */
    a: usize,
    /** number of wild hosts to start each simulation run */
    b: usize,
    /** length of the set of numbers associated with each host */
    c: usize,
    /** number of species of parasites */
    d: usize,
    /** number of individuals in each parasite species */
    e: usize,
    /** possible values of each number in the set of numbers associated with both parasites and hosts */
    f: usize,
    /** length of the set of numbers associated with each parasite */
    g: usize,
    /** number of parasite species a host is exposed to in each generation */
    h: usize,
    /** number of parasite species of additional exposure for some reservation individuals */
    i: usize,
    /** Match threshold */
    j: usize,
    /** mutation rate for parasites */
    k: f32,
    /** no additional exposure until after the Lth generation */
    l: i32,
    /** no additional exposure unless less than an M fraction of host (reservation and wild) individuals have been killed */
    m: f32,
    /** host killed if it has a match score with at least X parasite individuals that’s lower than N */
    n: usize,
    /** relatedness factor */
    o: f32,
    /** reservation factor */
    p: f32,
    /** the number of parasite species that are replaced each generation */
    q: usize,
    /** reserve_constant */
    r: usize,
    /** wild_constant */
    s: usize,
    /** If a host individual has a match score with at least X parasite individuals that is lower than N, then that host individual is considered killed. */
    x: usize,
    /** inefficiency factor */
    y: f32,
    /** T variables */
    z: f32,
    /** percentage of reservation hosts exposed to additional parasites */
    aa: f32,
    /** all parasite individuals with a total match score that is higher than BB% of all parasite individuals in the same parasite species become eliminated */
    bb: f32,
    /** on the conditions that it is after the Lth generation and if less than an M fraction of host individuals (a total of reservation and wild host individuals) have been killed if a reservation host individual has a match score (including this additional exposure) with at least CC parasite individuals that is lower than DD, then that host individual is considered killed */
    cc: usize,
    /** on the conditions that it is after the Lth generation and if less than an M fraction of host individuals (a total of reservation and wild host individuals) have been killed if a reservation host individual has a match score (including this additional exposure) with at least CC parasite individuals that is lower than DD, then that host individual is considered killed */
    dd: usize,
    /** mutation rate for hosts */
    ee: f32,
    /** number of generations per simulation run */
    ff: usize,
    /** number of simulation runs */
    gg: usize,
    /** HH-number of match scores below the threshold */
    hh: usize,

    #[serde(default = "default_usize")]
    oo: usize,

    #[serde(default = "default_f32")]
    jj: f32,

    #[serde(default = "default_usize")]
    pp: usize,

    #[serde(default = "default_usize")]
    ss: usize,

    #[serde(default = "default_usize")]
    tt: usize,

    // program configs, will be overridden if command line provided
    #[serde(default = "default_usize")]
    program_version: usize,

    #[serde(default = "default_usize")]
    death_rule: usize,

    #[serde(default = "default_f32")]
    vv: f32,

    #[serde(default = "default_f32")]
    ww: f32,
}

impl SimulationPref {

    pub(crate) fn ss(&self) -> usize {
        self.ss
    }

    pub(crate) fn tt(&self) -> usize {
        self.tt
    }
    /** number of reservation hosts to start each simulation run */
    pub fn a(&self) -> usize {
        self.a
    }

    /** number of wild hosts to start each simulation run */
    pub fn b(&self) -> usize {
        self.b
    }

    /** length of the set of numbers associated with each host */
    pub fn c(&self) -> usize {
        self.c
    }

    /** number of species of parasites */
    pub fn d(&self) -> usize {
        self.d
    }

    /** number of individuals in each parasite species */
    pub fn e(&self) -> usize {
        self.e
    }
    /** possible values of each number in the set of numbers associated with both parasites and hosts */
    pub fn f(&self) -> usize {
        self.f
    }

    /** length of the set of numbers associated with each parasite */
    pub fn g(&self) -> usize {
        self.g
    }

    /** number of parasite species a host is exposed to in each generation */
    pub fn h(&self) -> usize {
        self.h
    }

    /** number of parasite species of additional exposure for some reservation individuals */
    pub fn i(&self) -> usize {
        self.i
    }

    /** Match threshold */
    pub fn j(&self) -> usize {
        self.j
    }

    /** mutation rate for parasites */
    pub fn k(&self) -> f32 {
        self.k
    }

    /**
    no additional exposure until after the Lth generation
     */
    pub fn l(&self) -> i32 {
        self.l
    }

    /** no additional exposure unless less than an M fraction of host (reservation and wild) individuals have been killed */
    pub fn m(&self) -> f32 {
        self.m
    }

    /** host killed if it has a match score with at least X parasite individuals that’s lower than N */
    pub fn n(&self) -> usize {
        self.n
    }

    /** relatedness factor */
    pub fn o(&self) -> f32 {
        self.o
    }

    /** reservation factor */
    pub fn p(&self) -> f32 {
        self.p
    }

    /** the number of parasite species that are replaced each generation */
    pub fn q(&self) -> usize {
        self.q
    }

    /** reserve_constant */
    pub fn r(&self) -> usize {
        self.r
    }

    /** wild_constant */
    pub fn s(&self) -> usize {
        self.s
    }

    /** If a host individual has a match score with at least X parasite individuals that is lower than N, then that host individual is considered killed. */
    pub fn x(&self) -> usize {
        self.x
    }

    /** inefficiency factor */
    pub fn y(&self) -> f32 {
        self.y
    }

    /** T variables */
    pub fn z(&self) -> f32 {
        self.z
    }

    /** percentage of reservation hosts exposed to additional parasites */
    pub fn aa(&self) -> f32 {
        self.aa
    }

    /** all parasite individuals with a total match score that is higher than BB% of all parasite individuals in the same parasite species become eliminated */
    pub fn bb(&self) -> f32 {
        self.bb
    }

    /** On the conditions that it is after the *L*th generation and if less than an *M* fraction of host
             individuals (a total of reservation and wild host individuals) have been killed if a reservation
             host individual has a match score (including this additional exposure) with at least *CC* parasite
             individuals that is lower than *DD*, then that host individual is considered killed
     */
    pub fn cc(&self) -> usize {
        self.cc
    }

    /** on the conditions that it is after the Lth generation and if less than an M fraction of host individuals
          (a total of reservation and wild host individuals) have been killed if a reservation host
          individual has a match score (including this additional exposure) with at least CC
          parasite individuals that is lower than DD, then that host individual is considered killed */
    pub fn dd(&self) -> usize {
        self.dd
    }

    /** mutation rate for hosts */
    pub fn ee(&self) -> f32 {
        self.ee
    }

    /** number of generations per simulation run */
    pub fn ff(&self) -> usize {
        self.ff
    }

    /** number of simulation runs */
    pub fn gg(&self) -> usize {
        self.gg
    }

    /** HH-number of match scores below the threshold */
    pub fn hh(&self) -> usize {
        self.hh
    }

    pub fn oo(&self) -> usize {
        self.oo
    }

    pub fn jj(&self) -> f32 {
        self.jj
    }

    pub fn pp(&self) -> usize {
        self.pp
    }

    pub fn death_rule(&self) -> usize {
        self.death_rule
    }

    pub fn program_version(&self) -> usize {
        self.program_version
    }

    pub fn vv(&self) -> f32 {
        self.vv
    }
    pub fn ww(&self) -> f32 {
        self.ww
    }
}

fn default_usize() -> usize {
    0
}

fn default_f32() -> f32 {
    0.
}
