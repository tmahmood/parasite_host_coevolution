# Making the program executable
In Linux, program has to be made executable, So, after uploading to the server, you can make the program executable by running

`chmod +x parasite_host_coevolutiono`

# Running the programs

## V1
`./parasite_host_coevolution params.conf 1`

## V2
`./parasite_host_coevolution params.conf 2`

## V3
`./parasite_host_coevolution params.conf 3`

## V4
`./parasite_host_coevolution params.conf 4`

# Generated files
```
report
├── output.log
├── sim_0
│   ├── 0
│   │   ├── additional_exposure_cond_check.txt
│   │   ├── host_additional_exposure.txt
│   │   ├── host_dying_initial_exposure.txt
│   │   ├── host_exposed_to.txt
│   │   ├── hosts_at_end.txt
│   │   ├── hosts_birth.txt
│   │   ├── mutated_hosts.txt
│   │   ├── mutated_parasites.txt
│   │   ├── parasite_birth.txt
│   │   ├── parasites_at_end.txt
│   │   └── replaced_parasites.txt
│   ├── 1
│   ├── 2
│   ├── hosts
│   └── parasites
├── sim_1
└── sim_2

12 directories, 109 files
```

### additional_exposure_cond_check.txt
if additional exposure conditions were met in this run
### host_additional_exposure.txt
hosts that are selected for additional exposure
### host_dying_initial_exposure.txt
Which hosts are killed at initial exposure
### host_exposed_to.txt
Which hosts are exposed initially
### hosts_at_end.txt
Hosts survived at this run
### hosts_birth.txt
Chances and p values, hosts born in this run
### mutated_hosts.txt
hosts that are mutated, total mutation, and mutation rate
### mutated_parasites.txt
how many mutation occurred, and probablity
### parasite_birth.txt
Which individuals are killed and which individual used to birth new parasite individual
### parasites_at_end.txt
Parasites left after this run
### replaced_parasites.txt
Which parasites are replaced
### output.log
contains statistics at the end of each run