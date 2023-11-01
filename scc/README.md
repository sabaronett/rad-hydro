# SCC Workflow

## [Computing Resources Overview](https://wiki.flatironinstitute.org/SCC/Overview)

## Account

- [Status and resource usage](https://fido.flatironinstitute.org/home/self)
- [Usage over time (on FI network)](https://grafana.flatironinstitute.org/d/GM2HFVR7k/user-resource-usage?var-user=sbaronett)
- [Planned resource usage](https://wiki.flatironinstitute.org/SCC/AccountManagement/PlannedResourceUsage)


## Interactive Slurm Jobs
An example of requesting an [interactive node with srun](https://wiki.flatironinstitute.org/SCC/Software/Slurm#srun_Run_a_program_on_allocated_resources):
```bash
srun -N1 -p genx -C rome --pty bash -i
```


## [Data Transfer](https://wiki.flatironinstitute.org/SCC/Hardware/DataTransfer)

### For small amounts of data (< 25GB)

#### [scp](https://wiki.flatironinstitute.org/SCC/Hardware/DataTransfer#scp)

```bash
# Tranfer from a remote host
$ scp username@rusty.flatironinstitute.org:~/experiments/mycode.tar.gz .
```

#### [rsync](https://wiki.flatironinstitute.org/SCC/Hardware/DataTransfer#rsync)

```bash
# Syncing data from a remote host
rsync -a username@rusty.flatironinstitute.org:/home/username/code_folder ~/my_local_code
```


## [Athena++](https://github.com/PrincetonUniversity/athena/wiki)

### [Configuring](https://github.com/PrincetonUniversity/athena/wiki/Configuring)

#### AMD Rome Nodes (Rusty)

##### [GNU Compiler](https://www.nas.nasa.gov/hecc/support/kb/preparing-to-run-on-aitken-rome-nodes_657.html#:~:text=on%20Rome%20processors.-,GNU%20Compilers,-%3A)
In Athena++'s root,
```bash
./configure.py --prob=[PROBLEM] -implicit_radiation -mpi -hdf5 -h5double --cxx=gcc --cflag="-march=znver2"
```


##### AMD Optimizing C/C++ and Fortran Compilers ([AOCC](https://www.nas.nasa.gov/hecc/support/kb/preparing-to-run-on-aitken-rome-nodes_657.html))

1. `module load aocc`
2. In Athena++'s root,
```bash
./configure.py --prob=[PROBLEM] -implicit_radiation -mpi -hdf5 -h5double --cxx=clang++ --cflag="-march=znver2"
```


### [Compiling](https://github.com/PrincetonUniversity/athena/wiki/Compiling)

Parallel compilation (e.g., on a single Rome node):
1. Request an [interactive node with srun](https://wiki.flatironinstitute.org/SCC/Software/Slurm#srun_Run_a_program_on_allocated_resources):
   ```bash
   srun -N1 -p cca -C rome --exclusive -t 0:15:00 --pty bash -i
   ```
2. In Athena++'s root,
   ```bash
   make clean
   ```
   then
   ```bash
   make -j && cp bin/athena bin/athena.[problem_id]
   ```
Compilation should take less than a minute.


### [Running the Code](https://github.com/PrincetonUniversity/athena/wiki/Running-the-Code)

#### [sbatch](https://wiki.flatironinstitute.org/SCC/Software/Slurm#sbatch_Allocating_Resources)

See [`sample`](/scc/sample).
```bash
sbatch <script>
```

#### [scancel](https://slurm.schedmd.com/scancel.html)

```bash
scancel <job_id>
```

#### [Continuous Restarts](https://slurm.schedmd.com/sbatch.html#OPT_dependency)
```bash
sbatch <script> -d afterok:<job_id>
```
