# SCC Workflow

## [Computing Resources Overview](https://wiki.flatironinstitute.org/SCC/Overview)

## Account

- [Status and resource usage](https://fido.flatironinstitute.org/home/self)
- [Usage over time (on FI network)](https://grafana.flatironinstitute.org/d/GM2HFVR7k/user-resource-usage?var-user=sbaronett)
- [Planned resource usage](https://wiki.flatironinstitute.org/SCC/AccountManagement/PlannedResourceUsage)


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
   srun -N1 -p genx -C rome -t 0:10:00 --pty bash -i
   ```
2. In Athena++'s root,
   ```bash
   make clean
   make -j
   ```
Compilation should take only a few minutes.


### [Running the Code](https://github.com/PrincetonUniversity/athena/wiki/Running-the-Code)

#### [sbatch](https://wiki.flatironinstitute.org/SCC/Software/Slurm#sbatch_Allocating_Resources)

See [`sample`](/scc/sample).
```bash
sbatch <script>
```

#### [Continuous Restarts](https://slurm.schedmd.com/sbatch.html#OPT_dependency)
```bash
sbatch <script> -d afterok:<job_id>
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





# OLD

The following compilation command and flags at the end [generate code for any processor type](https://www.nas.nasa.gov/hecc/support/kb/recommended-compiler-options_99.html#:~:text=Generate%20Code%20for%20Any%20Processor%20Type):

1. ```bash
   pfeXX:~> ./configure.py --prob=streaming_instability -p --eos=isothermal --nghost=3 -mpi -hdf5 -h5double --cxx=icpc --mpiccmd="icpc -lmpi -lmpi++" --cflag="-O3 -axCORE-AVX512,CORE-AVX2 -xAVX"
   ```
2. Manually remove `-xhost` from `athena/Makefile`, under
   ```Makefile
   # General compiler specifications
   ...
   CXXFLAGS := ... -xhost ...
   ```


## [Archiving to Lou](https://www.nas.nasa.gov/hecc/support/kb/using-shift-for-transfers-and-tar-operations-between-two-nas-hosts_513.html)

- Check quota status:
```bash
lfs quota -h -u $USERNAME /nobackupp12
```
- Snapshotting:
```bash
lfeX:~$ mkdir snapshots/$(date +"%Y-%m-%d")
lfeX:~$ cd /nobackup/$USERNAME
lfeX:/nobackup/$USERNAME$ shiftc --hosts=8 --create-tar --index-tar github lfe:~/snapshots/$(date +"%Y-%m-%d")/github.tar
```
- Archiving:
```bash
lfeX:~$ mkdir archives/.../dir
lfeX:~$ cd /nobackup/$USERNAME
lfeX:/nobackup/$USERNAME$ shiftc --hosts=8 --create-tar --index-tar github/.../dir lfe:~/archives/.../dir/dir.$(date +"%Y-%m-%d").tar
```


### Interactive Jobs

Using [VNC Xterm](https://www.nas.nasa.gov/hecc/support/kb/vnc-a-faster-alternative-to-x11_257.html): 
1. Start VNC server on PFE:
```bash
pfeXX:~> vncserver -localhost
```
2. Copy-paste this line into PFE prompt: `~C`
3. In FE prompt: `-L 2222X:localhost:592X` (or any available local port)
4. Locally (same FE port above): 
```bash
$ vncviewer localhost:2222X
``` 
5. Submit job:
```bash
pfeXX:~> qsub -I -X -lselect=1:ncpus=28:mpiprocs=28:model=bro_ele,walltime=1:00:00 -q devel
```
6. When finished, before logging off, make sure to kill the VNC server:
```bash
pfeXX:~> vncserver -kill :XX
```
where `XX` is the ID of the server initiated earlier.
