# `Loki v1.4`
## C++ code to analize polarization properties of radiowaves propagating in the pulsar magnetosphere.

Theoretical basis for the program is given in the following papers <sup>[1](#myfootnote1)</sup>
<sup>[2](#myfootnote2)</sup> <sup>[3](#myfootnote3)</sup>.

[![Loki](Loki_lo.jpg)](#loki-v14)

### Libraries description

##### `lib/`:
* `constants.h`: stores all the necessary constants within the `constants::` namespace;
* `initialize.h`: some of the initialization routines used, such as finding the radiation point, reading input etc;
* `process_functions.h`: all the physical functions involved;
* `b_field.h`: functions to find the magnetic field for a given point;
* `RHS.h`: the right-hand side of the ODE (both in `float` and `double`);
* `pc_dens.h`: back integration routine to find the density on a polar cap for any given point on a ray trajectory.

##### `lib/aux/`:
* `read_write.h`: reading and writing routines;
* `functions.h`: custom functions for vector algebra;
* `integrator.h`: custom 1d function integrator;
* `NRutil.h`: all auxiliary structures and interpolation function;
* `dsolve.h`: some Runge-Kutta solvers (partly from "Numerical Recipes in C" <sup>[4](#myfootnote4)</sup> textbook and from [here](https://people.sc.fsu.edu/~jburkardt/c_src/)) to solve nonlinear ODE with rk4, rkf45 and rk8 routines;
* `rk4.h`, `rkf45.h`: rk routines.


### Compiling and running
0. Configure the compiler by running `python configure.py -<flag1> -<flag2> ...`. See below the available flags.
1. Compile with `make all` in the root directory, the executable is generated in the `bin/` directory.
2. The input parameters are passed with the `loki.input` file, from the same directory as the executable file.
* The input file name is specified with a flag when running the code `... -i loki.input`.
* Output path is specified with the `-o` flag, default value is "output/".
3. Run the code as follows:
```bash
$ ./loki -i loki.input
```
or for MPI version:
```bash
$ mpiexec -n <number of CPUs> loki -i loki.input
```
#### Available flags for `configure.py`
* `-mpi`: MPI version of the code;
* `-intback`: enable backwards integration to find the polar cap density;
* `-findnan`: exit when find NaN-s;
* if `-mpi` not specified, compiles with `g++`.


#### Example run:
```bash
$ python configure.py -intback
$ make all
$ cd bin
$ ./loki -i loki.input -o my_output
```
or
```bash
$ python configure.py -mpi -intback -findnan
$ make all
$ cd bin
$ mpiexec -n 16 loki -i loki.input -o my_mpi_output
```

### Visualization

For a non-MPI run the output is written in a standard `utf-8` encoding and can be directly read in python as follows:
```python
import numpy as np
data = np.loadtxt("<output_file>")
```

In MPI run the file is written in a binary format and should be treated differently. `vis/` folder in this repository contains a library to read the binary data. Import the `loki` module (`loki.py` file should be in the same directory as the python file or jupyter notebook) and read the data as follows:
```python
import loki
data = loki.loki_read("<output_file>")
```

### Things to do
- Make initial radiation point for any given magnetic field (right now only takes dipole field, which is not bad, since we start close to the star).
- Fix `Step size too small in routine RKQC` issue and `nan`-s for high density plasma.
- Read `hdf5` __b__-field data.

### Literature
<a name="myfootnote1">1</a>: _A.S. Andrianov, V.S. Beskin_, [Limiting polarization effect—a key link in investigating the mean profiles of radio pulsars](https://link.springer.com/article/10.1134/S1063773710040031) __(2010)__

<a name="myfootnote2">2</a>: _V.S. Beskin, A.A. Philippov_, [On the mean profiles of radio pulsars I: Theory of the propagation effects](https://arxiv.org/pdf/1107.3775.pdf) __(2012)__

<a name="myfootnote3">3</a>: _H. Hakobyan, V. Beskin, A. Philippov_, [On the mean profiles of radio pulsars II: Reconstruction of complex pulsar light-curves and other new propagation effects](https://arxiv.org/abs/1704.08743) __(2017)__

<a name="myfootnote4">4</a>: _Brian P. Flannery, Saul Teukolsky, William H. Press, and William T. Vetterling_, [Numerical Recipes in C: The Art of Scientific Computing](http://www2.units.it/ipl/students_area/imm2/files/Numerical_Recipes.pdf)

---

### Memo for git

- Pull updates every time before starting a work with:
```bash
$ git pull
```

- Switch between branches using:
```bash
$ git checkout <branchname>
```

> NOTE! Always `commit` or `stash` your changes on a branch before changing to another branch. On committing and stashing see below.

- Display all branches including the one you're using right now:
```bash
$ git branch
```

- Create new branch by doing
```bash
$ git checkout -b <newbranchname>
```

- Adding updated files, committing and pushing to the repo:
```bash
$ git add <files>
$ git commit -m 'Some comment here'
$ git push origin <branch>
```
> To add all files you can simply do `$ git add .`. To push to `master` branch do `$ git push origin master`.

- If you made some changes in `branch_A` and want to change to a `branch_B`, but you don't want to commit your changes on `branch_A` yet, you can (should!) `stash` them (keep them save until you get back). You can do that with:
```bash
$ git stash # <- stash files on branch_A
$ git checkout branch_B # <- switch to branch_B
$ ... # <- do some work on branch_B
$ git checkout branch_A # <- switch back to branch_A
$ git stash apply # <- get your stashed files back on branch_A
```

---

__BONUS.__ To see the active branch in your terminal with a highlighted color, you can add this to your `~/.bash_profile` file, save and then restart your terminal:
```bash
parse_git_branch() {
  git branch 2> /dev/null | sed -e '/^[^*]/d' -e 's/* \(.*\)/ (\1) /'
}
export PS1="\u \W \[\033[32m\]\$(parse_git_branch)\[\033[00m\]$ "
```
This will change your terminal appearance to:
```bash
<user> <directory> (<branch>) $
```
