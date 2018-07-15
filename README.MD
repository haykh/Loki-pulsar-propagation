# `Loki v1.4`
## C++ code to analize polarization properties of radiowaves propagating in the pulsar magnetosphere.

Theoretical basis for the program is given in the following papers <sup>[1](#myfootnote1)</sup>
<sup>[2](#myfootnote2)</sup> <sup>[3](#myfootnote3)</sup>.

### Libraries description

#####`lib/`:
* `constants.h`: stores all the necessary constants within the `constants::` namespace;
* `initialize.h`: some of the initialization routines used, such as finding the radiation point, reading input etc;
* `process_functions.h`: all the physical functions involved;
* `b_field.h`: functions to find the magnetic field for a given point;
* `RHS.h`: the right-hand side of the ODE (both in `float` and `double`);
* `pc_dens.h`: back integration routine to find the density on a polar cap for any given point on a ray trajectory.

#####`lib/aux/`:
* `read_write.h`: reading and writing routines;
* `functions.h`: custom functions for vector algebra;
* `integrator.h`: custom 1d function integrator;
* `NRutil.h`: all auxiliary structures and interpolation function;
* `dsolve.h`: some Runge-Kutta solvers (partly from "Numerical Recipes in C" <sup>[4](#myfootnote4)</sup> textbook) to solve nonlinear ODE with rk4, rkf45 and rk8 routines;
* `rk4.h`, `rkf45.h`: rk routines.


### Running
* Compile with `make all` in the root directory, the executable is generated in the `bin/` directory.
* The input parameters are passed with the `loki.input` file, from the same directory as the executable file.
* The input file name is specified with a flag when running the code `./loki -i loki.input`.
* Output path is specified with the `-o` flag, default value is "output/".

#### Example run:
```bash
$ make all
$ cd bin
$ ./loki -i loki.input -o my_output
```

### Things to do
- Write density routine for any given ___B___-field.
- Make initial radiation point for any given magnetic field (right now only takes dipole field, which is not bad, since we start close to the star).
- Fix `Step size too small in routine RKQC` issue for high density plasma.

### Literature
<a name="myfootnote1">1</a>: _A.S. Andrianov, V.S. Beskin_, [Limiting polarization effect—a key link in investigating the mean profiles of radio pulsars](https://link.springer.com/article/10.1134/S1063773710040031) __(2010)__

<a name="myfootnote2">2</a>: _V.S. Beskin, A.A. Philippov_, [On the mean profiles of radio pulsars I: Theory of the propagation effects](https://arxiv.org/pdf/1107.3775.pdf) __(2012)__

<a name="myfootnote3">3</a>: _H. Hakobyan, V. Beskin, A. Philippov_, [On the mean profiles of radio pulsars II: Reconstruction of complex pulsar light-curves and other new propagation effects](https://arxiv.org/abs/1704.08743) __(2017)__

<a name="myfootnote4">4</a>: _Brian P. Flannery, Saul Teukolsky, William H. Press, and William T. Vetterling_, [Numerical Recipes in C: The Art of Scientific Computing](http://www2.units.it/ipl/students_area/imm2/files/Numerical_Recipes.pdf)

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
$ git checkout -n <newbranchname>
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
