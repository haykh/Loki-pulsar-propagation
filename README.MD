# `Loki v1.4`
## C++ code to analize polarization properties of radiowaves propagating in the pulsar magnetosphere.

Theoretical basis for the program is given in the following papers <sup>[1](#myfootnote1)</sup>
<sup>[2](#myfootnote2)</sup> <sup>[3](#myfootnote3)</sup>.

### Libraries description
`lib` folder contains all the necessary libraries.
* `Nrutil.h` and `diffeqsolver.h` are from "Numerical Recipes in C" <sup>[4](#myfootnote4)</sup> textbook and contain everything to solve nonlinear ODE with rk4 routine.
* `integrator.h` is a custom one dimensional function integrator. `functions.h` contains custom functions for vector algebra.
* All the physical functions involved are in `process_functions.h` and the right-hand side of the ODE is written in `RHS.h`. `initialize.h` is to store some of the initialization routines used, such as finding the radiation point.
* `constants.h` is to store all the necessary constants within the `constants::` namespace.

### Running
* Compile with `make all` in the root directory, the executable is generated in the `bin/` directory.
* The input parameters are passed with the `loki.input` file, from the same directory as the executable file.
* The input file name is specified with a flag when running the code `./loki -i loki.input`.
* Output path is specified with the `-o` flag, default value is "output/".

#### Example run:
```
$ make all
$ cd bin
$ ./loki -i loki.input -o my_output
```

### Things to do
- Write density routine for any given ___B___-field.
- Make initial radiation point for any given magnetic field.

### Literature
<a name="myfootnote1">1</a>: _A.S. Andrianov, V.S. Beskin_, [Limiting polarization effect—a key link in investigating the mean profiles of radio pulsars](https://link.springer.com/article/10.1134/S1063773710040031) __(2010)__

<a name="myfootnote2">2</a>: _V.S. Beskin, A.A. Philippov_, [On the mean profiles of radio pulsars I: Theory of the propagation effects](https://arxiv.org/pdf/1107.3775.pdf) __(2012)__

<a name="myfootnote3">3</a>: _H. Hakobyan, V. Beskin, A. Philippov_, [On the mean profiles of radio pulsars II: Reconstruction of complex pulsar light-curves and other new propagation effects](https://arxiv.org/abs/1704.08743) __(2017)__

<a name="myfootnote4">4</a>: _Brian P. Flannery, Saul Teukolsky, William H. Press, and William T. Vetterling_, [Numerical Recipes in C: The Art of Scientific Computing](http://www2.units.it/ipl/students_area/imm2/files/Numerical_Recipes.pdf)
