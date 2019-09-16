# Code for Rapid Bayesian inference using preconditioning and Moment-Matching

This repository contains research code used to demonstrate the algorithms described
in the paper,

David J. Warne, Ruth E. Baker, Matthew J. Simpson.
Rapid Bayesian inference for expensive stochastic models. 2019.
Pre-print available on ArXiv [Link TBA](TBA).

## Developer
David J. Warne (david.warne@qut.edu.au), School of Mathematical Sciences, Science and Engineering Faculty, Queensland University of Technology.

Google Scholar: (https://scholar.google.com.au/citations?user=t8l-kuoAAAAJ&hl=en)

## License
This source code is licensed under the GNU General Public License Version 3. MCL: Monte Carlo Library
Copyright (C) 2019  David J. Warne

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

## Contents

This repositoty contains a number of programs, functions and scripts that are
designed to assist in reproducability of our work. 

```bash
The following are descriptions of the contents of each folder
|-- bin                  Containes C programs that implement the numerical examples presented in the paper
|-- data                 Synthetic data used for the weak Allee model and scratch assay model
|-- include              C header files required to build the examples programs
|-- src/
     |-- abc/              Contains C implementations of ABC rejection, SMC-ABC, PC-SMC-ABC and MM-SMC-ABC.
     |-- sim/              Contains C implementations of ODE/PDE numerical schemes with error control, and the stochastic lattice-based random-walk model.
     |-- util/             Various useful C functions used in the main programs and functions.
```
