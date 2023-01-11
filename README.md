## General Purpose Trajectory Analyser (GPTA)
Copyright (C) 2021 by Paolo Raiteri

GPTA (General Purpose Trajectory Analyser) is a program that I developed over the years to perform routine day-to-day tasks required for my research as a computational materials scientist.
In particular, my main research tool is Molecular Dynamics simulations and I use GPTA on a daily basis to prepare the input structure and analyse the output trajectories.
The tasks that GPTA is able to perform can be broadly divided into two categories; manipulation of the coordinates of individual structures and post processing of (large) trajectory files.
Although GPTA has been rewritten multiple times, there are some legacy routines and there are some inconsistencies in style and appearance of the routines.
It is provided here in the hope that it may be useful to others.

Although you are not required to do so, the authors would consider it a courtesy if you submit to them any bug fixs and changes you consider to be worthwhile.

### Installation

To install GPTA from source you need compatible FORTRAN and C compilers. 
For a basic compilation, it should be enough to type
```
make serial
```
For more detailed information about the installation options see the user manual in the doc folder. 
 
### Contacts

If you want to report a bug please email the command used and a small input file to reproduce to error to 
<gpta.help@gmail.com>
or
<p.raiteri@curtin.edu.au>.

Although you are not required to do so, the author would consider it a courtesy if you submit to them any changes you consider to be worthwhile adding to the code. 
The goal would be to keep the development of GDIS more or less centralized.

