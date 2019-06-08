# Information Theory Metrics for Q-State Potts Model

This is the companion code to the doctoral thesis `Information Theoretic Measures of Transitions to Collective Behaviour`, in which two key information theoretic metrics, Mutual Information (MI) and Global Transfer Entropy (GTE), are measured.

This code base deals with the Potts model portion of the thesis, and is capable of performing both the Glauber simulation code as well as the Density of States approach (Wang & Landau, 2001).

## Requirements

 - CMake 3.12 or above
 - C++11 for running simulations
 - MATLAB for visualing results

Libraries (included in `/external`):
 - Catch2, Clara, fmt

## Install

The code works on both Linux and Windows. Simply run cmake to configre and compile:

```
mkdir build; cd $_
cmake ..
cmake --build .
```

## Running the code

There are three main components available in the program, accessible via `./potts-entropy <component>`. 

 - `glauber` - Runs the straight-forward Glauber simulation approach, measuring a variety of metrics, included MI, TE, GTE, magnetisation, autocorrelation length and interface length. See sample script in `/tools/run_scripts` which will run a complete experiment for a given Q value
 - `dos` - Calculates the Density of States for given Q value and lattice size. The output of this is used as the input to other runs.
 - `dos_quantities` and `dos_gte_et` - Used for calculating thermodynamic MI and GTE, respectively, from the DoS files.

Use `./potts-entropy -h` for more information about all available components, or `./potts-entropy <component> -h` for help with a specific component.
