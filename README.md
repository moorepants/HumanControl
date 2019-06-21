Description
===========

This software implements a human operator control model for a person
controlling a bicycle in lateral flat ground maneuvers. The bicycle model is
based on the linear Whipple bicycle dynamic model. The controller is based on
the principles of the crossover model and includes basic visual preview. The
model is capable of tracking arbitrary a paths in the ground plane.

This software serves as the companion code to several publications. A tagged
release corresponds to each related publication.

- [[tag:ieee](https://github.com/moorepants/HumanControl/releases/tag/ieee)] Hess, R.; Moore, J.K.; Hubbard, M., "Modeling the Manually
  Controlled Bicycle," Systems, Man and Cybernetics, Part A: Systems and
  Humans, IEEE Transactions on, vol.42, no.3, pp.545,557, May 2012,
  http://dx.doi.org/10.1109/TSMCA.2011.2164244.
- [[tag:dissertation](https://github.com/moorepants/HumanControl/releases/tag/dissertation)] J. K. Moore, "Human Control of a Bicycle," Doctor of
  Philosophy, University of California, Davis, CA, August 2012,
  http://moorepants.github.io/dissertation/.
- [[tag:bmd2016](https://github.com/moorepants/HumanControl/releases/tag/bmd2016)]
  J. Moore, M. Hubbard, and R. A. Hess, "An Optimal Handling Bicycle," in
  Proceedings of the 2016 Bicycle and Motorcycle Dynamics Conference, 2016,
  http://moorepants.github.io/dissertation/.

Requirements
============

- Matlab 2010a (7.10.0)
- Matlab Control Systems Toolbox
- Simulink
- [pshacker](http://bicycle.tudelft.nl/schwab/pshacker/) (for the bicycle
  geometry plot)
- [cmaes.m 3.61.beta](http://cma.gforge.inria.fr/cmaes_sourcecode_page.html#matlab)

Installation
============

Download the [current development
version](https://github.com/moorepants/HumanControl/archive/master.zip),
[desired release](https://github.com/moorepants/HumanControl/releases), or
clone the git repository. Set the root of the repository or extracted directory
as the working directory in MATLAB. Download `cmaes.m` and place `cmaes.m` into
the working directory.

File Descriptions
===================

`benchmark_geometry_figure.m` : Generates a postscript drawing of the bicycle
dimensions. Requires [pshacker](http://bicycle.tudelft.nl/schwab/pshacker/).

`bicycle_state_space.m` : Function that returns the state space system of the
Whipple model linearized about the nominal configuration and the supplied
speed.

`bmd2016.m` : Script that plots the peak of the handling quality transfer
function magnitude versus change in trail.

`convert_variable.m` : Function that returns the name and order of the given
variable in the output type.

`create_ieee_paper_plots.m` : Generates most of the plots for our first journal
paper on the topic.

`fix_ps_linstyle.m` : [A file from the Matlab file
exchange](http://www.mathworks.com/matlabcentral/fileexchange/17928) for
improving plot lines.

`generate_data.m` : Generates the data (transfer functions, simulation results,
and handling quality metric) by simulating and perturbing the Simulink model.

`handling_vs_par.m` : Generates a plot of peak HQM magnitude versus changes in
a particular bicycle parameter.

`handling_vs_speed.m` : Generates a plot of the peak HQM versus change in speed
for all six bicycles in the IEEE paper.

`heading_track_analytic.py` : Generates the closed heading loop state space
description in symbolic form.

`ieee.m` : Loads the data needed for `create_ieee_paper_plots` and runs the
function. This is not part of `create_ieee_paper_plots` to facilitate the
separation of the data generation and the plotting mainly for debugging
purposes.

`lane_change.m` : Generates path data for a single or double lane change
maneuver at a particular speed.

`lateral_track_analytic.py` : Generates the closed lateral deviation loop state
space description in symbolic form.

`load_bikes.m` : Loads the data from `generate_data` for multiple bicycles and
speeds into one structure.

`lookup_gains.m` : Function that returns a guess for the gains based on
precomputed gains at various speeds using linear interpolation/extrapolation.

`objective.m` : Function that computes the optimization objective value, i.e.
the peak of the HQM.

`optimal_bicycle.m` : Script that runs the CMES optimizer to discover bicycles
with minimal peak HQM values.

`par_text_to_struct.m` : Loads in a parameter file to Matlab's structure data
type.

`plot_gains.m` : Function to plot the gains stored in a gain file.

`system_state_space.m` : Function that returns the closed loop system state
space of the Hess bicycle-rider system for heading or lateral deviation
tracking with bicycle steer torque as the controlled input.

`test_system_state_space.m` : Unit tests for system_state_space.m.

`varargin.m` : Function that returns a structure from a cell array of pairs.

`WhippleModel.mdl` : The Simulink model which is used to calculate the gains,
the various transfer function for the open and closed loops, and the handling
qualities metric for both steer torque and roll torque inputs.

`whipple_pull_force_abcd.m` : Generates the linearized Whipple model about the
upright constant velocity equilibrium point for various parameter sets. It
includes an additional lateral pull force input.

`write_gains.m` : Function that adds the provided gains to a gain file.

`zipforieee.py` : Script that zips up all the IEEE paper figures.

Example Use
===========

The core function in this package is `generate_data.m`. The simplest method of
using it is to supply it with a bicycle name and a speed. It will then find the
appropriate gains for a stable model, compute all of the loop transfer
functions and simulate a double lane change maneuver. The function outputs all
of this data for easy plotting and analysis.

Generate data for the Benchmark bicycle at 4.8 m/s.

```matlab
>> data = generate_data('Benchmark', 4.8);
```

You should see this basic output:

```
-------------------------------------------------------------------------------
Benchmark at 4.80 m/s.
-------------------------------------------------------------------------------
Parameters for the Benchmark bicycle and rider have been loaded.
Calculating the A, B, C, D matrices for 4.80 m/s
A, B, C, D calculated in 0.0291 seconds.
Finding the loop transfer function of the Delta loop with a start guess of 44.6406.
Delta loop gain set to 44.6406.
Finding the loop transfer function of the PhiDot loop with a start guess of -0.0529.
PhiDot loop gain set to -0.0529.
Finding the loop transfer function of the Phi loop with a start guess of 13.2253.
Phi loop gain set to 13.2253.
Finding the loop transfer function of the Psi loop with a start guess of 0.1720.
Psi loop gain set to 0.1720.
Finding the loop transfer function of the Y loop with a start guess of 0.1001.
Y loop gain set to 0.1001.
Gains are set to: kDelta = 44.641
                  kPhiDot = -0.053
                  kPhi = 13.225
                  kPsi = 0.172
                  kY = 0.100
Finding the Delta closed loop transfer function.
Finding the PhiDot closed loop transfer function.
Finding the Phi closed loop transfer function.
Finding the Psi closed loop transfer function.
Finding the Y closed loop transfer function.
Finding the Delta open loop transfer function.
Finding the PhiDot open loop transfer function.
Finding the Phi open loop transfer function.
Finding the Psi open loop transfer function.
Finding the Y open loop transfer function.
Finding the handling quality metric.
Gains written to gains/BenchmarkSteerGains.txt
Simulating the tracking task.
Simulation finished in 0.185 seconds.
Done.
```

`data` is now a structure that contains all of the output from the function.

```matlab
>> data

data =

             speed: 4.8000
               par: [1x1 struct]
          modelPar: [1x1 struct]
       closedLoops: [1x1 struct]
         openLoops: [1x1 struct]
    handlingMetric: [1x1 struct]
              time: [209x1 double]
           command: [209x5 double]
            inputs: [209x3 double]
           outputs: [209x18 double]
        outputsDot: [209x18 double]
              path: [209x1 double]
          gainMuls: [1 1 1 1 1]

```

You can check the eigenvalues of the bicycle model:

```matlab
>> eig(data.modelPar.A)
```

The loop transfer functions (open and closed) can be visualized with a Bode
plot:

```matlab
>> bode(tf(data.closedLoops.delta.num, data.closedLoops.delta.den))
```

The handling quality metric can be accessed by:

```matlab
>> tf(data.handlingMetric.num, data.handlingMetric.den)
```

The inputs and outputs from the simulation can be be visualized with:

```matlab
>> figure(1)
>> plot(data.time, data.inputs)
>> figure(2)
>> plot(data.time, data.outputs)
```

The model gains can be accessed with:

```matlab
>> data.modelPar.kDelta
>> data.modelPar.kPhiDot
>> data.modelPar.kPhi
>> data.modelPar.kPsi
>> data.modelPar.kY
```

`generate_data.m` can also take many optional arguments.

```matlab
% Generate the data set for the Benchmark bicycle at 5.0 m/s with roll as the
% input.
data = generate_data('Benchmark', 5.0, 'input', 'Roll');

% Generate the data set for the Fisher bicycle at 7.5 m/s with steer input
% and show the graphs.
data = generate_data('Fisher', 7.5, 'input', 'Steer', 'plot', 1);

% Generate the data set for the Browser bicycle at 2.5 m/s with steer as an
% input and multiply the five gains by various values and show the graphs.
data = generate_data('Browser', 2.5, 'plot', 1, 'gainMuls', [1.1, 1.1, 0.9, 1.0, 0.8])

% Generate the data set for the Bianchi Pista bicycle at 7.5 m/s with steer as the
% input and a single lane change as the manuever.
data = generate_data('Pista', 7.5, 'laneType', 'single');
```

We used `generate_data.m` to create the data and plots for the journal paper on
the subject. Run `ieee.m` to generate all of the plots. This code has examples
of how to extract and plot all of the data that is made available from
`generate_data.m`

```matlab
% Generate the plots for the paper on this model.
ieee
```
