Description
===========
This is a human operator control model for a person controlling a bicycle. The
bicycle model is based on the Whipple bicycle dynamic model. The control
portion is based on the crossover model and includes basic preview. The model
is capable of tracking a path.

File Descriptions
===================
whipple_pull_force_abcd.m : Generates the linearized Whipple model about the
upright constant velocity equilibrium point for various parameter sets. It
includes an additional lateral pull force input.

par_text_to_struct.m : Loads in a parameter file to Matlab's structure data
type.

WhippleModel.mdl : The simulink model which is used to calculate the
gains, the various transfer function for the open and closed loops, and the
handling qualities metric for both steer torque and roll torque inputs.

generate_data.m : Generates the data (transfer functions, simulation results,
and handling quality metric) by simulating and perturbing the Simulink model.

create_ieee_paper_plots.m : Generates most of the plots for our first journal
paper on the topic.

fix_ps_linstyle.m : [A file from the Matlab file
exchange](http://www.mathworks.com/matlabcentral/fileexchange/17928) for
improving plot lines.

load_bikes.m : Loads the data from generate_data for multiple bicycles and
speeds into one structure.

ieee.m : Loads the data needed for create_ieee_paper_plots and runs the
function. This is not part of create_ieee_paper_plots to faciltate the
separation of the data generation and the plotting mainly for debugging
purposes.

lane_change.m : Generates path data for a single or double lane change manuever
at a particular speed.

benchmark_geometry_figure.m : Generates a postscript drawing of the bicycle
dimenions. Requires [pshacker](http://bicycle.tudelft.nl/schwab/pshacker/).

Requirements
============
- Matlab 2010a (7.10.0)
- Matlab Control Systems Toolbox
- Simulink
- [pshacker](http://bicycle.tudelft.nl/schwab/pshacker/) (for the bicycle
  geometry plot)

Example Use
===========
The core function in this package is `generate_data.m`. The simplest method of
using it is to supply it with a bicycle name and a speed. It will then find the
appropriate gains for a stable model, compute all of the loop transfer
functions and simulate a double lane change manuever. The function outputs all
of this data for easy plotting and analysis.

Generate data for the Benchmark bicyle at 4.8 m/s and explore some of the
output.

```matlab
>> data = generate_data('Benchmark', 4.8);
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

% Generate the plots for the paper on this model.
ieee
```
