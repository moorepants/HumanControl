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
This will generate the transfer function and simulation data for the a bicycle
in a lane change maneuver and plot basic time histories of the inputs and
outputs along with bode plots of the various transfer functions.

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
