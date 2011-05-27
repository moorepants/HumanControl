Description
-----------
This is a human operator control model for a person controlling a bicycle. The
bicycle model is based on the Whipple bicycle dynamic model. The control
portion is based on the crossover model and includes a crude preview. The model
is capable of tracking a path.

Requirements
------------
- Matlab 2010a (7.10.0)
- Simulink

Example Use
-----------
This will generate the transfer function and simulation data for the Fisher
bicycle in a lane change maneuver and plot basic time histories of the inputs
and outputs along with bode plots of the various transfer functions.

```matlab
% generate the data set for the Fisher bicycle at 7.5 m/s with only steer input
and show the graphs.
data = generate_data('Fisher', 7.5, 'Steer', 1., 1);
```
