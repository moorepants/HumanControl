Description
-----------
This is a human operator control model for a
person controlling a bicycle. The bicycle model is based on the Whipple bicycle
dynamic model. The control portion is based on the crossover model.

Example Use
-----------
This will generate the transfer function and simulation data for the Fisher
bicycle in a lane change maneuver and plot basic time histories of the inputs
and outputs along with bode plots of the various transfer functions.

```matlab
data = generate_data('Fisher', 7.5, 1., 1);
```
