To Do List
==========

- Move the pull force location parameters into the parameter text files instead
  of whipple_pull_force_abcd.
- Move all of the data generation into ieee.m (phase portrait stuff in particular).
- Add the new values for the Browserins instead of the old.
- Make the varargin_to_structure take the default structure and add the
  varargins to the default stucture. This may need to be done recursively.
- Fix xlim to match all plots on the basic plots from generate_data.m
- Make the path filter frequency and time delay dependent on the
  characteristics of each bike at each speed, see the paper.
- Try out Ron's idea to find the gains for low speeds.
- Add a way for the user to graphically choose where to measure the PhiDot neuro
  peak height from.
- The simulink model could be significantly simplified by removing all switches
  related to the outputs and simply outputing a vector quantify for all output
  variables. It may be possible to do something similar with inputs too.
