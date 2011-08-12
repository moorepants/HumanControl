To Do List
==========

- Move the pull force location parameters into the parameter text files instead
  of whipple_pull_force_abcd.
- Move all of the data generation into ieee.m (phase portrait stuff in particular).
- Add the new values for the Browserins instead of the old.
- Make the varargin_to_structure take the default structure and add the
  varargins to the default stucture. This may need to be done recursively.
- Save correctly found gains to a file
- Add a way for the user to graphically choose where to measure the PhiDot neuro
  peak height from.
- Fix xlim to match all plots on the basic plots from generate_data.m
- Setup load_bikes to use Ron's gains instead of the one found by our autok
  routine
- Allow the user to set the neuro transfer function parameters as an option
  argument to generata_data
