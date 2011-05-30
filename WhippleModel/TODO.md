To Do List
==========

IEEE Graphs
-----------
- Add the roll torque input handling quality to the main HQ graph.
- benchmarkRollOpen.eps
    - replace grid with single line at 0 or remove vertical grid
    - add markers for crossover frequencies
    - decrease font size of tick markers
    - psi may look better with solid line (switch with phi)
    - remove title
    - xlabel is cut off
- benchmarkSteerOpen.eps
    - replace grid with single line at 0 or remove vertical grid
    - add markers for crossover frequencies
    - remove title
    - xlabel is cut off
- *Distance.eps
    - decrease font size in legend
    - change legend to numbers instead of bicycle names
    - enlarge size (add -loose)
    - could plot all on one graph instead of subplots (would preserve scale)
    - ylabels are too crowded, go with two lines or move speed in as anotation
    - start at same time as paths.eps
- handling.eps
    - increase the size and add -loose
    - add the metric for the roll bike
    - decrease font size in legend
    - move level 3 up to 15 on the y axis
    - make the legend have numbers instead of bike names
- openBode.eps
    - decrease legend font size
    - increase size (-loose)
    - numbers instead of names in legend
    - decrease tick fontsize
    - make y axes black instead of gray
    - add an xlabel
- paths.eps
    - increase line thickness
    - numbers instead of names in legend
    - smaller legend font
    - make starting point for lane change more distinct
    - change y axis tick marks to say 2 meters at each height
    - plot -y instead so that you are looking from the top of the bike (label
      ticks correcly)
- phasePortraits.eps
    - decrease legend fontsize
    - decrease tick mark font size
    - increase size (add -loose)
    - lower xlabels are cut off
- rollDistance.eps
    - make the path the solid line
    - show enlargement of the countersteer
- eigenvalues
    - add the riders

Other
-----
- Move the full force location parameters into the parameter text files.
- Move all of the data generation into ieee.m (phase portrait in particular).
- Make path generation from simulink file to m file.
- Generate the lane change at the same distance instead of at the same
  time.
- Add the new values for the Browserins instead of the old.
