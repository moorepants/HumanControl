clc; close all;

bikes = {'Benchmark', 'Browserins', 'Browser', 'Pista', ...
         'Fisher', 'Yellow', 'Yellowrev'};

data = load_bikes(bikes);

create_ieee_paper_plots(data)