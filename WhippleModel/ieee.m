% This script loads in the data used for the IEEE paper and creates the plots.
clc; close all;

bikes = {'Browserins', 'Browser', 'Pista', ...
         'Fisher', 'Yellow', 'Yellowrev'};

data = load_bikes(bikes, 'Steer');

data.Benchmark = generate_data('Benchmark', 5.0);

rollData = generate_data('Benchmark', 5.0, 'input', 'Roll');

create_ieee_paper_plots(data, rollData)
