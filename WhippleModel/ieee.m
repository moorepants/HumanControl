clc; close all;

bikes = {'Benchmark', 'Browserins', 'Browser', 'Pista', ...
         'Fisher', 'Yellow', 'Yellowrev'};

data = load_bikes(bikes, 'Steer');

rollData = generate_data('Benchmark', 5.0, 'Roll', 0)

create_ieee_paper_plots(data, rollData)
