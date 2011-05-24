% make the bode plots for the benchmark example at 5 m/s
data = generate_data('Benchmark', 5.0, 1.0, 0)

figure()
freq = {0.1, 20.0}
num = data.closedLoops.Delta.num
den = data.closedLoops.Delta.den
hold all
deltaBode = bodeplot(tf(num, den), freq)
setoptions(deltaBode,'PhaseMatching', 'on', 'PhaseMatchingValue', 0);
num = data.closedLoops.PhiDot.num
den = data.closedLoops.PhiDot.den
phiDotBode = bodeplot(tf(num, den), freq)
setoptions(phiDotBode,'PhaseMatching', 'on', 'PhaseMatchingValue', 0);
hold off
