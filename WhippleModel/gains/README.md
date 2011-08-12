File Descriptions
=================
These files contain the five gains for a particular bicycle and at particular
speeds for both steer torque and roll torque inputs. These values are a mixture
of values chosen manually by Ron Hess and ones computed with the software.

The files are named in this fashion: '[bicycle short name][input]Gains.txt'
where the '[bicycle short name]' is a captilized word that corresponds to a
particular bicycle parameter set and '[input]' is either 'Steer' or 'Roll' for
the roll input, e.g. 'BenchmarkSteerGains.txt'.

The file contains a comma seperated header with the speed and names of the five
gains for the control model is this order: 'speed,kDelta,kPhiDot,kPhi,kPsi,kY'.
The following lines should have the values for the speed and the associated
gains for that speed.

These gains are used as starting guesses for the algorithm that homes in on the
best gain for the bicycle at a speed.
