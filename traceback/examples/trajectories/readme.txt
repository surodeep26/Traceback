This example calculates either potential, or epicyclic, or linear orbits
of up to three stars over time and writes formatted Galactic coordinates
into the output file.

The usage is

	trajectories <configfile>

Besides the usual constants the following settings from the
config file are used:

Simulation.InFile:
	name of input file containing the star parameters

Simulation.Star1, Star2, Assoc:
	line number in InFile for the stars the orbits of
	which shall be calculated; at least Star1 must be
	different from 0

Simulation.Steps:
	number of points for each orbit

Simulation.Width:
	number of years between the steps

Simulation.StepSize:
	direction of orbit (forward or backward) and initial
	step size for numerical integration

Simulation.Type:
	either 'Potential', 'Epicycle', or 'Linear'
	Note that using 'Epicycle' is still adventurous

Simulation.OutFile:
	name of the output file
