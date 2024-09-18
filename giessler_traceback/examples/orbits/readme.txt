This example calculates potential, epicyclic, and linear orbits
over time and writes X-, Y-, and Z-components into the file
'orbit.out'. The orbits can then be plotted with the scripts
'Orbit3d.gplt' and 'OrbitXY.gplt'.

The usage is

	orbits <configfile>

Besides the usual constants the following settings from the
configfile are used:

Simulation.InFile:
	name of input file containing the star parameters

Simulation.Star1:
	line number in InFile for the star the orbits of
	which shall be calculated

Simulation.Steps:
	number of points for each orbit

Simulation.Width:
	number of years between the steps

Simulation.StepSize:
	direction of orbit (forward or backward) and initial
	step size for numerical integration

The scripts use gnuplot for plotting the curves.
