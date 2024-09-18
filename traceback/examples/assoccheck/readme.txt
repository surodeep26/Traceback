This example traces back the orbits of star/association pairs
and checks whether the star had been inside the association
at any time during the tracing period. In case of a match
the IDs (HIP) are written to stdout.

The usage is

	assoccheck <configfile>

Besides the usual constants the following settings from the
configfile are used:

Simulation.Star1:
	File name and line numbers for the stars the orbit of
	which shall be calculated

Simuation.Assoc:
	File name and line numbers for the associations the orbit of
	which shall be calculated

Simulation.Steps:
	Number of points for each orbit

Simulation.Width:
	Number of years between the steps

Simulation.StepSize:
	Direction of orbit (forward or backward) and initial
	step size for numerical integration

Simulation.Limit:
	Used as a factor to inflate the association