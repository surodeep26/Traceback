This will read a result file from examples/trace and produce a
2D histogram similar to the contour plots in Tetzlaff et al.,
MNRAS 402, 2010. It will also draw a single contour line at
68% from maximum. It will produce a new data file which will
contain only the lines falling into the 68% area.

The usage is

	histogram <datafile> [dmin dmax tmin tmax] [create frequency density?]

where
  <datafile>: name of the output file from traceback
  dmin      : lower distance limit of the 68% area [pc]
  dmax      : upper distance limit of the 68% area [pc]
  tmin      : lower τ limit of the 68% area [Myr]
  tmax      : upper τ limit of the 68% area [Myr]
  [create frequency density?]: either 'yes' and last argument, or omitted (meaning 'no').

An executable file with the extension .hist2d will be produced
which can be used to re-create the plot.
