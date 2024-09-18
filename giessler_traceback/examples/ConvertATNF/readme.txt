This example will convert a list of objects downloaded from
the ATNF database

	http://www.atnf.csiro.au/people/pulsar/psrcat/

and containing the predefined variables
"Name, JName, RaJ, DecJ, PMRA, PMDec, PX, Dist, Dist_DM, Age, Age_i"

to a valid input file for traceback. The usage is

	convert <atnf_file> <output file>

When downloading from ATNF the list must not contain "last digit errors".
Choose the output style "long with erros" instead (see pulsar.txt).
