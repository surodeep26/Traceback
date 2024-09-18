#include <QFile>
#include <QString>
#include <QTextStream>
#include <QList>
#include <Vector.h>
#include "starparam.h"

using namespace std;

	/* simulParams is a global variable defined elsewhere.
	 * It will (and must) be read by StarGalacticParam::readConfig() */
extern SimulParams simulParams;

static void printUsage( const char *progname )
{
	fprintf( stderr, "Usage: %s <config file>\n", progname );
}

int main(int argc, const char *argv[])
{
	QTextStream out( stdout );	// For terminal output
	out.setCodec( "UTF-8" );

	/* Check that the program is called with an argument (the name of the config file) */
	if ( argc < 2 ) {
		printUsage( argv[0] );
		exit( 1 );
	}

	/* Read the argument */
	const char *cfName = argv[1];

	/* Read the config file */
	StarGalacticParam::readConfig( cfName );

	/* Assign some values to shorter names for convenience */
	const QString fName   = StarGalacticParam::simStarToFname( simulParams.Star1 );
	const int line1       = StarGalacticParam::simStarToLineNr( simulParams.Star1 );
	/* Depending on StepSize width can be positive or negative */
	const int width       = simulParams.StepSize > 0. ? simulParams.Width : -simulParams.Width;
	const int steps       = simulParams.Steps;

	/* Create three instances for the same star but with different orbit types */
	StarGalacticParam starPotential( fName, line1, StarGalacticParam::Numeric ),
	                     starLinear( fName, line1, StarGalacticParam::Linear ),
	                   starEpicycle( fName, line1, StarGalacticParam::Epicycle );

	/* It is convenient (but not necessary) to use a list of
	 * the stars through which we can iterate later. The list
	 * members are pointers to StarGalacticParam */
	QList<StarGalacticParam*> myStars;
	myStars << &starPotential << &starEpicycle << &starLinear;

	/* Define the output file and open it */
	QFile outfile( "orbit.out" );
	outfile.open( QIODevice::ReadWrite | QIODevice::Truncate );
	QTextStream outStream( &outfile );
	outStream.setCodec("UTF-8");

	/* Write the start values into the output file for later reference.
	 * Prepend each line with a hash sign so gnuplot will ignore it. */
	starPotential.printStartValues( outStream, true );

#if ( QT_VERSION < QT_VERSION_CHECK(5,14,0) )
# define ENDL endl
#else
# define ENDL Qt::endl
#endif

	outStream << ENDL;

	/* Now iterate through steps and orbits */
	for ( int i=0; i<steps; i++ ) {
		Vector<mytype> coords(3);
		QString outLine;
		outLine.clear();
		foreach ( StarGalacticParam *star, myStars ) {
			star->progressOrbit( (mytype) width );
			coords = star->coordsInertial()[0];
			outLine += QString( "%1 %2 %3     " )
			           .arg( coords[0], 9, 'f', 2 )
			           .arg( coords[1], 9, 'f', 2 )
			           .arg( coords[2], 9, 'f', 2 );
		}
		outStream << outLine << ENDL;
	}

	outfile.close();

	return 0;
}

