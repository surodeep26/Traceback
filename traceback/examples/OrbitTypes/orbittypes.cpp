#include <QTextStream>
#include <QString>
#include "starparam.h"

using namespace std;

	// simulParams is a global variable defined elsewhere.
	// It will (and must) be read by StarGalacticParam::readConfig()
extern SimulParams simulParams;

static void printUsage( const char *progname )
{
	fprintf( stderr, "Usage: %s <config file> <epoch>\n", progname );
}

int main(int argc, const char *argv[])
{
	QTextStream out( stdout );	// For terminal output
	out.setCodec( "UTF-8" );

	// Check that the program is called with two arguments
	if ( argc < 3 ) {
		printUsage( argv[0] );
		exit( 1 );
	}

	// Read the arguments
	const char *cfName = argv[1];
	const mytype epoch = QString( argv[2] ).toDouble();

	// Read the config file
	StarGalacticParam::readConfig( cfName );

	// Assign some values to shorter names for convenience
	const QString fName = StarGalacticParam::simStarToFname( simulParams.Star1 );
	const int line1 = StarGalacticParam::simStarToLineNr( simulParams.Star1 );

	// Create three instances for the same star but with different orbit types
	StarGalacticParam starPotential( fName, line1, StarGalacticParam::Numeric ),
	                     starLinear( fName, line1, StarGalacticParam::Linear ),
	                   starEpicycle( fName, line1, StarGalacticParam::Epicycle );

	// It is convenient (but not necessary) to use a list of
	// the stars through which we can iterate later. The list
	// members are pointers to StarGalacticParam
	QList<StarGalacticParam*> myStars;
	myStars << &starPotential << &starEpicycle << &starLinear;

	// Print out what we've just read in
	starPotential.printStartValues( out, false );

	// Now we decide whether the calculated orbit will be
	// printed in heliocentric coordinates or in the
	// inertial system (giving this on the command
	// line is left as an exercise for the reader)
	const bool heliocentric = true;

#if ( QT_VERSION < QT_VERSION_CHECK(5,14,0) )
# define ENDL endl
#else
# define ENDL Qt::endl
#endif

	out << ENDL;

	// Print out the values for epoch 0
//	starPotential.printEquatorial( out, heliocentric );
//	starPotential.printGalactic( out, heliocentric );
	starPotential.printXYZ( out, heliocentric );
	out << QString( "  Distance from Sun:           %L1 pc\n" ).arg(starPotential.getDistance());
	starPotential.printUVW( out, heliocentric );
//	out << QString( "Heliocentric space velocity: %L1 km/s\n" ).arg(starPotential.getSpaceVelocity(heliocentric));
//	out << QString( "Radial velocity:             %L1 km/s\n" ).arg(starPotential.getRadialVelocity());
	out << ENDL;

	// We must assure that the initial step goes in the
	// correct direction. Otherwise the GSL will barf.
	// We simply set the width of the first step to
	// 100 years thus overwriting simulParams.StepSize
	if ( epoch >= 0.0 )
		starPotential.resetODE( 100 );
	else
		starPotential.resetODE( -100 );

	// Now we calculate the orbit parameters for the desired
	// epoch and print out what we are interested in.
	foreach (StarGalacticParam *star, myStars) {
		star->orbitAt( epoch );
		star->printInfo( out);

		star->printXYZ( out, heliocentric );
		out << QString( "  Distance from Sun:           %L1 pc\n" )
		       .arg(star->getDistance());

//		star->printUVW( out, heliocentric );
//		out << QString( "Heliocentric space velocity: %L1 km/s\n" )
//		       .arg(star->getSpaceVelocity(heliocentric));

//		star->printPM( out );
//		out << QString( "Radial velocity:             %L1 km/s\n" )
//		       .arg(star->getRadialVelocity(heliocentric));

		out << ENDL;
	} /* foreach */

	// That's it.
	return 0;
}

