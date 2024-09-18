#include <QString>
#include <QTextStream>
#include <QList>
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

	/* Assign some values to shorter names for convenience.
	 * Note: we misuse the Limit parameter as a factor to inflate the association
	 * in the hope that might somehow compensate for the fact that we do not vary
	 * any input parameters */

	const mytype inflation   = simulParams.Limit;
	const QString fNameStar  = StarGalacticParam::simStarToFname( simulParams.Star1 );
	const QString fNameAssoc = StarGalacticParam::simStarToFname( simulParams.Assoc );

	/* Depending on StepSize width can be positive or negative */
	const int width          = simulParams.StepSize > 0. ? simulParams.Width : -simulParams.Width;
	const int steps          = simulParams.Steps;

#if ( QT_VERSION < QT_VERSION_CHECK(5,14,0) )
# define ENDL endl
#else
# define ENDL Qt::endl
#endif

	out << "Inflation = " << inflation << ENDL;
	foreach ( int starNo, StarGalacticParam::simStarToLineList( simulParams.Star1 ) ) {
		StarGalacticParam *star = new StarGalacticParam( fNameStar, starNo, StarGalacticParam::Numeric );
		foreach ( int assocNo, StarGalacticParam::simStarToLineList( simulParams.Assoc ) ) {
			StarGalacticParam *assoc = new StarGalacticParam( fNameAssoc, assocNo, StarGalacticParam::Numeric );
			QString output;
			for ( int i=0; i<steps; i++ ) {
				star->progressOrbit( width );
				assoc->progressOrbit( width );
				/* Check whether the star is inside the inflated association */
				if ( (*star || *assoc) < assoc->getRadius()*inflation ) {
					output = QString( "Star '%1' crossed '%2'\n" )
					         .arg( star->getHIP() )
					         .arg( assoc->getHIP() );
					out << output;
					break;
				}
			}
			delete assoc;
		}
		delete star;
	}
	return 0;
}

