#include <QFile>
#include <QString>
#include <QTextStream>
#include <QList>
#include <Vector.h>
#include <starparam.h>

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
	/* Make sure the decimal separator is not a comma but a point */
	QLocale::setDefault( QLocale("en_US") );

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
	const QString star1 = simulParams.Star1;
	const QString star2 = simulParams.Star2;
	const QString assoc = simulParams.Assoc;

	const QString fName1 = StarGalacticParam::simStarToFname( star1 );
	const QString fName2 = StarGalacticParam::simStarToFname( star2 );
	const QString fName3 = StarGalacticParam::simStarToFname( assoc );

	const int line1       = StarGalacticParam::simStarToLineNr( star1 );
	const int line2       = StarGalacticParam::simStarToLineNr( star2 );
	const int line3       = StarGalacticParam::simStarToLineNr( assoc );

	/* Depending on StepSize width can be positive or negative */
	const int width       = simulParams.StepSize > 0. ? simulParams.Width : -simulParams.Width;
	const int steps       = simulParams.Steps;
	const StarGalacticParam::OrbitType
	        oType = simulParams.Type == "Potential" ? StarGalacticParam::Numeric
	                                                : simulParams.Type == "Epicycle"
	                                                  ? StarGalacticParam::Epicycle
	                                                  : StarGalacticParam::Linear;

	/* It is convenient (but not necessary) to use a list of
	 * the stars through which we can iterate later. The list
	 * members are pointers to StarGalacticParam */
	QList<StarGalacticParam*> myStars;

	/* Create instances for the stars */
	myStars << new StarGalacticParam( fName1, line1, oType );
	if ( line2 )
		myStars << new StarGalacticParam( fName2, line2, oType );
	if ( line3 )
		myStars << new StarGalacticParam( fName3, line3, oType );

	/* Define the output file and open it */
	// Get a name for the output file
	QString oName( simulParams.OutFile );
	if ( oName == "auto" ) {
		QString outName( cfName );

		outName.truncate( outName.lastIndexOf( '.' ) );
		oName = outName + ".out";
	}
	QFile outfile( oName );
	outfile.open( QIODevice::ReadWrite | QIODevice::Truncate );
	QTextStream outStream( &outfile );
	outStream.setCodec("UTF-8");

	/* Write the start values into the output file for later reference.
	 * Prepend each line with a hash sign so pyxplot will ignore it. */
	foreach ( StarGalacticParam *star, myStars ) {
		star->printStartValues( outStream, true );
	}

#if ( QT_VERSION < QT_VERSION_CHECK(5,14,0) )
# define ENDL endl
#else
# define ENDL Qt::endl
#endif

	outStream << ENDL;

//	outStream << "# Columns: Time" << ENDL;

#define XYZ 0

#if !XYZ
	outStream << "# Columns:     Time Long1 Lat1 Dist1  ..." << ENDL;
	outStream << "# ColumnUnits: year degree degree parsec ";
	if ( line2 )
		outStream << "degree degree parsec ";
	if ( line3 )
		outStream << "degree degree parsec ";
	if ( line2 )
		outStream<< "parsec";
#else
	outStream << "# Columns: Time    X1 Y1 Z1    X2 Y2 Z2   Separation" << ENDL;
	outStream << "# ColumnUnits: year parsec parsec parsec ";
	if ( line2 )
		outStream << "parsec parsec parsec ";
	if ( line3 )
		outStream << "parsec parsec parsec ";
	if ( line2 )
		outStream<< "parsec";
#endif
	outStream << ENDL << ENDL;

	/* Now iterate through steps and stars */
	for ( int i=0; i<steps; i++ ) {
		Vector<mytype> coords(3);
		QString outLine;
		outLine = QString( "%1   " )
		          .arg( i*width, 8 );
		foreach ( StarGalacticParam *star, myStars ) {
			star->orbitAt( (mytype) i*width );
#if !XYZ
			coords = star->getGalactic( false, true );
#else
			coords = star->coordsHeliocentric()[0];
#endif
			outLine += QString( "%1 %2 %3     " )
			           .arg( coords[0], 10, 'f', 5 )
			           .arg( coords[1], 10, 'f', 5 )
			           .arg( coords[2], 10, 'f', 5 );
		}
		if ( line2 ) {
			outLine += QString( "%1" )
			           .arg( *myStars[0] || *myStars[1], 7, 'f', 4 );
		}
		outStream << outLine << ENDL;
	}

	outfile.close();

	while ( !myStars.isEmpty() )
		delete myStars.takeLast();

	return 0;
}

