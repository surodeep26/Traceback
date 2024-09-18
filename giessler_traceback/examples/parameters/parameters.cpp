#include <QTextStream>
#include "starparam.h"

using namespace std;

extern SimulParams simulParams;

static void printUsage( const char *progname )
{
	fprintf( stderr, "Usage: %s <config file>\n", progname );
}

int main( int argc, const char *argv[] )
{
	QLocale::setDefault( QLocale( "en_US" ) );

	QTextStream out( stdout );
	out.setCodec( "UTF-8" );
	const char *cfName = argv[1];

	if ( argc != 2 ) {
		printUsage( argv[0] );
		exit( 1 );
	}

	StarGalacticParam::readConfig( cfName );

	const QString fName1 = StarGalacticParam::simStarToFname( simulParams.Star1 );
	const QString fName2 = StarGalacticParam::simStarToFname( simulParams.Star2 );
	const QString fNameAssoc = StarGalacticParam::simStarToFname( simulParams.Assoc );

	QList<StarGalacticParam *> stars;

	if ( QFile::exists( fName1 ) ) {
		foreach( int line, StarGalacticParam::simStarToLineList( simulParams.Star1 ) ) {
			if ( line ) stars.append( new StarGalacticParam( fName1, line, StarGalacticParam::Numeric ) );
		}
	}

	if ( QFile::exists( fName2 ) ) {
		foreach( int line, StarGalacticParam::simStarToLineList( simulParams.Star2 ) ) {
			if ( line ) stars.append( new StarGalacticParam( fName2, line, StarGalacticParam::Numeric ) );
		}
	}

	if ( QFile::exists( fNameAssoc ) ) {
		foreach( int line, StarGalacticParam::simStarToLineList( simulParams.Assoc ) ) {
			if ( line ) stars.append( new StarGalacticParam( fNameAssoc, line, StarGalacticParam::Numeric ) );
		}
	}

	foreach( StarGalacticParam *star, stars ) {
		star->printStartValues( out, false );
	}

	while ( !stars.isEmpty() )
		delete stars.takeLast();

	return 0;
}
