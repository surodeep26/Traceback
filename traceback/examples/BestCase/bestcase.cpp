#include <QFile>
#include <QFileInfo>
#include <QStringList>
#include <QTextStream>
#include <float.h>	// DBL_MAX

using namespace std;

static void printUsage( const char *progname )
{
	fprintf( stderr, "Usage: %s <data file>\n", progname );
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

	/* Create an object for our data file */
	QFile dataFile( cfName );
	dataFile.open( QIODevice::ReadOnly );
	QTextStream inStream ( &dataFile );
	inStream.setCodec( "UTF-8" );

	QString bestLine, readLine;
	double minimum = DBL_MAX, criterion, d, r1, r2;

	QStringList values;
	/* Now read the file line by line */
	while ( !inStream.atEnd() ) {
		readLine = inStream.readLine();

		/* Drop comments and empty lines */
		if ( readLine.startsWith('#') || readLine.isEmpty() ) continue;
		values = readLine.simplified().split(' ');

		/* Here we define the criterion for "best" */
		d  = values[0].toDouble();	// separation between star 1 and 2
		r1 = values[2].toDouble();	// separation between star 1 and assoc.
		r2 = values[3].toDouble();	// separation between star 2 and assoc.
		values.clear();
		criterion = /*2.**/d*d + r1*r1 + r2*r2;
//		criterion = d;
		if ( criterion < minimum ) {
			minimum = criterion;
			bestLine = readLine;
		}
	} /* while() */

#if ( QT_VERSION < QT_VERSION_CHECK(5,14,0) )
# define ENDL endl
#else
# define ENDL Qt::endl
#endif

	out << bestLine << ENDL << ENDL;
	out << "-------------------------------------------" << ENDL;
	out << QString( "Abst%1nde : %2 pc, %3 pc, %4 pc" )
	       .arg(QString::fromUtf8("ä"))
	       .arg(bestLine.simplified().split(' ')[0])
	       .arg(bestLine.simplified().split(' ')[2])
	       .arg(bestLine.simplified().split(' ')[3])
	        << ENDL;

	double zeit = bestLine.simplified().split(' ')[1].toDouble() / 1.e6;

	out << QString("Zeitpunkt: %1 Myr")
	       .arg( zeit )
	    << ENDL << ENDL;

	out << QString("%1 Oph:")
	       .arg(QString::fromUtf8("ζ"))
	    << ENDL;

	QString Long( bestLine.simplified().split(' ')[8] ),
	        lat( bestLine.simplified().split(' ')[9] ),
	        dist( bestLine.simplified().split(' ')[10] );

	out << QString("l=%1%4 b=%2%4 R=%3 pc")
	       .arg( Long )
	       .arg( lat )
	       .arg( dist )
	       .arg( QString::fromUtf8( "°" ) )
	    << ENDL << ENDL;

	QString pulsar( QFileInfo(dataFile).baseName().simplified().split('+')[2] + ':');

	out << pulsar << ENDL;

	Long = bestLine.simplified().split(' ')[15];
	lat = bestLine.simplified().split(' ')[16];
	dist = bestLine.simplified().split(' ')[17];
	out << QString("l=%1%4 b=%2%4 R=%3 pc")
	       .arg( Long )
	       .arg( lat )
	       .arg( dist )
	       .arg( QString::fromUtf8( "°" ) )
	    << ENDL;

	out << QString("RV=%1 km/s")
	       .arg(bestLine.simplified().split(' ')[12])
	        << ENDL << ENDL;

	QString nova( dataFile.fileName().simplified().split('+')[1] + ':');

	out << nova << ENDL;

	Long = bestLine.simplified().split(' ')[22];
	lat = bestLine.simplified().split(' ')[23];
	dist = bestLine.simplified().split(' ')[24];
	out << QString("l=%1%4 b=%2%4 R=%3 pc")
	       .arg( Long )
	       .arg( lat )
	       .arg( dist )
	       .arg( QString::fromUtf8( "°" ) )
	    << ENDL;

	out << "-------------------------------------------" << ENDL;
	dataFile.close();
	return 0;
}

