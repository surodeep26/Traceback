#include <math.h>
#include <QFile>
#include <QStringList>
#include <QTextStream>

using namespace std;

static void printUsage( const char *progname )
{
	fprintf( stderr, "Usage: %s <input file> <output file>\n", progname );
}

#if ( QT_VERSION < QT_VERSION_CHECK(5,14,0) )
# define ENDL endl
#else
# define ENDL Qt::endl
#endif

int main(int argc, const char *argv[])
{
	/* Make sure the decimal separator is not a comma but a point */
	QLocale::setDefault( QLocale("en_US") );

	QTextStream out( stdout );	// For terminal output
	QTextStream errout( stderr );
	out.setCodec( "UTF-8" );
	errout.setCodec( "UTF-8" );

	/* Check that the program is called with an argument (the name of the input file) */
	if ( argc != 3 ) {
		printUsage( argv[0] );
		exit( 1 );
	}

	/* Read the argument */
	const char *cfName = argv[1];
	const char *ofName = argv[2];

	/* Create objects for our input and output files */
	QFile inFile( cfName ), outFile( ofName );
	inFile.open( QIODevice::ReadOnly );
	QTextStream inStream ( &inFile );
	inStream.setCodec( "UTF-8" );

	outFile.open( QIODevice::WriteOnly | QIODevice::Truncate );
	QTextStream outStream ( &outFile );
	outStream.setCodec( "UTF-8" );

	outStream << "#    RA                  DEC               Parallax        RV Maxw.        PMRA              PMDEC      ID           Remarks";
	outStream << ENDL;
	outStream << "# ---------------------------------------------------------------------------------------------------------------------------";
	outStream << ENDL;

	QString readLine, writeLine;
	double plx, plxerr, pma, pmaerr, pmd, pmderr;
	bool ok;
	int line = 3;

	/* Now read the input file line by line */
	while ( !inStream.atEnd() ) {
		readLine = inStream.readLine();

		/* Drop comments and empty lines */
		if ( readLine.startsWith('#') || readLine.isEmpty() ) continue;
		QStringList values = readLine.simplified().split(' ');

		if ( values.count() != 25 ) {
			errout << "Zeile " << line << "Spaltenzahl: " << values.count() << ENDL;
			continue;
		}

		plx = values[17].toDouble( &ok );
		if ( !ok ) {
			plx = values[20].toDouble( &ok );
			if ( !ok ) {
				plx = values[21].toDouble( &ok );
				if ( !ok ) {
					errout << "Zeile " << line << ": keine Parallaxe!" << ENDL;
					continue;
				}
			}
			plx = 1./plx;
		}

		plxerr = values[18].toDouble( &ok );
		if ( !ok || plxerr == 0.0 ) {
			plxerr = fabs( 0.1 * plx );
		}

		pma = values[11].toDouble( &ok );
		if ( !ok ) {
			errout << "Zeile " << line << ": kein PMRA!" << ENDL;
			continue;
		}

		pmaerr = values[12].toDouble( &ok );
		if ( !ok || pmaerr == 0.0 ) {
			pmaerr = fabs( 0.1 * pma );
		}

		pmd = values[14].toDouble( &ok );
		if ( !ok ) {
			errout << "Zeile " << line << ": kein PMDEC!" << ENDL;
			continue;
		}

		pmderr = values[15].toDouble( &ok );
		if ( !ok || pmderr == 0.0 ) {
			pmderr = fabs( 0.1 * pmd );
		}

		writeLine = QString("1  %1 %2 %3 %4     0 0 1   %5 %7   %8 %9   %10 (aus %11, Zeile %12)")
		            .arg( values[5].replace( ':', ' ' ), -17 )
		        .arg( values[8].replace( ':', ' ' ), -17 )
		        .arg( plx,    7, 'g', 4 )
		        .arg( plxerr, 6, 'g', 3 )
		        .arg( pma,    9, 'g', 8 )
		        .arg( pmaerr, 6, 'g', 3 )
		        .arg( pmd,    9, 'g', 8 )
		        .arg( pmderr, 6, 'g', 3 )
		        .arg( values[1], -12 )
		        .arg( cfName )
		        .arg( line++ );

		outStream << writeLine << ENDL;
	} /* while() */

	inFile.close();
	return 0;
}

