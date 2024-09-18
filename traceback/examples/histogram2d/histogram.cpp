#include <QString>
#include <QFile>
#include <QTextStream>
#include <QStringList>
#include <QProcess>
#include <gsl/gsl_histogram2d.h>
#include <math.h>

using namespace std;

#if ( QT_VERSION < QT_VERSION_CHECK(5,14,0) )
# define ENDL endl
#else
# define ENDL Qt::endl
#endif

static void plot( QString message, QProcess &p, QTextStream &file )
{
	p.write( message.toUtf8().constData() );
	p.write( "\n" );
	file << message << ENDL;
}

static QString tau = QString::fromUtf8("τ");

static void stop( QString msg )
{
	QTextStream err(stderr);
	err << msg << ", stop." << ENDL;
	exit( 2 );
}

static void usage(const char *progname)
{
	QTextStream err( stderr );

	err.setCodec( "UTF-8" );
	err << QString("usage: %1 <datafile> [dmin dmax tmin tmax] [create frequency density?]")
	       .arg(progname) << ENDL << ENDL;
	err << QString("where\n"
	               "  <datafile>: name of the output file from traceback\n"
	               "  dmin      : lower distance limit of the 68% area [pc]\n"
	               "  dmax      : upper distance limit of the 68% area [pc]\n"
	               "  tmin      : lower %1 limit of the 68% area [Myr]\n"
	               "  tmax      : upper %1 limit of the 68% area [Myr]\n"
	               "  [create frequency density?]: either 'yes' and last argument, or omitted (meaning 'no')")
	       .arg( tau );
	err << ENDL;

	exit( 1 );
}

int main(int argc, const char *argv[])
{
	QLocale::setDefault( QLocale("en_US") );

	if ( argc < 2 )
		usage(argv[0]);

	bool valuesOk = true, ok = true;

	const bool withLimits = (argc > 5);

	const double dmin = withLimits ? QString(argv[2]).toDouble( &ok ) : 0.0;     valuesOk &= ok;
	const double dmax = withLimits ? QString(argv[3]).toDouble( &ok ) : 10000.0; valuesOk &= ok;
	const double tmin = withLimits ? QString(argv[4]).toDouble( &ok ) : 0.0;     valuesOk &= ok;
	const double tmax = withLimits ? QString(argv[5]).toDouble( &ok ) : 10000.0; valuesOk &= ok;

	if ( withLimits && !valuesOk )
		usage(argv[0]);

	QString fName( argv[1] );

	const int distCol = 1;
	const int tauCol  = 2;

	const bool density = (QString(argv[argc-1]).toUpper() == "YES");

	// Open data file
	QFile dataFile( fName );

	if ( !dataFile.open( QIODevice::ReadOnly ) ) {
		stop( QString("Unable to open file '%1'")
		      .arg(argv[1]) );
	}

	const double distBin = 2.0;		// [pc]
	const double tauBin  = 0.02;	// [Myr]

	// In case of PDF the unit for the bin size is kpc*Myr
	const double binsize = density ? distBin/1000.0 * tauBin : 1.;

	QTextStream dataStream( &dataFile );
	dataStream.setCodec( "UTF-8" );
	QString dataLine;

	/* Read the number of orbits from datafile */
	double orbits = 0.0;
	do {
		dataLine = dataStream.readLine().simplified();
	} while (!dataLine.startsWith( "# Orbits:" ) && !dataStream.atEnd() );

	if ( dataLine.startsWith( "# Orbits:" ) ) {
		orbits = dataLine.section(':', 1, 1).toDouble( &ok );
	}

	if ( !(dataLine.startsWith("# Orbits:") && ok) ) {
		orbits = 3000000.0;
		qWarning( "Cannot read number of orbits! Assuming %s\n",
		          QString("%L1").arg(orbits, 0, 'f', 0).toLatin1().constData() );
	}

	/* This will create a histogram structure with bin size 2pc x 0.02Myr
	 * as in Tetzlaff et al. 2010:
	 * First, allocate a histogram with 15 bins in X direction and
	 * 100 bins in Y direction */
	gsl_histogram2d *hist = gsl_histogram2d_alloc( lrint(30.0/distBin), lrint(2.0/tauBin) );

	/* Now we set the minimum and maximum values for both directions */
	gsl_histogram2d_set_ranges_uniform( hist, 0., 30., 0., 2. );

	/* We normalize the values to get some kind of probability
	 * densitiy per Megayear and parsec.
	 *
	 * By not dividing through the bin size we get a "probability
	 * per bin" just like the above mentioned paper.
	 */
	const double scale = density ?
	                         1. :		// unity
	                         10000.;	// 0.01 %

	const double accu = scale /
	                    binsize /
	                    orbits;

	// TODO: read the bin size from the results file

	QTextStream out( stdout );
	out << QString( "Reading file '%1'..." )
	       .arg( fName );

	dataStream.seek( 0 );
	int linenumber = 0;

	while ( !dataStream.atEnd() ) {
		double dist, tau;

		linenumber++;

		// Read line and remove any whitespace
		dataLine = dataStream.readLine().simplified();

		// Skip comment and empty lines
		if ( dataLine.isEmpty() || dataLine.startsWith( '#' ) )
			continue;

		// Convert columnX and columnY to doubles
		dist = dataLine.section(' ', distCol-1, distCol-1).toDouble( &ok ); valuesOk &= ok;
		tau  = dataLine.section(' ', tauCol-1,  tauCol-1 ).toDouble( &ok ); valuesOk &= ok;

		if ( !ok ) {
			gsl_histogram2d_free( hist );
			dataFile.close();
			stop( QString("Cannot read values at line #%1")
			      .arg(linenumber) );
		}

		// Convert to Myr and change sign
		tau /= -1e6;

		// Increment the histogram by the value accu
		gsl_histogram2d_accumulate( hist, dist, tau, accu );
	}
	dataFile.close();

	out << "done." << ENDL;

	/* Write the histogram to a file *.hist for later use.
	 * The file is made executable so just running it
	 * will recreate the histogram. It can also be loaded
	 * in gnuplot. Note that gnuplot version 5 or newer
	 * must be used. */
	QString outFilename( fName );
	outFilename.truncate( outFilename.lastIndexOf('.'));
	outFilename += ".hist2d";

	QFile outFile( outFilename );
	if ( !outFile.open( QIODevice::WriteOnly ) ) {
		gsl_histogram2d_free( hist );
		stop( QString("Cannot open output file '%1'")
		      .arg(outFilename) );
	}

	outFile.setPermissions( outFile.permissions() |
	                        QFile::ExeOwner |
	                        QFile::ExeGroup |
	                        QFile::ExeOther );

	QTextStream outStream( &outFile );
	outStream.setCodec( "UTF-8" );

	// Add some information
	outStream << "#!/usr/bin/gnuplot" << ENDL;
	outStream << "# Gnuplot version 5 is needed!" << ENDL;

	/* Draw the histogram using gnuplot */
	QProcess p;

	p.start( "/usr/bin/gnuplot", QStringList() );
	p.waitForStarted();
	plot( "reset", p, outStream );
	plot( "set encoding utf8", p, outStream );
	plot( "set term qt persist size 900,700 enhanced font \"CMU Serif,14\"", p, outStream );
	plot( "set size .95,1", p, outStream );
	plot( "set title \"{/:Bold"
	      + tau
	      + QString( " %1 {/:Italic d_{min} }} (%2)\"" )
	      .arg( QString::fromUtf8("–"))
	      .arg( fName ), p, outStream );
	plot( "set xlabel \"{/:Italic d_{min} } [pc]\"", p, outStream );
	plot( "set ylabel \"" + tau + " [Myr]\"", p, outStream );

	if ( !density ) {
		plot( "set cblabel \"probability per bin [10^{&{n}-2} %]\"", p, outStream );
	} else {
		plot( "set cblabel \"relative frequency density [ (kpc"
		      + QString::fromUtf8("∙")
		      + "Myr)^{&{n}-1} ]\"", p, outStream);
	}
	plot( "unset key", p, outStream );
	plot( "set pm3d map", p, outStream );
	plot( "set pm3d interpolate 0,0", p, outStream );

	/* Create the Jet colormap from Matlab */
	plot( "set palette defined (0  0.0 0.0 0.5, \\\n"
	      "                     1  0.0 0.0 1.0, \\\n"
	      "                     2  0.0 0.5 1.0, \\\n"
	      "                     3  0.0 1.0 1.0, \\\n"
	      "                     4  0.5 1.0 0.5, \\\n"
	      "                     5  1.0 1.0 0.0, \\\n"
	      "                     6  1.0 0.5 0.0, \\\n"
	      "                     7  1.0 0.0 0.0, \\\n"
	      "                     8  0.5 0.0 0.0 )", p, outStream );
	plot( "set contour base", p, outStream );
	plot( "set cntrparam bspline", p, outStream );

	/* Draw the 68% contour line */
	double cont68 = 0.68*gsl_histogram2d_max_val(hist);
	plot( QString("set cntrparam levels discrete %1")
	      .arg( cont68 ), p, outStream );

	plot( QString("set cbrange [0:%1]")
	      .arg(gsl_histogram2d_max_val(hist)), p, outStream );

	plot( "$MyHist << EOD", p, outStream );

	out << QString( "Writing file '%1'..." )
	       .arg( outFilename );

	for (uint i=0; i<gsl_histogram2d_nx(hist); i++) {
		double xlower, xupper, ylower, yupper;
		QString plotLine;

		gsl_histogram2d_get_xrange( hist, i, &xlower, &xupper );
		for (uint j=0; j<gsl_histogram2d_ny(hist); j++) {
			gsl_histogram2d_get_yrange( hist, j, &ylower, &yupper );
			double value = gsl_histogram2d_get( hist, i, j );
			plotLine = QString( "%1 %2 %3" )
			           .arg( (xlower+xupper)/2., 5 )
			           .arg( (ylower+yupper)/2., 5 )
			           .arg( value );
			plot( plotLine, p, outStream );
		}
		plot( "", p, outStream );
	}
	plot( "EOD", p, outStream );
	plot( "splot $MyHist w pm3d", p, outStream );

	outStream << ENDL;
	outFile.close();
	p.closeWriteChannel();
	p.waitForFinished();

	out << "done." << ENDL;

	/* Now we have to seek all the input lines that fall into the 68% area
	 * and create a new output file containing only those lines */
	QString fName68( fName );
	fName68.truncate( fName68.lastIndexOf('.'));
	fName68 += "-68.out";

	QFile file68( fName68 );
	if ( !file68.open( QIODevice::WriteOnly ) ) {
		gsl_histogram2d_free( hist );
		stop( QString("Cannot open output file '%1'")
		      .arg(fName68) );
	}

	QTextStream out68( &file68 );
	out68.setCodec( "UTF-8" );
	out68 << QString( "# This file contains only the orbits "
	                  "that end up in the 68% area of the %1 - d_{min} "
	                  "histogram, and are limited by")
	            .arg( QString::fromUtf8( "τ")) << ENDL;
	out68 << QString( "#   %1 %2 d_min [pc] %2 %3\n"
	                  "#   %4 %2   %5  [Myr] %2 %6\n" )
	         .arg(dmin, 4)
	         .arg(QString::fromUtf8("≤"))
	         .arg(dmax, -4)
	         .arg(tmin, 4)
	         .arg(QString::fromUtf8("τ"))
	         .arg(tmax, -4) << ENDL;

	out << QString( "Writing file '%1'..." )
	       .arg(fName68);

	/* Re-open the data file */
	dataFile.open( QIODevice::ReadOnly );
	dataStream.seek( 0 );

	while ( !dataStream.atEnd() ) {
		double dist, tau, val;
		size_t i,j;

		dataLine = dataStream.readLine();
		if ( dataLine.isEmpty() || dataLine.startsWith( '#' ) ) {
			out68 << dataLine << ENDL;
			continue;
		}

		dist = dataLine.section( ' ', distCol-1, distCol-1, QString::SectionSkipEmpty ).toDouble();
		tau  = dataLine.section( ' ', tauCol-1,  tauCol-1,  QString::SectionSkipEmpty ).toDouble();

		tau /= -1e6;

		gsl_histogram2d_find( hist, dist, tau, &i, &j );
		val = gsl_histogram2d_get( hist, i, j );

		if ( (val > cont68) &&
		     (dist >= dmin) &&
		     (dist <= dmax) &&
		     (tau  >= tmin) &&
		     (tau  <= tmax) ) {
			out68 << dataLine << ENDL;
		}
	}

	dataFile.close();
	file68.close();

	out << "done." << ENDL;

	gsl_histogram2d_free( hist );
	return 0;
}
