#include <cmath>		// M_PI
#include <QFile>
#include <QTextStream>
#include <QTime>
#include <QElapsedTimer>
#include "Matrix.h"
#include "starparam.h"
#include "rng.h"
#include "threads.h"
#include "VecNorm.h"

using namespace std;

extern SimulParams simulParams;

static void printUsage( const char *progname )
{
	fprintf( stderr, "Usage: %s <config file>\n", progname );
}

int main(int argc, const char *argv[])
{
	QLocale::setDefault( QLocale("en_US") );

	QTextStream out( stdout );	// For terminal output
	out.setCodec( "UTF-8" );
	QElapsedTimer timer;				// To measure computational time
	int duration;
	const char *cfName = argv[1];

	// Check that the program is called with an argument (the name of the config file)
	if ( argc != 2 ) {
		printUsage( argv[0] );
		exit( 1 );
	}

	// Read the config file
	StarGalacticParam::readConfig( cfName );

	// To check whether all values are read correctly
	// one could write them to a new file and compare
//	StarGalacticParam::writeConfig( "./traceback2.conf" ); return 0;

	// Assign some values to shorter names for convenience
	int numThreads = simulParams.Threads;

	// Find out the number of processors/cores
	if ( numThreads == 0 ) {
		numThreads = QThread::idealThreadCount();

		if ( numThreads < 1 ) numThreads = 1;
	}

	QString logName( cfName );
	logName.truncate( logName.lastIndexOf( '.' ) );
	logName = QString( "%1.StarsOnly-%2.log" )
	          .arg( logName )
	          .arg( QDateTime::currentDateTime().toString( "hhmmss") );
	QFile logFile( logName );
	logFile.open( QIODevice::WriteOnly );
	QTextStream log( &logFile );

	Output *threadOut;

	// Create a list of threads
	QList<TraceThread *> traceThreads;

	const unsigned int numOrbits = simulParams.Orbits;
	const QString fName1 = StarGalacticParam::simStarToFname( simulParams.Star1 );
	const QString fName2 = StarGalacticParam::simStarToFname( simulParams.Star2 );
	const QString fNameAssoc = StarGalacticParam::simStarToFname( simulParams.Assoc );
	const StarGalacticParam::OrbitType oType = simulParams.Type == "Potential" ? StarGalacticParam::Numeric :
	                                           simulParams.Type == "Epicycle"  ? StarGalacticParam::Epicycle :
	                                                                             StarGalacticParam::Linear;
	Pool *pool;
	bool abort;
	QList<int> assocList = StarGalacticParam::simStarToLineList( simulParams.Assoc );
	if ( assocList.isEmpty() )
		assocList.append( 0 );

	bool someOutput, firstAssoc;
	foreach( int line1, StarGalacticParam::simStarToLineList( simulParams.Star1 ) ) {
		foreach( int line2, StarGalacticParam::simStarToLineList( simulParams.Star2 ) ) {
			// Create the pool of orbits
			pool = new Pool( numOrbits );
			// Create an object that serializes output from the threads and stores the results
			threadOut = new Output;
			abort = false;

			// Create the threads themselves on the heap
			// and populate them with the stars.
			// Note, that the threads do not run, yet,
			// and associations will be added later.
			for ( int i=0; i<numThreads; i++ ) {
				TraceThread *t = new TraceThread( i,
				                                  pool,
				                                  &abort,
				                                  *threadOut );
				t->addStar( fName1, line1, oType );
				t->addStar( fName2, line2, oType );
				traceThreads << t;
			}
			// Now add the associations and one output file
			// per association to the threads
			someOutput = false;
			foreach( int lineAssoc, assocList ) {
				// Get a name for the output file
				QString oName( simulParams.OutFile );
				if ( oName == "auto" ) {
					QString outName( cfName );

					outName.truncate( outName.lastIndexOf( '.' ) );
					oName = outName + ".out";
				}
				QString ext( oName.section( '.', -1, -1, QString::SectionIncludeLeadingSep ) );
				if ( ext == oName ) ext.clear();

				QString file( oName.left( oName.lastIndexOf( '.' ) ) );

				oName = QString( "%1.%2" ).arg( file ).arg( line1, 3, 10, QChar('0') );
				if ( line2 ) oName += QString( "+%1" ).arg( line2, 3, 10, QChar('0') );
				if ( lineAssoc ) oName += QString( "+%1" ).arg( lineAssoc, 3, 10, QChar('0') );
				oName += ext;

				// Initialize the file for serialized output
				// If initialization fails, continue with the next case
				if ( threadOut->initFile( lineAssoc, oName ) ) {
					someOutput = true;
					if ( lineAssoc ) {
						foreach( TraceThread *t, traceThreads )
							t->addAssoc( fNameAssoc, lineAssoc, oType );
					}
				}
			}
			if ( !someOutput ) {
				// We were unable to initialize any output file so we destroy
				// the threads and continue with the next star2
				while ( !traceThreads.isEmpty() )
					delete traceThreads.takeLast();
				continue;
			}
			// Write some text to the output files
			foreach( int lineAssoc, assocList ) {
				threadOut->writeToFile( lineAssoc,
				                        QString("# Orbits: %1\n")
				                       .arg(numOrbits) );
				threadOut->writeToFile( lineAssoc,
				                        QString( "# Number of steps = %L1\n"
				                                 "# Width of steps  = %L2 years\n"
				                                 "# Criterion       = '%3'\n")
				                        .arg( simulParams.Steps )
				                        .arg( simulParams.Width )
				                        .arg( simulParams.Criterion ) );
				threadOut->writeToFile( lineAssoc,
				                        QString("# Calculating %L1 '%2' orbits for the following stars and association:\n")
				                       .arg( numOrbits )
				                       .arg( simulParams.Type ) );
				if ( threadOut->isValidOut( lineAssoc ) )
					traceThreads.first()->printStarStartValues( lineAssoc, *threadOut->getOut( lineAssoc ), true );

				QString msg = QString( "# column 1 = %1separation between star 1 and 2 [pc]\n" )
				              .arg( simulParams.Criterion == CRITVARLERI1 ? "" : "minimum " );
				threadOut->writeToFile( lineAssoc, msg);
				msg =         QString( "# column 2 = time since %1separation [yr]\n" )
				              .arg( simulParams.Criterion == CRITVARLERI1 ? "" : "minimum " );
				threadOut->writeToFile( lineAssoc, msg);

				if ( lineAssoc ) {
					threadOut->writeToFile( lineAssoc, "# column 3 = distance of star 1 to center of association [pc]\n");
					threadOut->writeToFile( lineAssoc, "# column 4 = distance of star 2 to center of association [pc]\n");
					threadOut->writeToFile( lineAssoc, "#\n# The following columns depend on the CoordType of the stars above:\n");
					threadOut->writeToFile( lineAssoc,
					              QString::fromUtf8(   "# columns 5, 12, 19  = π [mas] of star 1, 2, and assoc, resp.\n") );
					threadOut->writeToFile( lineAssoc, "# columns 6, 13, 20 = RV  or U   [km/s]                         of star 1, 2, and assoc, resp.\n");
					threadOut->writeToFile( lineAssoc,
					              QString::fromUtf8(   "# columns 7, 14, 21 = µ_α* or µ_l* or µ_λ* [mas/yr] or V [km/s] of star 1, 2, and assoc, resp.\n") );
					threadOut->writeToFile( lineAssoc,
					              QString::fromUtf8(   "# columns 8, 15, 22 = µ_δ or µ_b or µ_β [mas/yr] or W [km/s]    of star 1, 2, and assoc, resp.\n") );
					threadOut->writeToFile( lineAssoc, "# columns 9, 16, 23 = Galactic longitude [degrees] at minimum distance of star 1, 2, and assoc, resp.\n" );
					threadOut->writeToFile( lineAssoc, "# columns 10, 17, 24 = Galactic latitude [degrees] at minimum distance of star 1, 2, and assoc, resp.\n" );
					threadOut->writeToFile( lineAssoc, "# columns 11, 18, 25 = Distance from Sun [pc] of star 1, 2, and assoc, resp.\n\n" );
				} else {
					threadOut->writeToFile( 0, "# column 3 = always 0\n");
					threadOut->writeToFile( 0, "# column 4 = always 0\n");
					threadOut->writeToFile( 0, "#\n# The following columns depend on the CoordType of the stars above:\n");
					threadOut->writeToFile( 0,
					        QString::fromUtf8( "# columns 5, 12  = π [mas] of star 1 and 2, resp.\n" ) );
					threadOut->writeToFile( 0, "# columns 6, 13 = RV  or U   [km/s]                         of star 1 and 2, resp.\n");
					threadOut->writeToFile( 0,
					        QString::fromUtf8( "# columns 7, 14 = µ_α* or µ_l* or µ_λ* [mas/yr] or V [km/s] of star 1 and 2, resp.\n") );
					threadOut->writeToFile( 0,
					        QString::fromUtf8( "# columns 8, 15 = µ_δ or µ_b or µ_β [mas/yr] or W [km/s]    of star 1 and 2, resp.\n#\n") );
					threadOut->writeToFile( 0, "# columns 9, 16 = Galactic longitude [degrees] at minimum distance of star 1, 2, resp.\n");
					threadOut->writeToFile( 0, "# columns 10, 17 = Galactic latitude [degrees] at minimum distance of star 1 and 2, resp.\n");
					threadOut->writeToFile( 0, "# columns 11, 18 = Distance from Sun [pc] of star 1 and 2, resp.\n\n");
				}
				threadOut->writeToFile( lineAssoc, "\n" );
			}

#if ( QT_VERSION < QT_VERSION_CHECK(5,14,0) )
# define ENDL endl
#else
# define ENDL Qt::endl
#endif

			// Start timer and threads
			out << QString( "%1: starting %2 threads for %L3 '%4' orbits:" )
			       .arg( threadOut->getFileName( assocList[0] ) )
			       .arg( numThreads )
			       .arg( numOrbits )
			       .arg( simulParams.Type ) << ENDL << ENDL;

			timer.start();
			foreach ( TraceThread *t, traceThreads )
				t->start();		// Now the threads actually run

			// Wait until all threads have finished and memorize elapsed time
			foreach ( TraceThread *t, traceThreads )
				t->wait();

			duration = timer.elapsed();

			// Delete/destroy the threads
			while ( !traceThreads.isEmpty() )
				delete traceThreads.takeLast();

			QString outLine;
			out << ENDL;

			bool deleteFile;
			firstAssoc = true;
			foreach( int lineAssoc, assocList ) {
				deleteFile = false;
				// If there is no output file for this association we just jump to the next,
				// this will also keep the noise away from the terminal
				if ( !threadOut->isValidOut( lineAssoc ) ) continue;

				threadOut->writeToFile( lineAssoc, "\n" );
				QString logEntry( threadOut->getFileName( lineAssoc ) );
				logEntry.truncate( logEntry.lastIndexOf( '+' ) );

				// Get the results from the output object and write them to the file
				if ( simulParams.Criterion == CRITTETZLAFF ) {
					outLine = QString::fromUtf8( "In %L1 orbits both stars will get closer than 15 pc to the center of the association #%2 simultaneously.\n" )
					          .arg( threadOut->getCloseToAssoc( lineAssoc ) )
					          .arg( lineAssoc );
					out << outLine;
					threadOut->writeToFile( lineAssoc, "# " + outLine );

					outLine = QString::fromUtf8( "The minimum distance between the stars was %L1 pc %L2 years ago.\n" )
					          .arg( threadOut->getMinDist() )
					          .arg( -threadOut->getTimePoint() );
					out << outLine;
					threadOut->writeToFile( lineAssoc, "# " + outLine );
					if ( threadOut->getCloseToAssoc( lineAssoc ) == 0 )
						deleteFile = true;
				}

				if ( simulParams.Criterion == CRITHOOGERWERF ) {
					outLine = QString::fromUtf8( "In %L1 of %L2 orbits the stars get closer than %3 pc to each other.\n" )
					          .arg( threadOut->getHowMany() )
					          .arg( simulParams.Orbits )
					          .arg( simulParams.Limit );
					out << outLine;
					threadOut->writeToFile( lineAssoc, "# " + outLine );

					outLine = QString::fromUtf8( "The minimum separation of all orbits was %L1 pc %L2 years ago.\n" )
					          .arg( threadOut->getMinDist() )
					          .arg( -threadOut->getTimePoint() );
					out << outLine;
					threadOut->writeToFile( lineAssoc, "# " + outLine );

					const Vector<mytype> sep ( threadOut->getMeanSeparation() );
					const Vector<mytype> sigma( threadOut->getSeparationSigma() );

					outLine = QString::fromUtf8( "The mean minimum separation vector of the two stars (X,Y,Z) was:\n" );
					out << outLine;
					threadOut->writeToFile( lineAssoc, "# " + outLine );
					outLine = QString::fromUtf8( "(%L1 ± %L2, %L3 ± %L4, %L5 ± %L6) pc.\n")
					          .arg( sep[0],   0, 'f', 2 )
					          .arg( sigma[0], 0, 'f', 2 )
					          .arg( sep[1],   0, 'f', 2 )
					          .arg( sigma[1], 0, 'f', 2 )
					          .arg( sep[2],   0, 'f', 2 )
					          .arg( sigma[2], 0, 'f', 2 );
					out << outLine;
					threadOut->writeToFile( lineAssoc, "# " + outLine );

					outLine = QString::fromUtf8( "The mean minimum separation was %L1 pc (%L2 ± %L3) Myr ago.\n" )
					          .arg( norm2( sep ) )
					          .arg( -threadOut->getMeanMinTimes()*1e-6, 0, 'f', 3 )
					          .arg( threadOut->getMinTimesSigma()*1e-6, 0, 'f', 3 );
					out << outLine;
					threadOut->writeToFile( lineAssoc, "# " + outLine );

					if ( lineAssoc ) {
						outLine = QString::fromUtf8( "In %L1 of those orbits both stars were also closer than 10 pc to the center of the association #%2.\n" )
						          .arg( threadOut->getCloseToAssoc( lineAssoc ) )
						          .arg( lineAssoc );
						out << outLine;
						threadOut->writeToFile(lineAssoc, "# " + outLine );
					}
					if ( threadOut->getHowMany() == 0 )
						deleteFile = true;
				}

				if ( simulParams.Criterion == CRITRALPH1 ) {
					const mytype radius = StarGalacticParam( fNameAssoc, lineAssoc, oType ).getRadius();

					outLine = QString::fromUtf8( "In %L1 of %L2 orbits the stars get closer than %3 pc to each other.\n" )
					          .arg( threadOut->getHowMany() )
					          .arg( simulParams.Orbits )
					          .arg( simulParams.Limit );
					out << outLine;
					if ( firstAssoc && threadOut->getHowMany() ) {
						log << logEntry << ": " << outLine;
					}
					threadOut->writeToFile( lineAssoc, "# " + outLine );

					outLine = QString::fromUtf8( "The minimum distance was %L1 pc %L2 years ago.\n" )
					          .arg( threadOut->getMinDist() )
					          .arg( -threadOut->getTimePoint() );
					out << outLine;
					if ( firstAssoc && threadOut->getHowMany() ) {
						log << logEntry << ": " << outLine;
					}
					threadOut->writeToFile( lineAssoc, "# " + outLine );

					const Vector<mytype> sep ( threadOut->getMeanSeparation() );
					const Vector<mytype> sigma( threadOut->getSeparationSigma() );

					outLine = QString::fromUtf8( "The mean minimum separation vector of the two stars (X,Y,Z) was:\n" );
					out << outLine;
					if ( firstAssoc && threadOut->getHowMany() ) {
						log << logEntry << ": " << outLine;
					}
					threadOut->writeToFile( lineAssoc, "# " + outLine );
					outLine = QString::fromUtf8( "(%L1 ± %L2, %L3 ± %L4, %L5 ± %L6) pc.\n")
					          .arg( sep[0],   0, 'f', 2 )
					          .arg( sigma[0], 0, 'f', 2 )
					          .arg( sep[1],   0, 'f', 2 )
					          .arg( sigma[1], 0, 'f', 2 )
					          .arg( sep[2],   0, 'f', 2 )
					          .arg( sigma[2], 0, 'f', 2 );
					out << outLine;
					if ( firstAssoc && threadOut->getHowMany() ) {
						log << logEntry << ": " << outLine;
					}
					threadOut->writeToFile( lineAssoc, "# " + outLine );

					outLine = QString::fromUtf8( "The mean minimum separation was %L1 pc (%L2 ± %L3) Myr ago.\n" )
					          .arg( norm2( sep ) )
					          .arg( -threadOut->getMeanMinTimes()*1e-6, 0, 'f', 3 )
					          .arg( threadOut->getMinTimesSigma()*1e-6, 0, 'f', 3 );
					out << outLine;
					if ( firstAssoc && threadOut->getHowMany() ) {
						log << logEntry << ": " << outLine;
					}
					threadOut->writeToFile( lineAssoc, "# " + outLine );

					outLine = QString::fromUtf8( "In %L1 of those orbits both stars were also inside the association #%2.\n" )
					          .arg( threadOut->getCloseToAssoc( lineAssoc ) )
					          .arg( lineAssoc );
					out << outLine;
					threadOut->writeToFile( lineAssoc, "# " + outLine );

					outLine = QString::fromUtf8( "(Less than %L1 pc from the center)\n" )
					          .arg( radius );
					out << outLine;
					threadOut->writeToFile( lineAssoc, "# " + outLine );
					if ( threadOut->getCloseToAssoc( lineAssoc ) == 0 )
						deleteFile = true;
					if ( firstAssoc && !threadOut->getHowMany() ) {
						log << logEntry << ": no matches" << ENDL;
					}
					if ( firstAssoc )
						log << ENDL;
				}

				if ( simulParams.Criterion == CRITVARLERI1 ) {
					const mytype radius = StarGalacticParam( fNameAssoc, lineAssoc, oType ).getRadius();

					outLine = QString::fromUtf8( "In %L1 time points of orbits the stars get closer than %2 pc to each other.\n" )
					          .arg( threadOut->getHowMany() )
					          .arg( simulParams.Limit );
					out << outLine;
					threadOut->writeToFile( lineAssoc, "# " + outLine );

					outLine = QString::fromUtf8( "The minimum distance of the two stars was %L1 pc %L2 years ago, regardless of the association.\n" )
					          .arg( threadOut->getMinDist() )
					          .arg( -threadOut->getTimePoint() );
					out << outLine;
					threadOut->writeToFile( lineAssoc, "# " + outLine );

					outLine = QString::fromUtf8( "In %L1 of those time points of orbits (%L2 orbits) both stars were also inside the association #%3.\n" )
					          .arg( threadOut->getCloseToAssoc( lineAssoc ) )
					          .arg( threadOut->getEncounterOrbits( lineAssoc ) )
					          .arg( lineAssoc );
					out << outLine;
					threadOut->writeToFile( lineAssoc, "# " + outLine );

					outLine = QString::fromUtf8( "(Less than %L1 pc from the center)\n" )
					          .arg( radius );
					out << outLine;
					threadOut->writeToFile( lineAssoc, "# " + outLine );
					if ( threadOut->getHowMany() == 0 )
						deleteFile = true;
				}

				outLine = QString::fromUtf8( "Computational time: %L1 s\n" )
				          .arg(duration/1000.,0,'f',1);
				out << outLine << ENDL;
				threadOut->writeToFile( lineAssoc, "# " + outLine );

				if ( deleteFile && simulParams.RemoveFile ) {
					outLine = QString::fromUtf8( "%1 deleted\n\n" )
					          .arg( threadOut->getFileName( lineAssoc ) );
					out << outLine;
					threadOut->deleteFile( lineAssoc );
				}
				firstAssoc = false;
			}	// next Assoc
			// Spit out a bell character which, if configured, will ring
			// a beep or flash the terminal window
			out << "Delete this line to remove the beep ->\a<-" << ENDL << ENDL;
			delete threadOut;
			delete pool;
		}		// next Star2
	}			// next Star1

	// That's it!
	if ( logFile.isOpen() )
		logFile.close();

	if ( simulParams.Criterion != CRITRALPH1 )
		logFile.remove();

	return 0;
}

