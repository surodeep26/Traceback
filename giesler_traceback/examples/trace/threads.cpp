/***************************************************************************
 *   Copyright (C) 2012 by Frank Giessler                                  *
 *   ifg@astro.uni-jena.de                                                 *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
// $Id: threads.cpp 143 2022-07-28 10:52:27Z ifg $

#include "threads.h"

extern SimulParams simulParams;

/*!
 * \brief Return the parameters of the 'successful' orbits
 * \param assoc Pointer to the assoiation with which the orbits were successful
 * \return Formatted string
 *
 * Since the parameters are stored as unions, the returned values
 * are as follows:
 *
 * \li \a para [mas]
 * \li \a RV or \a U [km/s]
 * \li \a PMa or \a PMl [mas/yr] or \a V [km/s]
 * \li \a PMd or \a PMb [mas/yr] or \a W [km/s]
 * \li \a Galactic \a longitude at closest encounter [°]
 * \li \a Galactic \a latitude at closest encounter [°]
 * \li Distance from Sun at closest encounter [pc]
 *
 * Galactic coordinates are relative to the rotating heliocentric system..
 */
QString TraceThread::getSuccessful(StarGalacticParam *assoc) const
{
	QString remark;
	Vector<mytype> coords(3);

	QList<StarGalacticParam *> list(stars);
	if ( assoc )
		list.append(assoc);

	foreach ( StarGalacticParam *star, list ) {
//		raDec = star->getEquatorial();
		coords = star->getGalactic();
		remark += QString( "   %1 %2 %3 %4 %5 %6 %7" )
		          .arg( star->getInitialPara(), 11, 'f', 5 )	// para [mas]
		          .arg( star->getInitialU(), 10, 'f', 5 )		// RV or U [km/s]
		          .arg( star->getInitialV(), 10, 'f', 5 )		// PMa or PMl [mas/yr] or V [km/s]
		          .arg( star->getInitialW(), 10, 'f', 5 )		// PMd or PMb [mas/yr] or W [km/s]
		          .arg( coords[0], 10, 'f', 5 )					// l in [degrees]
		          .arg( coords[1], 10, 'f', 5 )					// b [degrees]
		          .arg( coords[2], 10, 'f', 5 );				// Dist [pc]
	}
	return remark;
}

/*!
 * \brief The function that is called when the thread is started
 *
 * It does the calculation that was previously in main.cpp
 */
void TraceThread::run()
{
	QTextStream out( stdout );

	foreach ( StarGalacticParam *star, stars + assocs ) {
		if ( star->getOrbitType() == StarGalacticParam::Numeric )
			star->zeroSigmas();
	}

//	number of steps for each orbit
	const int steps = simulParams.Steps;

//	width of steps in years
	const int width = simulParams.Width;

//	distance for criteria
	const mytype limit = simulParams.Limit;

//	things to remember
	mytype dist;
	mytype mindist;
	int timeMinDist;

//	Iterate over the number of orbits
	unsigned int variation;

	bool orbitFound;

	while ( (variation = pool->available()) > 0 ) {
		// Get random values for starting parameters for the two stars and all associations
		foreach ( StarGalacticParam *star, stars + assocs ) {
			// Processing the original values is a bad idea because of the unknown RV of a neutron star
			if ( !star->varyAllStartValues( rng ) || *abortAll ) {
				QString msg = QString( *abortAll ? "aborting on request" : "parameter variation failed" );
				*abortAll = true;
				QList<int> fileList;
				if ( !assocs.isEmpty() ) {
					foreach( StarGalacticParam *assoc, assocs )
						fileList << assoc->getLineNumber();
				} else
					fileList << 0;

				foreach( int lineAssoc, fileList ) {
					output.writeToFile( lineAssoc,
					                    QString("Thread #%1: %2, finished only %L3 orbits\n")
					                    .arg(ID, 2, 10, QLatin1Char('0'))
					                    .arg(msg)
					                    .arg(vars) );
				}

				return;
			}
			star->coordsFromStartValues( false );
		}

		mindist = dist = *stars[0] || *stars[1];
		timeMinDist = 0;

		// start the new orbit by first resetting the ODE driver in the GSL
		foreach ( StarGalacticParam *star, stars + assocs )
			star->resetODE( simulParams.StepSize );

		// now step through the orbits of stars 1 & 2 over time
		orbitFound = false;
		for ( int j=1; j<=steps; j++ ) {
			int t = -width*j;
			for ( int i=0; i<2; i++ ) {
				stars[i]->progressOrbit( -width );
			}
			// after each step calculate the separation between the stars
			// and remember if they are getting closer
			dist =  *stars[0] || *stars[1];
			if ( dist < mindist ) {
				mindist = dist;
				timeMinDist = t;
			}
			if ( simulParams.Criterion == CRITVARLERI1 ) {
				if ( assocs.isEmpty() ) qFatal( "Criterion \"Valeri1\" needs three objects!\n" );
				/*********************
				Auswahlkriterium nach VVH 10.01.2019:
				Die Distanz zwischen den Sternen (nicht minimale!) ist kleiner
				1 pc und beide befinden sich zu diesem Zeitpunkt innerhalb
				der Assoziation. Jeder solcher Zeitpunkt wird gezählt. Die
				Orbits, bei denen so etwas auftritt, werden auch gezählt.
				*/
				mytype dist0, dist1, radius;
				int lineAssoc;
				if ( dist < limit ) {
					output.incHowMany( stars[1]->coordsInertial()[0] - stars[0]->coordsInertial()[0], (mytype) t );
					foreach ( StarGalacticParam *assoc, assocs ) {
						assoc->resetODE( simulParams.StepSize );
						assoc->orbitAt( t );
						dist0 = *stars[0] || *assoc;
						dist1 = *stars[1] || *assoc;
						lineAssoc = assoc->getLineNumber();
						radius = assoc->getRadius();
						if ( dist0 < radius && dist1 < radius ) {
							output.updateFile( lineAssoc, dist, t, dist0, dist1, getSuccessful( assoc ) );
							output.incCloseToAssoc( lineAssoc );
							if ( !orbitFound ) {
								output.incEncounterOrbits( lineAssoc );
								orbitFound = true;
							}
						}
					}
				}
			} // end of Valeri1 criterion
		} // next time point for orbits

		// remember the minimum separation of all pairs of orbits
		if ( mindist < output.getMinDist() ) {
			output.updateMinDist( mindist, timeMinDist );
		}

		// decide whether the orbit was a success
		if ( simulParams.Criterion == CRITTETZLAFF ) {
			/**********************
			Auswahlkriterium nach Tetzlaff et al.:
			Zum Zeitpunkt der minimalen Distanz befinden sich beide Sterne höchstens
			15 pc von der Assoziation entfernt.
			*/
			foreach ( StarGalacticParam *star, stars + assocs ) {
				star->resetODE( simulParams.StepSize );
				star->orbitAt( timeMinDist );
			}

			mytype dist0, dist1;
			int lineAssoc;
			foreach( StarGalacticParam *assoc, assocs ) {
				dist0 = *stars[0] || *assoc;
				dist1 = *stars[1] || *assoc;
				lineAssoc = assoc->getLineNumber();
				// TODO: replace with limit
				if ( dist0 <= 15.0 && dist1 <= 15.0 ) {
					output.updateFile( lineAssoc, mindist, timeMinDist, dist0, dist1, getSuccessful(assoc) );
					output.incHowMany( stars[1]->coordsInertial()[0] - stars[0]->coordsInertial()[0], timeMinDist );
					output.incCloseToAssoc( lineAssoc );
				}
			}
		} // Tetzlaff criterion

		if ( simulParams.Criterion == CRITHOOGERWERF ) {
			/**********************
			Auswahlkriterium nach Hoogerwerf et al.:
			Die minimale Distanz ist kleiner als 10 pc und beide Sterne sind
			weniger als 10 pc von der Assoziation entfernt
			*/
			if ( mindist < limit ) {
				foreach ( StarGalacticParam *star, stars + assocs ) {
					star->resetODE( simulParams.StepSize );
					star->orbitAt( timeMinDist );
				}
				output.incHowMany( stars[1]->coordsInertial()[0] - stars[0]->coordsInertial()[0], timeMinDist );
				if ( !assocs.isEmpty() ) {
					mytype dist0, dist1;
					int lineAssoc;
					foreach( StarGalacticParam *assoc, assocs ) {
						dist0 = *stars[0] || *assoc;
						dist1 = *stars[1] || *assoc;
						lineAssoc = assoc->getLineNumber();
						if ( dist0 < 10.0 && dist1 < 10.0 ) {
							output.updateFile( lineAssoc, mindist, timeMinDist, dist0, dist1, getSuccessful( assoc ) );
							output.incCloseToAssoc( lineAssoc );
						}
					}
				} else {
					output.updateFile( 0, mindist, timeMinDist, 0, 0, getSuccessful( NULL ) );
				}
			}
		} // Hoogerwerf criterion

		if ( simulParams.Criterion == CRITRALPH1 ) {
			/**********************
			Auswahlkriterium nach RNE 23.11.2018:
			Die minimale Distanz ist kleiner als 10 pc und beide Sterne
			befinden sich innerhalb der Assoziation
			*/
			if ( assocs.isEmpty() ) qFatal( "Criterion \"Ralph1\" needs three objects!\n" );

			if ( mindist <= limit ) {
				foreach ( StarGalacticParam *star, stars + assocs ) {
					star->resetODE( simulParams.StepSize );
					star->orbitAt( timeMinDist );
				}
				output.incHowMany( stars[1]->coordsInertial()[0] - stars[0]->coordsInertial()[0], timeMinDist );

				mytype dist0, dist1, radius;
				int lineAssoc;
				foreach( StarGalacticParam *assoc, assocs ) {
					dist0 = *stars[0] || *assoc;
					dist1 = *stars[1] || *assoc;
					lineAssoc = assoc->getLineNumber();
					radius = assoc->getRadius();
					if ( dist0 < radius && dist1 < radius ) {
						output.updateFile( lineAssoc, mindist, timeMinDist, dist0, dist1, getSuccessful( assoc ) );
						output.incCloseToAssoc( lineAssoc );
					}
				}
			}
		} // Ralph1 criterion
		vars++;
	} // next orbit of stars 1 & 2

#if ( QT_VERSION < QT_VERSION_CHECK(5,14,0) )
# define ENDL endl
#else
# define ENDL Qt::endl
#endif

	out << QString( "Thread #%1: finished %L2 orbits." )
	       .arg( ID, 2, 10, QLatin1Char('0') )
	       .arg( vars ) << ENDL;
}

/*!
 * \brief Create a TraceThread instance
 * \param threadId Identifier for the thread
 * \param cnt Pointer to a pool of integers
 * \param abort Pointer to global variable that can
 * be set from other threads to abort calculation
 * \param out Where to store the results
 */
TraceThread::TraceThread(const int threadId, Pool *cnt, bool *abort, Output &out)
    : ID( threadId ),
      output( out ),
      pool( cnt ),
      vars( 0 ),
      abortAll( abort )
{
}

TraceThread::~TraceThread()
{
	while ( !stars.isEmpty() )
		delete stars.takeLast();

	while ( !assocs.isEmpty() )
		delete assocs.takeLast();
}

StarGalacticParam *TraceThread::getStar(int index)
{
	return index < getNumStars() ? stars.at( index ) : NULL;
}

StarGalacticParam *TraceThread::getAssoc(int index)
{
	return index < getNumAssocs() ? assocs.at( index ) : NULL;
}

/*!
 * \brief Add a star for which the orbit is calculated
 * \param fname File with the parameters of the star
 * \param lineNumber Line number in \a fname
 * \param type orbit calculation method to use
 */
void TraceThread::addStar(QString fname, int lineNumber, StarGalacticParam::OrbitType type)
{
	stars.append( new StarGalacticParam(fname, lineNumber, type) );
}

/*!
 * \brief Add an association against which the stars are checked
 * \param fname File with the parameters of the association
 * \param lineNumber Line number in \a fname
 * \param type orbit calculation method to use
 */
void TraceThread::addAssoc(QString fname, int lineNumber, StarGalacticParam::OrbitType type)
{
	assocs.append( new StarGalacticParam(fname, lineNumber, type) );
}

/*!
 * \brief Print the parameters of the stars and current association read
 * \param lineAssoc Identifier of the association
 * \param out Where the output goes
 * \param prependHash Whether a hash sign ('\#') should be put in front of the lines
 */
void TraceThread::printStarStartValues(int lineAssoc, QTextStream &out, bool prependHash)
{
	QList<StarGalacticParam *>list( stars );
	StarGalacticParam *assoc( NULL );

	if ( lineAssoc && !assocs.isEmpty() ) {
		assoc = new StarGalacticParam( assocs[0]->getInputFileName(), lineAssoc, assocs[0]->getOrbitType() );
		list.append( assoc );
	}
	foreach ( StarGalacticParam *star, list ) {
		star->printStartValues( out, prependHash );
		out << ENDL;
	}
	out.flush();
	if ( assoc )
		delete assoc;
}

/*!
 * \brief Print the parameters of the stars and all associations read
 * \param out Where the output goes
 * \param prependHash Whether a hash sign ('\#') should be put in front of the lines
 */
void TraceThread::printStarStartValues(QTextStream &out, bool prependHash)
{
	foreach ( StarGalacticParam *star, stars + assocs ) {
		star->printStartValues( out, prependHash );
		out << ENDL;
	}
	out.flush();
}

/*!
 * \brief Print information on all stars and associations
 * \param out Where output goes
 * \param heliocentric Whether coordinates are true Galactic ones
 */
void TraceThread::printCoordinates(QTextStream &out, bool heliocentric) const
{
	foreach ( StarGalacticParam *star, stars + assocs ) {
		star->printInfo( out );
		star->printEquatorial( out, heliocentric );
		star->printGalactic( out, heliocentric );
		star->printXYZ( out, heliocentric );
		star->printUVW( out, heliocentric );
//		star->printSigmas( out );
		out << ENDL;
	}
}
