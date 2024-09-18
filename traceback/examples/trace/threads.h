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
// $Id: threads.h 143 2022-07-28 10:52:27Z ifg $


#ifndef TRACETHREAD_H
#define TRACETHREAD_H

#include <QThread>
#include <QMutex>
#include <QTextStream>
#include <QFile>
#include <QString>
#include <QMap>
#include <float.h>			// DBL_MAX
#include "starparam.h"
#include "rng.h"
#include "Vector.h"

/*!
 * \brief A class where threads will store their results
 *
 * At program startup a single instance should be created that
 * will be used by all threads concurrently. The instance is meant to be
 * used once only, do not re-use it with other pairs of stars and
 * associations!
 */
class Output
{
private:
	QMutex mutexParam;					///< Protect private attributes of class
	QMutex mutexFile;					///< Serialize write access of output files
	int timePoint;						///< When #minDist occurred
	int howMany;						///< Count how often two objects get close to each other
	mytype minDist;						///< Minimum separation of all pairs of orbits (two objects)
	Vector<mytype> *sumSeparations;		///< Sum of all separation vectors between two objects
	Vector<mytype> *sumSepsSquared;		///< Sum of all squared (component-wise) separations between two objects
	mytype sumMinimumTimes;				///< Sum of all times since minimum separations
	mytype sumMinTimesSquared;			///< Sum of all times since minimum separations squared
	QMap<int, int> closeToAssoc;		///< Mapping between line number of assoc and number of close encounters of stars/assoc
	QMap<int, int> encounterOrbits;		///< Mapping between assoc and number of orbits with close encounters
	QMap<int, QTextStream *> out;		///< Mapping between line number of assoc and output stream

public:
	/*!
	 * \brief Constructor. Some counting variables are initialized here,
	 * so do not re-use the instance when done.
	 */
	Output()
	    : timePoint( 0 ),
	      howMany( 0 ),
	      minDist( DBL_MAX ),
	      sumSeparations( new Vector<mytype>(3, 0.0) ),
	      sumSepsSquared( new Vector<mytype>(3, 0.0) ),
	      sumMinimumTimes( 0.0 ),
	      sumMinTimesSquared( 0.0 )
	{}
	/*!
	 * \brief Destructor. All files are closed.
	 */
	~Output()
	{
		foreach ( QTextStream *out, out ) {
			delete out->device();
			delete out;
		}
		delete sumSepsSquared;
		delete sumSeparations;
	}
	/*!
	 * \brief Raise the number of close encounters of two objects
	 * \param relPosition Difference vector of the two objects
	 * \param when Time of \a relPosition
	 *
	 * This will also update #sumSeparations and #sumSepsSquared which are
	 * needed to get a mean separation %vector and standard deviations
	 * in \a X, \a Y, and \a Z, as well as #sumMinimumTimes an
	 * #sumMinTimesSquared.
	 *
	 * \sa getMeanSeparation()
	 * \sa getSeparationSigma()
	 */
	void incHowMany(Vector<mytype> relPosition, mytype when) {

		mutexParam.lock();
		*sumSeparations += relPosition;
		*sumSepsSquared += (relPosition*relPosition);
		sumMinimumTimes += when;
		sumMinTimesSquared += (when*when);
		howMany++;
		mutexParam.unlock();
	}
	/*!
	 * \brief Raise the number of close encounters of three objects
	 * \param assoc Identifier of the association (line number in input file)
	 */
	void incCloseToAssoc(int assoc) {
		mutexParam.lock();
		closeToAssoc[assoc]++;
		mutexParam.unlock();
	}

	void incEncounterOrbits(int assoc) {
		mutexParam.lock();
		encounterOrbits[assoc]++;
		mutexParam.unlock();
	}

	mytype getMinDist() const {return minDist;}
	void updateMinDist(mytype dist, int when) {
		mutexParam.lock();
		if ( dist < minDist ) {
			minDist = dist;
			timePoint = when;
		}
		mutexParam.unlock();
	}
	int getHowMany() const {return howMany;}
	int getCloseToAssoc(int assoc) const {return closeToAssoc.value(assoc);}
	int getEncounterOrbits(int assoc) const {return encounterOrbits.value(assoc);}
	int getTimePoint() const {return timePoint;}
	/*!
	 * \brief Initialize output per association
	 * \param assoc Identifier for association. If 0, no association is used (2 objects only). That's the idea.
	 * \param fName Name of the output file
	 * \return Success
	 */
	bool initFile( int assoc, QString fName )
	{
		if ( QFile::exists( fName ) ) {
			qDebug( "Error: File \"%s\" exists, will not overwrite!\n",
			        fName.toLatin1().constData() );
			return false;
		}
		QFile *outFile = new QFile( fName );
		if ( !outFile ) return false;
		bool status = outFile->open( QIODevice::WriteOnly );
		QTextStream *outStream = new QTextStream( outFile );
		outStream->setCodec("UTF-8");
		out.insert( assoc, outStream );
		closeToAssoc.insert( assoc, 0 );
		encounterOrbits.insert( assoc, 0 );
		return status;
	}
	void writeToFile(int assoc, QString string)
	{
		if ( !out.contains(assoc) ) return;
		mutexFile.lock();
		*out.value(assoc) << string;
		mutexFile.unlock();
	}

#if ( QT_VERSION < QT_VERSION_CHECK(5,14,0) )
# define ENDL endl
#else
# define ENDL Qt::endl
#endif

	void updateFile(int assoc, mytype dist, int time, mytype dist1, mytype dist2, QString remarks = QString())
	{
		if ( !out.contains(assoc) )
			return;
		else {
			QTextStream *out = this->out.value(assoc);
			mutexFile.lock();
			*out << QString( "%1 %2 %3 %4%5" )
			        .arg( dist, 8, 'f', 4 )
			        .arg( time, 9 )
			        .arg( dist1, 8, 'f', 4 )
			        .arg( dist2, 8, 'f', 4 )
			        .arg( remarks );
			*out << ENDL;
			mutexFile.unlock();
		}
	}
	void deleteFile(int assoc)
	{
		if ( !out.contains(assoc) ) return;
		mutexFile.lock();
		QIODevice *qDev = out.value(  assoc )->device();
		static_cast<QFile *>( qDev )->remove();
		delete qDev;
		out.remove( assoc );
		mutexFile.unlock();
	}
	QString getFileName( int assoc )
	{
		if ( !out.contains(assoc) ) return( "Error" );

		QIODevice *qDev = out.value( assoc )->device();
		return static_cast<QFile *>( qDev )->fileName();
	}
	// FIXME: Using 'out' is not thread safe
	QTextStream *getOut( int assoc )
	{
		if ( !out.contains( assoc ) )
			return NULL;
		return out.value(assoc);
	}
	bool isValidOut( int assoc )
	{
		return out.contains( assoc );
	}
	/*!
	 * \brief Get the mean minimum separation %vector of all "successful orbits"
	 * \param ok Pointer to a success indicator that will be set if present
	 * \return Vector with components in parsec.
	 */
	Vector<mytype> getMeanSeparation(bool *ok = NULL)
	{
		if ( !howMany ) {
			Vector<mytype> rc( 3, 0.0 );
			if ( ok ) *ok = false;
			return rc;
		}
		if ( ok ) *ok = true;

		return *sumSeparations/howMany;
	}
	/*!
	 * \brief Get the standard deviations in \a X, \a Y, and \a Z
	 *        of the minimum separation %vetors of all "successful orbits"
	 * \param ok Pointer to a success indicator that will be set if present
	 * \return Vector with \f$(\sigma_x, \sigma_y, \sigma_z)\f$ in parsec.
	 *
	 * The standard deviation is calculated according to the following formula:
	 *
	 * \f$
	 * \sigma_j = \displaystyle \sqrt{\frac{1}{n-1}\left[\left( \sum_{i=1}^{n}j_i^2 \right)
	 *            - \frac{1}{n} \left( \sum_{i=1}^{n}j_i \right)^2 \right]}
	 *            ~~, ~~~j = \{ x,y,z \}
	 * \f$
	 *
	 * No special care has been taken to avoid rounding erros, yet.
	 *
	 * \sa <a href="https://de.wikipedia.org/wiki/Stichprobenvarianz_(Sch%C3%A4tzfunktion)#Stichprobenstandardabweichung">
	 *     https://de.wikipedia.org/wiki/Stichprobenvarianz_(Schätzfunktion)#Stichprobenstandardabweichung</a>
	 */
	Vector<mytype> getSeparationSigma(bool *ok = NULL)
	{
		Vector<mytype> result( *sumSeparations );

		if ( howMany<2 ) {
			Vector<mytype> rc(3, 0.0);
			if ( ok ) *ok = false;
			return rc;
		}
		result *= result;
		result /= howMany;
		result = *sumSepsSquared - result;
		result /= (howMany-1);
		for ( int i=0; i<result.length(); i++) {
			assert(result[i]>=0);
			result[i] = sqrt( result[i] );
		}
		if ( ok ) *ok = true;

		return result;
	}
	mytype getMeanMinTimes(bool *ok = NULL)
	{
		if ( !howMany ) {
			if ( ok ) *ok = false;
			return nan("");
		}
		if ( ok ) *ok = true;

		return sumMinimumTimes/howMany;
	}
	mytype getMinTimesSigma(bool *ok = NULL)
	{
		mytype rc = sumMinimumTimes*sumMinimumTimes;

		if ( howMany<2 ) {
			if ( ok ) *ok = false;
			return nan("");
		}
		rc /= howMany;
		rc = sumMinTimesSquared - rc;
		rc /= (howMany-1);
		assert( rc>=0 );
		rc = sqrt( rc );
		if ( ok ) *ok = true;

		return rc;
	}
};

/*!
 * \brief Class to distribute the orbits over the threads. It will count
 * backwards starting from limit. It is supposed to be thread-safe.
 */
class Pool
{
private:
	const unsigned int limit;
	unsigned int current;
	QMutex mx;

public:
	/*!
	 * \brief Constructor
	 * \param variations Number of overall orbits
	 */
	Pool(const unsigned int variations)
	    : limit( variations ), current ( variations )
	{}
	/*!
	 * \brief Check whether the pool has been drained.
	 * Calling this function will consume an orbit.
	 * \return the number of available integers (orbits)
	 */
	unsigned int available()
	{
		unsigned int rc;
		mx.lock();
		rc = current;
		if ( rc ) current--;
		mx.unlock();

		return rc;
	}
	unsigned int getLimit() const { return limit; }
};

/*!
 * \brief Class that does calculations in a thread
 *
 * The number of threads that are used can be set by the
 * configuration setting Simulation/Threads
 */
class TraceThread : public QThread
{
private:
	const int ID;							///< Identification number of thread
	Rng rng;								///< Random number generator for the Gaussian deviates
	Output &output;							///< Where to store the output
	QList<StarGalacticParam *> stars;		///< List of stars the orbits of which is calculated
	QList<StarGalacticParam *> assocs;		///< List of associations the stars are checked against
	Pool *pool;								///< A guardian that keeps an eye over the number of orbits
	int vars;								///< Number of orbits that have been calculated
	bool *abortAll;							///< Whether the calculation should be aborted
	QString getSuccessful(StarGalacticParam *assoc) const;

protected:
	virtual void run();

public:
	TraceThread(const int threadId, Pool *cnt, bool *abort, Output &out);
	~TraceThread();
	int getNumStars() const { return stars.size(); }
	int getNumAssocs() const { return assocs.size(); }
	StarGalacticParam *getStar(int index);
	StarGalacticParam *getAssoc(int index);
	void addStar(QString fname, int lineNumber, StarGalacticParam::OrbitType type);
	void addAssoc(QString fname, int lineNumber, StarGalacticParam::OrbitType type);
	void printStarStartValues(int lineAssoc, QTextStream &out, bool prependHash);
	void printCoordinates(QTextStream &out, bool heliocentric) const;
	void printStarStartValues(QTextStream &out, bool prependHash);
};

#endif // TRACETHREAD_H
