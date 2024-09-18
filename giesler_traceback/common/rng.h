/***************************************************************************
 *   Copyright (C) 2013 by Frank Giessler                                  *
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
// $Id: rng.h 123 2019-02-07 12:15:28Z ifg $


#ifndef RNG_H
#define RNG_H

#define HAVE_INLINE

#include <gsl/gsl_randist.h>
#include <time.h>
#include <math.h>

#define RNG_TYPE gsl_rng_mt19937

/*!
 * \brief The random number generator for the simulation
 *
 * The class is just a wrapper for a random number generator
 * from the
 * <a href="http://www.gnu.org/software/gsl/">
 * GNU scientific library (GSL)</a>.
 * It is probably not thread safe and it is too
 * time consuming to use mutexes, so one should have
 * one instance per thread.
 * The class is used to produce numbers from
 * either a uniform distribution between 0 and 1 and/or, more
 * important, from a Gaussian distribution with the
 * probability density function
 *
 * @f$ \displaystyle
 * p(x) = \frac{1}{\sigma \sqrt{2 \pi}} \cdot \exp \left \{{-\frac{(x - \mu)^2}{2 \sigma^2}} \right \}
 * @f$
 *
 * with given \a µ and \a σ.
 *
 * It can also produce random numbers from a Maxwell distribution
 * with the probability density function
 *
 * @f$ \displaystyle
 * p(x) = \frac{1}{a^3} \sqrt{\frac{2}{\pi}}\, x^2 \exp \left \{{-\frac{x^2}{2 a^2}} \right \}
 * @f$
 *
 * with given parameter @f$a@f$.
 * \sa
 * <a href="http://www.gnu.org/software/gsl/manual/html_node/">
 * GNU scientific library</a>, sections "Random Number
 * Generation:" and "Random Number Distributions:".
 */
class Rng
{
private:
	gsl_rng *rng;	//!< Pointer to the random number generator from the GSL

public:
	/*!
	 * \brief Initialize the random number generator
	 *
	 * This will allocate a random number generator from
	 * the GSL and set its seed value to the
	 * nano second part of a \a clock_gettime() call.
	 */
	Rng() : rng( gsl_rng_alloc(RNG_TYPE) )
	{
		struct timespec tp;

		clock_gettime( CLOCK_REALTIME, &tp );
		gsl_rng_set( rng, tp.tv_nsec );
//		qDebug( "nanosec = %i\n",(unsigned int) tp.tv_nsec );
//		uniform();
	}
	/*!
	 * \brief Have the GSL release the resources for the generator
	 */
	~Rng()
	{
		gsl_rng_free( rng );
	}
	/*!
	 * \brief Produce a uniform deviate in the interval 0 ≤ \a x < 1.
	 * \return Random number
	 *
	 * The current implementation uses the "Mersenne Twister" generator.
	 * \sa
	 * Makoto Matsumoto and Takuji Nishimura (1998): “Mersenne Twister: A 623-dimensionally
	 * equidistributed %uniform pseudorandom number generator”. <em>ACM Trans.
	 * Model. Comput. Simul.</em> \b 8(1), 3–30
	 */
	double uniform() const
	{
		return gsl_rng_uniform( rng );
	}
private:
	/*!
	 * \brief Produce a gaussian deviate using the Marsaglia-Tsang ziggurat algorithm
	 * \param mean Mean value of the desired Gaussian distribution
	 * \param sigma Standard deviation of the Gaussian distribution
	 * \return Random number
	 *
	 * \sa
	 * \li Marsaglia, G. and Tsang, W. W. (2000): The ziggurat method for generating random variables. <em>J. Statis.
	 * Softw.</em> \b 5(8), 1–7
	 * \li <a href="http://www.gnu.org/software/gsl/manual/html_node/Random-number-generator-algorithms.html#Random-number-generator-algorithms">
	 * GNU scientific library</a> for possible algorithms.
	 */
	double gauss( double mean, double sigma ) const
	{
		return (mean + gsl_ran_gaussian_ziggurat( rng, sigma ));
	}
public:
	/*!
	 * \brief Produce a gaussian deviate in a specified interval
	 * \param mean Mean value of the desired Gaussian distribution
	 * \param sigma Standard deviation of the Gaussian distribution
	 * \param limit Interval half width @f$l@f$ in units of @f$\sigma@f$
	 * \return Random number
	 *
	 * The function will generate a number in the interval
	 * @f$(\mu \pm l \sigma)@f$.
	 * Please note that if the interval ist very narrow the function
	 * might fail because it will try only 100 times and then return NaN.
	 *
	 * Setting \a limit to zero (the default) will remove any limitation, thus
	 * setting the interval
	 * to @f$-\infty \dots \infty@f$.
	 */
	double gaussian( double mean, double sigma, double limit = 0.0 ) const
	{
		double tmp, diff;
		int counter = 0;

		if ( limit == 0.0 )
			return gauss( mean, sigma );

		do {
			tmp = gauss( mean, sigma );
			diff = fabs( tmp - mean );
			counter++;
		} while ( (diff > limit*sigma) && (counter<100) );

		if ( counter >= 100 ) tmp = nan("");

		return tmp;
	}
	/*!
	 * \brief Produce a maxwellian deviate wih a 1-D rms \a sigma
	 * \param sigma Standard deviation in one dimension
	 * \return Random number
	 *
	 * The Maxwell distribution is the distribution of the absolute value of
	 * the space velocities the three components of which being normally
	 * distributed with mean value @f$\mu=0@f$ and standard deviation
	 * @f$\sigma = a@f$.
	 *
	 * The algorithm used is the one given by Mohamed. Two other algorithms
	 * are also given but commented out because Mohamed's seems to be slightly
	 * faster.
	 *
	 * \sa
	 * \li <a href="http://adsabs.harvard.edu/abs/2011JSP...145.1653M">
	 * Nader M.A. Mohamed: Efficient Algorithm for Generating Maxwell Random Variables</a>
	 * \li Walck, C.: Handbook of Statistical Distributions for Experimentalists. University
	 * of Stockholm, Stockholm (2007)
	 */
	double maxwellian( double sigma ) const
	{
#if 0
		// Direct method
		const double z1 = gsl_ran_gaussian_ziggurat( rng, sigma );
		const double z2 = gsl_ran_gaussian_ziggurat( rng, sigma );
		const double z3 = gsl_ran_gaussian_ziggurat( rng, sigma );

		return sqrt( z1*z1 + z2*z2 + z3*z3 );
#endif
#if 1
		// Mohamed's algorithm
		const double g2 = 1.647*1.647;
		double y;
		bool reject;

		do {
			double r1 = gsl_rng_uniform_pos( rng );
			double r2 = gsl_rng_uniform( rng );

			y = -2. * log(r1);
			reject = (g2*y < r2*r2/r1/r1);
		} while ( reject );

		return sigma * sqrt( 2.*y );
#endif
#if 0
		// Johnk's algorithm
		const double xi1 = gsl_rng_uniform( rng );
		double r = -log(xi1);
		double w, w1;

		do {
			double xi2 = gsl_rng_uniform( rng );
			double xi3 = gsl_rng_uniform( rng );
			w1 = xi2*xi2;
			double w2 = xi3*xi3;
			w = w1 + w2;
		} while ( w>1 );

		r -= w1*log(gsl_rng_uniform_pos( rng ))/w;

		return sigma * sqrt( 2.*r );
#endif
	}
	/*!
	 * \brief This is an overloaded function.
	 *        It will produce a uniform deviate in the intervall \a lowerLimit < x < \a upperLimit
	 * \param lowerLimit
	 * \param upperLimit
	 * \return Random number
	 */
	double uniform(double lowerLimit, double upperLimit) const
	{
		const double result = gsl_rng_uniform_pos( rng )*(upperLimit-lowerLimit)+lowerLimit;

		return result;
	}
	/*!
	 * \brief Generate random sign on a 50/50 basis
	 * \return +1 or -1
	 */
	int randomsign() const
	{
		double z;

		do {
			z = gsl_rng_uniform_pos( rng );
		} while ( z == 0.5 );

		return (z>0.5 ? 1 : -1);
	}
};

#undef RNG_TYPE

#endif // RNG_H
