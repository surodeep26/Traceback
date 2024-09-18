// $Id: VecNorm.h 126 2019-02-14 12:49:42Z ifg $
// Norms of vectors

/*!
 * \file VecNorm.h
 *
 * Willi-Hans Steeb:
 * %Matrix Calculus and Kronecker Product with Applications and C++ Programs,
 * World Scientific Publishing Co. Pte. Ltd. 1997, Singapore, New Jersey, London, Hong Kong,
 * ISBN 9810232411
 */

#ifndef MVECNORM_H
#define MVECNORM_H

#include <iostream>
#include <math.h>

template <class T>
/*!
 * \brief One-norm of vector
 * \return Norm as sum of the absolute values of the elements of \p v
 *
 * \f$
 * \|\vec{x}\|_1 := |x_1| + |x_2| + \dots + |x_n|
 * \f$
 */
T norm1(const Vector<T> &v)
{
	T result( 0 );

	for ( int i=0; i<v.length(); i++ )
		result = result + abs( v[i] );

	return result;
}

/*!
 * \brief Override of norm1()
 * \return Norm1 of \p v
 *
 * Uses fabs() instead of abs()
 */
double norm1(const Vector<double> &v);

template <class T>
/*!
 * \brief Two-norm of vector
 * \return Norm as root of sum of squared elements of \p v
 *
 * \f$
 * \|\vec{x}\|_2 := \sqrt{x_1^2 + x_2^2 + \dots + x_n^2}
 * \f$
 */
double norm2(const Vector<T> &v)
{
	T result( 0 );

	for ( int i=0; i<v.length(); i++ )
		result = result + v[i]*v[i];

	return sqrt( double( result ) );
}

template <class T>
/*!
 * \brief Infinite-norm of vector
 * \return Norm as largest element of \p v
 *
 * \f$
 * \|\vec{x}\|_\infty := \max\left\{\, |x_1|, |x_2|, \dots, |x_n| \,\right\}
 * \f$
 */
T normI(const Vector<T> &v)
{
	T maxItem( abs( v[0] ) ), temp;

	for ( int i=1; i<v.length(); i++ ) {
		temp = abs( v[i] );
		if ( temp > maxItem ) maxItem = temp;
	}
	return maxItem;
}

/*!
 * \brief Override of normI()
 *
 * Uses fabs() instead of abs()
 */
double normI(const Vector<double> &v);

template <class T>
/*!
 * \brief Normalization of vector \p v
 *
 * Returns
 * \f$
 * \vec{v}/\|\vec{v}\|_2
 * \f$
 *
 */
Vector<T> normalize(const Vector<T> &v)
{
	Vector<T> result( v.length() );
	double length = norm2( v );

	for ( int i=0; i<v.length(); i++ )
		result[i] = v[i] / length;

	return result;
}

#endif // MVECNORM_H
