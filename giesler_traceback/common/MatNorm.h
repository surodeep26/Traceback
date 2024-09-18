// $Id: MatNorm.h 59 2016-08-29 09:12:25Z ifg $
// Norms of Matrices

/*!
 * \file MatNorm.h
 *
 * Willi-Hans Steeb:
 * %Matrix Calculus and Kronecker Product with Applications and C++ Programs,
 * World Scientific Publishing Co. Pte. Ltd. 1997, Singapore, New Jersey, London, Hong Kong,
 * ISBN 9810232411
 */

#ifndef MATNORM_H
#define MATNORM_H

#include <iostream>
#include <math.h>
#include "Vector.h"
#include "VecNorm.h"

template <class T>
/*!
 * \brief One-norm of the matrix
 *
 * Defined as the maximum value of the sum of the entries in column vectors,
 *
 * \f$
 * \| A \|_1 = \displaystyle\max_{1 \leq j \leq n} \left( \sum_{i=1}^n | a_{ij} | \right)
 * \f$
 */
T norm1(const Matrix<T> &m)
{
	T maxItem( 0 ), temp;
	int i, j;

	for ( i=0; i<m.rows(); i++ )
		maxItem += m[i][0];

	for ( i=0; i<m.cols(); i++ ) {
		temp = T( 0 );
		for ( j=0; j<m.rows(); j++ )
			temp += abs( m[j][i] );
		if ( temp > maxItem ) maxItem = temp;
	}
	return maxItem;
}

template <class T>
/*!
 * \brief Infinite-norm of the matrix
 *
 * Defined as the maximum value of the sum of the entries in row vectors,
 *
 * \f$
 * \| A \|_{\infty} = \displaystyle\max_{1 \leq i \leq n} \left( \sum_{j=1}^n | a_{ij} | \right)
 * \f$
 */
T normI(const Matrix<T> &m)
{
	T maxItem( norm1( m[0] ) );

	for ( int i=1; i<m.rows(); i++ )
		if ( norm1( m[i] ) > maxItem ) maxItem = norm1( m[i] );

	return maxItem;
}

template <class T>
/*!
 * \brief Hilbert-Schmidt norm of the matrix
 *
 * Defined as
 *
 * \f$
 * \| A \|_H = (\mbox{tr}(A^*A))^{1/2} = (\mbox{tr}(AA^*))^{1/2} =
 * \displaystyle\sqrt{\sum_{i=1}^n \sum_{j=1}^m |a_{ij}|^2}
 * \f$
 */
T normH(const Matrix<T> &m)
{
	return sqrt( (m*(m.transpose())).trace() );
}

#endif // MATNORM_H
