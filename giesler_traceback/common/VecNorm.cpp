//$Id: VecNorm.cpp 126 2019-02-14 12:49:42Z ifg $
// Some functions which are not templates mustn't be inlined

#include "Vector.h"

double norm1(const Vector<double> &v)
{
	double result( 0 );

	for ( int i=0; i<v.length(); i++ )
		result = result + fabs( v[i] );

	return result;
}

double normI(const Vector<double> &v)
{
	double maxItem( fabs( v[0] ) ), temp;

	for ( int i=1; i<v.length(); i++ ) {
		temp = fabs( v[i] );
		if ( temp > maxItem ) maxItem = temp;
	}
	return maxItem;
}
