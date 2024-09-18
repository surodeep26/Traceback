// $Id: Vector.h 79 2016-09-27 12:09:41Z ifg $
// Vector class

/*!
 * \file Vector.h
 *
 * Willi-Hans Steeb:
 * %Matrix Calculus and Kronecker Product with Applications and C++ Programs,
 * World Scientific Publishing Co. Pte. Ltd. 1997, Singapore, New Jersey, London, Hong Kong,
 * ISBN 9810232411
 */

#ifndef MVECTOR_H
#define MVECTOR_H

#include <iostream>
#include <string>
#include <math.h>
#include <assert.h>
#include <QTextStream>

using namespace std;

// forward declaration
template <class T> class Matrix;

/*!
 * \brief The Vector class
 */
template <class T> class Vector
{
private:
	/// @name Data fields
	/// @{
	int size;	///< Size of the vector (number of elements)
	T *data;	///< Pointer to the data

	/// @}

public:
	/// @name Constructors, destructor
	/// @{
	Vector();
	Vector(int n);
	Vector(int n, T value);
	Vector(int n, const T p[]);
	Vector(const Vector<T> &v);
	~Vector() { delete [] data; }
	/// @}

	/// @name Member functions
	/// @{
	/*!
	 * \brief Index operator (zero based)
	 * \param i index
	 * \return i-th element
	 */
	T& operator[] (int i) const { assert( i >= 0 && i < size ); return data[i]; }

	/*!
	 * \brief Size of Vector
	 * \return Number of elements
	 */
	int length() const { return size; }
	void swap( int i, int j);
	void resize( int length);
	void resize( int length, T value);
	void reset(int length);
	void reset(int length, T value);
	/// @}

	/// @name Arithmetic operators
	/// @{
	const Vector<T> &operator =(const Vector<T> &v);
	const Vector<T> &operator =(T value);
	/*!
	 * \brief (unary)+
	 * \return The vector itself
	 */
	Vector<T> &operator +() const { return *this; }
	/*!
	 * \brief (unary)-
	 * \return The "negative" of the vector
	 */
	Vector<T> &operator -() const { return *this * T(-1); }
	Vector<T> &operator +=(const Vector<T> &v);
	Vector<T> &operator -=(const Vector<T> &v);
	Vector<T> &operator *=(const Vector<T> &v);
	Vector<T> &operator /=(const Vector<T> &v);
	Vector<T> operator +(const Vector<T> &v) const;
	Vector<T> operator -(const Vector<T> &v) const;
	Vector<T> operator *(const Vector<T> &v) const;
	Vector<T> operator /(const Vector<T> &v) const;

	Vector<T> &operator +=(T x);
	Vector<T> &operator -=(T x);
	Vector<T> &operator *=(T x);
	Vector<T> &operator /=(T x);
	Vector<T> operator +(T x) const;
	Vector<T> operator -(T x) const;
	Vector<T> operator *(T x) const;
	Vector<T> operator /(T x) const;

	// Dot product / Inner product
	T operator |(const Vector<T> &v);

	// Cross Product
	Vector<T> operator %(const Vector<T> &v);
	/// @}

	/// @name I/O stream functions
	/// @{
	friend class Matrix<T>;
	friend ostream &operator<<(ostream &s, const Vector<T> &v)
	{
		int lastnum = v.length();
		for ( int i=0; i<lastnum; i++ )
			s << "[" << v[i] << "]" << endl;
		return s;
	}
	friend QTextStream &operator<<(QTextStream &s, const Vector<T> &v)
	{
		int lastnum = v.length();
		for ( int i=0; i<lastnum; i++ )
			s << "[" << v[i] << "]" << endl;
		return s;
	}
	friend istream &operator >>(istream &s, Vector<T> &v)
	{
		int i, num;
		s.clear();					// set stream state to good
		s >> num;					// read size of vector
		if ( !s.good() ) return s;	// can't get an integer, just return
		v.resize( num );			// resize Vector v
		for ( i=0; i<num; i++ ) {
			s >> v[i];				// read in entries
			if ( !s.good() ) {
				s.clear( s.rdstate() | ios::badbit );
				return s;
			}
		}
		return s;
	}
	friend QTextStream &operator >>(QTextStream &s, Vector<T> &v)
	{
		int i, num;
		s >> num;
		if ( !num ) return s;
		v.resize( num );
		for ( i=0; i<num; i++ ) {
			s >> v[i];
		}
		return s;
	}
	/// @}
};

// implementation of class Vector
template <class T>
/*!
 * \brief Default constructor
 *
 * Initialize data pointer and size to zero
 */
Vector<T>::Vector()
	:size( 0 ), data( NULL ) {}

template <class T>
/*!
 * \brief Overloaded constructor
 * \param n Number of elements of the new vector
 *
 * Create a new vector with \p n elements and initialize them to zero
 */
Vector<T>::Vector(int n)
	:size( n ), data( new T[n] )
{
	assert( data != NULL );
	for ( int i=0; i<size; i++ )
		data[i] = 0;
}

template <class T>
/*!
 * \brief Overloaded constructor
 * \param n Number of elements of the new vector
 * \param value Initial value for each element
 *
 * Create a new vector with \p n elements and initialize them to \p value
 */
Vector<T>::Vector(int n, T value)
	:size( n ), data( new T[n] )
{
	assert( data != NULL );
	for ( int i=0; i<n; i++ ) data[i] = value;
}

template <class T>
/*!
 * \brief Overloaded constructor
 * \param n Number of elements of the new vector
 * \param p Array with initial values
 *
 * No range checking is done, so handle with care!
 */
Vector<T>::Vector(int n, const T p[])
    :size( n ), data( new T[n] )
{
	for ( int i=0; i<size; i++ ) data[i] = p[i];
}

template <class T>
/*!
 * \brief Copy constructor
 * \param v Vector to copy
 *
 * Create a vector as a deep copy of another one
 */
Vector<T>::Vector(const Vector<T> &v)
	:size( v.size ), data( new T[v.size] )
{
	assert( data != NULL );
	for ( int i=0; i<v.size; i++ ) data[i] = v.data[i];
}

template <class T>
/*!
 * \brief Swap elements with indexes \p i and \p j
 */
void Vector<T>::swap(int i, int j)
{
	T tmp = data[i];
	data[i] = data[j];
	data[j] = tmp;
}

template <class T>
/*!
 * \brief Resize vector
 * \param length New size of the vector
 *
 * Resizes the vector to \p length and fills the rest of the
 * vector with value 0 typecasted to \a T if \p length is bigger than #size,
 * otherwise drops the odd ones out.
 */
void Vector<T>::resize(int length)
{
	int i;
	T zero(0);
	T *newData = new T[length]; assert( newData != NULL );
	if ( length <= size )
		for ( i=0; i<length; i++ ) newData[i] = data[i];
	else {
		for ( i=0; i<size; i++ )   newData[i] = data[i];
		for ( i=size; i<length; i++ ) newData[i] = zero;
	}
	delete [] data;
	size = length;
	data = newData;
}

template <class T>
/*!
 * \brief Overloaded function
 *
 * Resizes the vector to \p length and fills the rest of
 * the vector with \p value.
 */
void Vector<T>::resize(int length, T value)
{
	int i;
	T *newData = new T[length]; assert( newData != NULL );
	if ( length <= size )
		for ( i=0; i<length; i++ )   newData[i] = data[i];
	else {
		for ( i=0; i<size; i++)      newData[i] = data[i];
		for ( i=size; i<length; i++ ) newData[i] = value;
	}
	delete [] data;
	size = length;
	data = newData;
}

template <class T>
/*!
 * \brief Resets the vector to size \p length
 * \param length New size of vector
 *
 * This will delete all elements and create a new vector with
 * T(0) elements.
 */
void Vector<T>::reset(int length)
{
	T zero( 0 );
	delete [] data;
	data = new T[length]; assert( data != NULL );
	size = length;
	for ( int i=0; i<size; i++ ) data[i] = zero;
}

template <class T>
/*!
 * \brief Overloaded function
 *
 * Resets the vector to size \p length and initializes entries to \p value
 */
void Vector<T>::reset(int length, T value)
{
	delete [] data;
	data = new T[length]; assert( data != NULL );
	size = length;
	for ( int i=0; i<size; i++ ) data[i] = value;
}

template <class T>
/*!
 * \brief Assignment operator
 *
 * Assign one vector to another (deep copy)
 */
const Vector<T> &Vector<T>::operator =(const Vector<T> &v)
{
	if ( this == &v ) return *this;
	if ( size != v.size ) {
		delete [] data;
		data = new T[v.size]; assert( data != NULL );
		size = v.size;
	}
	for ( int i=0; i<size; i++ ) data[i] = v.data[i];

	return *this;
}

template <class T>
/*!
 * \brief Overloaded assignment operator
 *
 * Assign each element to \p value
 */
const Vector<T> &Vector<T>::operator =(T value)
{
	for ( int i=0; i< size; i++ ) data[i] = value;
	return *this;
}

template <class T>
/*!
 * \brief Add different values to each element
 * \param v Vector with elements to add
 */
Vector<T> &Vector<T>::operator +=(const Vector<T> &v)
{
	assert( size == v.size );
	for ( int i=0; i<size; i++ ) data[i] += v.data[i];
	return *this;
}

template <class T>
/*!
 * \brief Subtract different values from each element
 * \param v Vector with elements to subtract
 */
Vector<T> &Vector<T>::operator -=(const Vector<T> &v)
{
	assert( size == v.size );
	for ( int i=0; i<size; i++ ) data[i] -= v.data[i];
	return *this;
}

template <class T>
/*!
 * \brief Multiply different values to each element
 * \param v Vector with elements to multiply
 */
Vector<T> &Vector<T>::operator *=(const Vector<T> &v)
{
	assert( size == v.size );
	for ( int i=0; i<size; i++ ) data[i] *= v.data[i];
	return *this;
}

template <class T>
/*!
 * \brief Divide each element by different values
 * \param v Vector with elements for division
 */
Vector<T> &Vector<T>::operator /=(const Vector<T> &v)
{
	assert( size == v.size );
	for ( int i=0; i<size; i++ ) data[i] /= v.data[i];
	return *this;
}

template <class T>
/*!
 * \brief Add two vectors
 */
Vector<T> Vector<T>::operator +(const Vector<T> &v) const
{
	Vector<T> result( *this );
	return result += v;
}

template <class T>
/*!
 * \brief Subtract a vector
 */
Vector<T> Vector<T>::operator -(const Vector<T> &v) const
{
	Vector<T> result(*this);

	return result -= v;
}

template <class T>
/*!
 * \brief Multiply a vector element-wise
 */
Vector<T> Vector<T>::operator *(const Vector<T> &v) const
{
	Vector<T> result(*this);
	return result *= v;
}

template <class T>
/*!
 * \brief Divide by a vector element-wise
 */
Vector<T> Vector<T>::operator /(const Vector<T> &v) const
{
	Vector<T> result( *this );
	return result /= v;
}

template <class T>
/*!
 * \brief Add a single value to each element
 * \param x Value to add
 */
Vector<T> &Vector<T>::operator +=(T x)
{
	for ( int i=0; i<size; i++ ) data[i] += x;
	return *this;
}

template <class T>
/*!
 * \brief Subtract a single value from each element
 * \param x Value to subtract
 */
Vector<T> &Vector<T>::operator -=(T x)
{
	for ( int i=0; i<size; i++ ) data[i] -= x;
	return *this;
}

template <class T>
/*!
 * \brief Multiply a single value to each element
 * \param x Value to multiply
 *
 * Multiply a scalar to a vector i.e., scaling. Modifies original vector.
 */
Vector<T> &Vector<T>::operator *=(T x)
{
	for ( int i=0; i<size; i++ ) data[i] *= x;
	return *this;
}

template <class T>
/*!
 * \brief Divide each element by a value
 * \param x Value for division
 *
 * Another form of scaling. Modifies original vector.
 */
Vector<T> &Vector<T>::operator /=(T x)
{
	for ( int i=0; i<size; i++ ) data[i] /= x;
	return *this;
}

template <class T>
/*!
 * \brief Add a single value to each element
 * \param x Value to add
 * \return New vector with result
 *
 * Operates on a deep copy and does not modify original vector
 */
Vector<T> Vector<T>::operator +(T x) const
{
	Vector<T> result(*this);
	return result += x;
}

template <class T>
/*!
 * \brief Subtract a value from each element
 * \param x Value to subtract
 * \return New vector with result
 *
 * Operates on a deep copy and does not modify original vector
 */
Vector<T> Vector<T>::operator -(T x) const
{
	Vector<T> result(*this);
	return result -= x;
}

template <class T>
/*!
 * \brief Multiply each element with a value
 * \param x Scaling factor
 * \return New vector with the result
 *
 * Operates on a deep copy and does not modify the original vector
 */
Vector<T> Vector<T>::operator *(T x) const
{
	Vector<T> result(*this);
	return result *= x;
}

template <class T>
/*!
 * \brief Divide each element by a value
 * \param x Value for division
 * \return New vector with result
 *
 * Operates on a deep copy and does not modify the original vector
 */
Vector<T> Vector<T>::operator /(T x) const
{
	Vector<T> result(*this);
	return result /= x;
}

template <class T>
Vector<T> operator +(T x, const Vector<T> &v)
{
	return v+x;
}

template <class T>
Vector<T> operator -(T x, const Vector<T> &v)
{
	return -v+x;
}

template <class T>
Vector<T> operator *(T x, const Vector<T> &v)
{
	return v*x;
}

template <class T>
Vector<T> operator /(T x, const Vector<T> &v)
{
	Vector<T> result( v.length() );
	for ( int i=0; i<result.length(); i++ )
		result[i] = x / v[i];
	return result;
}

template <class T>
/*!
 * \brief Scalar (aka dot) product
 * \return Scalar of type \a T
 *
 * Will calculate scalar product of two vectors
 *
 * \f$
 * \vec{u} \cdot \vec{v} = \displaystyle\sum_{j=0}^{n-1} u_j v_j
 * \f$
 */
T Vector<T>::operator |(const Vector<T> &v)
{
	assert( size == v.size );
	T result( 0 );
	for ( int i=0; i<size; i++ ) result += data[i] * v.data[i];
	return result;
}

template <class T>
/*!
 * \brief %Vector (aka cross) product.
 *
 * Will calculate vector product of two vectors
 *
 * \f$
 * \vec{u} \times \vec{v} =
 * \left|
 * \begin{array}{ccc}
 * \vec{e_x} & \vec{e_y} & \vec{e_z} \\
 * u_x & u_y & u_z \\
 * v_x & v_y & v_z
 * \end{array}
 * \right|
 * \f$.
 *
 * Only implemented for #size = 3.
 */
Vector<T> Vector<T>::operator %(const Vector<T> &v)
{
	assert( size == 3 && v.size == 3);
	Vector<T> result(3);

	result[0] = data[1] * v.data[2] - v.data[1] * data[2];
	result[1] = v.data[0] * data[2] - data[0] * v.data[2];
	result[2] = data[0] * v.data[1] - v.data[0] * data[1];

	return result;
}

// Equality
template <class T>
/*!
 * \brief Equality of vectors
 * \param u
 * \param v
 * \return \a True if \p u and \p v contain the sme elements in the same order \a false otherwise
 */
int operator ==(const Vector<T> &u, const Vector<T> &v)
{
	if (u.length() != v.length() ) return 0;
	for ( int i=0; i<u.length(); i++ )
		if ( u[i] != v[i] ) return 0;
	return 1;
}

// Inequality
template <class T>
/*!
 * \brief Inequality of vectors
 * \param u
 * \param v
 * \return True or false
 *
 * Converse of equality
 * \sa operator ==()
 */
int operator !=(const Vector<T> &u, const Vector<T> &v)
{
	return !( u == v );
}

#endif // MVECTOR_H
