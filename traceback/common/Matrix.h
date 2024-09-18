// $Id: Matrix.h 133 2020-03-16 10:52:30Z ifg $
// Matrix class

/*!
 * \file Matrix.h
 *
 * Willi-Hans Steeb:
 * %Matrix Calculus and Kronecker Product with Applications and C++ Programs,
 * World Scientific Publishing Co. Pte. Ltd. 1997, Singapore, New Jersey, London, Hong Kong,
 * ISBN 9810232411
 */

#ifndef MATRIX_H
#define MATRIX_H

#define HAVE_INLINE

#include <iostream>
#include <math.h>
#include <assert.h>
#include <QTextStream>
#include "Vector.h"
#include <gsl/gsl_matrix.h>

template <class T>
/*!
 * \brief The Matrix class
 */
class Matrix
{
protected:
	/// @name Data fields
	/// @{
	int rowNum,			///< Number of rows
	    colNum;			///< Number of columns
	Vector<T> *mat;		///< Pointer to store the data item of the matrix

	/// @}

public:
	/// @name Constructors, destructor
	/// @{
	Matrix();
	Matrix(int nr, int nc);
	Matrix(int nr, int nc, T value);
	Matrix(const Vector<T> &v);
	Matrix(const Matrix<T> &m);
	Matrix(const gsl_matrix *gslm);
	~Matrix();
	/// @}

	/// @name Member functions
	/// @{
	Vector<T> &operator [](int index) const;
	Vector<T> operator ()(int index) const;

	Matrix<T> identity();
	Matrix<T> transpose() const;
	Matrix<T> inverse() const;
	T trace() const;
	T determinant() const;

	int rows() const;
	int cols() const;
	void resize(int r, int c);
	void resize(int r, int c, T value);
	void fill(T value);
	gsl_matrix *toGSL() const;
	/// @}

	/// @name Arithmetic operators
	/// @{
	const Matrix<T> &operator =(const Matrix<T> &m);
	const Matrix<T> &operator =(T value);

	Matrix<T> operator +() const;
	Matrix<T> operator -() const;
	Matrix<T> operator +=(const Matrix<T> &m);
	Matrix<T> operator -=(const Matrix<T> &m);
	Matrix<T> operator *=(const Matrix<T> &m);
	Matrix<T> operator +(const Matrix<T> &m) const;
	Matrix<T> operator -(const Matrix<T> &m) const;
	Matrix<T> operator *(const Matrix<T> &m) const;
	Vector<T> operator *(const Vector<T> &v) const;

	Matrix<T> operator +=(T c);
	Matrix<T> operator -=(T c);
	Matrix<T> operator *=(T c);
	Matrix<T> operator /=(T c);
	Matrix<T> operator +(T c) const;
	Matrix<T> operator -(T c) const;
	Matrix<T> operator *(T c) const;
	Matrix<T> operator /(T c) const;

	/*!
	 * \brief Vectorize operator
	 * \param m %Matrix to create a vector from
	 * \return %Vector that contains elements of \p m, one column after the other
	 */
	friend Vector<T> vec(const Matrix<T> &m)
	{
		int i=0, j, k, size = m.rowNum * m.colNum;
		Vector<T> result( size );

		for ( j=0; j<m.colNum; j++ )
			for ( k=0; k<m.rowNum; k++ ) result[i++] = m.mat[k][j];
		return result;
	}

	/*!
	 * \brief Kronecker product
	 * \param s First matrix for product
	 * \param m Second matrix for product
	 * \return Kronecker product
	 *
	 * Let \f$A\f$ be an \f$m \times n\f$ matrix and let \f$B\f$ be a
	 * \f$p \times q\f$ matrix. Then the \a Kronecker \a product of
	 * \f$A\f$ and \f$B\f$ is that \f$(mp) \times (nq)\f$ matrix defined by
	 * \f[
	 * A \otimes B := \left(
	 * \begin{array}{cccc}
	 * a_{11}B & a_{12}B & \cdots & a_{1n}B \\
	 * a_{21}B & a_{22}B & \cdots & a_{2n}B \\
	 * \vdots \\
	 * a_{m1}B & a_{m2}B & \cdots & a_{mn}B
	 * \end{array}
	 * \right)
	 * \f]
	 * Sometimes the Kronecker product is also called \a direct \a product
	 * or \a tensor \a product.
	 */
	friend Matrix<T> kron(const Matrix<T> &s, const Matrix<T> &m)
	{
		int size1 = s.rowNum * m.rowNum,
			size2 = s.colNum * m.colNum,
				i, j, k, p;
		Matrix<T> result( size1, size2 );

		for ( i=0; i<s.rowNum; i++ )
			for ( j=0; j<s.colNum; j++ )
				for ( k=0; k<m.rowNum; k++ )
					for ( p=0; p<m.colNum; p++ )
						result[k + i*m.rowNum][p + j*m.colNum] = s.mat[i][j] * m.mat[k][p];
		return result;
	}
	/// @}

	/// @name Stream operators
	/// @{
	friend ostream &operator << (ostream &s, const Matrix<T> &m)
	{
		int t = m.cols() - 1;

		for ( int i=0; i<m.rows(); i++ ) {
			s << "[";
			for ( int j=0; j<t; j++ ) s << m[i][j] << " ";
			s << m[i][t] << "]" << endl;
		}

		return s;
	}
	friend QTextStream &operator << (QTextStream &s, const Matrix<T> &m)
	{
		int t = m.cols() - 1;

		for ( int i=0; i<m.rows(); i++ ) {
			s << "[";
			for ( int j=0; j<t; j++ ) s << m[i][j] << " ";
			s << m[i][t] << "]" << endl;
		}

		return s;
	}
	friend istream &operator >> (istream &s, Matrix<T> &m)
	{
		int i, j, num1, num2;

		s.clear();					// set stream state to good
		s >> num1;					// read in row number
		if ( !s.good() ) return s;	// can't get an integer, just return

		s >> num2;					// read in column number
		if ( !s.good() ) return s;	// can't get an integer, just return

		m.resize( num1, num2 );		// resize to Matrix into right order

		for ( i=0; i<num1; i++ )
			for ( j=0; j<num2; j++ ) {
				s >> m[i][j];
				if ( !s.good() ) {
					s.clear( s.rdstate() | ios::badbit );
					return s;
				}
			}
		return s;
	}
	friend QTextStream &operator >> (QTextStream &s, Matrix<T> &m)
	{
		int i, j, num1, num2;

		s >> num1;
		if ( !num1 ) return s;

		s >> num2;
		if ( !num2 ) return s;

		m.resize( num1, num2 );

		for ( i=0; i<num1; i++ )
			for ( j=0; j<num2; j++ ) {
				s >> m[i][j];
			}
		return s;
	}
	/// @}
};

// implementation of class Matrix
template <class T>
/*!
 * \brief Default constructor
 *
 * Matrix() declares a matrix with no size specified.
 * Such a matrix is not usable. To activate the matrix, resize()
 * must be called.
 */
Matrix<T>::Matrix()
	:rowNum(0), colNum(0), mat(NULL) {}

template <class T>
/*!
 * \brief Overloaded constructor
 * \param nr Number of rows
 * \param nc Number of columns
 *
 * Declares an \p nr*nc matrix with the entries value undefined.
 */
Matrix<T>::Matrix(int nr, int nc)
	:rowNum(nr), colNum(nc), mat( new Vector<T>[nr] )
{
	assert( mat != NULL );
	for ( int i=0; i<nr; i++ ) mat[i].resize( nc );
}

template <class T>
/*!
 * \brief Overloaded constructor
 * \param nr Number of rows
 * \param nc Number of columns
 * \param value Initial value
 *
 * Declares an \p nr*nc matrix with all the entries initialized to \p value.
 */
Matrix<T>::Matrix(int nr, int nc, T value)
	:rowNum(nr), colNum(nc), mat(new Vector<T>[nr])
{
	assert( mat != NULL );
	for ( int i=0; i<nr; i++ ) mat[i].resize( nc, value );
}

template <class T>
/*!
 * \brief Overloaded constructor
 * \param v Vector
 *
 * Construct a matrix from a Vector \p v. It is understood
 * that the resultant matrix will contain only one column.
 */
Matrix<T>::Matrix(const Vector<T> &v)
	:rowNum(v.length()), colNum(1), mat(new Vector<T>[rowNum])
{
	assert( mat != NULL );
	for ( int i=0; i<rowNum; i++ ) mat[i].resize( 1, v[i] );
}

template <class T>
/*!
 * \brief Copy constructor
 * \param m Matrix to copy
 */
Matrix<T>::Matrix(const Matrix<T> &m)
	:rowNum(m.rowNum), colNum(m.colNum), mat(new Vector<T>[m.rowNum])
{
	assert( mat != NULL );
	for ( int i=0; i<m.rowNum; i++ ) mat[i] = m.mat[i];
}

template <class T>
/*!
 * \brief Overloaded constructor
 * \param gslm Pointer to a \e gsl_matrix struct
 *
 * Construct a Matrix and copies the values from a matrix struct used by the
 * GNU scientific library
 */
Matrix<T>::Matrix(const gsl_matrix *gslm)
	:rowNum(gslm->size1), colNum(gslm->size2), mat( new Vector<T>[rowNum] )
{
	assert( mat != NULL );
	for ( int i=0; i<rowNum; i++ ) {
		mat[i].resize( colNum );
		for ( int j=0; j<colNum; j++ )
			mat[i][j] = gsl_matrix_get( gslm, i, j );
	}
}

template <class T>
/*!
 * \brief Destructor
 *
 * Release memory back to the free memory pool.
 */
Matrix<T>::~Matrix()
{
	delete [] mat;
}

template <class T>
/*!
 * \brief Subscript operator to access a specific row vector of the matrix
 * \param index Number of row to access (zero based)
 * \return Vector containing values of row with number \a index
 */
Vector<T> &Matrix<T>::operator [] (int index) const
{
	assert( index >= 0 && index < rowNum );
	return mat[index];
}

template <class T>
/*!
 * \brief Parenthesis operator to access a specific column vector of the matrix
 * \param index Number of column to access (zero based)
 * \return Vector containing values of column with number \a index
 */
Vector<T> Matrix<T>::operator () (int index) const
{
	assert( index >= 0 && index < colNum );
	Vector<T> result( rowNum );

	for ( int i=0; i<rowNum; i++ ) result[i] = mat[i][index];

	return result;
}

template <class T>
/*!
 * \brief Create an identity matrix
 * \return %Matrix with principal diagonal set to 1, all other elements set to 0.
 *
 * The matrix must already have the correct size i.e.,
 * number of rows and columns.
 */
Matrix<T> Matrix<T>::identity()
{
	for ( int i=0; i<rowNum; i++ )
		for ( int j=0; j<colNum; j++)
			if ( i==j) mat[i][j] = T( 1 );
			else mat[i][j] = T( 0 );

	return *this;
}

template <class T>
/*!
 * \brief Transpose of matrix
 * \return Matrix with rows and columns swapped
 */
Matrix<T> Matrix<T>::transpose() const
{
	Matrix<T> result( colNum, rowNum );

	for ( int i=0; i<rowNum; i++ )
		for ( int j=0; j<colNum; j++ )
			result[j][i] = mat[i][j];

	return result;
}

// Symbolical Inverse using Leverrier's Method
template <class T>
/*!
 * \brief Symbolical Inverse using Leverrier's Method
 * \return Matrix inverted
 */
Matrix<T> Matrix<T>::inverse() const
{
	assert( rowNum == colNum );
	Matrix<T> B( *this ), D, I( rowNum, colNum );
	T c0(B.trace()), c1;
	int i;

	I.identity();

	for ( i=2; i<rowNum; i++ ) {
		B = *this * (B - c0*I);
		c0 = B.trace() / T( i );
	}
	D = *this * (B - c0*I);
	c1 = D.trace() / T( i );

	return (B - c0*I) / c1;
}

template <class T>
/*!
 * \brief Trace of the matrix
 * \return Trace (sum of elements in principal diagonal)
 */
T Matrix<T>::trace() const
{
	assert( rowNum == colNum );
	T result( 0 );

	for ( int i=0; i<rowNum; i++ ) result += mat[i][i];

	return result;
}

template <class T>
/*!
 * \brief Determinant of the matrix
 *
 */
T Matrix<T>::determinant() const
{
	assert( rowNum == colNum );
	Matrix<T> B( *this ), I( rowNum, colNum, T(0) );
	T c( B.trace() );
	int i;

	for ( i=0; i<rowNum; i++ ) I[i][i] = T( 1 );

	/// Note that determinant of int-type gives zero
	/// because of division by T(i)
	for ( i=2; i<=rowNum; i++ ) {
		B = *this * (B - c*I);
		c = B.trace() / T( i );
	}

	if ( rowNum%2 ) return c;
	return -c;
}

template <class T>
/*!
 * \brief Number of rows of the matrix
 */
int Matrix<T>::rows() const
{
	return rowNum;
}

template <class T>
/*!
 * \brief Number of columns of the matrix
 */
int Matrix<T>::cols() const
{
	return colNum;
}

template <class T>
/*!
 * \brief Resize the matrix
 * \param r New number of rows
 * \param c New number of columns
 *
 * Reallocate the number of rows and column according to the new specification
 * provided in the arguments \p r and \p c.
 */
void Matrix<T>::resize(int r, int c)
{
	int i;
	Vector<T> *newMat = new Vector<T>[r]; assert( newMat != NULL );

	if ( r <= rowNum ) {
		for ( i=0; i<r; i++ ) {
			(mat+i)->resize( c );
			newMat[i] = mat[i];
		}
	} else {
		for ( i=0; i<rowNum; i++ ) {
			(mat+i)->resize( c );
			newMat[i] = mat[i];
		}
		for ( i=rowNum; i<r; i++ ) newMat[i].resize( c );
	}
	delete [] mat;

	rowNum = r; colNum = c;
	mat = newMat;
}

template <class T>
/*!
 * \brief Resize the matrix and initialize the rest of the entries
 * \param r New number of rows
 * \param c New number of columns
 * \param value Value for additional elements
 */
void Matrix<T>::resize(int r, int c, T value)
{
	int i;
	Vector<T> *newMat = new Vector<T>[r]; assert( newMat != NULL );

	if ( r<= rowNum ) {
		for ( i=0; i<r; i++ ) {
			(mat+i)->resize( c, value );
			newMat[i] = mat[i];
		}
	} else {
		for ( i=0; i<rowNum; i++) {
			(mat+i)->resize( c, value );
			newMat[i] = mat[i];
		}
		for ( i=rowNum; i<r; i++ ) newMat[i].resize( c, value );
	}
	delete [] mat;

	rowNum = r, colNum = c;
	mat = newMat;
}

template <class T>
/*!
 * \brief Fill the matrix with \p value
 */
void Matrix<T>::fill(T value) {
	for ( int i=0; i<rowNum; i++ )
		for ( int j=0; j<colNum; j++ )
			mat[i][j] = value;
}

template <class T>
/*!
 * \brief Convert Matrix to a GSL type
 * \return Pointer to a newly allocated \e gsl_matrix.
 *
 * Remember to free the GSL matrix after use!
 */
gsl_matrix *Matrix<T>::toGSL() const
{
	gsl_matrix *m = gsl_matrix_alloc( rows(), cols() );

	for ( int i=0; i<rows(); i++ )
		for ( int j=0; j<cols(); j++ )
			gsl_matrix_set( m, i, j, mat[i][j] );

	return m;
}

template <class T>
/*!
 * \brief Assingment operator
 *
 * Assign one matrix to another (deep copy)
 */
const Matrix<T> &Matrix<T>::operator = (const Matrix<T> &m)
{
	if ( this == &m ) return *this;

	delete [] mat;
	rowNum = m.rowNum; colNum = m.colNum;

	mat = new Vector<T>[m.rowNum]; assert( mat != NULL );

	for ( int i=0; i<m.rowNum; i++ ) mat[i] = m.mat[i];

	return *this;
}

template <class T>
/*!
 * \brief Assignment operator
 *
 * Assign each element to \p value
 */
const Matrix<T> &Matrix<T>::operator = (T value)
{
	for ( int i=0; i<rowNum; i++ ) mat[i] = value;
	return *this;
}

template <class T>
/*!
 * \brief (unary)+
 * \return The matrix itself ("positive" of the matrix)
 */
Matrix<T> Matrix<T>::operator + () const
{
	return *this;
}

template <class T>
/*!
 * \brief (unary)-
 * \return The matrix with each element multiplied with -1
 * ("negative" of the matrix)
 */
Matrix<T> Matrix<T>::operator - () const
{
	return *this * T(-1);
}

template <class T>
/*!
 * \brief Add a matrix and assign the result
 * \param m %Matrix to add
 *
 * Addition according to its normal definition (element-by-element).
 * \sa operator+(const Matrix<T> &) const
 */
Matrix<T> Matrix<T>::operator += (const Matrix<T> &m)
{
	return *this = *this + m;
}

template <class T>
/*!
 * \brief Subtract a matrix and assign the result
 * \param m Matrix to subtract
 *
 * Subtraction according to its normal definition (element-by-element).
 * \sa operator-(const Matrix<T> &) const
 */
Matrix<T> Matrix<T>::operator -= (const Matrix<T> &m)
{
	return *this = *this - m;
}

template <class T>
/*!
 * \brief Multiply a matrix and assign the result
 * \param m %Matrix to multiply
 *
 * Multiplication according to its normal definition (row\f$\times\f$column)
 * \sa operator*(const Matrix<T>&) const
 */
Matrix<T> Matrix<T>::operator *= (const Matrix<T> &m)
{
	return *this = *this * m;
}

template <class T>
/*!
 * \brief Add two matrices and return the result
 * \param m %Matrix to add
 *
 * Addition according to its normal definition( element-by-element).
 * Does not modify the original matrices.
 * \sa operator+=(const Matrix<T> &)
 */
Matrix<T> Matrix<T>::operator + (const Matrix<T> &m) const
{
	assert( rowNum == m.rowNum && colNum == m.colNum );

	Matrix<T> result( *this );

	for ( int i=0; i<rowNum; i++ ) result[i] += m[i];
	return result;
}

template <class T>
/*!
 * \brief Subtract two matrices from each other and return the result
 * \param m Matrix to subtract from the first one
 *
 * Subtraction according to its normal definition (element-by-element).
 * Does not modify the original matrices
 * \sa operator-=(const Matrix<T> &)
 */
Matrix<T> Matrix<T>::operator - (const Matrix<T> &m) const
{
	assert( rowNum == m.rowNum && colNum == m.colNum );

	Matrix<T> result( *this );

	for ( int i=0; i<rowNum; i++ ) result[i] -= m[i];
	return result;
}

template <class T>
/*!
 * \brief Multiply two matrizes and return the result
 * \param m %Matrix to multiply
 *
 * Multiplication according to its normal definition (row\f$\times\f$column).
 * Does not modify the original matrices.
 * \sa operator*=(const Matrix<T>&)
 */
Matrix<T> Matrix<T>::operator * (const Matrix<T> &m) const
{
	assert( colNum == m.rowNum );

	Matrix<T> result( rowNum, m.colNum, T(0) );

	for ( int i=0; i<rowNum; i++ )
		for ( int j=0; j<m.colNum; j++ )
			for ( int k=0; k<colNum; k++ )
				result[i][j] += mat[i][k] * m[k][j];
	return result;
}

template <class T>
/*!
 * \brief Multiply matrix to a vector (apply a linear transformation)
 * \param v %Vector to transform
 * \return Transformed vector
 */
Vector<T> Matrix<T>::operator *(const Vector<T> &v) const
{
	assert( colNum == v.length() );

	Vector<T> result( rowNum );

	// Dot product | is used
	for ( int i=0; i<rowNum; i++) result[i] = (mat[i] | v);

	return result;
}

template <class T>
/*!
 * \brief Add the product of a scalar and the identity matrix \f$cI\f$ and assign the result
 * \param c Constant to add
 * \return \f$A+cI\f$
 *
 * This just means adding \p c to all elements
 * in the principal diagonal.
 */
Matrix<T> Matrix<T>::operator += (T c)
{
	assert( rowNum == colNum );
	for ( int i=0; i<rowNum; i++ ) mat[i][i] += c;
	return *this;
}

template <class T>
/*!
 * \brief Subtract the product of a scalar and the identity matrix \f$cI\f$ and assign the result
 * \param c Constant to subtract
 * \return \f$A-cI\f$
 *
 * This just means subtracting \p c from all elements
 * in the principal diagonal.
 */
Matrix<T> Matrix<T>::operator -= (T c)
{
	assert( rowNum == colNum );
	for ( int i=0; i<rowNum; i++ ) mat[i][i] -= c;
	return *this;
}

template <class T>
/*!
 * \brief Multiply with a "scaled" identity matrix
 * \param c Scaling factor for identity matrix
 * \return \f$A*cI\f$
 *
 * Each element of the matrix will be multiplied by \p c?
 */
Matrix<T> Matrix<T>::operator *= (T c)
{
	for ( int i=0; i<rowNum; i++ ) mat[i] *= c;
	return *this;
}

template <class T>
/*!
 * \brief Divide by a "scaled" identity matrix
 * \param c Scaling factor fo identity matrix
 * \return \f$A/cI\f$
 *
 * Each element of the matrix will be divided by \p c?
 */
Matrix<T> Matrix<T>::operator /= (T c)
{
	for ( int i=0; i<rowNum; i++ ) mat[i] /= c;
	return *this;
}

template <class T>
/*!
 * \brief Add the product of a scalar and the identity matrix \f$cI\f$ and return the result
 * \param c Constant to add
 * \return \f$A+cI\f$
 *
 * This just means adding \p c to all elements
 * in the principal diagonal.
 * \sa operator+=(T)
 */
Matrix<T> Matrix<T>::operator + (T c) const
{
	assert( rowNum == colNum );
	Matrix<T> result( *this );
	return result += c;
}

template <class T>
/*!
 * \brief Subtract the product of a scalar and the identity matrix \f$cI\f$ and return the result
 * \param c Constant to subtract
 * \return \f$A-cI\f$
 *
 * This just means subtracting \p c from all elements
 * in the principal diagonal.
 * \sa operator-=(T)
 */
Matrix<T> Matrix<T>::operator - (T c) const
{
	assert( rowNum == colNum );
	Matrix<T> result( *this );
	return result -= c;
}

template <class T>
/*!
 * \brief Multiply with a "scaled" identity matrix \f$cI\f$ and return the result
 * \param c Sacling factor for identity matrix
 * \return \f$A*cI\f$
 *
 * \sa operator*=(T)
 */
Matrix<T> Matrix<T>::operator * (T c) const
{
	Matrix<T> result( *this );
	return result *= c;
}

template <class T>
/*!
 * \brief Divide by a "scaled" identity matrix \f$cI\f$ and return the result
 * \param c Scaling factor for identity matrix
 * \return \f$A/cI\f$
 *
 * Each element of the matrix will be divided by \p c?
 */
Matrix<T> Matrix<T>::operator / (T c) const
{
	Matrix<T> result( *this );
	return result /= c;
}

template <class T>
Matrix<T> operator + (T value, const Matrix<T> &m)
{
	return m + value;
}

template <class T>
Matrix<T> operator - (T value, const Matrix<T> &m)
{
	return -m + value;
}

template <class T>
Matrix<T> operator * (T value, const Matrix<T> &m)
{
	return m * value;
}

template <class T>
Matrix<T> operator / (T value, const Matrix<T> &m)
{
	Matrix<T> result( m.rows(), m.cols() );

	for ( int i=0; i<result.rows(); i++ ) result[i] = value / m[i];
	return result;
}

// Vectorize operator
//template <class T>

// Kronecker product
//template <class T>

template <class T>
int operator == (const Matrix<T> &m1, const Matrix<T> &m2)
{
	if (m1.rows() != m2.rows() ) return 0;
	for ( int i=0; i<m1.rows(); i++ )
		if ( m1[i] != m2[i] ) return 0;

	return 1;
}

template <class T>
int operator != (const Matrix<T> &m1, const Matrix<T> &m2)
{
	return !( m1 == m2 );
}

//template <class T>

//template <class T>

#endif // MATRIX_H

