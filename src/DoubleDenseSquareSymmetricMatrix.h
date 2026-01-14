//--------------------------------------------------------------------------
// Copyright(c) 2024, Yoshiya Usui
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met :
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and /or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED.IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//--------------------------------------------------------------------------
#ifndef DBLDEF_DOUBLE_DENSE_SQUARE_SYMMETRIC_MATRIX
#define DBLDEF_DOUBLE_DENSE_SQUARE_SYMMETRIC_MATRIX

#include "DoubleDenseSquareMatrix.h"

#include <set>

class DoubleDenseSquareSymmetricMatrix : public DoubleDenseSquareMatrix {

public:

	// Default Constructer
	explicit DoubleDenseSquareSymmetricMatrix();

	// Destructer
	virtual ~DoubleDenseSquareSymmetricMatrix();

	// Factorize matrix
	virtual void factorizeMatrix();

	// Calculate determinant
	double determinant();

	// Factorize and solve a linear equation with a right-hand-side vector
	virtual void factorizeAndSolveLinearEquation( const double* const rhsVector, double* result ) const;

	// Set Degree of equation
	virtual void setDegreeOfEquation( const int nEq );
	
	// Set number of rows and columns
	virtual void setNumRowsAndColumns( const int nrows, const int ncols );

	// Add value to matrix
	virtual void addValue( const int row, const int col, const double val );

	// Get value
	virtual double getValue( const int row, const int col ) const;

	// Set value to matrix
	virtual void setValue( const int row, const int col, const double val );

	// Solve a linear equation with a right-hand-side vector
	virtual void solveLinearEquation( const double* const rhsVector, double* result ) const;

	// Zero clear matrix
	virtual void zeroClearMatrix();

private:

	// Copy constructer
	DoubleDenseSquareSymmetricMatrix(const DoubleDenseSquareSymmetricMatrix &matrix );

	// Assignment operator
	DoubleDenseSquareSymmetricMatrix& operator=(const DoubleDenseSquareSymmetricMatrix& rhs);

};

#endif
