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
#ifndef DBLDEF_DOUBLE_DENSE_SQUARE_MATRIX
#define DBLDEF_DOUBLE_DENSE_SQUARE_MATRIX

#include "DoubleDenseMatrix.h"

#include <string>

class DoubleDenseSquareMatrix : public DoubleDenseMatrix{

public:

	// Default Constructer
	explicit DoubleDenseSquareMatrix();

	// Destructer
	virtual ~DoubleDenseSquareMatrix();

	// Factorize matrix
	virtual void factorizeMatrix();

	// Calculate determinant
	double calculateDeterminant() const;

	// Set Degree of equation
	virtual void setDegreeOfEquation( const int nEq );

	// Get Degree of equation
	int getDegreeOfEquation() const;

	// Solve a linear equation with a right-hand-side vector
	virtual void solveLinearEquation( const double* const rhsVector, double* result ) const;

private:
	// Copy constructer
	DoubleDenseSquareMatrix(const DoubleDenseSquareMatrix &matrix );

	// Assignment operator
	DoubleDenseSquareMatrix& operator=(const DoubleDenseSquareMatrix& rhs);

};

#endif
