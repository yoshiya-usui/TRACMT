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
#include "DoubleDenseSquareMatrix.h"
#include "Util.h"
#include "OutputFiles.h"

#include <stddef.h> // For null pointer
#include <stdlib.h> // For exit
#include <iostream>
#include <assert.h>

// Default Constructer
DoubleDenseSquareMatrix::DoubleDenseSquareMatrix():
	DoubleDenseMatrix()
{}

// Destructer
DoubleDenseSquareMatrix::~DoubleDenseSquareMatrix(){
}

// Factorize matrix
void DoubleDenseSquareMatrix::factorizeMatrix(){

	if( m_factorizedMatrix == NULL ){
		m_factorizedMatrix = new double[m_numRows * m_numColumns];
	}

	if( m_ipiv == NULL ){
		m_ipiv = new int[m_numRows];
	}

	Util::factorizeRealSquareMatrix(m_numRows, m_matrix, m_factorizedMatrix, m_ipiv);

}

// Calculate determinant
double DoubleDenseSquareMatrix::calculateDeterminant() const{

	double* matrix = new double[m_numRows*m_numColumns];

	// Copy matrix
	for( int col = 0; col < m_numColumns; ++col ){
		for( int row = 0; row < m_numRows; ++row ){
			// Column major
			const int index = col * m_numRows + row;
			matrix[index] = getValue(row, col);
		}
	}
	const double determinant = Util::calculateDeterminantOfMatrix( m_numRows, matrix );
	delete [] matrix;

	return determinant;

}

// Set Degree of equation
void DoubleDenseSquareMatrix::setDegreeOfEquation( const int nEq ){

	assert( nEq > 0 );

	setNumRowsAndColumns( nEq, nEq );
}

// Get Degree of equation
int DoubleDenseSquareMatrix::getDegreeOfEquation() const{
	return m_numRows;
}

// Solve a linear equation with a right-hand-side vector
void DoubleDenseSquareMatrix::solveLinearEquation( const double* const rhsVector, double* result ) const{

	if( m_factorizedMatrix == NULL ){
		OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
		ptrOutputFiles->writeErrorMessage("Matrix has not been factorized");
	}

	Util::solveLinearEquationRealSquareMatrix( m_numRows, 1, m_ipiv, m_factorizedMatrix, rhsVector, result );

}

// Copy constructer
DoubleDenseSquareMatrix::DoubleDenseSquareMatrix(const DoubleDenseSquareMatrix &matrix ){
	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeErrorMessage("Copy constructer of the class DoubleDenseSquareMatrix is not implemented");
}

// Assignment operator
DoubleDenseSquareMatrix& DoubleDenseSquareMatrix::operator=(const DoubleDenseSquareMatrix& rhs){
	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeErrorMessage("Assignment operator of the class DoubleDenseSquareMatrix is not implemented");
	exit(1);
}
