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
#include "DoubleDenseSquareSymmetricMatrix.h"
#include "OutputFiles.h"
#include "Util.h"

#include <stddef.h> // For null pointer
#include <stdlib.h> // For exit
#include <iostream>
#include <assert.h>

// Default Constructer
DoubleDenseSquareSymmetricMatrix::DoubleDenseSquareSymmetricMatrix():
	DoubleDenseSquareMatrix()
{}

// Destructer
DoubleDenseSquareSymmetricMatrix::~DoubleDenseSquareSymmetricMatrix(){
}

// Factorize matrix
void DoubleDenseSquareSymmetricMatrix::factorizeMatrix(){

	if( m_factorizedMatrix == NULL ){
		m_factorizedMatrix = new double[m_numRows * m_numColumns];
	}

	if( m_ipiv == NULL ){
		m_ipiv = new int[m_numRows];
	}
	Util::factorizeRealSymmetricMatrix(m_numRows, m_matrix, m_factorizedMatrix, m_ipiv);

//#ifdef _DEBUG_WRITE
//	for( int row = 0; row < m_numRows; ++row ){
//		for( int col = 0; col < m_numColumns; ++col ){
//			// Column major
//			const int index = col * m_numRows + row;
//			std::cout << "row col val: " << row << " " << col << " " << m_factorizedMatrix[index] << std::endl;
//		}
//	}
//#endif

}

// Factorize and solve a linear equation with a right-hand-side vector
void DoubleDenseSquareSymmetricMatrix::factorizeAndSolveLinearEquation( const double* const rhsVector, double* result ) const{

	Util::factorizeAndSolveLinearEquationRealSymmetricMatrix(m_numRows, 1, m_matrix, rhsVector, result);

}

// Set Degree of equation
void DoubleDenseSquareSymmetricMatrix::setDegreeOfEquation( const int nEq ){

	assert( nEq > 0 );

	//Total number of rows
	m_numRows = nEq;

	//Total number of columns
	m_numColumns = nEq;

	// Total number of components
	m_numComponents = ( nEq * nEq + nEq ) / 2;

	m_matrix = new double [m_numComponents];

	// Zero clear
	zeroClearMatrix();

}

// Set number of rows and columns
void DoubleDenseSquareSymmetricMatrix::setNumRowsAndColumns( const int nrows, const int ncols ){

	assert( nrows == ncols );
	setDegreeOfEquation(nrows);

}

// Add value to matrix
void DoubleDenseSquareSymmetricMatrix::addValue( const int row, const int col, const double val ){

	// Allow lower triangle part only
	assert( row >= col );

	// Column major
	const int index = m_numColumns * col + ( col - col * col ) / 2 + ( row - col );
	m_matrix[index] += val;

}

// Get value
double DoubleDenseSquareSymmetricMatrix::getValue( const int row, const int col ) const{

	if( row >= col ){
		// Lower triangle part
		// Column major
		const int index = m_numColumns * col + ( col - col * col ) / 2 + ( row - col );
		return m_matrix[index];
	}else{
		// Upper triangle part
		return getValue( col, row );
	}

}

// Set value to matrix
void DoubleDenseSquareSymmetricMatrix::setValue( const int row, const int col, const double val ){

	// Allow lower triangle part only
	assert( row >= col );

	// Column major
	const int index = m_numColumns * col + ( col - col * col ) / 2 + ( row - col );
	m_matrix[index] = val;

}

// Solve a linear equation with a right-hand-side vector
void DoubleDenseSquareSymmetricMatrix::solveLinearEquation( const double* const rhsVector, double* result ) const{

	if( m_factorizedMatrix == NULL ){
		OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
		ptrOutputFiles->writeErrorMessage("Matrix has not been factorized");
	}

	Util::solveLinearEquationRealSymmetricMatrix( m_numRows, 1, m_ipiv, m_factorizedMatrix, rhsVector, result );

}

// Zero clear matrix
void DoubleDenseSquareSymmetricMatrix::zeroClearMatrix(){

	DoubleDenseMatrix::zeroClearMatrix();

}

//Copy constructer
DoubleDenseSquareSymmetricMatrix::DoubleDenseSquareSymmetricMatrix(const DoubleDenseSquareSymmetricMatrix &matrix ){
	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeErrorMessage("Copy constructer of the class DoubleDenseSquareSymmetricMatrix is not implemented");
}

// Assignment operator
DoubleDenseSquareSymmetricMatrix& DoubleDenseSquareSymmetricMatrix::operator=(const DoubleDenseSquareSymmetricMatrix& rhs){
	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeErrorMessage("Assignment operator of the class DoubleDenseSquareSymmetricMatrix is not implemented");
	exit(1);
}

