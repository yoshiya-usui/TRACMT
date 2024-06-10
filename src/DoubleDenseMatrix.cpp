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
#include "DoubleDenseMatrix.h"
#include "OutputFiles.h"
#include "Util.h"

#include <stddef.h> // For null pointer
#include <stdlib.h> // For exit
#include <iostream>
#include <assert.h>

//Default Constructer
DoubleDenseMatrix::DoubleDenseMatrix():
	m_numRows(0),
	m_numColumns(0),
	m_numComponents(0),
	m_matrix(NULL),
	m_factorizedMatrix(NULL),
	m_ipiv(NULL)
{}

// Destructer
DoubleDenseMatrix::~DoubleDenseMatrix(){

	if( m_matrix != NULL ){
		delete[] m_matrix;
		m_matrix = NULL;
	}

	if( m_factorizedMatrix != NULL ){
		delete[] m_factorizedMatrix;
		m_factorizedMatrix = NULL;
	}

	if( m_ipiv != NULL ){
		delete [] m_ipiv;
	}

}

// Set number of rows and columns
void DoubleDenseMatrix::setNumRowsAndColumns( const int nrows, const int ncols ){

	// Total number of rows
	m_numRows = nrows;

	// Total number of columns
	m_numColumns = ncols;

	// Total number of components
	m_numComponents = nrows * ncols;

	m_matrix = new double [m_numComponents];

	// Zero clear
	zeroClearMatrix();

}

// Add value to matrix
void DoubleDenseMatrix::addValue( const int row, const int col, const double val ){

	assert( row <= m_numRows - 1 );
	assert( row >= 0 );
	assert( col <= m_numColumns - 1 );
	assert( col >= 0  );

	// Column major
	const int index = col * m_numRows + row;
	m_matrix[index] += val;
		
}

// Multiply a scalar to the matrix
void DoubleDenseMatrix::multiplyScalarToMatrix( const double scalar ){

	for( int i = 0; i < m_numComponents; ++i ){
		m_matrix[i] *= scalar;
	}

}

// Set value to matrix
void DoubleDenseMatrix::setValue( const int row, const int col, const double val ){

	assert( row <= m_numRows - 1 );
	assert( row >= 0 );
	assert( col <= m_numColumns - 1 );
	assert( col >= 0  );

	// Column major
	const int index = col * m_numRows + row;
	m_matrix[index] = val;

}

//Zero clear matrix
void DoubleDenseMatrix::zeroClearMatrix(){

	for( int i = 0; i < m_numComponents; ++i ){
		m_matrix[i] = 0.0; // Zero clear
	}

	if( m_factorizedMatrix != NULL ){
		delete[] m_factorizedMatrix;
		m_factorizedMatrix = NULL;
	}

	if( m_ipiv != NULL ){
		delete [] m_ipiv;
		m_ipiv = NULL;
	}

}

// Get total number of rows
int DoubleDenseMatrix::getNumRows() const{
	return m_numRows;
}

// Get total number of columns
int DoubleDenseMatrix::getNumColumns() const{
	return m_numColumns;
}

// Get value
double DoubleDenseMatrix::getValue( const int row, const int col ) const{

	// Column major
	const int index = col * m_numRows + row;
	return m_matrix[index];

}

// Debug write the matrix componets
void DoubleDenseMatrix::debugWriteMatrix() const{

	std::cout << "[";
	for( int row = 0; row < m_numRows; ++row ){
		for( int col = 0; col < m_numColumns; ++col ){		
			std::cout << getValue(row, col) << " ";
		}
		if( row+1 < m_numRows ){
			std::cout << ";";
		}
	}
	std::cout << "]" << std::endl;

}

//Copy constructer
DoubleDenseMatrix::DoubleDenseMatrix(const DoubleDenseMatrix &matrix ){
	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeErrorMessage("Copy constructer of the class DoubleDenseMatrix is not implemented");
}

// Assignment operator
DoubleDenseMatrix& DoubleDenseMatrix::operator=(const DoubleDenseMatrix& rhs){
	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeErrorMessage("Assignment operator of the class DoubleDenseMatrix is not implemented");
	exit(1);
}
