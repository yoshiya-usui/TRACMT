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
#ifndef DBLDEF_DOUBLE_DENSE_MATRIX
#define DBLDEF_DOUBLE_DENSE_MATRIX

class DoubleDenseMatrix{

public:

	//Default Constructer
	explicit DoubleDenseMatrix();

	//Constructer
	explicit DoubleDenseMatrix( const int nrows, const int ncols );

	//Destructer
	virtual ~DoubleDenseMatrix();

	// Set number of rows and columns
	virtual void setNumRowsAndColumns( const int nrows, const int ncols );

	// Add value to matrix
	virtual void addValue( const int row, const int col, const double val );

	// Multiply a scalar to the matrix
	void multiplyScalarToMatrix( const double scalar );

	// Set value to matrix
	virtual void setValue( const int row, const int col, const double val );

	// Zero clear matrix
	virtual void zeroClearMatrix();

	// Get number of rows
	int getNumRows() const;

	// Get number of columns
	int getNumColumns() const;

	// Get value
	virtual double getValue( const int row, const int col ) const;

	//Debug write the matrix componets
	void debugWriteMatrix() const;

protected:

	//Total number of rows
	int m_numRows;

	//Total number of columns
	int m_numColumns;

	//Total number of components
	int m_numComponents;

	//Values of non-zero compnents
	double* m_matrix;

	// Factorized matrix
	double* m_factorizedMatrix;

	// Array for LAPACK
	int* m_ipiv;
	
private:
	//Copy constructer
	DoubleDenseMatrix(const DoubleDenseMatrix &matrix );

	// Assignment operator
	DoubleDenseMatrix& operator=(const DoubleDenseMatrix& rhs);

};

#endif
