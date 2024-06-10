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
#include "LapackInterface.h"
#include "OutputFiles.h"
#include "Util.h"

extern "C"
{
#include "f2c.h"
#include "clapack.h"
}

// Calcualte all eigenvalues and eigenvectors of a real symmetric matrix
// @note Matrix and vectors are overwritten
void LapackInterface::calculateEigenValuesAndVectorsOfRealSymmetricMatrix( const int dimension, double* matrix, double* vectors ){

	integer n = static_cast<integer>(dimension);
	integer lda = n;
	integer lwork = n * n;
	double* work = new double[lwork];
	integer info = 0;

	dsyev_("V", "L", &n, matrix, &lda, vectors, work, &lwork, &info);
	delete [] work;

	if( info < 0 ){
		OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
		ptrOutputFiles->writeErrorMessage("An argument had an illegal value : info="+Util::toString(info) );
	}else if( info > 0 ){
		OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
		ptrOutputFiles->writeErrorMessage("Eigenvalue calculation is not converged : info="+Util::toString(info) );
	}

}

// Factrize and solve linear equation with real coefficents matrix
// @note Matrix and vectors are overwritten
void LapackInterface::factorizeAndSolveLinearEquationRealMatrix( const int dimension, const int nRhs, double* matrix, double* vectors ){

	integer n = static_cast<integer>( dimension );
	integer nb = static_cast<integer>( nRhs );
	integer lda = n;
	integer* ipiv = new integer[dimension];
	integer ldb = n;
	integer info(0);

	dgesv_(&n, &nb, matrix, &lda, ipiv, vectors, &ldb, &info);
	delete [] ipiv;

	if( info < 0 ){
		OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
		ptrOutputFiles->writeErrorMessage("An argument had an illegal value : info="+Util::toString(info) );
	}else if( info > 0 ){
		OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
		ptrOutputFiles->writeErrorMessage("Singular matrix : info="+Util::toString(info) );
	}

}

// Factrize and solve linear equation with real symmetric matrix
// @note Matrix and vectors are overwritten
void LapackInterface::factorizeAndSolveLinearEquationRealSymmetricMatrix( const int dimension, const int nRhs, double* matrix, double* vectors ){

	integer n = static_cast<integer>( dimension );
	integer nb = static_cast<integer>( nRhs );
	integer lda = n;
	integer* ipiv = new integer[dimension];
	integer ldb = n;
	integer lwork = n * nb;
	double* work = new double[lwork];
	integer info(0);

	dsysv_("L", &n, &nb, matrix, &lda, ipiv, vectors, &ldb, work, &lwork, &info);
	delete [] ipiv;
	delete [] work;

	if( info < 0 ){
		OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
		ptrOutputFiles->writeErrorMessage("An argument had an illegal value : info="+Util::toString(info) );
	}else if( info > 0 ){
		OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
		ptrOutputFiles->writeErrorMessage("D(i,i) is exactly zero : info="+Util::toString(info) );
	}

}

// Factrize and solve linear equation with real symmetric positive definite matrix
// @note Matrix and vectors are overwritten
void LapackInterface::factorizeAndSolveLinearEquationRealSymmetricPositiveDefiniteMatrix( const int dimension, const int nRhs, double* matrix, double* vectors  ){

	integer n = static_cast<integer>( dimension );
	integer nb = static_cast<integer>( nRhs );
	integer lda = n;
	integer ldb = n;
	integer info(0);

	dposv_("L", &n, &nb, matrix, &lda, vectors, &ldb, &info);

	if( info < 0 ){
		OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
		ptrOutputFiles->writeErrorMessage("An argument had an illegal value : info="+Util::toString(info) );
	}else if( info > 0 ){
		OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
		ptrOutputFiles->writeErrorMessage("D(i,i) is exactly zero : info="+Util::toString(info) );
	}

}

// Factorize a real square matrix
// @note Matrix and vectors are overwritten
void LapackInterface::factorizeRealSquareMatrix( const int dimension, double* matrix, int* ipivInt ){

	integer n = static_cast<integer>(dimension);
	integer m = n;
	integer lda = m;
	integer info(0);
	integer* ipiv = new integer[dimension];

    dgetrf_( &m, &n, matrix, &lda, ipiv, &info );
	for( int i = 0; i < dimension; ++i ){
		ipivInt[i] = ipiv[i];
	}
	delete [] ipiv;

	if( info < 0 ){
		OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
		ptrOutputFiles->writeErrorMessage("An argument had an illegal value : info="+Util::toString(info) );
	}else if( info > 0 ){
		OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
		ptrOutputFiles->writeErrorMessage(" D(i,i) is exactly zero : info="+Util::toString(info) );
	}

}

// Factorize a real symmetric matrix
// @note Matrix and vectors are overwritten
void LapackInterface::factorizeRealSymmetricMatrix( const int dimension, double* matrix, int* ipivInt ){

	integer n = static_cast<integer>(dimension);
	integer dum1 = 1;
	integer dum2 = -1;
	integer nb = ilaenv_( &dum1, "DSYTRF", "L", &n, &dum2, &dum2, &dum2 );
	integer lda = n;
	integer* ipiv = new integer[dimension];
	integer lwork = n * nb;
	double* work = new double[lwork];
	integer info(0);

	dsytrf_("L", &n, matrix, &lda, ipiv, work, &lwork, &info);
	for( int i = 0; i < dimension; ++i ){
		ipivInt[i] = ipiv[i];
	}
	delete [] ipiv;
	delete [] work;

	if( info < 0 ){
		OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
		ptrOutputFiles->writeErrorMessage("An argument had an illegal value : info="+Util::toString(info) );
	}else if( info > 0 ){
		OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
		ptrOutputFiles->writeErrorMessage(" D(i,i) is exactly zero : info="+Util::toString(info) );
	}

}

// Factorize a real symmetric positive definite matrix
// @note Matrix and vectors are overwritten
void LapackInterface::factorizeRealSymmetricPositiveDefiniteMatrix( const int dimension, double* matrix ){

	integer n = static_cast<integer>(dimension);
	integer lda = n;
	integer info(0);

	dpotrf_("L", &n, matrix, &lda, &info);

	if( info < 0 ){
		OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
		ptrOutputFiles->writeErrorMessage("An argument had an illegal value : info="+Util::toString(info) );
	}else if( info > 0 ){
		OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
		ptrOutputFiles->writeErrorMessage(" D(i,i) is exactly zero : info="+Util::toString(info) );
	}

}

// Solve a linear equation with real square matrix
void LapackInterface::solveLinearEquationRealSquareMatrix( const int dimension, const int nRhs, const int* const ipivInt, double* factorizedMatrix, double* vectors ){

	integer n = static_cast<integer>(dimension);
	integer nb = static_cast<integer>(nRhs);
	integer lda = n;
	integer* ipiv = new integer[dimension];
	for( int i = 0; i < dimension; ++i ){
		ipiv[i] = ipivInt[i];
	}
	integer ldb = n;
	integer info(0);

	dgetrs_("N", &n, &nb, factorizedMatrix, &lda, ipiv, vectors, &ldb, &info);
	delete [] ipiv;

	if( info < 0 ){
		OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
		ptrOutputFiles->writeErrorMessage("An argument had an illegal value : info="+Util::toString(info) );
	}

}

// Solve a linear equation with real symmetric matrix
void LapackInterface::solveLinearEquationRealSymmetricMatrix( const int dimension, const int nRhs, const int* const ipivInt, double* factorizedMatrix, double* vectors ){

	integer n = static_cast<integer>(dimension);
	integer nb = static_cast<integer>(nRhs);
	integer lda = n;
	integer* ipiv = new integer[dimension];
	for( int i = 0; i < dimension; ++i ){
		ipiv[i] = ipivInt[i];
	}
	integer ldb = n;
	integer info(0);

	dsytrs_("L", &n, &nb, factorizedMatrix, &lda, ipiv, vectors, &ldb, &info);
	delete [] ipiv;

	if( info < 0 ){
		OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
		ptrOutputFiles->writeErrorMessage("An argument had an illegal value : info="+Util::toString(info) );
	}

}

// Solve a linear equation with real symmetric positive definite matrix
void LapackInterface::solveLinearEquationRealSymmetricPositiveDefiniteMatrix( const int dimension, const int nRhs, double* factorizedMatrix, double* vectors ){

	integer n = static_cast<integer>(dimension);
	integer nb = static_cast<integer>(nRhs);
	integer lda = n;
	integer* ipiv = new integer[dimension];
	integer ldb = n;
	integer info(0);

	dpotrs_("L", &n, &nb, factorizedMatrix, &lda, vectors, &ldb, &info);

	if( info < 0 ){
		OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
		ptrOutputFiles->writeErrorMessage("An argument had an illegal value : info="+Util::toString(info) );
	}

}
