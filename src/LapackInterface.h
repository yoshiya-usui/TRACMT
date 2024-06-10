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
#ifndef DBLDEF_LAPACK_INTERFACE
#define DBLDEF_LAPACK_INTERFACE

namespace LapackInterface
{

// Calcualte all eigenvalues and eigenvectors of a real symmetric matrix
// @note Matrix and vectors are overwritten
void calculateEigenValuesAndVectorsOfRealSymmetricMatrix( const int dimension, double* matrix, double* vectors );

// Factorize and solve a linear equation with real coefficents matrix
// @note Matrix and vectors are overwritten
void factorizeAndSolveLinearEquationRealMatrix( const int dimension, const int nRhs, double* matrix, double* vectors );

// Factorize and solve a linear equation with real symmetric matrix
// @note Matrix and vectors are overwritten
void factorizeAndSolveLinearEquationRealSymmetricMatrix( const int dimension, const int nRhs, double* matrix, double* vectors );

// Factorize and solve a linear equation with real symmetric positive definite matrix
// @note Matrix and vectors are overwritten
void factorizeAndSolveLinearEquationRealSymmetricPositiveDefiniteMatrix( const int dimension, const int nRhs, double* matrix, double* vectors );

// Factorize a real square matrix
// @note Matrix and vectors are overwritten
void factorizeRealSquareMatrix( const int dimension, double* matrix, int* ipivInt );

// Factorize a real symmetric matrix
// @note Matrix and vectors are overwritten
void factorizeRealSymmetricMatrix( const int dimension, double* matrix, int* ipivInt );

// Factorize a real symmetric positive definite matrix
// @note Matrix and vectors are overwritten
void factorizeRealSymmetricPositiveDefiniteMatrix( const int dimension, double* matrix );

// Solve a linear equation with real square matrix
void solveLinearEquationRealSquareMatrix( const int dimension, const int nRhs, const int* const ipivInt, double* factorizedMatrix, double* vectors );

// Solve a linear equation with real symmetric matrix
void solveLinearEquationRealSymmetricMatrix( const int dimension, const int nRhs, const int* const ipivInt, double* factorizedMatrix, double* vectors );

// Solve a linear equation with real symmetric positive definite matrix
void solveLinearEquationRealSymmetricPositiveDefiniteMatrix( const int dimension, const int nRhs, double* factorizedMatrix, double* vectors );

}

#endif
