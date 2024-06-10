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
#ifndef DBLDEF_ROBUST_PREWHITENING
#define DBLDEF_ROBUST_PREWHITENING

#include "CommonParameters.h"
#include "DoubleDenseSquareSymmetricPositiveDefiniteMatrix.h"

#include <vector>

// Class of robust prewhitening 
class RobustPrewhitening{

public:

	// Default constructer
	RobustPrewhitening();

	// Destructer
	~RobustPrewhitening();

	// Return the the instance of the class
    static RobustPrewhitening* getInstance();

	// Perform robust prewhitening
	void robustPrewhitening( std::vector<CommonParameters::DataFileSet>& dataFileSets, std::vector<double>* coeffsAROutput ) const;

	// Perform prewhitening using user-defined AR coefficients
	void prewhiteningUsingUserDefinedARCoeffs( std::vector<CommonParameters::DataFileSet>& dataFileSets, std::vector<double>* coeffsAROutput ) const;

private:

	// Copy constructer
	RobustPrewhitening(const RobustPrewhitening& rhs);

	// Assignment operator
	RobustPrewhitening& operator=(const RobustPrewhitening& rhs);

	// Calculate robust-filtered value
	void calculateRobustFilteredValue( const int iChan, const int numOfData, const int degreesOfAR, 	const double* const coeffsOfAR,
		const double sigma, const double* const yOrg, const double* const autoCovariance, double* yMod ) const;

	// Calculate robust auto-covariance matrix
	void calculateRobustAutoCovarianceMatrix( const int degreesOfAR, const int numOfDataSets, const int* const numOfData,
		double** data, DoubleDenseSquareSymmetricMatrix& covarianceMatrix ) const;

	// Calculate Mahalanobis distances
	// @note covariance matrix is factorized in this function
	//void calculateMD( const int degreesOfAR, const int numOfDataAll, const double* const dataAll,
	//	DoubleDenseSquareSymmetricMatrix& covarianceMatrix, double* residualVector, double* MD ) const;
	void calculateMD( const int degreesOfAR, const int numOfDataSets, const int* const numOfData, double** data,
		DoubleDenseSquareSymmetricMatrix& covarianceMatrix, double* MD ) const;

	// Calculate weight for state vector in robust filter
	double calculateWeightForStateVectorOfRobustFilter( const double val, const double paramA, const double paramB, const bool residualAssumedToBeZero ) const;

	// Calculate weight for covariance matrix in robust filter
	double calculateWeightForCovarianceMatrixOfRobustFilter( const double val, const double paramA, const double paramB, const bool residualAssumedToBeZero  ) const;
};


#endif
