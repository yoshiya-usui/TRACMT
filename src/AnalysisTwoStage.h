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
#ifndef DBLDEF_ANALYSIS_TWO_STAGE
#define DBLDEF_ANALYSIS_TWO_STAGE

#include "Analysis.h"

// Class of analysis for two stage remote reference approach
class AnalysisTwoStage : public Analysis {

public:

	// Default constructer
	AnalysisTwoStage();

	// Destructer
	virtual ~AnalysisTwoStage();

private:

	// Calculate the maximum value of hat matrix diagonals
	double calculateMaximumValueOfHatMatrixDiagonals(  ) const;

	// Calculate response functions
	virtual void calculateResponseFunctions( const int iSegLen,const int freqDegree, const double timeLength, const double freq,
		const int numSegmentsTotal, std::complex<double>** ftval, const std::vector< std::pair<std::string, std::string> >& times,
		std::ofstream& ofsResp, std::ofstream& ofsRhoaPhs );

	void calculateResponseFunctionsAux( const int iSegLen,const int freqDegree, const double timeLength, const double freq,
		const int numSegmentsTotal, std::complex<double>** ftval, 	const std::vector< std::pair<std::string, std::string> >& times, 
		std::ofstream& ofsResp, std::ofstream& ofsRhoaPhs, const bool forJackknife, std::complex<double>* respOut0, std::complex<double>* respOut1, double** hatDiagsOut ) const;

	// Calculate response functions by iteratively reweighted least squares for the first stage
	void calculateResponseFunctionsByIRWLSForFirstStage( const int iRobustWeight, const int inputVariableIndex, 
		const int numSegments, const bool fixScale, double& scale, std::complex<double>** data, const double* const leverageWeights,
		double* weights, std::complex<double>* residuals, std::complex<double>& resp0, std::complex<double>& resp1, 
		double& coherence, std::vector<std::string>& titles, std::vector<double>* outputValues ) const;

	// Calculate response functions by iteratively reweighted least squares for the second stage
	void calculateResponseFunctionsByIRWLSForSecondStage( const int iRobustWeight, const int outputVariableIndex, 
		const int numSegments, const bool fixScale, double& scale, std::complex<double>** data, std::complex<double>** dataSyn, const double* const leverageWeights,
		double* weights, std::complex<double>* residuals, std::complex<double>& resp0, std::complex<double>& resp1,
		double& coherence, std::vector<std::string>& titles, std::vector<double>* outputValues ) const;
	
	// Calculate response function by the weighted least square method for the first stage
	void calculateResponseFunctionsByWLSForFirstStage( const int inputVariableIndex, const int numSegments, std::complex<double>** data, 
		 double* weights, std::complex<double>* residuals, std::complex<double>& resp0, std::complex<double>& resp1, double& coherence ) const;

	// Calculate response function by the weighted least square method for the second stage
	void calculateResponseFunctionsByWLSForSecondStage( const int outputVariableIndex, const int numSegments, std::complex<double>** data, 
		std::complex<double>** dataSyn, double* weights, std::complex<double>* residuals, std::complex<double>& resp0, std::complex<double>& resp1,
		double& coherence ) const;

	// Perform fixed-weights jackknife
	void fixedWeightsJackknife( const double freq, const int numSegmentsTotal, double** weightsOrg, std::complex<double>** ftval,
		std::complex<double>** dataSyn, std::ofstream& ofsResp, std::ofstream& ofsRhoaPhs, 
		const std::complex<double>* const resp0, const std::complex<double>* const resp1, double** hatDiagonals ) const;

	// Perform fixed-weights bootstrap
	void fixedWeightsBootstrap( const double freq, const int numSegmentsTotal, double** weightsOrg, std::complex<double>** ftval,
		std::complex<double>** dataSyn, std::ofstream& ofsResp, std::ofstream& ofsRhoaPhs, 
		const std::complex<double>* const resp0, const std::complex<double>* const resp1 ) const;

	// Output synthetic input data
	void outputSyntheticInputData( const int iSegLen, const int freqDegree, const int numSegmentsTotal, std::complex<double>** dataSyn ) const;

	// Estimate error by strict bootstrap
	void strictBootstrap( const int iSegLen,const int freqDegree, const double timeLength, const double freq,
		const int numSegmentsTotal, std::complex<double>** ftval, const std::vector< std::pair<std::string, std::string> >& times,
		std::ofstream& ofsResp, std::ofstream& ofsRhoaPhs, const std::complex<double>* const resp0, const std::complex<double>* const resp1 ) const;

	// Perform subset deletion jackknife
	void subsetDeletionJackknife( const int iSegLen,const int freqDegree, const double timeLength, const double freq,
		const int numSegmentsTotal, std::complex<double>** ftval, const std::vector< std::pair<std::string, std::string> >& times,
		std::ofstream& ofsResp, std::ofstream& ofsRhoaPhs, const std::complex<double>* const resp0, const std::complex<double>* const resp1, double** hatDiagonals ) const;

	// Write residuals
	void writeResiduals( const std::string& fileName, const int numSegmentsTotal,
		const std::vector< std::pair<std::string, std::string> >& times, const std::vector<std::string>& titles, const std::vector<double>* outputValues ) const;

};

#endif
