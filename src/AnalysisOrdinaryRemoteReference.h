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
#ifndef DBLDEF_ANALYSIS_ORDINARY_REMOTE_REFERENCE
#define DBLDEF_ANALYSIS_ORDINARY_REMOTE_REFERENCE

#include "Analysis.h"

// Class of ordinary remote reference method
class AnalysisOrdinaryRemoteReference : public Analysis {

public:

	// Default constructer
	AnalysisOrdinaryRemoteReference();

	// Destructer
	virtual ~AnalysisOrdinaryRemoteReference();

private:

	// Calculate partial derivatives of responses for robust bootstrap
	void calculatePartialDerivativesOfResponses(const int numSegments, const double threshould, std::complex<double>** ftval,
		const std::complex<double> resp0, const std::complex<double> resp1, const double scale, const std::complex<double>* const residuals,
		const double* const weights, const int iOut, double** derivativesRegardingResps, double* derivativesRegardingScale) const;

	// Calculate partial derivatives of scale
	void calculatePartialDerivativesOfScale(const int numSegments, const double threshould, const double paramB, std::complex<double>** ftval, 
		const std::complex<double> resp0, const std::complex<double> resp1, const double scale, const std::complex<double>* const residuals,
		const double* const weights, const int iOut, double* derivativesRegardingResps, double& derivativesRegardingScale) const;

	// Calculate response functions
	virtual void calculateResponseFunctions( const int iSegLen,const int freqDegree, const double timeLength, const double freq,
		const int numSegmentsTotal, std::complex<double>** ftval, const std::vector< std::pair<std::string, std::string> >& times,
		std::ofstream& ofsResp, std::ofstream& ofsRhoaPhs );

	void calculateResponseFunctionsAux(const int iSegLen, const int freqDegree, const double timeLength, const double freq,
		const int numSegmentsTotal, std::complex<double>** ftval, const std::vector< std::pair<std::string, std::string> >& times,
		std::ofstream& ofsResp, std::ofstream& ofsRhoaPhs, const bool forJackknife, std::complex<double>* respOut0, std::complex<double>* respOut1, 
		double** hatDiagsOut) const;

	// Calculate a vector for partial derivatives
	void calculateVectorForPartialDerivatives(const int iSeg, std::complex<double>** ftval, const std::complex<double>* const residuals,
		std::complex<double>* hSigmaMatrix, double* vector) const;

	// Perform fixed-weights jackknife
	void fixedWeightsJackknife(const double freq, const int numSegmentsTotal, double** weightsOrg, std::complex<double>** ftval,
		std::ofstream& ofsResp, std::ofstream& ofsRhoaPhs, const std::complex<double>* const resp0, const std::complex<double>* const resp1, double** hatDiagonals) const;

	// Perform fixed-weights bootstrap
	void fixedWeightsBootstrap(const double freq, const int numSegmentsTotal, double** weights, std::complex<double>** ftval,
		std::ofstream& ofsResp, std::ofstream& ofsRhoaPhs, const std::complex<double>* const resp0, const std::complex<double>* const resp1) const;

	// Estimate errors by parametric method
	void parametricErrorEstimation(const int numSegmentsTotal, const double* const weights, std::complex<double>** ftval,
		const std::complex<double>* const residuals, const double scale, const bool noRobust, double& error0, double& error1) const;

	// Estimate error by fixed-weights bootstrap
	void fixedWeightsBootstrap(const double freq, const int numSegmentsTotal, double** weightsOrg, double* scalesOrg, std::complex<double>** ftval,
		std::ofstream& ofsResp, std::ofstream& ofsRhoaPhs, const std::complex<double>* const resp0Org, const std::complex<double>* const resp1Org) const;

	// Estimate errors by strict bootstrap
	void strictBootstrap(const int iSegLen, const int freqDegree, const double timeLength, const double freq,
		const int numSegmentsTotal, std::complex<double>** ftval, const std::vector< std::pair<std::string, std::string> >& times,
		std::ofstream& ofsResp, std::ofstream& ofsRhoaPhs, const std::complex<double>* const resp0, const std::complex<double>* const resp1 ) const;

	// Perform subset deletion jackknife
	void subsetDeletionJackknife(const int iSegLen, const int freqDegree, const double timeLength, const double freq,
		const int numSegmentsTotal, std::complex<double>** ftval, const std::vector< std::pair<std::string, std::string> >& times,
		std::ofstream& ofsResp, std::ofstream& ofsRhoaPhs, const std::complex<double>* const resp0, const std::complex<double>* const resp1, double** hatDiagonals) const;

	// Write residuals
	void writeResiduals(const std::string& fileName, const int numSegmentsTotal,
		const std::vector< std::pair<std::string, std::string> >& times, const std::vector<std::string>& titles, const std::vector<double>* outputValues) const;

};
#endif
