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
#ifndef DBLDEF_ANALYSIS_MULTIVARIATE_REGRESSION
#define DBLDEF_ANALYSIS_MULTIVARIATE_REGRESSION

#include "Analysis.h"

#include <vector>
#include <complex>

// Class of analysis using multivariate regression
class AnalysisMultivariateRegression : public Analysis {

public:

	// Default constructer
	AnalysisMultivariateRegression();

	// Destructer
	virtual ~AnalysisMultivariateRegression();

private:

	// Number of input and output variables
	int m_numOfOutputAndInputVariables;

	// Response functions between output/input variables and reference variables
	std::vector< std::complex<double> >** m_responseFunctions;

	// Frequencies
	std::vector<double> m_frequencies;

	// Calculate complex residuals
	void calculateComplexResiduals( const int numSegmentsTotal, std::complex<double>** ftval, 
		std::complex<double>** resp, std::complex<double>** complexResiduals ) const;

	// Calculate partial derivatives of responses for robust bootstrap
	void calculatePartialDerivativesOfResponses( const int numSegments, const double paramC,
		std::complex<double>** ftval, std::complex<double>** resp, const double* const variancesWithoutScale, const double scale, 
		std::complex<double>** complexResiduals, const double* const MD, const double* const weights,
		double** derivativesRegardingResps, double** derivativesRegardingVariancesWithoutScale, double* derivativesRegardingScale ) const;

	// Calculate partial derivatives of scale
	void calculatePartialDerivativesOfScale( const int numSegments, const double paramB, const double paramC,
		std::complex<double>** ftval, std::complex<double>** resp, const double* const variancesWithoutScale, const double scale, 
		std::complex<double>** complexResiduals, const double* const MD, const double* const weights,
		double* derivativesRegardingResps, double* derivativesRegardingVariancesWithoutScale, double& derivativesRegardingScale ) const;

	// Calculate partial derivatives of variances without scale for robust bootstrap
	void calculatePartialDerivativesOfVariancesWithoutScale( const int numSegments,	const double paramC,
		std::complex<double>** ftval, std::complex<double>** resp, const double* const variancesWithoutScale, const double scale, const double determinant,
		std::complex<double>** complexResiduals, const double* const MD, const double* const weights,
		double** derivativesRegardingResps, double** derivativesRegardingVariancesWithoutScale, double* derivativesRegardingScale ) const;

	// Calculate response functions
	virtual void calculateResponseFunctions( const int iSegLen, const int freqDegree, const double timeLength, const double freq,
		const int numSegmentsTotal, std::complex<double>** ftval, const std::vector< std::pair<std::string, std::string> >& times, 
		std::ofstream& ofsResp, std::ofstream& ofsRhoaPhs );

	// Calculate response functions from two samples
	bool calculateResponseFunctionsFromTwoSamples( 
		const std::complex<double>& out1,const std::complex<double>& x1, const std::complex<double>& y1, 
		const std::complex<double>& out2, const std::complex<double>& x2, const std::complex<double>& y2,
		std::complex<double>& resp1, std::complex<double>& resp2 ) const;

	// Calculate a vector for partial derivatives
	void calculateVectorForPartialDerivatives( const int iSeg, std::complex<double>** ftval, 
		const double* const variancesWithoutScale, std::complex<double>** complexResiduals,
		std::complex<double>** hSigmaMatrix, double* vector ) const;

	// Determine candidates
	void determineCandidates( const double freq, const int numSegmentsTotal, std::complex<double>** ftval, 
		const std::vector< std::pair<std::string, std::string> >& times, std::vector< std::pair<int,int> >& candidatesOfPairs, 
		std::vector< std::complex<double>** >& resp ) const;

	// Determine candidates for later frequencies
	void determineCandidatesForLaterFrequencies( const double freq, const int numSegmentsTotal, std::vector< std::complex<double>** >& resp ) const;

	// Determine candidates by random sampling
	void determineCandidatesByRandomSampling( const double freq, const int numSegmentsTotal, std::complex<double>** ftval, 
		const std::vector< std::pair<std::string, std::string> >& times, std::vector< std::pair<int,int> >& candidatesOfPairs, 
		std::vector< std::complex<double>** >& resp ) const;

	// Improve candidate by I-steps
	void improveCandidate( const int numSegmentsTotal, const int numOfMaxIterations, 
		const double convergenceCriteria, const bool initialCalculation, const double paramB, const double paramC,
		std::complex<double>** ftval, std::complex<double>** resp, double* variancesWithoutScale, double& scale, double& determinant,
		double* coherencesMin ) const;

	// Calculate Mahalanobis distances
	void calculateMD( const int numSegmentsTotal, 
		const int numOfOutputAndInputVariables, std::complex<double>** complexResiduals, const double* const variancesWithoutScale, 
		double* MD ) const;

	// Estimate error by fixed-weights bootstrap
	void estimateErrorByFixedWeightsBootstrap(const int numSegmentsTotal, const double paramB, const double paramC,
		std::complex<double>** ftval, std::complex<double>** respOrg, const double* const variancesWithoutScaleOrg, const double scaleOrg,
		const double determinantOrg, double** respErr) const;

	// Estimate error by strict bootstrap
	void estimateErrorByStrictBootstrap( const int numSegmentsTotal, const double paramB, const double paramC, 
		std::complex<double>** ftvalOrg, std::complex<double>** respOrg, const double* const variancesWithoutScaleOrg, const double scaleOrg,
		const double determinantOrg, double** respErr ) const;

	// Estimate error by fixed-weights jackknife
	void estimateErrorByFixedWeightsJackknife( const int numSegmentsTotal, const double paramB, const double paramC, 
		std::complex<double>** ftvalOrg, const std::complex<double>* const respOrg0, const std::complex<double>* const respOrg1,
		std::complex<double>** respOrg, const double* const variancesWithoutScaleOrg, const double scaleOrg, double** respErr ) const;

	// Estimate error by subset deletion jackknife
	void estimateErrorBySubsetDeletionJackknife( const int numSegmentsTotal, std::complex<double>** ftvalOrg, std::complex<double>** respOrg, 
		const double* const variancesWithoutScaleOrg, const double scaleOrg,const double determinantOrg,	double** respErr ) const;

	// Estimate error by a parametric approach
	void estimateErrorParametric( const int numSegmentsTotal, const double paramB, const double paramC,
		std::complex<double>** ftval, const std::complex<double>* const respOrg0, const std::complex<double>* const respOrg1,
		std::complex<double>** respOrg, const double* const variancesWithoutScaleOrg, const double scaleOrg, double** respErr ) const;

	// Write residuals
	void writeResiduals( const std::string& fileName, const int numSegmentsTotal, const int numOfOutputAndInputVariables,
		const std::vector< std::pair<std::string, std::string> >& times, std::complex<double>** complexResiduals, 
		const double* const MD, const double* const weights ) const;

};

#endif
