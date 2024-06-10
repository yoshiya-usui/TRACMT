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
#ifndef DBLDEF_UTIL_ROBUST
#define DBLDEF_UTIL_ROBUST

namespace UtilRobust
{

	// Calculate bisquare weights
	double calculateBisquareWeights( const double val, const double c );

	// Calculate loss function of bisquare weights
	double calculateLossFunctionOfBisquareWeights( const double val, const double c );

	// Calculate derivative of the loss function of bisquare weights
	double calculateDerivativeOfLossFunctionOfBisquareWeights(const double val, const double c);

	// Calculate M-estimator of scale
	double calculateMEstimatorOfScale( const int numOfData, const double* const residuals );

	// Calculate robust scale with bisquare weights
	double calculateRobustScaleWithBisquareWeights( const int numOfData, const double* residuals, const double scale, 
		const double paramb, const double paramc );

	// Calculate second order derivative of the loss function of bisquare weights
	double calculateSecondDerivativeOfLossFunctionOfBisquareWeights( const double val, const double c );

	// Compute least square estimator for univariate linear regression
	void computeLSEstimatorForUnivariateLinearRegression( const int numOfData,	const double* const y, const double* const x,
		const double coeffMin, const double coeffMax, double& scaleOut, double& coeffOut );

	// Compute S-estimator for univariate linear regression
	void computeSEstimatorForUnivariateLinearRegression( const int numOfData, const int numOfCandidates,
		const double* const y, const double* const x, const double* const coeffCandidates, const double coeffMin, const double coeffMax,
		double& scaleOut, double& coeffOut, double& AICS );

	// Compute S-estimator for univariate linear regression
	void computeSEstimatorForUnivariateLinearRegressionWithIntercept(const int numOfData, const int numOfCandidates,
		const double* const y, const double* const x, const double* const slopeCandidates, const double* const interceptCandidates,
		double& scaleOut, double& slopeOut, double& interceptOut);

	// Perform I-steps for calculating S-estimator for univariate linear regression
	void performISteps( const int numOfData, const double* const y, const double* const x,
		const double paramb, const double paramc, const int numOfISteps, const double coeffMin, const double coeffMax,
		double& coeff, double& scale, double* residuals );

	// Perform I-steps for calculating S-estimator for univariate linear regression with intercept
	void performIStepsWithIntercept(const int numOfData, const double* const y, const double* const x,
		const double paramb, const double paramc, const int numOfISteps, double& slope, double& intercept, double& scale, double* residuals);

}
#endif
