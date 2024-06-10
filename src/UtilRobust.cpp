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
#include "UtilRobust.h"
#include "Util.h"
#include "OutputFiles.h"
#include "CommonParameters.h"
#include <algorithm>

#include <assert.h>
#include <iomanip>
#include <iostream>

// Calculate bisquare weights
double UtilRobust::calculateBisquareWeights( const double val, const double c ){

	assert( c > 0.0 );
	if( fabs(val) > c ){
		return 0.0;
	}else{
		const double temp = 1.0 - pow(val/c, 2);
		return pow(temp, 2);
	}

}

// Calculate loss function of bisquare weights
double UtilRobust::calculateLossFunctionOfBisquareWeights( const double val, const double c ){

	assert( c > 0.0 );
	if( fabs(val) > c ){
		return 1.0;
	}else{
		const double temp = 1.0 - pow(val/c, 2);
		return 1.0 - pow(temp, 3);
	}

}

// Calculate derivative of the loss function of bisquare weights
double UtilRobust::calculateDerivativeOfLossFunctionOfBisquareWeights( const double val, const double c ){

	assert( c > 0.0 );
	if( fabs(val) > c ){
		return 0.0;
	}else{
		const double temp = 1.0 - pow(val/c, 2);
		return 6.0 * val / pow(c, 2) * (temp, 2);
	}

}

// Calculate M-estimator of scale
double UtilRobust::calculateMEstimatorOfScale( const int numOfData, const double* const residuals ){

	double scale = Util::calculateMADN(numOfData, residuals);
	const double paramb = 0.5;
	const double paramc = 1.56;
	const int numOfISteps = 20;
	const double convergenceCriteria = 0.05;
	for( int iter = 0; iter < numOfISteps; ++iter ){
		const double scaleNew = calculateRobustScaleWithBisquareWeights(numOfData, residuals, scale, paramb, paramc);
		if( fabs(scale - scaleNew)/fabs(scale) < convergenceCriteria ){
			// Converged
			return scaleNew;
		}
		scale = scaleNew;
	}

	return scale;

}

// Calculate robust scale with bisquare weights
double UtilRobust::calculateRobustScaleWithBisquareWeights( const int numOfData, const double* residuals, const double scale, 
	const double paramb, const double paramc ){

	double squareScale(0.0);
	for( int iData = 0; iData < numOfData; ++iData ){
		const double u = residuals[iData] / scale;
		if( fabs(u) < CommonParameters::EPS ){
			// rho''(0) / 2 = 6 / c^2 / 2 = 3 / c^2
			squareScale += 3.0 / pow(paramc, 2) * pow(residuals[iData], 2);
		}else{
			squareScale += calculateLossFunctionOfBisquareWeights(residuals[iData]/scale, paramc) * pow(scale, 2);
		}
	}
	squareScale /= static_cast<double>(numOfData);
	squareScale /= paramb;
	return sqrt(squareScale);

}

// Calculate second derivative of the loss function of bisquare weights
double UtilRobust::calculateSecondDerivativeOfLossFunctionOfBisquareWeights( const double val, const double c ){

	assert( c > 0.0 );
	if( fabs(val) > c ){
		return 0.0;
	}else{
		const double temp = 1.0 - pow(val/c, 2);
		const double term1 = 6.0 / pow(c, 2) * pow(temp, 2);
		const double term2 = 24.0 * pow(val, 2) / pow(c, 4) * temp;
		return term1 - term2;
	}

}

// Compute least square estimator for univariate linear regression
void UtilRobust::computeLSEstimatorForUnivariateLinearRegression( const int numOfData,	const double* const y, const double* const x,
	const double coeffMin, const double coeffMax, double& scaleOut, double& coeffOut ){

	double xx = 0.0;
	double xy = 0.0;
	for( int iData = 0; iData < numOfData; ++iData ){
		xx += x[iData] * x[iData]; 
		xy += x[iData] * y[iData]; 
	}

	if( fabs(xx) < CommonParameters::EPS ){
		OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
		ptrOutputFiles->writeWarningMessage("Denominator is too small ("+ Util::toString(xx) + ") in the weighted least square method");
		coeffOut = xy / CommonParameters::EPS;
	}else{
		coeffOut = xy / xx;
	}

	if( coeffOut < coeffMin ){
		coeffOut = coeffMin;
	}
	if( coeffOut > coeffMax ){
		coeffOut = coeffMax;
	}

	double* residuals = new double[numOfData];
	for( int iData = 0; iData < numOfData; ++iData ){
		residuals[iData] = y[iData] - coeffOut * x[iData];
	}
	const double variance = Util::calculateVariance(numOfData, 0.0, residuals);
	scaleOut= sqrt(variance);
	delete [] residuals;

}

// Compute S-estimator for univariate linear regression
// This function is based on the algorithm of Salibian-Barrera & Yohai (2006)
void UtilRobust::computeSEstimatorForUnivariateLinearRegression( const int numOfData, const int numOfCandidates,
	const double* const y, const double* const x, const double* const coeffCandidates, const double coeffMin, const double coeffMax,
	double& scaleOut, double& coeffOut, double& AICS ){

	const double paramb = 0.5;
	const double paramc = 1.56;
	const int numOfISteps = 20;
	// Initialize
	scaleOut = 1.0e20;
	coeffOut = 0.0;
	double* residuals = new double[numOfData];
	for( int iCan = 0; iCan < numOfCandidates; ++iCan ){
		double coeff = coeffCandidates[iCan];
		for( int iData = 0; iData < numOfData; ++iData ){
			residuals[iData] = y[iData] - coeff * x[iData];
		}
		double scale = Util::calculateMADN(numOfData, residuals);
		performISteps( numOfData, y, x, paramb, paramc, numOfISteps, coeffMin, coeffMax, coeff, scale, residuals );
		if( scale < scaleOut ){
			scaleOut = scale;
			coeffOut = coeff;
		}
	}

	double Jsn(0.0);
	double Ksn(0.0);
	for( int iData = 0; iData < numOfData; ++iData ){
		const double normalizedResidual = residuals[iData] / scaleOut;
		Jsn += calculateSecondDerivativeOfLossFunctionOfBisquareWeights(normalizedResidual, paramc) * pow(x[iData]/scaleOut,2);
		Ksn += pow(calculateDerivativeOfLossFunctionOfBisquareWeights(normalizedResidual, paramc), 2) * pow(x[iData]/scaleOut,2);
	}
	// 1/n is omiteed because it is canceled out

	const double term1 = 2.0 * static_cast<double>(numOfData) * log(scaleOut);
	double term2(0.0);
	if( fabs(Jsn) > CommonParameters::EPS ){
		const double term2 = 2.0 * Ksn / Jsn;
	}

	AICS = term1 + term2;

	delete [] residuals;

}

// Compute S-estimator for univariate linear regression with intercept
void UtilRobust::computeSEstimatorForUnivariateLinearRegressionWithIntercept(const int numOfData, const int numOfCandidates,
	const double* const y, const double* const x, const double* const slopeCandidates, const double* const interceptCandidates,
	double& scaleOut, double& slopeOut, double& interceptOut) {

	const double paramb = 0.5;
	const double paramc = 1.56;
	const int numOfISteps = 20;
	// Initialize
	scaleOut = 1.0e20;
	slopeOut = 0.0;
	interceptOut = 0.0;
	double* residuals = new double[numOfData];
	for (int iCan = 0; iCan < numOfCandidates; ++iCan) {
		double slope = slopeCandidates[iCan];
		double intercept = interceptCandidates[iCan];
		for (int iData = 0; iData < numOfData; ++iData) {
			residuals[iData] = y[iData] - slope * x[iData] - intercept;
		}
		double scale = Util::calculateMADN(numOfData, residuals);
		performIStepsWithIntercept(numOfData, y, x, paramb, paramc, numOfISteps, slope, intercept, scale, residuals);
#ifdef _DEBUG_WRITE
		std::cout << "iCan " << iCan << " " << slope << " " << intercept << " " << scale << std::endl;
#endif
		if (scale < scaleOut) {
			slopeOut = slope;
			interceptOut = intercept;
			scaleOut = scale;
		}
	}

	delete[] residuals;
}

// Perform I-steps for calculating S-estimator for univariate linear regression
void UtilRobust::performISteps( const int numOfData, const double* const y, const double* const x,
	const double paramb, const double paramc, const int numOfISteps, const double coeffMin, const double coeffMax, 
	double& coeff, double& scale, double* residuals ){

	double* weights = new double[numOfData];
	const double convergenceCriteriaCoeff = 0.01;
	const double convergenceCriteriaScale = 0.05;
	double coeffPre = coeff;
	double scalePre = scale;
	for( int iter = 0; iter < numOfISteps; ++iter ){
		scale = calculateRobustScaleWithBisquareWeights(numOfData, residuals, scale, paramb, paramc);
		for( int iData = 0; iData < numOfData; ++iData ){
			weights[iData] = calculateBisquareWeights( residuals[iData]/scale, paramc );
		}
		coeff = Util::calculateRegressionCoefficientsByWLS( numOfData, y, x, weights );
		if( coeff < coeffMin ){
			coeff = coeffMin;
		}
		if( coeff > coeffMax ){
			coeff = coeffMax;
		}
		for( int iData = 0; iData < numOfData; ++iData ){
			residuals[iData] = y[iData] - coeff * x[iData];
		}
		if( fabs(coeff - coeffPre)/fabs(coeffPre) < convergenceCriteriaCoeff && 
			fabs(scale - scalePre)/fabs(scalePre) < convergenceCriteriaScale ){
			// Converged
			break;
		}
		coeffPre = coeff;
		scalePre = scale;
	}
	delete [] weights;

}

// Perform I-steps for calculating S-estimator for univariate linear regression with intercept
void UtilRobust::performIStepsWithIntercept(const int numOfData, const double* const y, const double* const x,
	const double paramb, const double paramc, const int numOfISteps, double& slope, double& intercept, double& scale, double* residuals) {

	double* weights = new double[numOfData];
	const double convergenceCriteria = 0.01;
	const double convergenceCriteriaScale = 0.05;
	double slopePre = slope;
	double interceptPre = intercept;
	double scalePre = scale;
	for (int iter = 0; iter < numOfISteps; ++iter) {
		scale = calculateRobustScaleWithBisquareWeights(numOfData, residuals, scale, paramb, paramc);
		for (int iData = 0; iData < numOfData; ++iData) {
			weights[iData] = calculateBisquareWeights(residuals[iData] / scale, paramc);
		}
		Util::calculateRegressionCoefficientsByWLSWithIntercept(numOfData, y, x, weights, slope, intercept);
		for (int iData = 0; iData < numOfData; ++iData) {
			residuals[iData] = y[iData] - slope * x[iData] - intercept;
		}
		if (fabs(slope - slopePre) / fabs(slopePre) < convergenceCriteria &&
			fabs(intercept - interceptPre) / fabs(interceptPre) < convergenceCriteria &&
			fabs(scale - scalePre) / fabs(scalePre) < convergenceCriteriaScale) {
			// Converged
			break;
		}
		slopePre = slope;
		interceptPre = intercept;
		scalePre = scale;
	}
	delete[] weights;

}
