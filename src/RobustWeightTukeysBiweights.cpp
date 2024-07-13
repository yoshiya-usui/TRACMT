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
#include "RobustWeightTukeysBiweights.h"
#include "OutputFiles.h"
#include "Util.h"
#include "TableOfTukeysBiweightParameters.h"

#include <assert.h>

// Default constructer
RobustWeightTukeysBiweights::RobustWeightTukeysBiweights():
	m_degreesOfFreedom(2)
{
	setNameOfRobustWeight("Tukey's biweights");
}

// Calculate sum of the second order derivative of loss function
double RobustWeightTukeysBiweights::calculateSumOf2ndOrderDerivativeOfLossFunction(const int numSegments, const std::complex<double>* const residuals, const double scale) const {

	double paramB(0.0);
	double paramC(0.0);
	calculateParams(getDegreesOfFreedom(), numSegments, paramB, paramC);

	double sumOfDerivatives(0.0);
	for (int iSeg = 0; iSeg < numSegments; ++iSeg){
		const double x = std::abs(residuals[iSeg]) / scale;
		sumOfDerivatives += calculateSecondDerivativeOfLossFunction(x, paramC);
	}
	return sumOfDerivatives;

}

// Calculate weights
double RobustWeightTukeysBiweights::calculateWeights(const int numSegments, const std::complex<double>* const residuals, const double scale,
	const double* const weightsPrior, double* weights) const {

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeCvgMessage("	Scale factor: " + Util::toString(scale));
	double paramB(0.0);
	double paramC(0.0);
	calculateParams(getDegreesOfFreedom(), numSegments, paramB, paramC);
	ptrOutputFiles->writeCvgMessage("	Parameter c: " + Util::toString(paramC));
	double sumOfWeights(0.0);
	for (int iSeg = 0; iSeg < numSegments; ++iSeg) {
		weights[iSeg] = weightsPrior[iSeg];
		const double x = std::abs(residuals[iSeg]) / scale;
		weights[iSeg] *= calculateWeights(x, paramC);
		sumOfWeights += weights[iSeg];
	}

	return sumOfWeights;

}
// Destructer
RobustWeightTukeysBiweights::~RobustWeightTukeysBiweights()
{
}

// Calculate scale of absolute values of complex residuals
double RobustWeightTukeysBiweights::calculateScale(const int numSegments, const std::complex<double>* const residuals, const double scalePre) const {

	double paramB(0.0);
	double paramC(0.0);
	calculateParams(getDegreesOfFreedom(), numSegments, paramB, paramC);

	double* realResiduals = new double[numSegments];
	for (int iSeg = 0; iSeg < numSegments; ++iSeg) {
		realResiduals[iSeg] = std::abs(residuals[iSeg]);
	}
	const double scaleNew =  calculateRobustScale(numSegments, realResiduals, scalePre, paramB, paramC);
	delete[] realResiduals;

	return scaleNew;

}

// Get degrees of freedom
int RobustWeightTukeysBiweights::getDegreesOfFreedom() const{
	return m_degreesOfFreedom;
}

// Set degrees of freedom
void RobustWeightTukeysBiweights::setDegreesOfFreedom( const int dof ){
	m_degreesOfFreedom = dof;
}

// Calculate derivative of the loss function of Tukey's biweights
double RobustWeightTukeysBiweights::calculateDerivativeOfLossFunction(const double val, const double c) {

	assert(c > 0.0);
	if (fabs(val) > c) {
		return 0.0;
	}
	else {
		const double temp = 1.0 - pow(val / c, 2);
		return val * pow(temp, 2);
	}

}

// Calculate loss function of Tukey's biweight function
double RobustWeightTukeysBiweights::calculateLossFunction(const double val, const double c) {

	assert(c > 0.0);
	if (fabs(val) > c) {
		// c^2 / 6
		return pow(c, 2) / 6;
	}
	else {
		// x^2 / 2 - x^4 / c^2 / 2 + x^6 / c^4 / 6
		return pow(val, 2) / 2.0 - pow(val, 4) / pow(c, 2) / 2.0 + pow(val, 6) / pow(c, 4) / 6.0;
	}

}

// Calculate expectation of Tukey's biweight function
double RobustWeightTukeysBiweights::calculateExpectation(const int dimension, const double c) {

	assert(dimension > 0);

	double* arguments = new double[dimension];
	double* factors = new double[dimension];
	double sum(0.0);
	calculateExpectationAux(dimension, c, 0, arguments, factors, sum);
	delete[] arguments;
	delete[] factors;
	return sum / pow(sqrt(2.0 * CommonParameters::PI), dimension);

}

// Auxiliary function for calculating the expectation of Tukey's biweight function
void RobustWeightTukeysBiweights::calculateExpectationAux(const int dimension, const double c, const int idim,
	double* arguments, double* factors, double& sum) {

	const double lowerLimit = -5.0;
	const double upperLimit = 5.0;
	const int numInterval = 20;
	const double interval = (upperLimit - lowerLimit) / static_cast<double>(numInterval);
	for (int iInterval = 0; iInterval < numInterval + 1; ++iInterval) {
		arguments[idim] = lowerLimit + static_cast<int>(iInterval) * interval;
		if (iInterval == 0 || iInterval == numInterval) {
			factors[idim] = 0.5 * interval;
		}
		else {
			factors[idim] = interval;
		}
		if (idim == dimension - 1) {
			double squareSumOfArguments(0.0);
			for (int i = 0; i < dimension; ++i) {
				squareSumOfArguments += pow(arguments[i], 2);
			}
			double temp = calculateLossFunction(sqrt(squareSumOfArguments), c) * exp(-0.5 * squareSumOfArguments);
			for (int i = 0; i < dimension; ++i) {
				temp *= factors[i];
			}
			sum += temp;
		}
		else {
			// Go to next inner loop
			calculateExpectationAux(dimension, c, idim + 1, arguments, factors, sum);
		}
	}

}

// Calculate parameters of Tukey's biweight function
void RobustWeightTukeysBiweights::calculateParams(const int dimension, const int numOfData, double& paramB, double& paramC) {

	assert(dimension > 0);
	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();

	if (dimension > TableOfTukeysBiweightParameters::numDimensions) {
		ptrOutputFiles->writeErrorMessage("Number of dimensions for Tukey's biweight function exceeds the upper limit ("
			+ Util::toString(TableOfTukeysBiweightParameters::numDimensions) + ")");
	}

	const int idim = dimension - 1;

	const double rhs = 0.5 * (1.0 - static_cast<double>(dimension) / static_cast<double>(numOfData));// (n-p)/2n
	for (int ic = 1; ic < TableOfTukeysBiweightParameters::numOfParamCs; ++ic) {
		const double c0 = TableOfTukeysBiweightParameters::paramCsOfTukeyBiweight[ic - 1];
		const double b0 = TableOfTukeysBiweightParameters::paramBsOfTukeyBiweight[idim][ic - 1];
		const double lhs0 = ic == 1 ? 1.0 : 6.0 * b0 / pow(c0, 2);// b / (c^2/6)
		const double diff0 = lhs0 - rhs;
		const double c1 = TableOfTukeysBiweightParameters::paramCsOfTukeyBiweight[ic];
		const double b1 = TableOfTukeysBiweightParameters::paramBsOfTukeyBiweight[idim][ic];
		const double lhs1 = 6.0 * b1 / pow(c1, 2);// b / (c^2/6)
		const double diff1 = lhs1 - rhs;
		// lhs monotonically decreases
		if (diff0 >= 0.0 && diff1 <= 0.0) {
			const double weight = fabs(diff0) / fabs(lhs1 - lhs0);
			paramC = (1.0 - weight) * c0 + weight * c1;
			paramB = (1.0 - weight) * b0 + weight * b1;
			return;
		}
	}

	ptrOutputFiles->writeErrorMessage("Parameters of Tukey's biweight function was not able to be determined");

}

// Calculate robust scale with Tukey's biweights
double RobustWeightTukeysBiweights::calculateRobustScale(const int numOfData, const double* residuals, const double scalePre,
	const double paramb, const double paramc) {

	double squareScale(0.0);
	for (int iData = 0; iData < numOfData; ++iData) {
		const double val = residuals[iData] / scalePre;
		if (fabs(val) < CommonParameters::EPS) {
			// rho''(0) / 2 = 1 / 2
			squareScale += 0.5 * pow(residuals[iData], 2);
		}
		else {
			squareScale += calculateLossFunction(val, paramc) * pow(scalePre, 2);
		}
	}
	squareScale /= static_cast<double>(numOfData);
	squareScale /= paramb;
	return sqrt(squareScale);

}

// Calculate second order derivative of the loss function of Tukey's biweights
double RobustWeightTukeysBiweights::calculateSecondDerivativeOfLossFunction(const double val, const double c) {

	assert(c > 0.0);
	if (fabs(val) > c) {
		return 0.0;
	}
	else {
		// 1 - 6 * x^2 / c^2 + 5 * x^4 / c^4
		return 1.0 - 6.0 * pow(val / c, 2) + 5.0 * pow(val / c, 4);
	}

}

// Calculate the term in the denominator of robust covariance matrix using Tukey's biweights
double RobustWeightTukeysBiweights::calculateTermInDenominatorOfRobustCovariance(const double val, const double c) {

	return calculateWeights(val, c) * pow(val, 2);

}

// Calculate Tukey's biweights
double RobustWeightTukeysBiweights::calculateWeights(const double val, const double c) {

	assert(c > 0.0);
	if (fabs(val) > c) {
		return 0.0;
	}
	else {
		// 1 - 2 * x^2 / c^2 + x^4 / c^4
		return 1.0 - 2.0 * pow(val / c, 2) + pow(val / c, 4);
	}

}