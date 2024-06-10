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
#ifndef DBLDEF_ROBUST_TUKEYS_BIWEIGHTS
#define DBLDEF_ROBUST_TUKEYS_BIWEIGHTS

#include "RobustWeight.h"
#include "CommonParameters.h"
#include <complex>

// Class of M-estimator using Tukey's biweights
class RobustWeightTukeysBiweights : public RobustWeight{

public:

	// Default constructer
	RobustWeightTukeysBiweights();

	// Destructer
	virtual ~RobustWeightTukeysBiweights();

	// Calculate sum of the second order derivative of loss function
	virtual double calculateSumOf2ndOrderDerivativeOfLossFunction(const int numSegments, const std::complex<double>* const residuals, const double scale) const;

	// Calculate weights
	virtual double calculateWeights(const int numSegments, const std::complex<double>* const residuals, const double scale, const double* const weightsPrior,
		double* weights) const;

	// Calculate scale of absolute values of complex residuals
	virtual double calculateScale(const int numSegments, const std::complex<double>* const residuals, const double scalePre) const;

	// Get degrees of freedom
	int getDegreesOfFreedom() const;

	// Set degrees of freedom
	void setDegreesOfFreedom( const int dof );

	// Calculate derivative of the loss function of Tukey's biweights
	static double calculateDerivativeOfLossFunction(const double val, const double c);

	// Calculate the expectation of Tukey's biweight function
	static double calculateExpectation(const int dimension, const double c);

	// Auxiliary function for calculating the expectation of Tukey's biweight function
	static void calculateExpectationAux(const int dimension, const double c, const int idim,
		double* arguments, double* factors, double& sum);

	// Calculate loss function of Tukey's biweight function
	static double calculateLossFunction(const double val, const double c);

	// Calculate parameters of Tukey's biweight function
	static void calculateParams(const int dimension, const int numOfData, double& paramB, double& paramC);

	// Calculate robust scale with Tukey's biweights
	static double calculateRobustScale(const int numOfData, const double* residuals, const double scalePre,
		const double paramb, const double paramc);

	// Calculate second order derivative of the loss function of Tukey's biweights
	static double calculateSecondDerivativeOfLossFunction(const double val, const double c);

	// Calculate the term in the denominator of robust covariance matrix using Tukey's biweights
	static double calculateTermInDenominatorOfRobustCovariance(const double val, const double c);

	// Calculate Tukey's biweights
	static double calculateWeights(const double val, const double c);

protected:

private:

	// Degrees of freedom
	int m_degreesOfFreedom;

};

#endif
