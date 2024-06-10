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
#ifndef DBLDEF_ROBUST_WEIGHT
#define DBLDEF_ROBUST_WEIGHT

#include "CommonParameters.h"
#include <complex>

// Class of Robust estimator
class RobustWeight{

public:

	// Default constructer
	RobustWeight();

	// Destructer
	virtual ~RobustWeight();

	// Calculate sum of the second order derivative of loss function
	virtual double calculateSumOf2ndOrderDerivativeOfLossFunction(const int numSegments, const std::complex<double>* const residuals, const double scale) const = 0;

	// Calculate scale of absolute values of complex residuals by MADN
	static double calculateScaleByMADN(const int numSegments, const std::complex<double>* const residuals);

	// Calculate scale of absolute values of complex residuals
	virtual double calculateScale(const int numSegments, const std::complex<double>* const residuals, const double scalePre) const;

	// Calculate weights
	virtual double calculateWeights(const int numSegments, const std::complex<double>* const residuals, const double scale, const double* const weightsPrior,
		double* weights) const = 0;

	// Get convergence criteria
	double getConvergenceCriteria() const;

	// Get name of M-estimator
	std::string getNameOfRobustWeight() const;

	// Get maximum number of iteration
	int getNumIterationMax() const;

	// Set convergence criteria
	void setConvergenceCriteria( const double criteria );

	// Set name of M-estimator
	void setNameOfRobustWeight( const std::string& name );

	// Set maximum number of iteration
	void setNumIterationMax( const int numIterMax );

protected:

private:
	// Convergence criteria
	double m_convergenceCriteria;

	// Name of M-estimator
	std::string m_name;

	// Maximum number of iteration
	int m_numIterationMax;

};

#endif
