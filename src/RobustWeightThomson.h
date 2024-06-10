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
#ifndef DBLDEF_ROBUST_WEIGHT_THOMSON
#define DBLDEF_ROBUST_WEIGHT_THOMSON

#include "RobustWeight.h"
#include "CommonParameters.h"
#include <complex>

// Class of M-estimator using Thomson weights
class RobustWeightThomson : public RobustWeight{

public:

	// Default constructer
	RobustWeightThomson();

	// Destructer
	virtual ~RobustWeightThomson();

	// Calculate sum of the second order derivative of loss function
	virtual double calculateSumOf2ndOrderDerivativeOfLossFunction(const int numSegments, const std::complex<double>* const residuals, const double scale) const;

	// Calculate weights
	virtual double calculateWeights(const int numSegments, const std::complex<double>* const residuals, const double scale, const double* const weightsPrior,
		double* weights) const;

	// Get parameter for determing the probability used for outlier rejection
	double getParameterForDetermingProbability() const;

	// Set parameter for determing the probability used for outlier rejection
	void setParameterForDetermingProbability( const double prob );

protected:

private:

	// Parameter for determing the probability used for outlier rejection
	double m_parameterForDetermingProbability;

};

#endif
