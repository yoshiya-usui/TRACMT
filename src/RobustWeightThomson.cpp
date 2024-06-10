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
#include "RobustWeightThomson.h"
#include "OutputFiles.h"
#include "Util.h"

#include <iostream>

// Default constructer
RobustWeightThomson::RobustWeightThomson():
	m_parameterForDetermingProbability(0.0)
{
	setNameOfRobustWeight("Thomson");
}

// Destructer
RobustWeightThomson::~RobustWeightThomson()
{
}

// Calculate sum of the second order derivative of loss function
double RobustWeightThomson::calculateSumOf2ndOrderDerivativeOfLossFunction(const int numSegments, const std::complex<double>* const residuals, const double scale) const {

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeErrorMessage("Parametric error estimation cannot be used for Thomson weight");
	return -1;

}

// Calculate weights
double RobustWeightThomson::calculateWeights(const int numSegments, const std::complex<double>* const residuals, const double scale,
	const double* const weightsPrior, double* weights) const {

	const double param = getParameterForDetermingProbability();
	double probability = (static_cast<double>(numSegments) - param - 0.5) / static_cast<double>(numSegments);
	if (param < 0.0 - CommonParameters::EPS) {
		probability = -param;
	}
	const double Q = sqrt(-2.0 * log(1.0 - probability));
	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeCvgMessage("	Probability for determining the point at which downweighting starts: " + Util::toString(probability));
	ptrOutputFiles->writeCvgMessage("	Scale factor: " + Util::toString(scale));
#ifdef _DEBUG_WRITE
	std::cout << "probability: " << probability << std::endl;
	std::cout << "Q: " << Q << std::endl;
	std::cout << "scaleFactor: " << scale << std::endl;
#endif
	double sumOfWeights(0.0);
	for (int iSeg = 0; iSeg < numSegments; ++iSeg) {
		weights[iSeg] = weightsPrior[iSeg];
		const double x = std::abs(residuals[iSeg]) / scale;
		weights[iSeg] *= exp(exp(-Q * Q)) * exp(-exp(Q * (fabs(x) - Q)));
		sumOfWeights += weights[iSeg];
#ifdef _DEBUG_WRITE
		std::cout << "iSeg absres x weightsPrior weights: "
			<< iSeg << " " << std::abs(residuals[iSeg]) << " " << x << " " << weightsPrior[iSeg] << " " << weights[iSeg] << std::endl;
#endif
	}

	return sumOfWeights;

}

// Get parameter for determing the probability used for outlier rejection
double RobustWeightThomson::getParameterForDetermingProbability() const{
	return m_parameterForDetermingProbability;
}

// Set parameter for determing the probability used for outlier rejection
void RobustWeightThomson::setParameterForDetermingProbability( const double prob ){
	m_parameterForDetermingProbability = prob;
}
