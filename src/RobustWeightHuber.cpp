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
#include "RobustWeightHuber.h"
#include "OutputFiles.h"
#include "Util.h"

#include <iostream>

// Default constructer
RobustWeightHuber::RobustWeightHuber():
	m_threshould(3.0)
{
	setNameOfRobustWeight("Huber");
}

//// Copy constructer
//RobustWeightHuber::RobustWeightHuber(const RobustWeightHuber& obj){
//	m_threshould = obj.m_threshould;
//}

// Destructer
RobustWeightHuber::~RobustWeightHuber()
{
}

// Calculate sum of the second order derivative of loss function
double RobustWeightHuber::calculateSumOf2ndOrderDerivativeOfLossFunction(const int numSegments, const std::complex<double>* const residuals, const double scale) const {

	const double threshould = getThreshould();
	double sumOfDerivatives(0.0);
	for (int iSeg = 0; iSeg < numSegments; ++iSeg) {
		const double x = std::abs(residuals[iSeg]) / scale;
		if (x <= threshould) {
			sumOfDerivatives += 1.0;
		}
	}
	return sumOfDerivatives;

}

// Calculate weights
double RobustWeightHuber::calculateWeights(const int numSegments, const std::complex<double>* const residuals, const double scale,
	const double* const weightsPrior, double* weights) const {

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeCvgMessage("	Scale factor: " + Util::toString(scale));
#ifdef _DEBUG_WRITE
	std::cout << "scaleFactor: " << scale << std::endl;
#endif
	const double threshould = getThreshould();
	double sumOfWeights(0.0);
	for (int iSeg = 0; iSeg < numSegments; ++iSeg) {
		weights[iSeg] = weightsPrior[iSeg];
		const double x = std::abs(residuals[iSeg]) / scale;
		if (x <= threshould) {
			// Huber weight is unity
		}
		else {
			weights[iSeg] *= threshould / x;
		}
		sumOfWeights += weights[iSeg];
#ifdef _DEBUG_WRITE
		std::cout << "iSeg absres x weightsPrior weights: "
			<< iSeg << " " << std::abs(residuals[iSeg]) << " " << x << " " << weightsPrior[iSeg] << " " << weights[iSeg] << std::endl;
#endif
	}

	return sumOfWeights;

}

// Get threshould value for downweighting
double RobustWeightHuber::getThreshould() const{
	return m_threshould;
}

// Set threshould value for downweighting
void RobustWeightHuber::setThreshould( const double threshould ){
	m_threshould = threshould;
}
