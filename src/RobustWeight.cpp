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
#include "RobustWeight.h"
#include "OutputFiles.h"
#include "Util.h"

// Default constructer
RobustWeight::RobustWeight():
	m_convergenceCriteria(0.01),
	m_name("no name"),
	m_numIterationMax(10)
{
}

// Destructer
RobustWeight::~RobustWeight()
{
}

// Get convergence criteria
double RobustWeight::getConvergenceCriteria() const{
	return m_convergenceCriteria;
}

// Get name of M-estimator
std::string RobustWeight::getNameOfRobustWeight() const{
	return m_name;
}

// Get maximum number of iteration
int RobustWeight::getNumIterationMax() const{
	return m_numIterationMax;
}

// Set convergence criteria
void RobustWeight::setConvergenceCriteria( const double criteria ){
	m_convergenceCriteria = criteria;
}

// Set name of M-estimator
void RobustWeight::setNameOfRobustWeight( const std::string& name ){
	m_name = name;
}

// Set maximum number of iteration
void RobustWeight::setNumIterationMax( const int numIterMax ){
	m_numIterationMax = numIterMax;
}

// Calculate scale of absolute values of complex residuals by MADN
double RobustWeight::calculateScaleByMADN(const int numSegments, const std::complex<double>* const residuals) {

	std::vector<double> absoluteResiduals;
	for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
		absoluteResiduals.push_back(std::abs(residuals[iSeg]));
	}

	const double sigmaMAD = Util::calculateMAD(absoluteResiduals);
	
	return sigmaMAD / 0.448453;

}

// Calculate scale of absolute values of complex residuals
double RobustWeight::calculateScale(const int numSegments, const std::complex<double>* const residuals, const double scalePre) const{
	return calculateScaleByMADN(numSegments, residuals);
}
