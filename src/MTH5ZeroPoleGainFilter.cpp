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
#include "MTH5ZeroPoleGainFilter.h"
#include "CommonParameters.h"

// Constructer
MTH5ZeroPoleGainFilter::MTH5ZeroPoleGainFilter() :
	MTH5Filter(),
	m_normalizationFactor(0.0)
{
}

// Destructer
MTH5ZeroPoleGainFilter::~MTH5ZeroPoleGainFilter() {
}

// Set normalization factor
void MTH5ZeroPoleGainFilter::setNormalizationFactor(const double normalizationFactor) {
	m_normalizationFactor = normalizationFactor;
}

// Set poles
void MTH5ZeroPoleGainFilter::setPoles(const std::vector< std::complex<double> >& poles) {
	m_poles = poles;
}

// Set zeros
void MTH5ZeroPoleGainFilter::setZeros(const std::vector< std::complex<double> >& zeros) {
	m_zeros = zeros;
}

// Get frequency response functions using the requency response functions of filter
// @note under construction
std::complex<double> MTH5ZeroPoleGainFilter::getFrequencyResponse(const double freq) const {

	// Analog ZPK filter in the Laplace domain: H(s) = K * prod(s - z_i) / prod(s - p_i)
	// Evaluate on the imaginary axis: s = j*omega, omega = 2*pi*freq
	const std::complex<double> s(0.0, 2.0 * CommonParameters::PI * freq);
	std::complex<double> numerator(1.0, 0.0);
	for (std::vector< std::complex<double> >::const_iterator itr = m_zeros.begin(); itr != m_zeros.end(); ++itr) {
		numerator *= (s - *itr);
	}
	std::complex<double> denominator(1.0, 0.0);
	for (std::vector< std::complex<double> >::const_iterator itr = m_poles.begin(); itr != m_poles.end(); ++itr) {
		denominator *= (s - *itr);
	}
	return m_normalizationFactor * numerator / denominator;

}
