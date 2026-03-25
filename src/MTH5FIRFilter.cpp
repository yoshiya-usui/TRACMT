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
#include "MTH5FIRFilter.h"
#include "CommonParameters.h"

// Constructer
MTH5FIRFilter::MTH5FIRFilter() :
	MTH5Filter(),
	m_decimationInputSampleRate(0),
	m_gain(0.0),
	m_typeOfSymmetry(MTH5FIRFilter::EVEN),
	m_decimationFactor(0)
{
}

// Destructer
MTH5FIRFilter::~MTH5FIRFilter() {
}

// Set input samping rate befor decimation
void MTH5FIRFilter::setDecimationInputSampleRate(const int decimationInputSampleRate){
	m_decimationInputSampleRate = decimationInputSampleRate;
}

// Set gain
void MTH5FIRFilter::setGain(const double gain) {
	m_gain = gain;
}

// Set type of symmetry
void MTH5FIRFilter::setTypeOfSymmetry(const int typeOfSymmetry) {
	m_typeOfSymmetry = typeOfSymmetry;
}

// Set decimation factor
void MTH5FIRFilter::setDecimationFactor(const int decimationFactor) {
	m_decimationFactor = decimationFactor;
}

// Set FIR coefficients
void MTH5FIRFilter::setFIRCoefficients(const std::vector<double>& coefficients) {
	m_FIRCoefficients = coefficients;
}

// Get frequency response functions using the requency response functions of filter
// @note under construction
std::complex<double> MTH5FIRFilter::getFrequencyResponse(const double freq) const {

	const double omega = 2.0 * CommonParameters::PI * freq / m_decimationInputSampleRate;
	const int numCoeffs = static_cast<int>(m_FIRCoefficients.size());
	const double adjust = (1.0 - static_cast<double>(numCoeffs)) / 2.0;
	std::complex<double> response(0.0, 0.0);
	switch (m_typeOfSymmetry) {
	case MTH5FIRFilter::EVEN:
		// Even-symmetric FIR: h[n] = h[N-1-n]. The general DTFT sum is used;
		// imaginary parts cancel due to symmetry, yielding a real-valued response.
		for (int i = 0; i < numCoeffs; ++i) {
			const double arg = (static_cast<double>(i + 1 - numCoeffs) - adjust) * omega;
			response += m_FIRCoefficients[i] * std::complex<double>(cos(arg), sin(arg));
		}
		response *= m_gain;
		return response;
		break;
	case MTH5FIRFilter::ADD:
		// Odd-symmetric (anti-symmetric) FIR: h[n] = -h[N-1-n]. The general DTFT sum
		// is used; real parts cancel due to anti-symmetry, yielding a purely imaginary response.
		for (int i = 0; i < numCoeffs; ++i) {
			const double arg = (static_cast<double>(i + 1 - numCoeffs) - adjust) * omega;
			response += m_FIRCoefficients[i] * std::complex<double>(cos(arg), sin(arg));
		}
		response *= m_gain;
		return response;
		break;
	case MTH5FIRFilter::ASYMMETRIC:
		for (int i = 0; i < numCoeffs; ++i) {
			const double arg = (static_cast<double>(i + 1 - numCoeffs) - adjust) * omega;
			response += m_FIRCoefficients[i] * std::complex<double>(cos(arg), sin(arg));
		}
		response *= m_gain;
		return response;
		break;
	default:
		return 1.0;
		break;
	}

}
