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
#include "MTH5FrequencyResponseTableFilter.h"
#include "Util.h"

// Constructer
MTH5FrequencyResponseTableFilter::MTH5FrequencyResponseTableFilter():
	MTH5Filter(),
	m_log10Frequencies(NULL),
	m_log10Amplitudes(NULL),
	m_phaseArray(NULL)
{
}

// Destructer
MTH5FrequencyResponseTableFilter::~MTH5FrequencyResponseTableFilter() {

	if (m_log10Frequencies != NULL) {
		delete [] m_log10Frequencies;
	}
	if (m_log10Amplitudes != NULL) {
		delete[] m_log10Amplitudes;
	}
	if (m_phaseArray != NULL) {
		delete[] m_phaseArray;
	}

}

// Set frequencies
void MTH5FrequencyResponseTableFilter::setFrequencies(const std::vector<double>& frequencies) {

	m_frequencies = frequencies;

	if (m_log10Frequencies != NULL) {
		delete[] m_log10Frequencies;
	}
	if (frequencies.size() > 0){
		m_log10Frequencies = new double[frequencies.size()];
		int icount(0);
		for (std::vector<double>::const_iterator itr = frequencies.begin(); itr != frequencies.end(); ++itr, ++icount){
			m_log10Frequencies[icount] = log10(*itr);
		}
	}

}

// Set amplitudes
void MTH5FrequencyResponseTableFilter::setAmplitude(const std::vector<double>& amplitudes) {

	m_amplitudes = amplitudes;

	if (m_log10Amplitudes != NULL) {
		delete[] m_log10Amplitudes;
	}
	if (amplitudes.size() > 0) {
		m_log10Amplitudes = new double[amplitudes.size()];
		int icount(0);
		for (std::vector<double>::const_iterator itr = amplitudes.begin(); itr != amplitudes.end(); ++itr, ++icount) {
			m_log10Amplitudes[icount] = log10(*itr);
		}
	}

}

// Set phases
void MTH5FrequencyResponseTableFilter::setPhases(const std::vector<double>& phases) {
	
	m_phases = phases;

	if (m_phaseArray != NULL) {
		delete[] m_phaseArray;
	}
	if (phases.size() > 0) {
		m_phaseArray = new double[phases.size()];
		int icount(0);
		for (std::vector<double>::const_iterator itr = phases.begin(); itr != phases.end(); ++itr, ++icount) {
			m_phaseArray[icount] = *itr;
		}
	}

}

// Get frequency response functions using the requency response functions of filter
std::complex<double> MTH5FrequencyResponseTableFilter::getFrequencyResponse(const double freq) const {

	if (m_frequencies.empty()) {
		return std::complex<double>(1.0, 0.0);
	}

	const int numFreqs = static_cast<int>(m_frequencies.size());
	const double log10Freq = log10(freq);
	const double logAmp = Util::interpolationAkima(numFreqs, m_log10Frequencies, m_log10Amplitudes, log10Freq);
	const double phs = Util::interpolationAkima(numFreqs, m_log10Frequencies, m_phaseArray, log10Freq);
	return pow(10.0, logAmp) * std::complex<double>(cos(phs), sin(phs));

}
