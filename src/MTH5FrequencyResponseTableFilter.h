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
#ifndef DBLDEF_MTH5_FREQ_RESP_TABLE_FILTER
#define DBLDEF_MTH5_FREQ_RESP_TABLE_FILTER

#include "MTH5Filter.h"

// Class to hold filter information for MTH5 frequency responsee table files
class MTH5FrequencyResponseTableFilter : public MTH5Filter {

public:

	// Constructer
	MTH5FrequencyResponseTableFilter();

	// Destructer
	~MTH5FrequencyResponseTableFilter();

	// Set frequencies
	void setFrequencies( const std::vector<double>& frequencies );

	// Set amplitudes
	void setAmplitude( const std::vector<double>& amplitudes );

	// Set phases
	void setPhases( const std::vector<double>& phases );

	// Get frequency response functions using the requency response functions of filter
	virtual std::complex<double> getFrequencyResponse(const double freq) const;

private:

	// Frequencies
	std::vector<double> m_frequencies;

	// Amplitudes
	std::vector<double> m_amplitudes;

	// Phases
	std::vector<double> m_phases;

	// Array of log10(frequency)
	double* m_log10Frequencies;

	// Array of log10(amplitude)
	double* m_log10Amplitudes;

	// Array of Phases
	double* m_phaseArray;

	// Calculate frequency response functions of filter
	std::complex<double> calcResponse(const double freq) const;

	// Copy constructer
	MTH5FrequencyResponseTableFilter(const MTH5FrequencyResponseTableFilter& rhs);

	// Assignment operator
	MTH5FrequencyResponseTableFilter& operator=(const MTH5FrequencyResponseTableFilter& rhs);

};

#endif
