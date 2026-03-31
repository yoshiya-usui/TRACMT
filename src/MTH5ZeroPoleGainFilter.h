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
#ifndef DBLDEF_MTH5_POLE_ZERO_FILTER
#define DBLDEF_MTH5_POLE_ZERO_FILTER

#include "MTH5Filter.h"

// Class to hold filter information for MTH5 pole-zero files
class MTH5ZeroPoleGainFilter : public MTH5Filter {

public:

	// Constructer
	MTH5ZeroPoleGainFilter();

	// Destructer
	~MTH5ZeroPoleGainFilter();

	// Set normalization factor
	void setNormalizationFactor(const double normalizationFactor);

	// Set poles
	void setPoles( const std::vector< std::complex<double> >& poles );

	// Set zeros
	void setZeros(const std::vector< std::complex<double> >& zeros);

	// Get frequency response functions using the requency response functions of filter
	virtual std::complex<double> getFrequencyResponse(const double freq) const;

private:

	// Normalization factor
	double m_normalizationFactor;

	// Poles
	std::vector< std::complex<double> > m_poles;

	// Zeros
	std::vector< std::complex<double> > m_zeros;

	// Calculate frequency response functions of filter
	std::complex<double> calcResponse(const double freq, const double samplingFreq, const std::complex<double>& respones) const;

	// Copy constructer
	MTH5ZeroPoleGainFilter(const MTH5ZeroPoleGainFilter& rhs);

	// Assignment operator
	MTH5ZeroPoleGainFilter& operator=(const MTH5ZeroPoleGainFilter& rhs);

};

#endif
