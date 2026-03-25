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
#ifndef DBLDEF_MTH5_FIR_FILTER
#define DBLDEF_MTH5_FIR_FILTER

#include "MTH5Filter.h"

// Class to hold filter information for MTH5 FIR files
class MTH5FIRFilter : public MTH5Filter {

public:

	// Constructer
	MTH5FIRFilter();

	// Destructer
	~MTH5FIRFilter();

	// Type of symmetry
	enum {
		EVEN,
		ADD,
		ASYMMETRIC,
	};

	// Set input samping rate befor decimation
	void setDecimationInputSampleRate( const int decimationInputSampleRate );

	// Set gain
	void setGain( const double gain );

	// Set type of symmetry
	void setTypeOfSymmetry( const int typeOfSymmetry );

	// Set decimation factor
	void setDecimationFactor( const int decimationFactor );

	// Set FIR coefficients
	void setFIRCoefficients ( const std::vector<double>& coefficients );

	// Get frequency response functions using the requency response functions of filter
	virtual  std::complex<double> getFrequencyResponse(const double freq) const;

private:

	// input samping rate befor decimation
	int m_decimationInputSampleRate;

	// Gain
	double m_gain;

	// Type of symmetry
	int m_typeOfSymmetry;

	// Decimation factor
	int m_decimationFactor;

	// FIR coefficients
	std::vector<double> m_FIRCoefficients;

	// Calculate frequency response functions of filter
	std::complex<double> calcResponse(const double freq, const double samplingFreq) const;

	// Copy constructer
	MTH5FIRFilter(const MTH5FIRFilter& rhs);

	// Assignment operator
	MTH5FIRFilter& operator=(const MTH5FIRFilter& rhs);

};
#endif