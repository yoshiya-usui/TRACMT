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
#ifndef DBLDEF_CALIBRATION_FUNCTION
#define DBLDEF_CALIBRATION_FUNCTION

#include <complex>
#include <string>

// Class of calibration function
class CalibrationFunction{

public:

	// Constructer
	CalibrationFunction( const std::string& calibrationFile );

	// Default constructer
	CalibrationFunction();

	// Destructer
	~CalibrationFunction();

	// Read calibration function
	void readCalibrationFunction( const std::string& calibrationFile );

	// Calculate calibration function for a input frequency
	std::complex<double> calculateCalibrationFunction( const double freq ) const;

private:
	// Constant factor
	double m_factor;

	// Number of frequencies in logarithmic scale
	int m_numFreqs;

	// Arrays of frequencies in logarithmic scale
	double* m_logFreqs;

	// Arrays of amplitudes in logarithmic scale
	double* m_logAmps;

	// Arrays of phase (radian)
	double* m_phases;

	// Copy constructer
	CalibrationFunction(const CalibrationFunction& CalibrationFunction);

	// Assignment operator
	CalibrationFunction& operator=(const CalibrationFunction& CalibrationFunction);

};

#endif
