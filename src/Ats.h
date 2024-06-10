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
#ifndef DBLDEF_ATS
#define DBLDEF_ATS

#include <complex>
#include <vector>
#include <string>

// Class of ats file
class Ats{

public:

	enum CoilType{
		MFS06 = 0,
		MFS07
	};

	// Return the the instance of the class
    static Ats* getInstance();

	// Read ats file
	void readAtsFile( const std::string& fileName, const int numSkipData, const int numDataPoints, double* data ) const;

	// Get calibration file name
	std::string getCalibrationFileName( const int channelIndex ) const;

	// Make calibration file
	void makeCalibrationFile( const std::string& inputString, const int channelIndex, const std::vector<double>& freq ) const;

	// Get flag specifing whether the calibration function for ADU is calculated
	bool isCalADUResp() const;

	// Set flag specifing whether the calibration function for ADU is calculated
	void setIsCalADUResp( const bool isCalADUResp );

private:

	short int bytesToInt16(unsigned char* ptr) const;

	int bytesToInt32(unsigned char* ptr) const;

	float bytesToFloat(unsigned char* ptr) const;

	double bytesToDouble(unsigned char* ptr) const;

	// Constructer
	Ats();

	// Destructer
	~Ats();

	// Copy constructer
	Ats(const Ats& rhs);

	// Assignment operator
	Ats& operator=(const Ats& rhs);

	// Calculate calibration function for coil
	void calculateCalibrationFunctionForCoil( const std::string& calibrationFile, const int coilType,
		const std::vector<double>& freq, std::complex<double>* calibrationFunction ) const;

	// Calculate calibration function for ADU
	std::complex<double> calculateCalibrationFunctionForADULFhannel( const double freq ) const;

	// Flag specifing whether the calibration function for ADU is calculated
	bool m_isCalADUResp;

};

#endif
