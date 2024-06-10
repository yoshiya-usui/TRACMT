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
#ifndef DBLDEF_ELOGDual
#define DBLDEF_ELOGDual

#include <stdio.h>
#include <stdint.h>
#include <string>
#include <vector>
#include <complex>

// Class of ELOG-Dual
class ElogDual {

public:

	struct recdata_head{
		uint8_t hour;
		uint8_t min;
		uint8_t sec;
	};

	const static int AD_CH = 2;
	const static int AD_BYTES = 3;
	const static int REC_DATALEN_HEAD = sizeof(recdata_head);

	// Return the the instance of the class
	static ElogDual* getInstance();

	// Read elog binary file
	void readElogBinaryFile( const std::string& fileName, const int numSkipData, const int numDataPoints, int& counter, double* ex, double* ey) const;

	// Read elog binary file under a directory
	void readElogBinaryFilesUnderADirectory(const std::string& directoryName, const int numSkipData, const int numDataPoints, double* ex, double* ey) const;

	// Get calibration file name
	std::string getCalibrationFileName( const int channelIndex ) const;

	// Make calibration file
	void makeCalibrationFile( const std::string& fileName, const double unitGroupDelay, const int channelIndexX, const int channelIndexY,
		const double dipoleLengthX, const double dipoleLengthY, const std::vector<double>& freq ) const;

private:

	// Calculate calibration function for analog filter
	void calculateCalibrationFunctionForAnalogFilter( const std::string& fileName, const double freq,
		std::complex<double>& calElogX, std::complex<double>& calElogY ) const;

	int  bytesToInt32(unsigned char* ptr) const;

};
#endif
