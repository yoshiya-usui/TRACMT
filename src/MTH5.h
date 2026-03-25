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
#ifndef DBLDEF_MTH5
#define DBLDEF_MTH5

#include <string>
#include <vector>
#include <H5Cpp.h>

#include "MTH5ChannelResponse.h"

// Class of MTH5 file
class MTH5{

public:

	// Return the the instance of the class
    static MTH5* getInstance();

	// Read MTH5 file
	void readMTH5File( const std::string& fileName, const std::string groupName, const int numSkipData, const int numDataPoints, double* data ) const;

	// Get name of the calibration file name made from the channel responses 
	static std::string getCalibrationFileName(const int channelIndex);

	// Get all filters and combine filters into a complete channel response for each channel
	void createChannelResponses(const std::vector< std::pair< std::string, std::string> >& fileNameAndPath);

	//// Correct frequency response functions using the requency response functions of all filter
	//void correctResponse(const int channelIndex, const double freq, const double samplingFreq, std::complex<double>& response) const;

	// Make calibration files using the requency response functions of all filter 
	void makeCalibrationFile(const int channelIndex, const std::vector<double>& freqs) const;

private:

	// Number of channel responses
	int m_numOfChannelRespones;

	// List of channel response (combination of all filters)
	MTH5ChannelResponse* m_channelResponses;

	// Constructer
	MTH5();

	// Destructer
	~MTH5();

	// Copy constructer
	MTH5(const MTH5& rhs);

	// Assignment operator
	MTH5& operator=(const MTH5& rhs);

};

#endif
