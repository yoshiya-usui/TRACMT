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
#ifndef DBLDEF_MTH5_FILTER
#define DBLDEF_MTH5_FILTER

#include <complex>
#include <string>
#include <vector>

// Class to hold filter information of MTH5 files
class MTH5Filter {

public:

	// Constructer
	MTH5Filter();

	// Destructer
	~MTH5Filter();

	// Copy constructer
	MTH5Filter(const MTH5Filter& rhs);

	// Assignment operator
	MTH5Filter& operator=(const MTH5Filter& rhs);

	// Get filter name
	std::string getName() const;

	// Get filter type (zpk, coefficient, time_delay, fap, fir)
	std::string getType() const;

	// Input units
	std::string getUnitsIn() const;

	// Output units
	std::string getUnitsOut() const;

	// Get order in filter chain
	int getSequenceNumber() const;

	// Set filter name
	void setName(const std::string& name);

	// Set filter type (zpk, coefficient, time_delay, fap, fir)
	void setType(const std::string& type);

	// Set input units
	void setUnitsIn(const std::string& unitsIn);

	// Set output units
	void setUnitsOut(const std::string& unitsOut);

	// Set order in filter chain
	void setSequenceNumber(const int sequenceNumber);

	// Get frequency response functions using the requency response functions of filter
	virtual  std::complex<double> getFrequencyResponse(const double freq) const;

protected:
	
	// Filter name
	std::string m_name;

	// Filter type (zpk, coefficient, time_delay, fap, fir)
	std::string m_type;

	// Input units
	std::string m_unitsIn;

	// Output units
	std::string m_unitsOut;

	// Order in filter chain
	int m_sequenceNumber;

	// Filter-specific data (usage depends on filter type)
	std::vector<double> m_data;

};

#endif
