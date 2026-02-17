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

// Structure to hold filter information
struct FilterInfo {
	std::string name;           // Filter name
	std::string type;           // Filter type (zpk, coefficient, time_delay, fap, fir)
	std::string units_in;       // Input units
	std::string units_out;      // Output units
	int sequence_number;        // Order in filter chain
	
	// Filter-specific data (usage depends on filter type)
	std::vector<double> data;   // Generic data array
	
	FilterInfo() : sequence_number(0) {}
};

// Class of MTH5 file
class MTH5{

public:

	// Return the the instance of the class
    static MTH5* getInstance();

	// Read MTH5 file
	void readMTH5File( const std::string& fileName, const std::string groupName, const int numSkipData, const int numDataPoints, double* data ) const;

	// Read filter names from channel dataset attribute
	std::vector<std::string> readFilterNamesFromChannel( const std::string& fileName, const std::string& channelPath ) const;
	
	// Read a single filter from the Filters group
	FilterInfo readFilter( const std::string& fileName, const std::string& filterPath, const std::string& filterName, int sequenceNumber ) const;
	
	// Get all filters for a channel (similar to Python's channel_response property)
	std::vector<FilterInfo> getChannelFilters( const std::string& fileName, const std::string& channelPath ) const;

private:
	
	// Helper: Read string attribute from HDF5 object
	std::string readStringAttribute( hid_t obj_id, const std::string& attrName ) const;
	
	// Helper: Read string array attribute from HDF5 object  
	std::vector<std::string> readStringArrayAttribute( hid_t obj_id, const std::string& attrName ) const;
	
	// Helper: Navigate channel path to Filters group path
	std::string getFiltersGroupPath( const std::string& channelPath ) const;

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
