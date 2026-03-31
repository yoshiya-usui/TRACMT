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
#ifndef DBLDEF_MTH5_CHANNEL_RESPONSE
#define DBLDEF_MTH5_CHANNEL_RESPONSE

#include <string>
#include <vector>
#include <H5Cpp.h>

#include "MTH5Filter.h"

// Class to hold complete channel response (combination of all filters)
class MTH5ChannelResponse {

public:

	// Default constructer
	MTH5ChannelResponse();

	// Destructer
	~MTH5ChannelResponse();

	// Get list of filter names
	std::vector<std::string> getFilterNames() const;

	// Get total number of filters
	int getFilterCount() const;

	// Check if response contains time delay filters
	bool hasDelayFilters() const;

	// Verify filters are properly ordered by sequence number
	bool isProperlyOrdered() const;

	// Get input units from the first filter
	std::string getUnitsIn() const;

	// Get output units from the last filter
	std::string getUnitsOut() const;

	// Get normalization frequency in Hz
	double getNormalizationFrequency() const;

	// When going from physical units to digital counts we are applying filters: true
	// When going from counts to physical units we are unapplying the filters: false
	bool isAppliedForwardly() const;

	// Whether response is valid
	bool isValid() const;

	// Get error message if invalid
	std::string getErrorMessage() const;

	// Get all filters for a channel and combine filters into a complete channel response
	// This method validates the filter chain and combines them into a single ChannelResponse object,
	// similar to the ChannelResponse class in Python. The response includes:
	// - All filters in the proper sequence
	// - Input units from the first filter
	// - Output units from the last filter
	// - Validation status and any error messages
	void createChannelResponse(const std::string& fileName, const std::string& channelPath);

	// Calculate frequency response functions using the requency response functions of all filter
	std::complex<double> calcResponse(const double freq) const;

private:

	// List of filters in order
	std::vector<MTH5Filter*> m_filtersList;

	// When going from physical units to digital counts we are applying filters: m_isAppliedForwardly = true
	// When going from counts to physical units we are unapplying the filters: m_isAppliedForwardly = false
	bool m_isAppliedForwardly;

	// Normalization frequency in Hz
	double m_normalizationFrequency;     

	// Whether response is valid
	bool m_isValid;                  

	// Error message if invalid
	std::string m_errorMessage;

	// Get all filters for a channel
	// Example usage:
	//   getChannelFilters("test.h5", "/Survey/Stations/MT001/MT001a/Ex", filters);
	void getChannelFilters(const std::string& fileName, const std::string& channelPath, std::vector<MTH5Filter*>& filters) const;

	// Read a single filter from the Filters group
	MTH5Filter* readFilter(const std::string& fileName, const std::string& filterPath, const std::string& filterName, int sequenceNumber) const;

	// Read filter names from channel dataset attribute
	std::vector<std::string> readFilterNamesFromChannel(const std::string& fileName, const std::string& channelPath) const;

	// Combine filters for a channel into a complete channel response
	std::string combineFilters();

	// Validate filter list for unit consistency
	// This ensures that the output units of each filter match the input units of the next filter,
	// similar to the _check_consistency_of_units method in Python's ChannelResponse class.
	bool validateFilterUnits(const std::vector<MTH5Filter*>& filters, std::string& errorMessage) const;

	// Helper: Navigate channel path to Filters group path
	std::string getFiltersGroupPath(const std::string& channelPath) const;

	// Helper: Read string attribute from HDF5 object
	std::string readStringAttribute(hid_t obj_id, const std::string& attrName) const;

	// Helper: Read string array attribute from HDF5 object
	std::vector<std::string> readStringArrayAttribute(hid_t obj_id, const std::string& attrName) const;

	// Sort filters by ascending order using bubble sort
	void sortFilters(std::vector<MTH5Filter*>& filters) const;

	// Copy constructer
	MTH5ChannelResponse(const MTH5ChannelResponse& rhs);

	// Assignment operator
	MTH5ChannelResponse& operator=(const MTH5ChannelResponse& rhs);

};

#endif
