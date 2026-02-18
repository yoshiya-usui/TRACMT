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
#include "MTH5.h"
#include "Util.h"
#include "OutputFiles.h"

#include <H5Cpp.h>
#include <algorithm>
#include <sstream>

// Default constructer
MTH5::MTH5()
{
}

// Destructer
MTH5::~MTH5(){
}

// Return the instance of the class
MTH5* MTH5::getInstance(){
   	static MTH5 instance;// The only instance
  	return &instance;
}

// Read MTH5 file
void MTH5::readMTH5File(const std::string& fileName, const std::string groupName, const int numSkipData, const int numDataPoints, double* data) const{

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeLogMessage("Read data from " + groupName + " of "  + fileName);

	hid_t file_id = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	hid_t dataset_id = H5Dopen(file_id, groupName.c_str(), H5P_DEFAULT);
	hid_t filespace_id = H5Dget_space(dataset_id);
	int ndims = H5Sget_simple_extent_ndims(filespace_id);
	if (ndims < 1) {
		ptrOutputFiles->writeErrorMessage("Number of dimensions is less than one (" + Util::toString(ndims) + ")");
	}
	hsize_t* dims = new hsize_t[ndims];
	H5Sget_simple_extent_dims(filespace_id, dims, NULL);
	const int numDataInGroup = dims[0];
	delete[] dims;
	if (numDataInGroup < numDataPoints)
	{
		ptrOutputFiles->writeErrorMessage("Number of data in this group (" + Util::toString(numDataInGroup) + ") is smaller than " + Util::toString(numDataPoints));
	}
	hid_t type = H5Dget_type(dataset_id);
	hsize_t dimsmr[1] = { numDataPoints };
	hid_t memspace_id = H5Screate_simple(1, dimsmr, NULL);
	hsize_t roffset[1] = { numSkipData };
	hsize_t rcount[1] = { numDataPoints };
	H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, roffset, NULL, rcount, NULL);
	H5Dread(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, filespace_id, H5P_DEFAULT, data);
	
	H5Sclose(memspace_id);
	H5Sclose(filespace_id);
	H5Dclose(dataset_id);
	H5Fclose(file_id);
}

// Helper: Read string attribute from HDF5 object
std::string MTH5::readStringAttribute(hid_t obj_id, const std::string& attrName) const {
	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	
	if (!H5Aexists(obj_id, attrName.c_str())) {
		return "";
	}
	
	hid_t attr_id = H5Aopen(obj_id, attrName.c_str(), H5P_DEFAULT);
	if (attr_id < 0) {
		ptrOutputFiles->writeWarningMessage("Could not open attribute: " + attrName);
		return "";
	}
	
	hid_t type_id = H5Aget_type(attr_id);
	size_t size = H5Tget_size(type_id);
	
	char* buffer = new char[size + 1];
	H5Aread(attr_id, type_id, buffer);
	buffer[size] = '\0';
	
	std::string result(buffer);
	delete[] buffer;
	
	H5Tclose(type_id);
	H5Aclose(attr_id);
	
	return result;
}

// Helper: Read string array attribute from HDF5 object
std::vector<std::string> MTH5::readStringArrayAttribute(hid_t obj_id, const std::string& attrName) const {
	std::vector<std::string> result;
	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	
	if (!H5Aexists(obj_id, attrName.c_str())) {
		return result;
	}
	
	hid_t attr_id = H5Aopen(obj_id, attrName.c_str(), H5P_DEFAULT);
	if (attr_id < 0) {
		ptrOutputFiles->writeWarningMessage("Could not open attribute: " + attrName);
		return result;
	}
	
	hid_t space_id = H5Aget_space(attr_id);
	hssize_t npoints = H5Sget_simple_extent_npoints(space_id);
	
	if (npoints > 0) {
		hid_t type_id = H5Aget_type(attr_id);
		
		// Handle variable length strings
		if (H5Tis_variable_str(type_id)) {
			char** rdata = new char*[npoints];
			H5Aread(attr_id, type_id, rdata);
			
			for (int i = 0; i < npoints; i++) {
				if (rdata[i] != NULL) {
					result.push_back(std::string(rdata[i]));
				}
			}
			
			H5Dvlen_reclaim(type_id, space_id, H5P_DEFAULT, rdata);
			delete[] rdata;
		}
		// Handle fixed length strings
		else {
			size_t size = H5Tget_size(type_id);
			char* buffer = new char[npoints * size];
			H5Aread(attr_id, type_id, buffer);
			
			for (int i = 0; i < npoints; i++) {
				std::string str(&buffer[i * size]);
				if (!str.empty()) {
					result.push_back(str);
				}
			}
			
			delete[] buffer;
		}
		
		H5Tclose(type_id);
	}
	
	H5Sclose(space_id);
	H5Aclose(attr_id);
	
	return result;
}

// Helper: Navigate channel path to Filters group path
std::string MTH5::getFiltersGroupPath(const std::string& channelPath) const {
	// Channel path format: /Survey/Stations/{station}/Grouped by {run}/{channel}
	// or: /Experiment/Surveys/{survey}/Stations/{station}/{run}/{channel}
	// Filters are at: /Survey/Filters or /Experiment/Surveys/{survey}/Filters
	
	size_t pos = channelPath.find("/Stations/");
	if (pos != std::string::npos) {
		std::string prefix = channelPath.substr(0, pos);
		return prefix + "/Filters";
	}
	
	// Default fallback
	return "/Survey/Filters";
}

// Read filter names from channel dataset attribute
std::vector<std::string> MTH5::readFilterNamesFromChannel(const std::string& fileName, const std::string& channelPath) const {
	std::vector<std::string> filterNames;
	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	
	ptrOutputFiles->writeLogMessage("Reading filter names from channel: " + channelPath);
	
	hid_t file_id = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	if (file_id < 0) {
		ptrOutputFiles->writeErrorMessage("Could not open file: " + fileName);
		return filterNames;
	}
	
	hid_t dataset_id = H5Dopen(file_id, channelPath.c_str(), H5P_DEFAULT);
	if (dataset_id < 0) {
		ptrOutputFiles->writeErrorMessage("Could not open dataset: " + channelPath);
		H5Fclose(file_id);
		return filterNames;
	}
	
	// Try new format first: "filter.name" attribute
	filterNames = readStringArrayAttribute(dataset_id, "filter.name");
	
	// If empty, try reading individual filter objects (newer format)
	// This would require parsing JSON-like attribute structure
	// For simplicity, we'll focus on filter.name array
	
	if (filterNames.empty()) {
		ptrOutputFiles->writeLogMessage("No filters found in channel attributes");
	} else {
		ptrOutputFiles->writeLogMessage("Found " + Util::toString(static_cast<int>(filterNames.size())) + " filter(s)");
	}
	
	H5Dclose(dataset_id);
	H5Fclose(file_id);
	
	return filterNames;
}

// Read a single filter from the Filters group
FilterInfo MTH5::readFilter(const std::string& fileName, const std::string& filterPath, 
                            const std::string& filterName, int sequenceNumber) const {
	FilterInfo filter;
	filter.sequence_number = sequenceNumber;
	filter.name = filterName;
	
	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	
	// Normalize filter name (replace "/" with " per ", convert to lowercase)
	std::string normalizedName = filterName;
	size_t pos = 0;
	while ((pos = normalizedName.find("/", pos)) != std::string::npos) {
		normalizedName.replace(pos, 1, " per ");
		pos += 5;
	}
	std::transform(normalizedName.begin(), normalizedName.end(), normalizedName.begin(), ::tolower);
	
	ptrOutputFiles->writeLogMessage("Reading filter: " + normalizedName + " from " + filterPath);
	
	hid_t file_id = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	if (file_id < 0) {
		ptrOutputFiles->writeErrorMessage("Could not open file: " + fileName);
		return filter;
	}
	
	// Check all filter type subdirectories
	std::vector<std::string> filterTypes = {"zpk", "coefficient", "time_delay", "fap", "fir"};
	
	bool found = false;
	for (const auto& ftype : filterTypes) {
		std::string fullPath = filterPath + "/" + ftype + "/" + normalizedName;
		
		// Check if this path exists
		H5E_BEGIN_TRY {
			hid_t filter_id = H5Gopen(file_id, fullPath.c_str(), H5P_DEFAULT);
			if (filter_id >= 0) {
				found = true;
				filter.type = ftype;
				
				// Read filter metadata attributes
				filter.units_in = readStringAttribute(filter_id, "units_in");
				filter.units_out = readStringAttribute(filter_id, "units_out");
				
				// Read filter data (if it exists as a dataset)
				// Different filter types have different data structures
				if (ftype == "zpk") {
					// ZPK filters have zeros, poles, and gain datasets
					// Read these if needed for calibration
				} else if (ftype == "coefficient" || ftype == "fir") {
					// Coefficient/FIR filters have coefficient array
					// Read if dataset exists
				} else if (ftype == "fap") {
					// FAP filters have frequency, amplitude, phase arrays
				}
				
				H5Gclose(filter_id);
				ptrOutputFiles->writeLogMessage("Found filter in " + ftype + " group");
				break;
			}
		} H5E_END_TRY;
	}
	
	if (!found) {
		ptrOutputFiles->writeWarningMessage("Could not find filter: " + normalizedName + " in any filter type group");
	}
	
	H5Fclose(file_id);
	
	return filter;
}

// Get all filters for a channel (similar to Python's channel_response property)
std::vector<FilterInfo> MTH5::getChannelFilters(const std::string& fileName, const std::string& channelPath) const {
	std::vector<FilterInfo> filters;
	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	
	ptrOutputFiles->writeLogMessage("Getting channel filters for: " + channelPath);
	
	// Step 1: Read filter names from channel attributes
	std::vector<std::string> filterNames = readFilterNamesFromChannel(fileName, channelPath);
	
	if (filterNames.empty()) {
		ptrOutputFiles->writeLogMessage("No filters to load");
		return filters;
	}
	
	// Step 2: Determine the Filters group path
	std::string filtersPath = getFiltersGroupPath(channelPath);
	ptrOutputFiles->writeLogMessage("Filters group path: " + filtersPath);
	
	// Step 3: Read each filter
	for (size_t i = 0; i < filterNames.size(); i++) {
		FilterInfo filter = readFilter(fileName, filtersPath, filterNames[i], i + 1);
		if (!filter.type.empty()) {
			filters.push_back(filter);
		}
	}
	
	// Step 4: Sort filters by sequence_number to ensure correct order
	// This is critical because filter order determines the signal processing chain
	if (filters.size() > 1) {
		std::sort(filters.begin(), filters.end(), 
			[](const FilterInfo& a, const FilterInfo& b) {
				return a.sequence_number < b.sequence_number;
			});
		ptrOutputFiles->writeLogMessage("Sorted filters by sequence number");
	}
	
	ptrOutputFiles->writeLogMessage("Loaded " + Util::toString(static_cast<int>(filters.size())) + " filter(s)");
	
	return filters;

}

// Validate filter list for unit consistency
// This ensures that the output units of each filter match the input units of the next filter,
// similar to the _check_consistency_of_units method in Python's ChannelResponse class.
bool MTH5::validateFilterUnits(const std::vector<FilterInfo>& filters, std::string& errorMessage) const {
	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	
	if (filters.empty()) {
		return true; // Empty filter list is valid
	}
	
	if (filters.size() == 1) {
		return true; // Single filter is always valid
	}
	
	// Check that output units of one filter match input units of the next
	// Also verify that sequence numbers are sequential
	for (size_t i = 0; i < filters.size() - 1; i++) {
		const FilterInfo& currentFilter = filters[i];
		const FilterInfo& nextFilter = filters[i + 1];
		
		// Verify sequence numbers are in order
		if (currentFilter.sequence_number >= nextFilter.sequence_number) {
			std::ostringstream oss;
			oss << "Filter sequence error: Filter '" << currentFilter.name 
			    << "' has sequence " << currentFilter.sequence_number
			    << " but next filter '" << nextFilter.name 
			    << "' has sequence " << nextFilter.sequence_number;
			errorMessage = oss.str();
			ptrOutputFiles->writeErrorMessage(errorMessage);
			return false;
		}
		
		if (currentFilter.units_out != nextFilter.units_in) {
			std::ostringstream oss;
			oss << "Unit consistency error: Filter '" << currentFilter.name 
			    << "' (seq " << currentFilter.sequence_number << ") outputs " << currentFilter.units_out 
			    << " but filter '" << nextFilter.name 
			    << "' (seq " << nextFilter.sequence_number << ") expects " << nextFilter.units_in;
			errorMessage = oss.str();
			ptrOutputFiles->writeErrorMessage(errorMessage);
			return false;
		}
	}
	
	ptrOutputFiles->writeLogMessage("Filter units are consistent");
	return true;
}

// Combine filters into a complete channel response
// This method validates the filter chain and combines them into a single ChannelResponse object,
// similar to the ChannelResponse class in Python. The response includes:
// - All filters in the proper sequence
// - Input units from the first filter
// - Output units from the last filter
// - Validation status and any error messages
//
// Example usage:
//   std::vector<FilterInfo> filters = mth5->getChannelFilters(fileName, "/Survey/Stations/MT001/MT001a/Ex");
//   ChannelResponse response = mth5->createChannelResponse(filters);
//   if (response.is_valid) {
//       // Use the response for calibration
//   }
ChannelResponse MTH5::createChannelResponse(const std::vector<FilterInfo>& filters) const {
	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	
	ChannelResponse response;
	response.filters_list = filters;
	
	// Validate the filter list
	if (filters.empty()) {
		ptrOutputFiles->writeWarningMessage("No filters provided for channel response");
		response.is_valid = false;
		response.error_message = "No filters provided";
		return response;
	}
	
	// Check unit consistency
	std::string errorMessage;
	if (!validateFilterUnits(filters, errorMessage)) {
		response.is_valid = false;
		response.error_message = errorMessage;
		return response;
	}
	
	// Set input and output units from first and last filters
	response.units_in = filters.front().units_in;
	response.units_out = filters.back().units_out;
	
	// Set normalization frequency (default to 1.0 Hz if not specified)
	// In a more complete implementation, this could be computed from pass band
	response.normalization_frequency = 1.0;
	
	response.is_valid = true;
	
	ptrOutputFiles->writeLogMessage("Created channel response with " + 
	                                Util::toString(static_cast<int>(filters.size())) + 
	                                " filter(s)");
	ptrOutputFiles->writeLogMessage("  Input units:  " + response.units_in);
	ptrOutputFiles->writeLogMessage("  Output units: " + response.units_out);
	
	// Log the filter sequence for verification
	for (const auto& filter : filters) {
		ptrOutputFiles->writeLogMessage("  [" + Util::toString(filter.sequence_number) + "] " + 
		                                filter.name + " (" + filter.type + "): " + 
		                                filter.units_in + " -> " + filter.units_out);
	}
	
	return response;
}

// Combine filters for a channel into a complete channel response
// This is a convenience method that combines getChannelFilters() and createChannelResponse()
// into a single call, similar to the channel_response property in Python's ChannelDataset.
//
// Example usage:
//   ChannelResponse response = mth5->getChannelResponse("test.h5", "/Survey/Stations/MT001/MT001a/Ex");
//   if (response.is_valid) {
//       std::cout << "Channel response has " << response.getFilterCount() << " filter(s)" << std::endl;
//       std::cout << "Units: " << response.units_in << " -> " << response.units_out << std::endl;
//   }
ChannelResponse MTH5::getChannelResponse(const std::string& fileName, const std::string& channelPath) const {
	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	
	ptrOutputFiles->writeLogMessage("Getting complete channel response for: " + channelPath);
	
	// Step 1: Get all filters for the channel
	std::vector<FilterInfo> filters = getChannelFilters(fileName, channelPath);
	
	// Step 2: Combine them into a channel response
	ChannelResponse response = createChannelResponse(filters);
	
	if (!response.is_valid) {
		ptrOutputFiles->writeErrorMessage("Failed to create valid channel response: " + response.error_message);
	} else {
		ptrOutputFiles->writeLogMessage("Successfully created channel response");
	}
	
	return response;
}
