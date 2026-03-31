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
#include "MTH5ChannelResponse.h"
#include "Util.h"
#include "OutputFiles.h"
#include "MTH5CoefficientFilter.h"
#include "MTH5ZeroPoleGainFilter.h"
#include "MTH5TimeDelayFilter.h"
#include "MTH5FrequencyResponseTableFilter.h"
#include "MTH5FIRFilter.h"

#include <algorithm>
#include <sstream>
#include <iomanip>

// Default constructer
MTH5ChannelResponse::MTH5ChannelResponse():
	m_isAppliedForwardly(false),
	m_normalizationFrequency(0.0),
	m_isValid(false)
{
}

// Destructer
MTH5ChannelResponse::~MTH5ChannelResponse(){
	for (std::vector<MTH5Filter*>::iterator itr = m_filtersList.begin(); itr != m_filtersList.end(); ++itr) {
		if (*itr != NULL) {
			delete* itr;
		}
	}
}

// Get list of filter names
std::vector<std::string> MTH5ChannelResponse::getFilterNames() const {
	std::vector<std::string> names;
	for(std::vector<MTH5Filter*>::const_iterator itr = m_filtersList.begin(); itr != m_filtersList.end(); ++itr){
		names.push_back((*itr)->getName());
	}
	return names;
}

// Get total number of filters
int MTH5ChannelResponse::getFilterCount() const {
	return static_cast<int>(m_filtersList.size());
}

// Check if response contains time delay filters
bool MTH5ChannelResponse::hasDelayFilters() const {
	for(std::vector<MTH5Filter*>::const_iterator itr = m_filtersList.begin(); itr != m_filtersList.end(); ++itr){
		if ((*itr)->getType() == "time_delay") {
			return true;
		}
	}
	return false;
}

// Verify filters are properly ordered by sequence number
bool MTH5ChannelResponse::isProperlyOrdered() const {
	const int numOfFilters = static_cast<int>(m_filtersList.size());
	for (int i = 1; i < numOfFilters; ++i) {
		if (m_filtersList[i]->getSequenceNumber() <= m_filtersList[i - 1]->getSequenceNumber()) {
			return false;
		}
	}
	return true;
}

// Get input units from the first filter
std::string MTH5ChannelResponse::getUnitsIn() const {
	return m_filtersList.front()->getUnitsIn();
}

// Get output units from the last filter
std::string MTH5ChannelResponse::getUnitsOut() const {
	return m_filtersList.back()->getUnitsOut();
}

// Get normalization frequency in Hz
double MTH5ChannelResponse::getNormalizationFrequency() const {
	return m_normalizationFrequency;
}

// When going from physical units to digital counts we are applying filters: true
// When going from counts to physical units we are unapplying the filters: false
bool MTH5ChannelResponse::isAppliedForwardly() const {
	return m_isAppliedForwardly;
}

// Whether response is valid
bool MTH5ChannelResponse::isValid() const {
	return m_isValid;
}

// Get error message if invalid
std::string MTH5ChannelResponse::getErrorMessage() const {
	return m_errorMessage;
}

// Get all filters for a channel and combine filters into a complete channel response
// This method validates the filter chain and combines them into a single ChannelResponse object,
// similar to the ChannelResponse class in Python. The response includes:
// - All filters in the proper sequence
// - Input units from the first filter
// - Output units from the last filter
// - Validation status and any error messages
void MTH5ChannelResponse::createChannelResponse(const std::string& fileName, const std::string& channelPath){

	std::string fileNameLower = fileName;
	transform(fileNameLower.begin(), fileNameLower.end(), fileNameLower.begin(), tolower);
	std::string channelPathLower = channelPath;
	transform(channelPathLower.begin(), channelPathLower.end(), channelPathLower.begin(), tolower);
	if (fileNameLower == "null" || fileNameLower == "none") {
		return;
	}
	if (channelPathLower == "null" || channelPathLower == "none") {
		return;
	}

	(OutputFiles::getInstance())->writeLogMessage("Getting complete channel response for: " + channelPath);

	// Step 1: Get all filters for the channel
	getChannelFilters(fileName, channelPath, m_filtersList);

	// Step 2: Combine them into a channel response
	combineFilters();

}

// Calculate frequency response functions using the requency response functions of all filter
std::complex<double> MTH5ChannelResponse::calcResponse(const double freq) const {

	std::complex<double> resp(1.0, 0.0);
	for (std::vector<MTH5Filter*>::const_iterator itr = m_filtersList.begin(); itr != m_filtersList.end(); ++itr) {
		resp *= (*itr)->getFrequencyResponse(freq);
	}
	return resp;

}

// Combine filters for a channel into a complete channel response
std::string MTH5ChannelResponse::combineFilters() {

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();

	std::string errorMessage;
	// Validate the filter list
	if (m_filtersList.empty()) {
		ptrOutputFiles->writeWarningMessage("No filters provided for channel response");
		m_isValid = false;
		errorMessage ="No filters provided";
		return errorMessage;
	}

	// Check unit consistency
	if (!validateFilterUnits(m_filtersList, errorMessage)) {
		m_isValid = false;
		return errorMessage;
	}

	// Set normalization frequency (default to 1.0 Hz if not specified)
	// In a more complete implementation, this could be computed from pass band
	m_normalizationFrequency = 1.0;
	m_isValid = true;

	ptrOutputFiles->writeLogMessage("Created channel response with " +
		Util::toString(static_cast<int>(m_filtersList.size())) +
		" filter(s)");
	ptrOutputFiles->writeLogMessage("  Input units:  " + getUnitsIn());
	ptrOutputFiles->writeLogMessage("  Output units: " + getUnitsOut());

	// Log the filter sequence for verification
	for (std::vector<MTH5Filter*>::const_iterator itr = m_filtersList.begin(); itr != m_filtersList.end(); ++itr) {
		const MTH5Filter* ptr = *itr;
		ptrOutputFiles->writeLogMessage("  [" + Util::toString(ptr->getSequenceNumber()) + "] " +
			ptr->getName() + " (" + ptr->getType() + "): " +
			ptr->getUnitsIn() + " -> " + ptr->getUnitsOut());
	}

	return "";

}

// Get all filters for a channel
// Example usage:
//   getChannelFilters("test.h5", "/Survey/Stations/MT001/MT001a/Ex", filters);
void MTH5ChannelResponse::getChannelFilters(const std::string& fileName, const std::string& channelPath, std::vector<MTH5Filter*>& filters) const {

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeLogMessage("Getting channel filters for: " + channelPath);

	// Step 1: Read filter names from channel attributes
	std::vector<std::string> filterNames = readFilterNamesFromChannel(fileName, channelPath);

	if (filterNames.empty()) {
		ptrOutputFiles->writeLogMessage("No filters to load");
		return;
	}

	// Step 2: Determine the Filters group path
	const std::string filtersPath = getFiltersGroupPath(channelPath);
	ptrOutputFiles->writeLogMessage("Filters group path: " + filtersPath);

	// Step 3: Read each filter
	const int numOfFilters = static_cast<int>(filterNames.size());
	for (int i = 0; i < numOfFilters; ++i) {
		MTH5Filter* ptrFilter = readFilter(fileName, filtersPath, filterNames[i], i + 1);
		if (!ptrFilter->getType().empty()) {
			filters.push_back(ptrFilter);
		}
	}

	// Step 4: Sort filters by sequence_number to ensure correct order
	// This is critical because filter order determines the signal processing chain
	if (filters.size() > 1) {
		sortFilters(filters);
		ptrOutputFiles->writeLogMessage("Sorted filters by sequence number");
	}

	ptrOutputFiles->writeLogMessage("Loaded " + Util::toString(static_cast<int>(filters.size())) + " filter(s)");

}

// Read a single filter from the Filters group
MTH5Filter* MTH5ChannelResponse::readFilter(const std::string& fileName, const std::string& filterPath, const std::string& filterName, int sequenceNumber) const {

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
		return NULL;
	}

	// Check all filter type subdirectories
	std::vector<std::string> filterTypes;
	filterTypes.push_back("zpk");
	filterTypes.push_back("coefficient");
	filterTypes.push_back("time_delay");
	filterTypes.push_back("fap");
	filterTypes.push_back("fir");

	bool found(false);
	MTH5Filter* ptrFilter(NULL);
	for (std::vector<std::string>::const_iterator itr = filterTypes.begin(); itr != filterTypes.end(); ++itr) {
		const std::string ftype = *itr;
		std::string fullPath = filterPath + "/" + ftype + "/" + normalizedName;

		// Check if this path exists
		H5E_BEGIN_TRY{
			hid_t filter_id = H5Gopen(file_id, fullPath.c_str(), H5P_DEFAULT);
			if (filter_id >= 0) {
				found = true;
				// Read filter data (if it exists as a dataset)
				// Different filter types have different data structures
				if (ftype == "zpk") {
					ptrFilter = new MTH5ZeroPoleGainFilter();
					MTH5ZeroPoleGainFilter* ptrZPK = dynamic_cast<MTH5ZeroPoleGainFilter*>(ptrFilter);
					// @note under construction
					// Read normalization factor (gain)
					double normFactor = 1.0;
					if (H5Aexists(filter_id, "gain") > 0) {
						hid_t attr_id = H5Aopen(filter_id, "gain", H5P_DEFAULT);
						if (attr_id >= 0) {
							H5Aread(attr_id, H5T_NATIVE_DOUBLE, &normFactor);
							H5Aclose(attr_id);
						}
					}
					ptrZPK->setNormalizationFactor(normFactor);
					// Compound type matching numpy complex128 stored by h5py ("r" and "i" fields)
					hid_t ctype = H5Tcreate(H5T_COMPOUND, 2 * sizeof(double));
					H5Tinsert(ctype, "r", 0, H5T_NATIVE_DOUBLE);
					H5Tinsert(ctype, "i", sizeof(double), H5T_NATIVE_DOUBLE);
					// Read zeros
					if (H5Lexists(filter_id, "zeros", H5P_DEFAULT) > 0) {
						hid_t ds_id = H5Dopen(filter_id, "zeros", H5P_DEFAULT);
						if (ds_id >= 0) {
							hid_t sp_id = H5Dget_space(ds_id);
							hsize_t dims[1];
							dims[0] = 0;
							H5Sget_simple_extent_dims(sp_id, dims, NULL);
							int n = static_cast<int>(dims[0]);
							if (n > 0) {
								double* buf = new double[2 * n];
								H5Dread(ds_id, ctype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
								std::vector< std::complex<double> > zeros;
								for (int k = 0; k < n; ++k) {
									zeros.push_back(std::complex<double>(buf[2 * k], buf[2 * k + 1]));
								}
								delete[] buf;
								ptrZPK->setZeros(zeros);
							}
							H5Sclose(sp_id);
							H5Dclose(ds_id);
						}
					}
					// Read poles
					if (H5Lexists(filter_id, "poles", H5P_DEFAULT) > 0) {
						hid_t ds_id = H5Dopen(filter_id, "poles", H5P_DEFAULT);
						if (ds_id >= 0) {
							hid_t sp_id = H5Dget_space(ds_id);
							hsize_t dims[1];
							dims[0] = 0;
							H5Sget_simple_extent_dims(sp_id, dims, NULL);
							int n = static_cast<int>(dims[0]);
							if (n > 0) {
								double* buf = new double[2 * n];
								H5Dread(ds_id, ctype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
								std::vector< std::complex<double> > poles;
								for (int k = 0; k < n; ++k) {
									poles.push_back(std::complex<double>(buf[2 * k], buf[2 * k + 1]));
								}
								delete[] buf;
								ptrZPK->setPoles(poles);
							}
							H5Sclose(sp_id);
							H5Dclose(ds_id);
						}
					}
					H5Tclose(ctype);
				}
				else if (ftype == "coefficient") {
					ptrFilter = new MTH5CoefficientFilter();
					MTH5CoefficientFilter* ptrCoef = dynamic_cast<MTH5CoefficientFilter*>(ptrFilter);
					// Read scalar gain attribute
					double gain = 1.0;
					if (H5Aexists(filter_id, "gain") > 0) {
						hid_t attr_id = H5Aopen(filter_id, "gain", H5P_DEFAULT);
						if (attr_id >= 0) {
							H5Aread(attr_id, H5T_NATIVE_DOUBLE, &gain);
							H5Aclose(attr_id);
						}
					}
					ptrCoef->setGain(gain);
				}
				else if (ftype == "fir") {
					ptrFilter = new MTH5FIRFilter();
					MTH5FIRFilter* ptrFIR = dynamic_cast<MTH5FIRFilter*>(ptrFilter);
					// @note under construction
					// Read gain
					double gain = 1.0;
					if (H5Aexists(filter_id, "gain") > 0) {
						hid_t attr_id = H5Aopen(filter_id, "gain", H5P_DEFAULT);
						if (attr_id >= 0) {
							H5Aread(attr_id, H5T_NATIVE_DOUBLE, &gain);
							H5Aclose(attr_id);
						}
					}
					ptrFIR->setGain(gain);
					// Read symmetry type: "even" -> EVEN, "odd" -> ADD, otherwise ASYMMETRIC
					std::string symmetry = readStringAttribute(filter_id, "symmetry");
					if (symmetry == "odd") {
						ptrFIR->setTypeOfSymmetry(MTH5FIRFilter::ADD);
					} else if (symmetry == "none") {
						ptrFIR->setTypeOfSymmetry(MTH5FIRFilter::ASYMMETRIC);
					} else {
						ptrFIR->setTypeOfSymmetry(MTH5FIRFilter::EVEN);
					}
					// Read decimation input sample rate
					double decimInputRate = 1.0;
					if (H5Aexists(filter_id, "decimation_input_sample_rate") > 0) {
						hid_t attr_id = H5Aopen(filter_id, "decimation_input_sample_rate", H5P_DEFAULT);
						if (attr_id >= 0) {
							H5Aread(attr_id, H5T_NATIVE_DOUBLE, &decimInputRate);
							H5Aclose(attr_id);
						}
					}
					ptrFIR->setDecimationInputSampleRate(static_cast<int>(decimInputRate));
					// Read decimation factor
					int decimFactor = 1;
					if (H5Aexists(filter_id, "decimation_factor") > 0) {
						hid_t attr_id = H5Aopen(filter_id, "decimation_factor", H5P_DEFAULT);
						if (attr_id >= 0) {
							H5Aread(attr_id, H5T_NATIVE_INT, &decimFactor);
							H5Aclose(attr_id);
						}
					}
					ptrFIR->setDecimationFactor(decimFactor);
					// Read FIR coefficients (try "numerator_coefficients" then "coefficients")
					const char* coeffNames[2] = { "numerator_coefficients", "coefficients" };
					for (int ci = 0; ci < 2; ++ci) {
						if (H5Lexists(filter_id, coeffNames[ci], H5P_DEFAULT) > 0) {
							hid_t ds_id = H5Dopen(filter_id, coeffNames[ci], H5P_DEFAULT);
							if (ds_id >= 0) {
								hid_t sp_id = H5Dget_space(ds_id);
								hsize_t dims[1];
								dims[0] = 0;
								H5Sget_simple_extent_dims(sp_id, dims, NULL);
								int n = static_cast<int>(dims[0]);
								if (n > 0) {
									double* buf = new double[n];
									H5Dread(ds_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
									std::vector<double> coefficients;
									for (int k = 0; k < n; ++k) {
										coefficients.push_back(buf[k]);
									}
									delete[] buf;
									ptrFIR->setFIRCoefficients(coefficients);
								}
								H5Sclose(sp_id);
								H5Dclose(ds_id);
							}
							break;
						}
					}
				}
				else if (ftype == "time_delay") {
					ptrFilter = new MTH5TimeDelayFilter();
					MTH5TimeDelayFilter* ptrTD = dynamic_cast<MTH5TimeDelayFilter*>(ptrFilter);
					// @note under construction
					// Read delay in seconds
					double delay = 0.0;
					if (H5Aexists(filter_id, "delay") > 0) {
						hid_t attr_id = H5Aopen(filter_id, "delay", H5P_DEFAULT);
						if (attr_id >= 0) {
							H5Aread(attr_id, H5T_NATIVE_DOUBLE, &delay);
							H5Aclose(attr_id);
						}
					}
					ptrTD->setDelay(delay);
				}
				else if (ftype == "fap") {
					ptrFilter = new MTH5FrequencyResponseTableFilter();
					MTH5FrequencyResponseTableFilter* ptrFAP = dynamic_cast<MTH5FrequencyResponseTableFilter*>(ptrFilter);
					// FAP data is stored as a single compound dataset "fap_table" with
					// fields "frequency", "amplitude", "phase" (all double).
					if (H5Lexists(filter_id, "fap_table", H5P_DEFAULT) > 0) {
						hid_t ds_id = H5Dopen(filter_id, "fap_table", H5P_DEFAULT);
						if (ds_id >= 0) {
							hid_t sp_id = H5Dget_space(ds_id);
							hsize_t dims[1];
							dims[0] = 0;
							H5Sget_simple_extent_dims(sp_id, dims, NULL);
							int n = static_cast<int>(dims[0]);
							if (n > 0) {
								// Memory compound type: three doubles packed as {freq, amp, phase}
								hid_t mem_type = H5Tcreate(H5T_COMPOUND, 3 * sizeof(double));
								H5Tinsert(mem_type, "frequency",                0, H5T_NATIVE_DOUBLE);
								H5Tinsert(mem_type, "amplitude",    sizeof(double), H5T_NATIVE_DOUBLE);
								H5Tinsert(mem_type, "phase",    2 * sizeof(double), H5T_NATIVE_DOUBLE);
								double* buf = new double[3 * n];
								H5Dread(ds_id, mem_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
								std::vector<double> freqs, amps, phases;
								for (int k = 0; k < n; ++k) {
									freqs.push_back(buf[3 * k]);
									amps.push_back(buf[3 * k + 1]);
									phases.push_back(buf[3 * k + 2]);
								}
								delete[] buf;
								H5Tclose(mem_type);
								ptrFAP->setFrequencies(freqs);
								ptrFAP->setAmplitude(amps);
								ptrFAP->setPhases(phases);
							}
							H5Sclose(sp_id);
							H5Dclose(ds_id);
						}
					}
				}
				ptrFilter->setType(ftype);
				// Read filter metadata attributes
				ptrFilter->setUnitsIn(readStringAttribute(filter_id, "units_in"));
				ptrFilter->setUnitsOut(readStringAttribute(filter_id, "units_out"));
				ptrFilter->setSequenceNumber(sequenceNumber);
				ptrFilter->setName(filterName);
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

	return ptrFilter;
}

// Validate filter list for unit consistency
// This ensures that the output units of each filter match the input units of the next filter,
// similar to the _check_consistency_of_units method in Python's ChannelResponse class.
bool MTH5ChannelResponse::validateFilterUnits(const std::vector<MTH5Filter*>& filters, std::string& errorMessage) const {
	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();

	if (filters.empty()) {
		return true; // Empty filter list is valid
	}

	if (filters.size() == 1) {
		return true; // Single filter is always valid
	}

	// Check that output units of one filter match input units of the next
	// Also verify that sequence numbers are sequential
	const int numFilters = static_cast<int>(filters.size());
	for (int i = 0; i < numFilters - 1; ++i) {
		const MTH5Filter* currentFilter = filters[i];
		const MTH5Filter* nextFilter = filters[i + 1];

		// Verify sequence numbers are in order
		if (currentFilter->getSequenceNumber() >= nextFilter->getSequenceNumber()) {
			std::ostringstream oss;
			oss << "Filter sequence error: Filter '" << currentFilter->getName()
				<< "' has sequence " << currentFilter->getSequenceNumber()
				<< " but next filter '" << nextFilter->getName()
				<< "' has sequence " << nextFilter->getSequenceNumber();
			errorMessage = oss.str();
			ptrOutputFiles->writeErrorMessage(errorMessage);
			return false;
		}

		if (currentFilter->getUnitsOut() != nextFilter->getUnitsIn()) {
			std::ostringstream oss;
			oss << "Unit consistency error: Filter '" << currentFilter->getName()
				<< "' (seq " << currentFilter->getSequenceNumber() << ") outputs " << currentFilter->getUnitsOut()
				<< " but filter '" << nextFilter->getName()
				<< "' (seq " << nextFilter->getSequenceNumber() << ") expects " << nextFilter->getUnitsIn();
			errorMessage = oss.str();
			ptrOutputFiles->writeErrorMessage(errorMessage);
			return false;
		}
	}

	ptrOutputFiles->writeLogMessage("Filter units are consistent");
	return true;
}

// Helper: Navigate channel path to Filters group path
std::string MTH5ChannelResponse::getFiltersGroupPath(const std::string& channelPath) const {
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
std::vector<std::string> MTH5ChannelResponse::readFilterNamesFromChannel(const std::string& fileName, const std::string& channelPath) const {
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
	if (filterNames.empty()) {
		ptrOutputFiles->writeLogMessage("No filters found in channel attributes");
	}
	else {
		ptrOutputFiles->writeLogMessage("Found " + Util::toString(static_cast<int>(filterNames.size())) + " filter(s)");
	}

	H5Dclose(dataset_id);
	H5Fclose(file_id);

	return filterNames;
}

// Helper: Read string attribute from HDF5 object
std::string MTH5ChannelResponse::readStringAttribute(hid_t obj_id, const std::string& attrName) const {
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
	std::string result;

	if (H5Tis_variable_str(type_id)) {
		// Variable-length string: H5Tget_size() returns pointer size, not string length.
		// Copy the file's own type so encoding (UTF-8 etc.) is preserved and the read succeeds.
		hid_t mem_type = H5Tcopy(type_id);
		char* vlen_str = NULL;
		if (H5Aread(attr_id, mem_type, &vlen_str) >= 0 && vlen_str != NULL) {
			result = std::string(vlen_str);
			hid_t space_id = H5Aget_space(attr_id);
#if (H5_VERS_MAJOR <= 1 && H5_VERS_MINOR < 12)
			H5Dvlen_reclaim(mem_type, space_id, H5P_DEFAULT, &vlen_str);
#else
			H5Treclaim(mem_type, space_id, H5P_DEFAULT, &vlen_str);
#endif
			H5Sclose(space_id);
		}
		H5Tclose(mem_type);
	}
	else {
		// Fixed-length string: H5Tget_size() gives the actual byte length.
		size_t size = H5Tget_size(type_id);
		char* buffer = new char[size + 1];
		buffer[size] = '\0';
		H5Aread(attr_id, type_id, buffer);
		result = std::string(buffer, size);
		// Strip null padding
		size_t nullPos = result.find('\0');
		if (nullPos != std::string::npos) {
			result = result.substr(0, nullPos);
		}
		delete[] buffer;
	}

	H5Tclose(type_id);
	H5Aclose(attr_id);

	return result;
}

// Helper: Read string array attribute from HDF5 object
std::vector<std::string> MTH5ChannelResponse::readStringArrayAttribute(hid_t obj_id, const std::string& attrName) const {
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
			char** rdata = new char* [npoints];
			H5Aread(attr_id, type_id, rdata);

			for (int i = 0; i < npoints; i++) {
				if (rdata[i] != NULL) {
					result.push_back(std::string(rdata[i]));
				}
			}
#if (H5_VERS_MAJOR <= 1 && H5_VERS_MINOR < 12)
			H5Dvlen_reclaim(type_id, space_id, H5P_DEFAULT, rdata);
#else
			H5Treclaim(type_id, space_id, H5P_DEFAULT, rdata);
#endif
			delete[] rdata;
		}
		// Handle fixed length strings
		else {
			size_t size = H5Tget_size(type_id);
			char* buffer = new char[npoints * size];
			H5Aread(attr_id, type_id, buffer);

			for (int i = 0; i < npoints; i++) {
				// Limit to exactly `size` bytes to avoid reading into the next element,
				// then strip any trailing null padding.
				std::string str(&buffer[i * size], size);
				size_t nullPos = str.find('\0');
				if (nullPos != std::string::npos) {
					str = str.substr(0, nullPos);
				}
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

// Sort filters by ascending order using bubble sort
void MTH5ChannelResponse::sortFilters(std::vector<MTH5Filter*>& filters) const {

	const int numFilters = static_cast<int>(filters.size());
	for (int i = 0; i < numFilters - 1; ++i) {
		for (int j = 0; j < numFilters - i - 1; ++j) {
			if (filters[j]->getSequenceNumber() > filters[j + 1]->getSequenceNumber()) {
				std::swap(filters[j], filters[j + 1]);
			}
		}
	}

}