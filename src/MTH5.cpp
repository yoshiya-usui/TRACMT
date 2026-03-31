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
#include <iomanip>
#include <assert.h>

// Default constructer
MTH5::MTH5():
	m_numOfChannelRespones(0)
{
}

// Destructer
MTH5::~MTH5(){
	
	for (std::vector<MTH5ChannelResponse*>::iterator itr = m_channelResponses.begin(); itr != m_channelResponses.end(); ++itr) {
		if ( *itr != NULL) {
			delete[] * itr;
		}
	}

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

// Get name of the calibration file name made from the channel responses 
std::string MTH5::getCalibrationFileName(const int channelIndex){

	std::ostringstream fileName;
	fileName << "channel" << channelIndex << ".cal";
	return fileName.str();

}

// Calculate frequency response functions using the frequency response functions of all filter
std::complex<double> MTH5::calcResponse(const int sectionIndex, const int channelIndex, const double freq) const {

	assert(sectionIndex >= 0 && sectionIndex < m_channelResponses.size());
	assert(channelIndex >= 0);

	return m_channelResponses[sectionIndex][channelIndex].calcResponse(freq);

}

// Read filters for indivial sections and channels
void MTH5::readFiltersAll(const int numChannels, const std::vector<CommonParameters::DataFileSet>& dataFileSets) {

	m_numOfChannelRespones = numChannels;

	for(std::vector<CommonParameters::DataFileSet>::const_iterator itr = dataFileSets.begin(); itr != dataFileSets.end(); ++itr){
		const std::vector<CommonParameters::DataFile>& dataFileList = itr->dataFile;
		MTH5ChannelResponse* channelResponses = new MTH5ChannelResponse[numChannels];
		for (int iChan = 0; iChan < numChannels; ++iChan) {
			const std::string fileName = dataFileList[iChan].fileName;
			const std::string channelPath = dataFileList[iChan].mth5GroupName;
			channelResponses[iChan].createChannelResponse(fileName, channelPath);
		}
		m_channelResponses.push_back(channelResponses);
	}

}
