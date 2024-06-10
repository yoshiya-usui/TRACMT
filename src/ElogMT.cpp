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
#include "ElogMT.h"
#include "Util.h"
#include "OutputFiles.h"
#include "CommonParameters.h"
#include "Control.h"

#include <iomanip>
#include <iostream>
#include <assert.h>
#include <complex>
#ifdef _USE_FILESYSTEM
#include <filesystem>
#endif

const static double MV_LSB_E = 2500. / 8388608.;
const static double MV_LSB_H = 10000. / 8388608.;

// Return the instance of the class
ElogMT* ElogMT::getInstance(){
   	static ElogMT instance;// The only instance
  	return &instance;
}

// Read ELOG-MT binary file
void ElogMT::readElogBinaryFile(const std::string& fileName, const int numSkipData, const int numDataPoints, int& counter, double* ex, double* ey, double* hz, double* hx, double* hy) const {

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeLogMessage("Read Ex, Ey, Hz, Hx, and Hy data from " + fileName);

	FILE* fp = fopen(fileName.c_str(), "rb");
	if (fp == NULL) {
		ptrOutputFiles->writeErrorMessage("File open error : " + fileName);
	}

	// Determine file type
	const Control* const ptrControl = Control::getInstance();
	const int freq = static_cast<int>(ptrControl->getSamplingFrequency());

	const int data_size = AD_CH * freq * AD_BYTES;
	unsigned char* data = new unsigned char[data_size];

	while (1) {
		// Read header of each block
		recdata_head	d;
		if (fread(&d, 1, REC_DATALEN_HEAD, fp) < 1) {
			break;
		}
		// Read AD data of each block
		if (fread(data, 1, data_size, fp) < 1) {
			break;
		}
		if (counter + freq <= numSkipData) {
			counter += freq;
			// Go to next block
			continue;
		}
		// Determine end position prior to the start position because counter can be modified in determing start position
		int endPos(freq);
		if (numSkipData + numDataPoints > counter && numSkipData + numDataPoints <= counter + freq) {
			// End position locates in this block
			endPos = numSkipData + numDataPoints - counter;
		}
		int startPos(0);
		if (numSkipData > counter && numSkipData <= counter + freq) {
			// Start position locates in this block
			startPos = numSkipData - counter;
			counter = numSkipData;
		}
		for (int j = startPos; j < endPos; ++j, ++counter) {
			long lbuf(0);
			// Ex
			lbuf = b3_to_long32(&data[j * AD_BYTES * AD_CH]);
			ex[counter - numSkipData] = lbuf * MV_LSB_E;
			// Ey
			lbuf = b3_to_long32(&data[j * AD_BYTES * AD_CH + AD_BYTES]);
			ey[counter - numSkipData] = lbuf * MV_LSB_E;
			// Hx
			lbuf = b3_to_long32(&data[j * AD_BYTES * AD_CH + AD_BYTES * 2]);
			hx[counter - numSkipData] = lbuf * MV_LSB_H;
			// Hy
			lbuf = b3_to_long32(&data[j * AD_BYTES * AD_CH + AD_BYTES * 3]);
			hy[counter - numSkipData] = lbuf * MV_LSB_H;
			// Hz
			lbuf = b3_to_long32(&data[j * AD_BYTES * AD_CH + AD_BYTES * 4]);
			hz[counter - numSkipData] = lbuf * MV_LSB_H;
		}
		if (counter - numSkipData >= numDataPoints) {
			break;
		}
	}

	delete[] data;

	fclose(fp);

}

// Read ELOG-MT binary file under a directory
void ElogMT::readElogBinaryFilesUnderADirectory(const std::string& directoryName, const int numSkipData, const int numDataPoints, double* ex, double* ey, double* hz, double* hx, double* hy) const {

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
#ifdef _USE_FILESYSTEM
	const Control* const ptrControl = Control::getInstance();
	const int freq = static_cast<int>(ptrControl->getSamplingFrequency());
	std::ostringstream oss;
	oss << "_" << freq << "Hz";
	const std::string stringCompared = oss.str();
	if (std::filesystem::is_directory(directoryName)) {
		auto dirItr = std::filesystem::directory_iterator(directoryName);
		int counter(0);
		for (auto& p : dirItr) {
			const std::string fileName = p.path().string();
			if (fileName.find(stringCompared) != std::string::npos && Util::extractExtensionOfFileName(fileName).find("dat") != std::string::npos ) {
				readElogBinaryFile(fileName, numSkipData, numDataPoints, counter, ex, ey, hz, hx, hy);
				if (counter >= numSkipData + numDataPoints) {
					break;
				}
			}
		}
	}else{
		ptrOutputFiles->writeErrorMessage("There is no directory of ELOG-MT data.");
	}
#else
	ptrOutputFiles->writeErrorMessage("Multiple inputs of ELOG data under a directory is not supported by this version.");
#endif

}

// Read ELOG-MT binary file (Ex and Ey)
void ElogMT::readElogBinaryFileExEyOnly(const std::string& fileName, const int numSkipData, const int numDataPoints, int& counter, double* ex, double* ey) const {

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeLogMessage("Read Ex and Ey data from "+fileName);

	FILE	* fp = fopen(fileName.c_str(), "rb");
	if(fp == NULL) {
		ptrOutputFiles->writeErrorMessage("File open error : " + fileName);
	}

	// Determine file type
	const Control* const ptrControl = Control::getInstance();
	const int freq = static_cast<int>(ptrControl->getSamplingFrequency());
	const int data_size = AD_CH * freq * AD_BYTES;
	unsigned char* data = new unsigned char[data_size];

	while(1) {
		// Read header of each block
		recdata_head	d;
		if(fread(&d, 1, REC_DATALEN_HEAD, fp) < 1){
			break;
		}
		// Read AD data of each block
		if(fread(data, 1, data_size, fp) < 1){
			break;
		}
		if( counter + freq <= numSkipData ){
			counter += freq;
			// Go to next block
			continue;
		}
		// Determine end position prior to the start position because counter can be modified in determing start position
		int endPos(freq);
		if( numSkipData + numDataPoints > counter && numSkipData + numDataPoints <= counter + freq ){
			// End position locates in this block
			endPos = numSkipData + numDataPoints - counter;
		}
		int startPos(0);
		if( numSkipData > counter && numSkipData <= counter + freq ){
			// Start position locates in this block
			startPos = numSkipData - counter;
			counter = numSkipData;
		}
		for( int j = startPos; j < endPos; ++j, ++counter ){
			// Ex
			long lbufx = b3_to_long32( &data[j*AD_BYTES*AD_CH] );
			ex[counter - numSkipData] = lbufx * MV_LSB_E;
			// Ey
			long lbufy = b3_to_long32( &data[j*AD_BYTES*AD_CH + AD_BYTES] );
			ey[counter - numSkipData] = lbufy * MV_LSB_E;
		}
		if( counter - numSkipData >= numDataPoints ){
			break;
		}
	}

	delete [] data;

	fclose(fp);

}

// Read ELOG-MT binary file under a directory (two electric field data only)
void ElogMT::readElogBinaryFilesUnderADirectoryExEyOnly(const std::string& directoryName, const int numSkipData, const int numDataPoints, double* ex, double* ey) const {

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
#ifdef _USE_FILESYSTEM
	const Control* const ptrControl = Control::getInstance();
	const int freq = static_cast<int>(ptrControl->getSamplingFrequency());
	std::ostringstream oss;
	oss << "_" << freq << "Hz";
	const std::string stringCompared = oss.str();
	if (std::filesystem::is_directory(directoryName)) {
		auto dirItr = std::filesystem::directory_iterator(directoryName);
		int counter(0);
		for (auto& p : dirItr) {
			const std::string fileName = p.path().string();
			if (fileName.find(stringCompared) != std::string::npos && Util::extractExtensionOfFileName(fileName).find("dat") != std::string::npos) {
				readElogBinaryFileExEyOnly(fileName, numSkipData, numDataPoints, counter, ex, ey);
				if (counter >= numSkipData + numDataPoints) {
					break;
				}
			}
		}
	} else {
		ptrOutputFiles->writeErrorMessage("There is no directory of ELOG-MT data.");
	}
#else
	ptrOutputFiles->writeErrorMessage("Multiple inputs of ELOG data under a directory is not supported by this version.");
#endif

}

// Read ELOG-MT binary file (Hx and Hy)
void ElogMT::readElogBinaryFileHxHyOnly(const std::string& fileName, const int numSkipData, const int numDataPoints, int& counter, double* hx, double* hy) const {

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeLogMessage("Read Hx and Hy data from " + fileName);

	FILE* fp = fopen(fileName.c_str(), "rb");
	if (fp == NULL) {
		ptrOutputFiles->writeErrorMessage("File open error : " + fileName);
	}

	// Determine file type
	const Control* const ptrControl = Control::getInstance();
	const int freq = static_cast<int>(ptrControl->getSamplingFrequency());

	const int data_size = AD_CH * freq * AD_BYTES;
	unsigned char* data = new unsigned char[data_size];

	while (1) {
		// Read header of each block
		recdata_head	d;
		if (fread(&d, 1, REC_DATALEN_HEAD, fp) < 1) {
			break;
		}
		// Read AD data of each block
		if (fread(data, 1, data_size, fp) < 1) {
			break;
		}
		if (counter + freq <= numSkipData) {
			counter += freq;
			// Go to next block
			continue;
		}
		// Determine end position prior to the start position because counter can be modified in determing start position
		int endPos(freq);
		if (numSkipData + numDataPoints > counter && numSkipData + numDataPoints <= counter + freq) {
			// End position locates in this block
			endPos = numSkipData + numDataPoints - counter;
		}
		int startPos(0);
		if (numSkipData > counter && numSkipData <= counter + freq) {
			// Start position locates in this block
			startPos = numSkipData - counter;
			counter = numSkipData;
		}
		for (int j = startPos; j < endPos; ++j, ++counter) {
			long lbuf(0);
			// Hx
			lbuf = b3_to_long32(&data[j * AD_BYTES * AD_CH + AD_BYTES * 2]);
			hx[counter - numSkipData] = lbuf * MV_LSB_H;
			// Hy
			lbuf = b3_to_long32(&data[j * AD_BYTES * AD_CH + AD_BYTES * 3]);
			hy[counter - numSkipData] = lbuf * MV_LSB_H;
		}
		if (counter - numSkipData >= numDataPoints) {
			break;
		}
	}

	delete[] data;

	fclose(fp);

}

// Read ELOG-MT binary file under a directory(two horizontal magnetic field data only)
void ElogMT::readElogBinaryFilesUnderADirectoryHxHyOnly(const std::string& directoryName, const int numSkipData, const int numDataPoints, double* hx, double* hy) const {

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
#ifdef _USE_FILESYSTEM
	const Control* const ptrControl = Control::getInstance();
	const int freq = static_cast<int>(ptrControl->getSamplingFrequency());
	std::ostringstream oss;
	oss << "_" << freq << "Hz";
	const std::string stringCompared = oss.str();
	if (std::filesystem::is_directory(directoryName)) {
		auto dirItr = std::filesystem::directory_iterator(directoryName);
		int counter(0);
		for (auto& p : dirItr) {
			const std::string fileName = p.path().string();
			if (fileName.find(stringCompared) != std::string::npos && Util::extractExtensionOfFileName(fileName).find("dat") != std::string::npos) {
				readElogBinaryFileHxHyOnly(fileName, numSkipData, numDataPoints, counter, hx, hy);
				if (counter >= numSkipData + numDataPoints) {
					break;
				}
			}
		}
	} else {
		ptrOutputFiles->writeErrorMessage("There is no directory of ELOG-MT data.");
	}
#else
	ptrOutputFiles->writeErrorMessage("Multiple inputs of ELOG data under a directory is not supported by this version.");
#endif

}

// Read ELOG-MT binary file (Ex, Ey, Hx, and Hy only)
void ElogMT::readElogBinaryFileExEyHxHyOnly(const std::string& fileName, const int numSkipData, const int numDataPoints, int& counter, double* ex, double* ey, double* hx, double* hy) const {

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeLogMessage("Read Ex, Ey, Hx, and Hy data from " + fileName);

	FILE* fp = fopen(fileName.c_str(), "rb");
	if (fp == NULL) {
		ptrOutputFiles->writeErrorMessage("File open error : " + fileName);
	}

	// Determine file type
	const Control* const ptrControl = Control::getInstance();
	const int freq = static_cast<int>(ptrControl->getSamplingFrequency());

	const int data_size = AD_CH * freq * AD_BYTES;
	unsigned char* data = new unsigned char[data_size];

	while (1) {
		// Read header of each block
		recdata_head	d;
		if (fread(&d, 1, REC_DATALEN_HEAD, fp) < 1) {
			break;
		}
		// Read AD data of each block
		if (fread(data, 1, data_size, fp) < 1) {
			break;
		}
		if (counter + freq <= numSkipData) {
			counter += freq;
			// Go to next block
			continue;
		}
		// Determine end position prior to the start position because counter can be modified in determing start position
		int endPos(freq);
		if (numSkipData + numDataPoints > counter && numSkipData + numDataPoints <= counter + freq) {
			// End position locates in this block
			endPos = numSkipData + numDataPoints - counter;
		}
		int startPos(0);
		if (numSkipData > counter && numSkipData <= counter + freq) {
			// Start position locates in this block
			startPos = numSkipData - counter;
			counter = numSkipData;
		}
		for (int j = startPos; j < endPos; ++j, ++counter) {
			long lbuf(0);
			// Ex
			lbuf = b3_to_long32(&data[j * AD_BYTES * AD_CH]);
			ex[counter - numSkipData] = lbuf * MV_LSB_E;
			// Ey
			lbuf = b3_to_long32(&data[j * AD_BYTES * AD_CH + AD_BYTES]);
			ey[counter - numSkipData] = lbuf * MV_LSB_E;
			// Hx
			lbuf = b3_to_long32(&data[j * AD_BYTES * AD_CH + AD_BYTES * 2]);
			hx[counter - numSkipData] = lbuf * MV_LSB_H;
			// Hy
			lbuf = b3_to_long32(&data[j * AD_BYTES * AD_CH + AD_BYTES * 3]);
			hy[counter - numSkipData] = lbuf * MV_LSB_H;
		}
		if (counter - numSkipData >= numDataPoints) {
			break;
		}
	}

	delete[] data;

	fclose(fp);

}

// Read ELOG-MT binary file under a directory (Ex, Ey, Hx, and Hy only)
void ElogMT::readElogBinaryFilesUnderADirectoryExEyHxHyOnly(const std::string& directoryName, const int numSkipData, const int numDataPoints, double* ex, double* ey, double* hx, double* hy) const {

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
#ifdef _USE_FILESYSTEM
	const Control* const ptrControl = Control::getInstance();
	const int freq = static_cast<int>(ptrControl->getSamplingFrequency());
	std::ostringstream oss;
	oss << "_" << freq << "Hz";
	const std::string stringCompared = oss.str();
	if (std::filesystem::is_directory(directoryName)) {
		auto dirItr = std::filesystem::directory_iterator(directoryName);
		int counter(0);
		for (auto& p : dirItr) {
			const std::string fileName = p.path().string();
			if (fileName.find(stringCompared) != std::string::npos && Util::extractExtensionOfFileName(fileName).find("dat") != std::string::npos) {
				readElogBinaryFileExEyHxHyOnly(fileName, numSkipData, numDataPoints, counter, ex, ey, hx, hy);
				if (counter >= numSkipData + numDataPoints) {
					break;
				}
			}
		}
	} else {
		ptrOutputFiles->writeErrorMessage("There is no directory of ELOG-MT data.");
	}
#else
	ptrOutputFiles->writeErrorMessage("Multiple inputs of ELOG data under a directory is not supported by this version.");
#endif

}

// Read ELOG-MT binary file (Hz, Hx, and Hy only)
void ElogMT::readElogBinaryFileHzHxHyOnly(const std::string& fileName, const int numSkipData, const int numDataPoints, int& counter, double* hz, double* hx, double* hy) const {

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeLogMessage("Read Hz, Hx, and Hy data from " + fileName);

	FILE* fp = fopen(fileName.c_str(), "rb");
	if (fp == NULL) {
		ptrOutputFiles->writeErrorMessage("File open error : " + fileName);
	}

	// Determine file type
	const Control* const ptrControl = Control::getInstance();
	const int freq = static_cast<int>(ptrControl->getSamplingFrequency());

	const int data_size = AD_CH * freq * AD_BYTES;
	unsigned char* data = new unsigned char[data_size];

	while (1) {
		// Read header of each block
		recdata_head	d;
		if (fread(&d, 1, REC_DATALEN_HEAD, fp) < 1) {
			break;
		}
		// Read AD data of each block
		if (fread(data, 1, data_size, fp) < 1) {
			break;
		}
		if (counter + freq <= numSkipData) {
			counter += freq;
			// Go to next block
			continue;
		}
		// Determine end position prior to the start position because counter can be modified in determing start position
		int endPos(freq);
		if (numSkipData + numDataPoints > counter && numSkipData + numDataPoints <= counter + freq) {
			// End position locates in this block
			endPos = numSkipData + numDataPoints - counter;
		}
		int startPos(0);
		if (numSkipData > counter && numSkipData <= counter + freq) {
			// Start position locates in this block
			startPos = numSkipData - counter;
			counter = numSkipData;
		}
		for (int j = startPos; j < endPos; ++j, ++counter) {
			long lbuf(0);
			// Hx
			lbuf = b3_to_long32(&data[j * AD_BYTES * AD_CH + AD_BYTES * 2]);
			hx[counter - numSkipData] = lbuf * MV_LSB_H;
			// Hy
			lbuf = b3_to_long32(&data[j * AD_BYTES * AD_CH + AD_BYTES * 3]);
			hy[counter - numSkipData] = lbuf * MV_LSB_H;
			// Hz
			lbuf = b3_to_long32(&data[j * AD_BYTES * AD_CH + AD_BYTES * 4]);
			hz[counter - numSkipData] = lbuf * MV_LSB_H;
		}
		if (counter - numSkipData >= numDataPoints) {
			break;
		}
	}

	delete[] data;

	fclose(fp);

}

// Read ELOG-MT binary file under a directory (Hz, Hx, and Hy only)
void ElogMT::readElogBinaryFilesUnderADirectoryHzHxHyOnly(const std::string& directoryName, const int numSkipData, const int numDataPoints, double* hz, double* hx, double* hy) const {

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
#ifdef _USE_FILESYSTEM
	const Control* const ptrControl = Control::getInstance();
	const int freq = static_cast<int>(ptrControl->getSamplingFrequency());
	std::ostringstream oss;
	oss << "_" << freq << "Hz";
	const std::string stringCompared = oss.str();
	if (std::filesystem::is_directory(directoryName)) {
		auto dirItr = std::filesystem::directory_iterator(directoryName);
		int counter(0);
		for (auto& p : dirItr) {
			const std::string fileName = p.path().string();
			if (fileName.find(stringCompared) != std::string::npos && Util::extractExtensionOfFileName(fileName).find("dat") != std::string::npos) {
				readElogBinaryFileHzHxHyOnly(fileName, numSkipData, numDataPoints, counter, hz, hx, hy);
				if (counter >= numSkipData + numDataPoints) {
					break;
				}
			}
		}
	}
	else {
		ptrOutputFiles->writeErrorMessage("There is no directory of ELOG-MT data.");
	}
#else
	ptrOutputFiles->writeErrorMessage("Multiple inputs of ELOG data under a directory is not supported by this version.");
#endif

}

// Get calibration file name
std::string ElogMT::getCalibrationFileName( const int channelIndex ) const{

	std::ostringstream fileName;
	fileName << "channel" << channelIndex << ".cal";

	return fileName.str();

}

// Make calibration file
void ElogMT::makeCalibrationFile( const std::string& fileName, const double unitGroupDelay, 	const std::vector<int>& channelIndexes,
	const double dipoleLengthX, const double dipoleLengthY, const std::vector<double>& freq ) const{

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();

	const int numFreq = static_cast<int>(freq.size());
	if( numFreq < 1 ){
		ptrOutputFiles->writeErrorMessage( "Number of the frequencies for which calibrations are estimated is less than 1 : " + Util::toString(numFreq) );
	}

	std::ofstream ofs[AD_CH];
	for( std::vector<int>::const_iterator itr = channelIndexes.begin(); itr != channelIndexes.end(); ++itr ){
		const int iCh = *itr;
		const std::string outputFileName = getCalibrationFileName(iCh);
		ofs[iCh].open( outputFileName.c_str(), std::ios::out );
		if( ofs[iCh].fail() ){
			ptrOutputFiles->writeErrorMessage( "File open error : " + outputFileName );
		}
		if( iCh == 0 ){
			const double dipoleLength = dipoleLengthX / -1000.0;// [m] -> [km] and invert sign
			ofs[iCh] << std::setw(20) << std::scientific << std::setprecision(9) << 1.0/dipoleLength << std::endl;
		}else if( iCh == 1 ){
			const double dipoleLength = dipoleLengthY / -1000.0;// [m] -> [km] and invert sign
			ofs[iCh] << std::setw(20) << std::scientific << std::setprecision(9) << 1.0/dipoleLength << std::endl;
		}else{
			ofs[iCh] << std::setw(20) << std::scientific << std::setprecision(9) << 1.0 << std::endl;
		}
		ofs[iCh] << std::setw(10) << numFreq << std::endl;
	}

	const Control* const ptrControl = Control::getInstance();
	for( int iFreq = 0; iFreq < numFreq; ++iFreq ){
		std::complex<double> calElog[AD_CH];
		for( int iCh = 0; iCh < AD_CH; ++iCh ){
			calElog[iCh] = std::complex<double>(1.0, 0.0);
		}
		calculateCalibrationFunctionForAnalogFilter(fileName, freq[iFreq], calElog);
		{
			//------- 1st fir filer from measured calibration table
			std::string path = "";
			if (!ptrControl->getDirectoryOfLoggerCalibrationFiles().empty()) {
				path = ptrControl->getDirectoryOfLoggerCalibrationFiles() + "\\";
			}
			std::string fileName = "firh.txt";
			if (ptrControl->getTypeOfElogMT() == Control::ELOGMT_ADU_MODE) {
				fileName = "firh_adu.txt";
			}
			else if (ptrControl->getTypeOfElogMT() == Control::ELOGMT_PHX_MODE) {
				fileName = "firh_phx.txt";
			}
			const std::complex<double> firh = Util::calculateCalibrationForFIRFilterType1(path + fileName, 9, 14336.0, freq[iFreq], true);
			for( int iCh = 0; iCh < AD_CH; ++iCh ){
				calElog[iCh] *= firh;
			}
		}
		{
			//------ 2nd fir filter ---
			// 4900 is described in elog cal file (4.9 Vpp input)
			for( int iCh = 0; iCh < AD_CH; ++iCh ){
				calElog[iCh] /= 4900.0;
			}
		}
		{
			//------ Gain correction ---
			const double tsGroupDelay = 1.0 / 14336.0 * unitGroupDelay;
			const double angle = 2.0 * CommonParameters::PI * freq[iFreq] * tsGroupDelay;
			const std::complex<double> groupDelay = std::complex<double>(cos(angle), sin(angle));
			for( int iCh = 0; iCh < AD_CH; ++iCh ){
				calElog[iCh] /= groupDelay;
			}
		}
		const double samplingFreq = ptrControl->getSamplingFrequency();
		if( fabs(samplingFreq - 32.0) < CommonParameters::EPS ) {
			//------ Group delay correction ---
			std::string path = "";
			if (!ptrControl->getDirectoryOfLoggerCalibrationFiles().empty()) {
				path = ptrControl->getDirectoryOfLoggerCalibrationFiles() + "\\";
			}
			std::string fileName = "firl.txt";
			if (ptrControl->getTypeOfElogMT() == Control::ELOGMT_ADU_MODE) {
				fileName = "firl_adu.txt";
			}
			else if (ptrControl->getTypeOfElogMT() == Control::ELOGMT_PHX_MODE) {
				fileName = "firl_phx.txt";
			}
			const std::complex<double> firl = Util::calculateCalibrationForFIRFilterType2(path + fileName, 9, 1024.0, freq[iFreq], -396, 0);
			for( int iCh = 0; iCh < AD_CH; ++iCh ){
				calElog[iCh] *= firl;
			}
		}
		for( std::vector<int>::const_iterator itr = channelIndexes.begin(); itr != channelIndexes.end(); ++itr ){
			const int iCh = *itr;
			const std::complex<double> invCal = 1.0 / calElog[iCh];
			ofs[iCh] << std::setw(20) << std::scientific << std::setprecision(9) << freq[iFreq];
			ofs[iCh] << std::setw(20) << std::scientific << std::setprecision(9) << invCal.real();
			ofs[iCh] << std::setw(20) << std::scientific << std::setprecision(9) << invCal.imag() << std::endl;
		}
	}

	for( std::vector<int>::const_iterator itr = channelIndexes.begin(); itr != channelIndexes.end(); ++itr ){
		const int iCh = * itr;
		ofs[iCh].close();
	}

}

// Calculate calibration function for analog filter
void ElogMT::calculateCalibrationFunctionForAnalogFilter( const std::string& fileName, const double freq, std::complex<double>* calElog ) const{

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();

	std::ifstream ifs( fileName.c_str(), std::ios::in );
	if( ifs.fail() ){
		ptrOutputFiles->writeErrorMessage( "File open error: " + fileName );
	}

	std::string line;
	std::getline(ifs, line);
	std::vector<std::string> dataLines;
	while(std::getline(ifs, line)){
		dataLines.push_back(line);
	}
	ifs.close();

	const int numFreqsInCalFile = static_cast<int>(dataLines.size());

	double* freqsTemp = new double[numFreqsInCalFile];
	double** amp = new double*[AD_CH];
	double** phs = new double*[AD_CH];
	for( int iCh = 0; iCh < AD_CH; ++iCh ){
		amp[iCh] = new double[numFreqsInCalFile];
		phs[iCh] = new double[numFreqsInCalFile];
		for( int i = 0; i < numFreqsInCalFile; ++i ){
			amp[iCh][i] = 0.0;
			phs[iCh][i] = 0.0;
		}
	}

	int iFreq(0);
	for( std::vector<std::string>::const_iterator itr = dataLines.begin(); itr != dataLines.end(); ++itr, ++iFreq ){
		std::istringstream iss(*itr);
		std::string sbuf;
		double dbuf(0.0);
		{// Freq(Hz)
			std::getline(iss, sbuf, ',');
			std::istringstream iss(sbuf);
			iss >> dbuf;
			freqsTemp[iFreq] = dbuf;
		}
		for( int iCh = 0; iCh < AD_CH; ++iCh ){
			{// Amp(Vpp)
				std::getline(iss, sbuf, ',');
				std::istringstream iss(sbuf);
				iss >> dbuf;
				amp[iCh][iFreq] = dbuf;
			}
			{// Phase(deg)
				std::getline(iss, sbuf, ',');
				std::istringstream iss(sbuf);
				iss >> dbuf;
				phs[iCh][iFreq] = dbuf * CommonParameters::DEG2RAD;
			}
		}
	}

	// Arrange the data in ascending order
	int* ids = new int[numFreqsInCalFile];
	for( int iFreq = 0; iFreq < numFreqsInCalFile; ++iFreq ){
		ids[iFreq] = iFreq;
	}
	Util::quickSort(numFreqsInCalFile, ids, freqsTemp);
	double* freqsInCalFile = new double[numFreqsInCalFile];
	std::complex<double>** cal = new std::complex<double>*[AD_CH];
	for( int iCh = 0; iCh < AD_CH; ++iCh ){
		cal[iCh] = new std::complex<double>[numFreqsInCalFile];
	}
	for( int iFreq = 0; iFreq < numFreqsInCalFile; ++iFreq ){
		const int index = ids[iFreq];
		if( iFreq >= 1 ){
			const int indexPre = ids[iFreq-1];
			assert( freqsTemp[index] > freqsTemp[indexPre] ); 
		}
		freqsInCalFile[iFreq] = freqsTemp[index];
		for( int iCh = 0; iCh < AD_CH; ++iCh ){
			cal[iCh][iFreq] = amp[iCh][index] * std::complex<double>(cos(phs[iCh][index]), sin(phs[iCh][index]));
		}
	}
	delete [] ids;
	delete [] freqsTemp;
	for( int i = 0; i < AD_CH; ++i ){
		delete [] amp[i];
		delete [] phs[i];
	}
	delete [] amp ;
	delete [] phs;

	if( freq < freqsInCalFile[0] ){
		const double rdDelay = -156.7 * 1.e-6 * freq;
		for( int iCh = 0; iCh < AD_CH; ++iCh ){
			calElog[iCh] = std::abs(cal[iCh][0]) * std::complex<double>(cos(rdDelay),sin(rdDelay));
		}
	}
	else if( freq > freqsInCalFile[numFreqsInCalFile-1] ){
		for( int iCh = 0; iCh < AD_CH; ++iCh ){
			calElog[iCh] = cal[iCh][numFreqsInCalFile-1];
		}
	}
	else{
		for( int iFreq = 0; iFreq < numFreqsInCalFile - 1; ++iFreq ){
			if( freq >= freqsInCalFile[iFreq] && freq <= freqsInCalFile[iFreq+1] ){
				const double factor = (freq - freqsInCalFile[iFreq]) / (freqsInCalFile[iFreq+1] - freqsInCalFile[iFreq]);
				for( int iCh = 0; iCh < AD_CH; ++iCh ){
					calElog[iCh] = cal[iCh][iFreq] + (cal[iCh][iFreq+1] - cal[iCh][iFreq]) * factor;
				}
			}
		}
	}

	delete [] freqsInCalFile;
	for( int iCh = 0; iCh < AD_CH; ++iCh ){
		delete [] cal[iCh];
	}
	delete [] cal;

}

int32_t ElogMT::b3_to_long32( unsigned char *ptr ) const{

	union {
		unsigned char c[4];
		int32_t i;
	} tmp;
	
	tmp.c[0] = *ptr++;
	tmp.c[1] = *ptr++;
	tmp.c[2] = *ptr;

	if(*ptr & 0x80){// 0x80 <-> 1000 0000
		// Sign bit is 1 => negative value
		tmp.c[3] = 0xFF;// 0xFF	<-> 1111 1111
	}else {
		// Sign bit is 0 => positive value
		tmp.c[3] = 0;
	}

#ifdef _DEBUG_WRITE
	printf("tmp.i     : 0x%08X\n", tmp.i);
	for( int i = 0; i < 4; i++){
        	printf("tmp.c[%d]  : 0x%02X\n", i, tmp.c[i]);
	}
	std::cout << tmp.i << std::endl;
#endif
	return tmp.i;

}
