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
#include "ElogDual.h"
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

// Return the instance of the class
ElogDual* ElogDual::getInstance() {
	static ElogDual instance;// The only instance
	return &instance;
}

// Read ELOG-Dual binary file
void ElogDual::readElogBinaryFile(const std::string& fileName, const int numSkipData, const int numDataPoints, int& counter, double* ex, double* ey) const {

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeLogMessage("Read Ex and Ey data from " + fileName);
	FILE	* fp = fopen(fileName.c_str(), "rb");
	if(fp == NULL) {
		ptrOutputFiles->writeErrorMessage("File open error : " + fileName);
	}

	// Determine file type
	const Control* const ptrControl = Control::getInstance();
	const int freq = static_cast<int>(ptrControl->getSamplingFrequency());
	const int data_size = AD_CH * freq * AD_BYTES;
	unsigned char* data = new unsigned char[data_size];
	const double factor = 2500.0 / static_cast<double>(pow(2,23));// Dynamic range +/- 2500 mV
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
			const int lbufx = bytesToInt32(&data[j * AD_BYTES * AD_CH]);
			ex[counter - numSkipData] = static_cast<double>(lbufx) * factor;
			// Ey
			const int lbufy = bytesToInt32(&data[j * AD_BYTES * AD_CH + AD_BYTES]);
			ey[counter - numSkipData] = static_cast<double>(lbufy) * factor;
#ifdef _DEBUG_WRITE
			std::cout << std::setw(10) << counter - numSkipData << std::endl;
			std::cout << std::setw(10) << lbufx << std::endl;
			std::cout << std::setw(10) << lbufy << std::endl;
			std::cout << std::setw(10) << ex[counter - numSkipData] << std::endl;
			std::cout << std::setw(10) << ey[counter - numSkipData] << std::endl;
#endif
		}
		if( counter - numSkipData >= numDataPoints ){
			break;
		}
	}

	delete [] data;

	fclose(fp);

}

// Read ELOG-Dual binary file under a directory
void ElogDual::readElogBinaryFilesUnderADirectory(const std::string& directory, const int numSkipData, const int numDataPoints, double* ex, double* ey) const {

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
#ifdef _USE_FILESYSTEM
	const Control* const ptrControl = Control::getInstance();
	const int freq = static_cast<int>(ptrControl->getSamplingFrequency());
	std::ostringstream oss;
	oss << "_" << freq << "Hz";
	const std::string stringCompared = oss.str();
	if (std::filesystem::is_directory(directory)) {
		auto dirItr = std::filesystem::directory_iterator(directory);
		int counter(0);
		for (auto& p : dirItr) {
			const std::string fileName = p.path().string();
			if (fileName.find(stringCompared) != std::string::npos) {
				readElogBinaryFile(fileName, numSkipData, numDataPoints, counter, ex, ey);
				if (counter >= numSkipData + numDataPoints){
					break;
				}
			}
		}
	}
	else
	{
		ptrOutputFiles->writeErrorMessage("There is no directory of ELOG data.");
	}
#else
	ptrOutputFiles->writeErrorMessage("Multiple inputs of ELOG data under a directory is not supported by this version.");
#endif
}

// Get calibration file name
std::string ElogDual::getCalibrationFileName(const int channelIndex) const {

	std::ostringstream fileName;
	fileName << "channel" << channelIndex << ".cal";

	return fileName.str();

}

// Make calibration file
void ElogDual::makeCalibrationFile( const std::string& fileName, const double unitGroupDelay, const int channelIndexX, const int channelIndexY,
	const double dipoleLengthX, const double dipoleLengthY, const std::vector<double>& freq) const {

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();

	const int numFreq = static_cast<int>(freq.size());
	if( numFreq < 1 ){
		ptrOutputFiles->writeErrorMessage( "Number of the frequencies for which calibrations are estimated is less than 1 : " + Util::toString(numFreq) );
	}

	const std::string outputFileNameX = getCalibrationFileName(channelIndexX);
	std::ofstream ofsX;
	ofsX.open( outputFileNameX.c_str(), std::ios::out );
	if( ofsX.fail() ){
		ptrOutputFiles->writeErrorMessage( "File open error : " + outputFileNameX );
	}
	const double dipoleLengthXkm = dipoleLengthX / -1000.0;// [m] -> [km] and invert sign
	ofsX << std::setw(20) << std::scientific << std::setprecision(9) << 1.0/dipoleLengthXkm << std::endl;
	ofsX << std::setw(10) << numFreq << std::endl;

	const std::string outputFileNameY = getCalibrationFileName(channelIndexY);
	std::ofstream ofsY;
	ofsY.open( outputFileNameY.c_str(), std::ios::out );
	if( ofsY.fail() ){
		ptrOutputFiles->writeErrorMessage( "File open error : " + outputFileNameY );
	}
	const double dipoleLengthYkm = dipoleLengthY / -1000.0;// [m] -> [km] and invert sign
	ofsY << std::setw(20) << std::scientific << std::setprecision(9) << 1.0/dipoleLengthYkm << std::endl;
	ofsY << std::setw(10) << numFreq << std::endl;

	const Control* const ptrControl = Control::getInstance();
	for( int iFreq = 0; iFreq < numFreq; ++iFreq ){
		std::complex<double> calElogX(1.0,0.0);
		std::complex<double> calElogY(1.0,0.0);
		calculateCalibrationFunctionForAnalogFilter(fileName, freq[iFreq], calElogX, calElogY);
#ifdef _DEBUG_WRITE
		std::cout << "elog analog filter" << std::endl;
		std::cout << std::setw(20) << std::setprecision(12) << 1.0/freq[iFreq];
		std::cout << std::setw(20) << std::setprecision(12) << std::abs(calElogX);
		std::cout << std::setw(20) << std::setprecision(12) << atan2(calElogX.imag(),calElogX.real()) * CommonParameters::RAD2DEG;
		std::cout << std::endl;
		std::cout << std::setw(20) << std::setprecision(12) << 1.0/freq[iFreq];
		std::cout << std::setw(20) << std::setprecision(12) << std::abs(calElogY);
		std::cout << std::setw(20) << std::setprecision(12) << atan2(calElogY.imag(),calElogY.real()) * CommonParameters::RAD2DEG;
		std::cout << std::endl;
#endif
		{
			//------- 1st fir filer from measured calibration table
			std::string path = "";
			if (!ptrControl->getDirectoryOfLoggerCalibrationFiles().empty()) {
#ifdef _LINUX
				path = ptrControl->getDirectoryOfLoggerCalibrationFiles() + "\/";
#else
				path = ptrControl->getDirectoryOfLoggerCalibrationFiles() + "\\";
#endif
			}
			std::string fileName = "firh.txt";
			if (ptrControl->getTypeOfElogDual() == Control::ELOGDUAL_ADU_MODE) {
				fileName = "firh_adu.txt";
			}
			else if (ptrControl->getTypeOfElogDual() == Control::ELOGDUAL_PHX_MODE) {
				fileName = "firh_phx.txt";
			}
			const std::complex<double> firh = Util::calculateCalibrationForFIRFilterType1(path + fileName, 9, 14336.0, freq[iFreq], true);
			calElogX *= firh;
			calElogY *= firh;
#ifdef _DEBUG_WRITE
			std::cout << "elog firh" << std::endl;
			std::cout << std::setw(20) << std::setprecision(12) << 1.0/freq[iFreq];
			std::cout << std::setw(20) << std::setprecision(12) << std::abs(calElogX);
			std::cout << std::setw(20) << std::setprecision(12) << atan2(calElogX.imag(),calElogX.real()) * CommonParameters::RAD2DEG;
			std::cout << std::endl;
			std::cout << std::setw(20) << std::setprecision(12) << 1.0/freq[iFreq];
			std::cout << std::setw(20) << std::setprecision(12) << std::abs(calElogY);
			std::cout << std::setw(20) << std::setprecision(12) << atan2(calElogY.imag(),calElogY.real()) * CommonParameters::RAD2DEG;
			std::cout << std::endl;
#endif
		}
		{
			//------ 2nd fir filter ---
			// 4900 is described in elog cal file (4.9 Vpp input)
			calElogX /= 4900.0;
			calElogY /= 4900.0;
#ifdef _DEBUG_WRITE
			std::cout << "elog 2nd fir filter" << std::endl;
			std::cout << std::setw(20) << std::setprecision(12) << 1.0/freq[iFreq];
			std::cout << std::setw(20) << std::setprecision(12) << std::abs(calElogX);
			std::cout << std::setw(20) << std::setprecision(12) << atan2(calElogX.imag(),calElogX.real()) * CommonParameters::RAD2DEG;
			std::cout << std::endl;
			std::cout << std::setw(20) << std::setprecision(12) << 1.0/freq[iFreq];
			std::cout << std::setw(20) << std::setprecision(12) << std::abs(calElogY);
			std::cout << std::setw(20) << std::setprecision(12) << atan2(calElogY.imag(),calElogY.real()) * CommonParameters::RAD2DEG;
			std::cout << std::endl;
#endif
		}
		{
			//------ Gain correction ---
			const double tsGroupDelay = 1.0 / 14336.0 * unitGroupDelay;
			const double angle = 2.0 * CommonParameters::PI * freq[iFreq] * tsGroupDelay;
			const std::complex<double> groupDelay = std::complex<double>(cos(angle), sin(angle));
			calElogX /= groupDelay;
			calElogY /= groupDelay;
#ifdef _DEBUG_WRITE
			std::cout << "elog gain correction" << std::endl;
			std::cout << std::setw(20) << std::setprecision(12) << 1.0/freq[iFreq];
			std::cout << std::setw(20) << std::setprecision(12) << std::abs(calElogX);
			std::cout << std::setw(20) << std::setprecision(12) << atan2(calElogX.imag(),calElogX.real()) * CommonParameters::RAD2DEG;
			std::cout << std::endl;
			std::cout << std::setw(20) << std::setprecision(12) << 1.0/freq[iFreq];
			std::cout << std::setw(20) << std::setprecision(12) << std::abs(calElogY);
			std::cout << std::setw(20) << std::setprecision(12) << atan2(calElogY.imag(),calElogY.real()) * CommonParameters::RAD2DEG;
			std::cout << std::endl;
#endif
		}
		if( fabs(ptrControl->getSamplingFrequencyOrg() - 32.0) < CommonParameters::EPS ) {
			//------ Group delay correction ---
			std::string path = "";
			if (!ptrControl->getDirectoryOfLoggerCalibrationFiles().empty()) {
#ifdef _LINUX
				path = ptrControl->getDirectoryOfLoggerCalibrationFiles() + "\/";
#else
				path = ptrControl->getDirectoryOfLoggerCalibrationFiles() + "\\";
#endif
			}
			std::string fileName = "firl.txt";
			if (ptrControl->getTypeOfElogDual() == Control::ELOGDUAL_ADU_MODE) {
				fileName = "firl_adu.txt";
			}
			else if (ptrControl->getTypeOfElogDual() == Control::ELOGDUAL_PHX_MODE) {
				fileName = "firl_phx.txt";
			}
			const std::complex<double> firl = Util::calculateCalibrationForFIRFilterType2(path + fileName, 9, 1024.0, freq[iFreq], -361, 0);
			calElogX *= firl;
			calElogY *= firl;
#ifdef _DEBUG_WRITE
			std::cout << "elog roup delay correction" << std::endl;
			std::cout << std::setw(20) << std::setprecision(12) << 1.0/freq[iFreq];
			std::cout << std::setw(20) << std::setprecision(12) << std::abs(calElogX);
			std::cout << std::setw(20) << std::setprecision(12) << atan2(calElogX.imag(),calElogX.real()) * CommonParameters::RAD2DEG;
			std::cout << std::endl;
			std::cout << std::setw(20) << std::setprecision(12) << 1.0/freq[iFreq];
			std::cout << std::setw(20) << std::setprecision(12) << std::abs(calElogY);
			std::cout << std::setw(20) << std::setprecision(12) << atan2(calElogY.imag(),calElogY.real()) * CommonParameters::RAD2DEG;
			std::cout << std::endl;
#endif
		}
		const std::complex<double> invCalX = 1.0 / calElogX;
		ofsX << std::setw(20) << std::scientific << std::setprecision(9) << freq[iFreq];
		ofsX << std::setw(20) << std::scientific << std::setprecision(9) << invCalX.real();
		ofsX << std::setw(20) << std::scientific << std::setprecision(9) << invCalX.imag() << std::endl;
		const std::complex<double> invCalY = 1.0 / calElogY;
		ofsY << std::setw(20) << std::scientific << std::setprecision(9) << freq[iFreq];
		ofsY << std::setw(20) << std::scientific << std::setprecision(9) << invCalY.real();
		ofsY << std::setw(20) << std::scientific << std::setprecision(9) << invCalY.imag() << std::endl;
	}

	ofsX.close();
	ofsY.close();

}

// Calculate calibration function for analog filter
void ElogDual::calculateCalibrationFunctionForAnalogFilter(const std::string& fileName, const double freq,
	std::complex<double>& calElogX, std::complex<double>& calElogY) const {

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
	//ptrOutputFiles->writeLogMessage("Number of frequencies in the ELOG cal-file: " + Util::toString(numFreqsInCalFile));

	double* freqsTemp = new double[numFreqsInCalFile];
	double* ampX = new double[numFreqsInCalFile];
	double* ampY = new double[numFreqsInCalFile];
	double* phsX = new double[numFreqsInCalFile];
	double* phsY = new double[numFreqsInCalFile];

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
		{// CH1 Amp(Vpp)
			std::getline(iss, sbuf, ',');
			std::istringstream iss(sbuf);
			iss >> dbuf;
			ampX[iFreq] = dbuf;
		}
		{// CH2 Amp(Vpp)
			std::getline(iss, sbuf, ',');
			std::istringstream iss(sbuf);
			iss >> dbuf;
			ampY[iFreq] = dbuf;
		}
		{// CH1 Phase(deg)
			std::getline(iss, sbuf, ',');
			std::istringstream iss(sbuf);
			iss >> dbuf;
			phsX[iFreq] = dbuf * CommonParameters::DEG2RAD;
		}
		{// CH2 Phase(deg)
			std::getline(iss, sbuf, ',');
			std::istringstream iss(sbuf);
			iss >> dbuf;
			phsY[iFreq] = dbuf * CommonParameters::DEG2RAD;
		}
	}

	// Arrange the data in ascending order
	int* ids = new int[numFreqsInCalFile];
	for( int iFreq = 0; iFreq < numFreqsInCalFile; ++iFreq ){
		ids[iFreq] = iFreq;
	}
	Util::quickSort(numFreqsInCalFile, ids, freqsTemp);
	double* freqsInCalFile = new double[numFreqsInCalFile];
	std::complex<double>* calX = new std::complex<double>[numFreqsInCalFile];
	std::complex<double>* calY = new std::complex<double>[numFreqsInCalFile];
	for( int iFreq = 0; iFreq < numFreqsInCalFile; ++iFreq ){
		const int index = ids[iFreq];
		if( iFreq >= 1 ){
			const int indexPre = ids[iFreq-1];
			assert( freqsTemp[index] > freqsTemp[indexPre] ); 
		}
		freqsInCalFile[iFreq] = freqsTemp[index];
		calX[iFreq] = ampX[index] * std::complex<double>(cos(phsX[index]), sin(phsX[index]));
		calY[iFreq] = ampY[index] * std::complex<double>(cos(phsY[index]), sin(phsY[index]));
#ifdef _DEBUG_WRITE
		std::cout << std::setw(20) << std::setprecision(12) << 1.0/freqsInCalFile[iFreq];
		std::cout << std::setw(20) << std::setprecision(12) << std::abs(calX[iFreq]);
		std::cout << std::setw(20) << std::setprecision(12) << atan2(calX[iFreq].imag(),calX[iFreq].real()) * CommonParameters::RAD2DEG;
		std::cout << std::endl;
		std::cout << std::setw(20) << std::setprecision(12) << 1.0/freqsInCalFile[iFreq];
		std::cout << std::setw(20) << std::setprecision(12) << std::abs(calY[iFreq]);
		std::cout << std::setw(20) << std::setprecision(12) << atan2(calY[iFreq].imag(),calY[iFreq].real()) * CommonParameters::RAD2DEG;
		std::cout << std::endl;
#endif
	}
	delete [] ids;
	delete [] freqsTemp;
	delete [] ampX;
	delete [] ampY;
	delete [] phsX;
	delete [] phsY;

	if( freq < freqsInCalFile[0] ){
#ifdef _DEBUG_WRITE
		std::cout << "cho-shuuki" << std::endl;
		std::cout << freq << " " << freqsInCalFile[0] << " " << freqsInCalFile[numFreqsInCalFile-1] << std::endl;
#endif
		const double rdDelay = -156.7 * 1.e-6 * freq;
		calElogX = std::abs(calX[0]) * std::complex<double>(cos(rdDelay),sin(rdDelay));
		calElogY = std::abs(calY[0]) * std::complex<double>(cos(rdDelay),sin(rdDelay));
	}
	else if( freq > freqsInCalFile[numFreqsInCalFile-1] ){
#ifdef _DEBUG_WRITE
		std::cout << "tan-shuuki" << std::endl;
		std::cout << freq << " " << freqsInCalFile[0] << " " << freqsInCalFile[numFreqsInCalFile-1] << std::endl;
#endif
		calElogX = calX[numFreqsInCalFile-1];
		calElogY = calY[numFreqsInCalFile-1];
	}
	else{
#ifdef _DEBUG_WRITE
		std::cout << "aida" << std::endl;
		std::cout << freq << " " << freqsInCalFile[0] << " " << freqsInCalFile[numFreqsInCalFile-1] << std::endl;
#endif
		for( int iFreq = 0; iFreq < numFreqsInCalFile - 1; ++iFreq ){
			if( freq >= freqsInCalFile[iFreq] && freq <= freqsInCalFile[iFreq+1] ){
#ifdef _DEBUG_WRITE
				std::cout << "haitta" << std::endl;
				std::cout << freq << " " << freqsInCalFile[iFreq] << " " << freqsInCalFile[iFreq+1] << std::endl;
				std::cout << calX[iFreq] << " " << calY[iFreq] << std::endl;
				std::cout << calX[iFreq+1] << " " << calY[iFreq+1] << std::endl;
#endif
				const double factor = (freq - freqsInCalFile[iFreq]) / (freqsInCalFile[iFreq+1] - freqsInCalFile[iFreq]);
				calElogX = calX[iFreq] + (calX[iFreq+1] - calX[iFreq]) * factor;
				calElogY = calY[iFreq] + (calY[iFreq+1] - calY[iFreq]) * factor;
			}
		}
	}

	delete [] freqsInCalFile;
	delete [] calX;
	delete [] calY;

}

int ElogDual::bytesToInt32(unsigned char* ptr) const {

	union {
		unsigned char c[4];
		int i;
	} tmp;

	tmp.c[0] = *ptr++;
	tmp.c[1] = *ptr++;
	tmp.c[2] = *ptr;

	if (*ptr & 0x80) {// 0x80 <-> 1000 0000
		// Sign bit is 1 => negative value
		tmp.c[3] = 0xFF;// 0xFF	<-> 1111 1111
	}
	else {
		// Sign bit is 0 => positive value
		tmp.c[3] = 0x00;
	}

	return tmp.i;

}