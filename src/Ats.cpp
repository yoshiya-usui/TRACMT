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
#include <vector>
#include <fstream>
#include <ctime>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <bitset>
#include <iomanip>
#include <iostream>
#define snprintf _snprintf

#include "Ats.h"
#include "Util.h"
#include "OutputFiles.h"
#include "CommonParameters.h"
#include "Control.h"

// Default constructer
Ats::Ats():
	m_isCalADUResp(false)
{
}

// Destructer
Ats::~Ats(){
}

// Return the instance of the class
Ats* Ats::getInstance(){
   	static Ats instance;// The only instance
  	return &instance;
}

// Read ats file
void Ats::readAtsFile( const std::string& fileName, const int numSkipData, const int numDataPoints, double* data ) const{

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeLogMessage("Read data from "+fileName);

	FILE* fp = fopen(fileName.c_str(), "rb");
	if (fp == NULL) {
		ptrOutputFiles->writeErrorMessage("File open error : " + fileName);
	}
	unsigned char* buff2bytes = new unsigned char[2];
	unsigned char* buff4bytes = new unsigned char[4];
	unsigned char* buff8bytes = new unsigned char[8];
	fread(buff2bytes, 2, 1, fp);
	const int headerLength = bytesToInt16(buff2bytes);
	fread(buff2bytes, 2, 1, fp);
	const int headerVersion = bytesToInt16(buff2bytes);
	fread(buff4bytes, 4, 1, fp);
	const int numDataInFile = bytesToInt32(buff4bytes);
	fread(buff4bytes, 4, 1, fp);
	const float samplingFreq = bytesToFloat(buff4bytes);
	fread(buff4bytes, 4, 1, fp);
	const int startDate = bytesToInt32(buff4bytes);
	fread(buff8bytes, 8, 1, fp);
	const double lsbmv = bytesToDouble(buff8bytes);
	unsigned char* buffRemaining = new unsigned char[headerLength - 24];
	fread(buffRemaining, headerLength - 24, 1, fp);
	delete[] buffRemaining;

	if (numDataPoints + numSkipData > numDataInFile) {
		ptrOutputFiles->writeErrorMessage("The number of data in the ATS file (" + Util::toString(numDataInFile) +
			") is lower than NSkip(" + Util::toString(numSkipData) + ") + NData(" + Util::toString(numDataPoints) + ")");
	};
	for (int i = 0; i < numSkipData; ++i) {
		fread(buff4bytes, 4, 1, fp);
	}
	int icount(0);
	for (int i = numSkipData; i < numSkipData + numDataPoints; ++i, ++icount) {
		fread(buff4bytes, 4, 1, fp);
		int32_t ibuf = bytesToInt32(buff4bytes);
		data[icount] = static_cast<double>(ibuf) * lsbmv;
	}
	fclose(fp);

	delete[] buff2bytes;
	delete[] buff4bytes;
	delete[] buff8bytes;

}

// Get calibration file name
std::string Ats::getCalibrationFileName( const int channelIndex ) const{

	std::ostringstream fileName;
	fileName << "channel" << channelIndex << ".cal";

	return fileName.str();

}

// Make calibration file
void Ats::makeCalibrationFile( const std::string& inputString, const int channelIndex, const std::vector<double>& freq ) const{

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();

	const std::string outputFileName = getCalibrationFileName(channelIndex);
	std::ofstream ofs;
	ofs.open( outputFileName.c_str(), std::ios::out );
	if( ofs.fail() ){
		ptrOutputFiles->writeErrorMessage( "File open error: " + outputFileName );
	}

	const int numFreq = static_cast<int>(freq.size());
	if( numFreq < 1 ){
		ptrOutputFiles->writeErrorMessage( "Number of the frequencies for which calibrations are estimated is less than 1 : " + Util::toString(numFreq) );
	}
	if( inputString.find("MFS06") != std::string::npos ){
		std::complex<double>* calibFuncCoil = new std::complex<double>[numFreq];
		calculateCalibrationFunctionForCoil( inputString, MFS06, freq, calibFuncCoil );
		ofs << std::setw(20) << std::scientific << std::setprecision(9) << 1.0 << std::endl;
		ofs << std::setw(10) << numFreq << std::endl;
		for( int iFreq = 0; iFreq < numFreq; ++iFreq ){
			std::complex<double> invCal = 1.0 / calibFuncCoil[iFreq];
			if( isCalADUResp() ){
				invCal /= calculateCalibrationFunctionForADULFhannel(freq[iFreq]);
			}
			ofs << std::setw(20) << std::scientific << std::setprecision(9) << freq[iFreq];
			ofs << std::setw(20) << std::scientific << std::setprecision(9) << invCal.real();
			ofs << std::setw(20) << std::scientific << std::setprecision(9) << invCal.imag() << std::endl;
		}
		delete [] calibFuncCoil;
	}else if( inputString.find("MFS07") != std::string::npos ){
		std::complex<double>* calibFuncCoil = new std::complex<double>[numFreq];
		calculateCalibrationFunctionForCoil( inputString, MFS07, freq, calibFuncCoil );
		ofs << std::setw(20) << std::scientific << std::setprecision(9) << 1.0 << std::endl;
		ofs << std::setw(10) << numFreq << std::endl;
		for( int iFreq = 0; iFreq < numFreq; ++iFreq ){
			std::complex<double> invCal = 1.0 / calibFuncCoil[iFreq];
			if( isCalADUResp() ){
				invCal /= calculateCalibrationFunctionForADULFhannel(freq[iFreq]);
			}
			ofs << std::setw(20) << std::scientific << std::setprecision(9) << freq[iFreq];
			ofs << std::setw(20) << std::scientific << std::setprecision(9) << invCal.real();
			ofs << std::setw(20) << std::scientific << std::setprecision(9) << invCal.imag() << std::endl;
		}
		delete [] calibFuncCoil;
	}else{
		// Channel for the electric field [mV/km]
		// The input should be dipole length [m]
		std::istringstream iss(inputString);
		double dipoleLength(0.0);
		iss >> dipoleLength;
		dipoleLength /= -1000.0;// [m] -> [km] and invert sign
		ofs << std::setw(20) << std::scientific << std::setprecision(9) << 1.0/dipoleLength << std::endl;
		if( m_isCalADUResp ){
			ofs << std::setw(10) << numFreq << std::endl;
			for( int iFreq = 0; iFreq < numFreq; ++iFreq ){
				const std::complex<double> invCal = 1.0 / calculateCalibrationFunctionForADULFhannel(freq[iFreq]);
				ofs << std::setw(20) << std::scientific << std::setprecision(9) << freq[iFreq];
				ofs << std::setw(20) << std::scientific << std::setprecision(9) << invCal.real();
				ofs << std::setw(20) << std::scientific << std::setprecision(9) << invCal.imag() << std::endl;
			}
		}else{
			ofs << std::setw(10) << 0 << std::endl;
		}
	}

	ofs.close();

}

// Get flag specifing whether the calibration function for ADU is calculated
void Ats::makeCalibrationFileOnlyWithDipoleLength(const std::string& inputString, const int channelIndex) const {

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();

	const std::string outputFileName = getCalibrationFileName(channelIndex);
	std::ofstream ofs;
	ofs.open(outputFileName.c_str(), std::ios::out);
	if (ofs.fail()) {
		ptrOutputFiles->writeErrorMessage("File open error: " + outputFileName);
	}

	// Channel for the electric field [mV/km]
	// The input should be dipole length [m]
	std::istringstream iss(inputString);
	double dipoleLength(0.0);
	iss >> dipoleLength;
	dipoleLength /= -1000.0;// [m] -> [km] and invert sign
	ofs << std::setw(20) << std::scientific << std::setprecision(9) << 1.0 / dipoleLength << std::endl;
	ofs << std::setw(10) << 0 << std::endl;

	ofs.close();

}

// Get flag specifing whether the calibration function for ADU is calculated
bool Ats::isCalADUResp() const{
	return m_isCalADUResp;
}

// Set flag specifing whether the calibration function for ADU is calculated
void Ats::setIsCalADUResp( const bool isCalADUResp ){
	m_isCalADUResp = isCalADUResp;
}

// Calculate calibration function for coil
void Ats::calculateCalibrationFunctionForCoil( const std::string& calibrationFile, const int coilType,
	const std::vector<double>& freq, std::complex<double>* calibrationFunction ) const{

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();

	std::ifstream ifs( calibrationFile.c_str(), std::ios::in );
	if( ifs.fail() ){
		ptrOutputFiles->writeErrorMessage("File open error : " + calibrationFile);
	}

	ptrOutputFiles->writeLogMessage("Read calibration function from " + calibrationFile , false);

#if 0
	int numFreqInCoilCalFile(0);
	ifs >> numFreqInCoilCalFile;
	double* freqInCoilCalFile = new double[numFreqInCoilCalFile];
	double* logFreqInCoilCalFile = new double[numFreqInCoilCalFile];
	double* amp = new double[numFreqInCoilCalFile];
	double* phs = new double[numFreqInCoilCalFile];
	for( int iFreq = 0; iFreq < numFreqInCoilCalFile; ++iFreq ){
		double dbuf(0.0);
		ifs >> dbuf;
		freqInCoilCalFile[iFreq] = dbuf;
		logFreqInCoilCalFile[iFreq] = log10(dbuf);
		ifs >> dbuf;
		// [V/nT/Hz] => [mV/nT]
		amp[iFreq] = dbuf * freqInCoilCalFile[iFreq] * 1000.0;
		ifs >> dbuf;
		// [deg] => [rad]
		phs[iFreq] = dbuf * CommonParameters::DEG2RAD;
		if( ifs.eof() ){
			ptrOutputFiles->writeErrorMessage("Number of frequencies in " + calibrationFile + " is less than " + Util::toString(numFreqInCoilCalFile) );
			break;
		}
	}
#else
	std::string sbuf;
	std::vector<std::string> lines;
	while(getline(ifs,sbuf)){
		if( sbuf.substr(0,1) == "+" ){
			lines.push_back(sbuf);
		}else if( !lines.empty() ){
			break;
		}
		if( ifs.eof() ){
			break;
		}
	}
	const int numFreqInCoilCalFile = static_cast<int>(lines.size());
	double* freqsTemp = new double[numFreqInCoilCalFile];
	double* amp = new double[numFreqInCoilCalFile];
	double* phs = new double[numFreqInCoilCalFile];
	for( int iFreq = 0; iFreq < numFreqInCoilCalFile; ++iFreq ){
		std::istringstream iss(lines[iFreq]);
		double dbuf(0.0);
		iss >> dbuf;
		freqsTemp[iFreq] = dbuf;
		iss >> dbuf;
		// [V/nT/Hz] => [mV/nT]
		amp[iFreq] = dbuf * freqsTemp[iFreq] * 1000.0;
		iss >> dbuf;
		// [deg] => [rad]
		phs[iFreq] = dbuf * CommonParameters::DEG2RAD;
	}
#endif
	ifs.close();

	// Arrange the data in ascending order
	int* ids = new int[numFreqInCoilCalFile];
	for( int iFreq = 0; iFreq < numFreqInCoilCalFile; ++iFreq ){
		ids[iFreq] = iFreq;
	}
	Util::quickSort(numFreqInCoilCalFile, ids, freqsTemp);
	double* freqInCoilCalFile = new double[numFreqInCoilCalFile];
	double* logFreqInCoilCalFile = new double[numFreqInCoilCalFile];
	double* coilCalAmp = new double[numFreqInCoilCalFile];
	double* coilCalLogAmp = new double[numFreqInCoilCalFile];
	double* coilCalPhs = new double[numFreqInCoilCalFile];
	for( int iFreq = 0; iFreq < numFreqInCoilCalFile; ++iFreq ){
		const int index = ids[iFreq];
		if( iFreq >= 1 ){
			const int indexPre = ids[iFreq-1];
			assert( freqsTemp[index] > freqsTemp[indexPre] ); 
		}
		freqInCoilCalFile[iFreq] = freqsTemp[index];
		logFreqInCoilCalFile[iFreq] = log10(freqsTemp[index]);
		coilCalAmp[iFreq] = amp[index];
		coilCalLogAmp[iFreq] = log10(amp[index]);
		coilCalPhs[iFreq] = phs[index];
	}	
	delete [] ids;
	delete [] freqsTemp;
	delete [] amp;
	delete [] phs;

	double cutoffFrequency(0.0);
	double fc1(0.0);
	double fc2(0.0);
	double AFactor(0.0);
	switch (coilType){
	case Ats::MFS06:
		cutoffFrequency = 0.5;
		fc1 = 4.0;
		fc2 = 8192.0;
		AFactor = 0.8;
		break;
	case Ats::MFS07:
		cutoffFrequency = 5.0;
		fc1 = 32.0;
		fc2 = 40000.0;
		AFactor = 0.64;
		break;
	default:
		ptrOutputFiles->writeErrorMessage("Unknown coil type : " + Util::toString(coilType) );
		break;
	}

	if( freqInCoilCalFile[0] < cutoffFrequency - CommonParameters::EPS ){
		int numFreqLower(0);
		double fc1Avg = 0.0;
		for( int iFreq = 0; iFreq < numFreqInCoilCalFile; ++iFreq ){
			const double freq = freqInCoilCalFile[iFreq];
			if( freq > cutoffFrequency ){
				continue;
			}
			const double phaseP2 = atan2( freq, fc2 );
			fc1Avg += freq * tan( coilCalPhs[iFreq] + phaseP2 );
			++numFreqLower;
		}
		fc1 = fc1Avg / static_cast<double>(numFreqLower);

		double AAvg = 0.0;
		for( int iFreq = 0; iFreq < numFreqInCoilCalFile; ++iFreq ){
			const double freq = freqInCoilCalFile[iFreq];
			if( freq > cutoffFrequency ){
				continue;
			}
			const double term1 = pow(fc1/freq, 2) + 1.0;// (f1/f)^2 + 1
			const double term2 = pow(freq/fc2, 2) + 1.0;// 1 + (f/f2)^2 
			AAvg += coilCalAmp[iFreq] * sqrt(term1 * term2);
		}
		AFactor = AAvg / static_cast<double>(numFreqLower);
	}

	int iFreq(0);
	for( std::vector<double>::const_iterator itrFreq = freq.begin(); itrFreq != freq.end(); ++itrFreq, ++iFreq ){
		const double freq = *itrFreq;
		if( freq < freqInCoilCalFile[0] ){
			// Calculate by estimated transfer function
			const double ffc1 = freq / fc1;
			const double ffc2 = freq / fc2;
			calibrationFunction[iFreq] =
				AFactor * std::complex<double>(0.0, ffc1) / std::complex<double>(1.0, ffc1) / std::complex<double>(1.0, ffc2);
		}else{
			// Interpolation in logarithmic axes
			const double logAmp =	Util::interpolationAkima( numFreqInCoilCalFile, logFreqInCoilCalFile, coilCalLogAmp, log10(freq) );
			const double phs = Util::interpolationAkima( numFreqInCoilCalFile, logFreqInCoilCalFile, coilCalPhs, log10(freq) );
			const double amp = pow(10.0, logAmp);
			calibrationFunction[iFreq] = amp * std::complex<double>( cos(phs), sin(phs) );
		}
	}

	delete [] freqInCoilCalFile;
	delete [] logFreqInCoilCalFile;
	delete [] coilCalAmp;
	delete [] coilCalLogAmp;
	delete [] coilCalPhs;

#ifdef _DEBUG_WRITE
	iFreq = 0;
	for( std::vector<double>::const_iterator itrFreq = freq.begin(); itrFreq != freq.end(); ++itrFreq, ++iFreq ){
		const double freq = *itrFreq;
		std::cout << std::setw(20) << std::setprecision(12) << 1.0/freq;
		std::cout << std::setw(20) << calibrationFunction[iFreq] << std::endl;
	}
#endif

}

// Calculate calibration function for LF-channel of ADU 
std::complex<double> Ats::calculateCalibrationFunctionForADULFhannel( const double freq ) const{

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeErrorMessage(static_cast<std::string>(__FUNCTION__) + "cannot be used in this version");
	return std::complex<double>(0.0, 0.0);

}

short int Ats::bytesToInt16(unsigned char* ptr) const {

	union {
		unsigned char c[2];
		short int i;
	} tmp;

	tmp.c[0] = *ptr++;
	tmp.c[1] = *ptr;

	return tmp.i;

}

int Ats::bytesToInt32(unsigned char* ptr) const {

	union {
		unsigned char c[4];
		int i;
	} tmp;

	tmp.c[0] = *ptr++;
	tmp.c[1] = *ptr++;
	tmp.c[2] = *ptr++;
	tmp.c[3] = *ptr;

	return tmp.i;

}

float Ats::bytesToFloat(unsigned char* ptr) const {

	union {
		unsigned char c[4];
		float f;
	} tmp;

	tmp.c[0] = *ptr++;
	tmp.c[1] = *ptr++;
	tmp.c[2] = *ptr++;
	tmp.c[3] = *ptr;

	return tmp.f;

}

double Ats::bytesToDouble(unsigned char* ptr) const {

	union {
		unsigned char c[8];
		double d;
	} tmp;

	tmp.c[0] = *ptr++;
	tmp.c[1] = *ptr++;
	tmp.c[2] = *ptr++;
	tmp.c[3] = *ptr++;
	tmp.c[4] = *ptr++;
	tmp.c[5] = *ptr++;
	tmp.c[6] = *ptr++;
	tmp.c[7] = *ptr;

	return tmp.d;

}