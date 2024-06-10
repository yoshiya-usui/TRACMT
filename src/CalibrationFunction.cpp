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
#include "CalibrationFunction.h"
#include "Util.h"
#include "OutputFiles.h"
#include "CommonParameters.h"

#include <assert.h>
#include <iostream>

// Constructer
CalibrationFunction::CalibrationFunction( const std::string& calibrationFile ):
	m_factor(1.0),
	m_numFreqs(0),
	m_logFreqs(NULL),
	m_logAmps(NULL),
	m_phases(NULL)
{
	readCalibrationFunction(calibrationFile);
}

// Default constructer
CalibrationFunction::CalibrationFunction():
	m_factor(1.0),
	m_numFreqs(0),
	m_logFreqs(NULL),
	m_logAmps(NULL),
	m_phases(NULL)
{
}

// Destructer
CalibrationFunction::~CalibrationFunction(){

	if( m_logFreqs != NULL ){
		delete [] m_logFreqs;
	};

	if( m_logAmps != NULL ){
		delete [] m_logAmps;
	};

	if( m_phases != NULL ){
		delete [] m_phases;
	};

}

// Read calibration function
void CalibrationFunction::readCalibrationFunction( const std::string& calibrationFile ){

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();

	std::ifstream ifs( calibrationFile.c_str(), std::ios::in );
	if( ifs.fail() ){
		ptrOutputFiles->writeErrorMessage("File open error : " + calibrationFile);
	}
	ptrOutputFiles->writeLogMessage("Read calibration function from " + calibrationFile , false);

	ifs >> m_factor;
	ifs >> m_numFreqs;
#ifdef _DEBUG_WRITE
	std::cout << calibrationFile.c_str() << std::endl;
	std::cout << "m_factor, m_numFreqs: " << m_factor << " " << m_numFreqs << std::endl;
#endif

	if( m_numFreqs < 1 ){
		return;
	}

	double* freqs = new double[m_numFreqs];
	double* calReal = new double[m_numFreqs];
	double* calImag = new double[m_numFreqs];
	for( int iFreq = 0; iFreq < m_numFreqs; ++iFreq ){
		double dbuf(0.0);
		ifs >> dbuf;
		freqs[iFreq] = dbuf;
		ifs >> dbuf;
		calReal[iFreq] = dbuf;
		ifs >> dbuf;
		calImag[iFreq] = dbuf;
		if( ifs.eof() ){
			ptrOutputFiles->writeErrorMessage("Number of frequencies in " + calibrationFile + " is less than " + Util::toString(m_numFreqs) );
			break;
		}
	}
	ifs.close();

	// Arrange the data in ascending order
	int* ids = new int[m_numFreqs];
	for( int iFreq = 0; iFreq < m_numFreqs; ++iFreq ){
		ids[iFreq] = iFreq;
	}
	Util::quickSort(m_numFreqs, ids, freqs);
	if( m_logFreqs != NULL ){
		delete [] m_logFreqs;
	};
	if( m_logAmps != NULL ){
		delete [] m_logAmps;
	};
	if( m_phases != NULL ){
		delete [] m_phases;
	};
	m_logFreqs = new double[m_numFreqs];
	m_logAmps = new double[m_numFreqs];
	m_phases = new double[m_numFreqs];
	for( int iFreq = 0; iFreq < m_numFreqs; ++iFreq ){
		const int index = ids[iFreq];
		if( iFreq >= 1 ){
			const int indexPre = ids[iFreq-1];
			assert( freqs[index] > freqs[indexPre] ); 
		}
		m_logFreqs[iFreq] = log10( freqs[index] );
		const double amplitude = hypot( calReal[index], calImag[index] );
		m_logAmps[iFreq] = log10(amplitude);
		m_phases[iFreq] = atan2( calImag[index], calReal[index] );
#ifdef _DEBUG_WRITE
		std::cout << "m_logFreqs, amplitude, m_logAmps, m_phases: "
			<< m_logFreqs[iFreq] << " " << amplitude << " " << m_logAmps[iFreq] << " " << m_phases[iFreq] << std::endl;
#endif
	}	

	delete [] freqs;
	delete [] calReal;
	delete [] calImag;
	delete [] ids;

}

// Calculate calibration function for a input log frequency
std::complex<double> CalibrationFunction::calculateCalibrationFunction( const double freq ) const{

	if( m_numFreqs < 1 ){
		return std::complex<double>(m_factor, 0.0);
	}

	const double log10Freq = log10(freq);

	const double logAmp =	Util::interpolationAkima( m_numFreqs, m_logFreqs, m_logAmps, log10Freq );
	const double phs = Util::interpolationAkima( m_numFreqs, m_logFreqs, m_phases, log10Freq );
	const double amp = m_factor * pow(10.0, logAmp);
	return amp * std::complex<double>( cos(phs), sin(phs) );

}
