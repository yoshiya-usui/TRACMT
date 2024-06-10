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
#include "Analysis.h"
#include "Control.h"
#include "OutputFiles.h"
#include "Util.h"
#include "UtilRobust.h"
#include "TableOfTukeysBiweightParameters.h"
#include "Ats.h"
#include "ElogDual.h"
#include "ElogMT.h"
#include "RobustWeightTukeysBiweights.h"
#include "RobustWeightHuber.h"
#include "RobustWeightThomson.h"
#include "RobustPrewhitening.h"

#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#ifdef _MERSENNE_TWISTER_ORIGINAL
#else
#include <random>
#endif

#ifdef _MERSENNE_TWISTER_ORIGINAL
#include "mt64.h"
#endif

#ifdef _USE_OMP
#include <omp.h>
#endif

// Default constructer
Analysis::Analysis():
	m_calibrationFunctions(NULL),
	m_coefficientsOfARModel(NULL)
{
	for( int i = 0; i < 2; ++i ){
		m_robustWeight[i] = NULL;
	}
}

// Destructer
Analysis::~Analysis()
{
	if( m_calibrationFunctions != NULL ){
		delete [] m_calibrationFunctions;
	}

	if( m_coefficientsOfARModel != NULL ){
		delete [] m_coefficientsOfARModel;
	}

	for( int i = 0; i < 2; ++i ){
		if( m_robustWeight[i] != NULL ){
			delete m_robustWeight[i];
		}
	}
}

// Run analysis
void Analysis::run( std::vector<CommonParameters::DataFileSet>& dataFileSets ){
	
	const Control* const ptrControl = Control::getInstance();
	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();

	std::vector<double> freqAll = ptrControl->getFrequenciesAllWithoutDuplications();
	readTimeSeriesData(dataFileSets);

	// Merge sections
	mergeSections(dataFileSets);

	// Apply decimation
	decimation(dataFileSets);

	// Read calibration files for main analysis
	readCalibrationFiles(freqAll);

	// Prior evaluation before preprocessing
	if( ptrControl->doesOutputTimeSeriesToCsv() ){
		outputTimeSeriesData( dataFileSets, false );
	}

	if( ptrControl->getParamsForTimeDomainEvaluation().doEvaluation ){
		const Control::ParamsForTimeDomainEvaluation params = ptrControl->getParamsForTimeDomainEvaluation();
		priorEvaluationOfTimeSeriesData( params.timeSeriesInterval, dataFileSets, false );
	}
	if( ptrControl->getParamsForFreqDomainEvaluation().doEvaluation ){
		priorEvaluationOfFrequencyDataFromAllData( dataFileSets, false );
	}

	// Pre-processings
	preprocessing(dataFileSets);

	// Prior evaluation after preprocessing
	if( ptrControl->getParamsForTimeDomainEvaluation().doEvaluation ){
		const Control::ParamsForTimeDomainEvaluation params = ptrControl->getParamsForTimeDomainEvaluation();
		priorEvaluationOfTimeSeriesData( params.timeSeriesInterval, dataFileSets, true );
	}
	if( ptrControl->getParamsForFreqDomainEvaluation().doEvaluation ){
		priorEvaluationOfFrequencyDataFromAllData( dataFileSets, true );
	}

	std::ofstream ofsResp;
	const std::string fileNameResp = "response_functions.csv";
	ofsResp.open( fileNameResp.c_str(), std::ios::out );
	if( ofsResp.fail() ){
		ptrOutputFiles->writeErrorMessage( "File open error : " + fileNameResp );
	}
	writeHeaderToOutputFileForResponseFunctions(ofsResp);

	std::ofstream ofsRhoaPhs;
	if( ptrControl->doesOutputApparentResistivityAndPhase() ){
		const std::string fileNameRhoaPhs = "apparent_resistivity_and_phase.csv";
		ofsRhoaPhs.open( fileNameRhoaPhs.c_str(), std::ios::out );
		if( ofsRhoaPhs.fail() ){
			ptrOutputFiles->writeErrorMessage( "File open error : " + fileNameRhoaPhs );
		}
		writeHeaderToOutputFileForApparentResistivityAndPhase(ofsRhoaPhs);
	}

	const int numChannels = ptrControl->getNumberOfChannels();
	const double samplingFreq = ptrControl->getSamplingFrequency();
	const int numSegmentLengths = ptrControl->getNumSegmentLengths();
	for( int iSegLen = 0; iSegLen < numSegmentLengths; ++iSegLen ){
		const int segmentLength = ptrControl->getSegmentLength(iSegLen);
		ptrOutputFiles->writeLogMessage("===============================================================================",false);
		ptrOutputFiles->writeLogMessage("Perform analysis for segment length : " + Util::toString(segmentLength));

		// Fourier transform of data
		std::complex<double>*** cdata = new std::complex<double>**[numChannels];
		int numSegmentsTotal(0);
		std::vector< std::pair<std::string, std::string> > timesOrg;
		std::vector<double>* meanSquares = new std::vector<double>[numChannels];
		convertToFrequencyData(segmentLength, dataFileSets, numSegmentsTotal, cdata, timesOrg, meanSquares);

		// Select segments to be excluded
		std::vector<bool> remainingSegments;
		remainingSegments.reserve(numSegmentsTotal);
		for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
			// Initially all segments are to be used
			remainingSegments.push_back(true);
		}

		std::ostringstream fileName;
		selectSegmentsToBeExcludedByMeanSquareCriteria(numSegmentsTotal, meanSquares, timesOrg, remainingSegments);

		// Loop of frequency at which response function is estimated
		const int numTargetFreqsInSegment = ptrControl->getNumTargetFrequencyInSegment(iSegLen);
		for( int iFreq = 0; iFreq < numTargetFreqsInSegment; ++iFreq ){
			ptrOutputFiles->writeLogMessage("-------------------------------------------------------------------------------",false);
			const int freqDegree = ptrControl->getTargetFrequencyDegreesInSegment(iSegLen, iFreq);
			const double freq = samplingFreq * static_cast<double>(freqDegree) / static_cast<double>(segmentLength); 
			ptrOutputFiles->writeLogMessage("Now Frequency(Hz): " + Util::toString(freq) + ", Period(s): " + Util::toString(1.0/freq));
			const double timeLength = static_cast<double>(segmentLength) / samplingFreq;

			// Copy the array because it will be overwrited in the selection using square coherence criteria
			std::vector<bool> remainingSegmentsForThisFrequency = remainingSegments;

			// Copy Fourier transformed values
			std::complex<double>** ftvalOrg = new std::complex<double>*[numChannels];
			for( int iChan = 0; iChan < numChannels; ++iChan ){
				ftvalOrg[iChan] = new std::complex<double>[numSegmentsTotal];
				for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
					assert(freqDegree > 0);
					ftvalOrg[iChan][iSeg] = cdata[iChan][iSeg][freqDegree];
				}
			}

			calibrationCorrectionAllChannels(numChannels, numSegmentsTotal, freq, ftvalOrg );
			calculateRotatedFields(numSegmentsTotal, ftvalOrg);

			// Prior evaluation of data segments
			if( ptrControl->getParamsForDataSegmentEvaluation().doEvaluation ){
				std::ostringstream fileName;
				fileName << "segm" << iSegLen << "_index" << freqDegree << ".csv";
				priorEvaluationOfDataSegments( numSegmentsTotal, ftvalOrg, timesOrg, timeLength, fileName.str() );
			}

			// Select segments to be excluded by square coherence criteria
			if( ptrControl->getParamsForSquareCoherenceCriteria().applySquareCoherenceCriteria ){
				selectSegmentsToBeExcludedBySquareCoherenceCriteria( numSegmentsTotal, ftvalOrg, remainingSegmentsForThisFrequency );
			}

			// Select segments to be excluded by square coherence criteria with random sampling
			if( ptrControl->getParamsForSquareCoherenceCriteriaWithRandomSampling().applySquareCoherenceCriteria ){
				selectSegmentsToBeExcludedBySquareCoherenceCriteriaWithRandomSampling( numSegmentsTotal, ftvalOrg, remainingSegmentsForThisFrequency );
			}

			// Select segments to be excluded by magnetic polarization diretion criteria
			if( ptrControl->getParamsForMagneticPolarizatitonDirectionCriteria().applyMagneticPolarizatitonDirectionCriteria ){
				selectSegmentsToBeExcludedByMagneticPolarizatitonDirectionCriteria( numSegmentsTotal, ftvalOrg, remainingSegmentsForThisFrequency );
			}

			// Select segments to be excluded by degree of magnetic polarization criteria
			if( ptrControl->getParamsForDegreeOfMagneticPolarizationCriteria().applyDegreeOfMagneticPolarizationCriteria ){
				selectSegmentsToBeExcludedByDegreeOfMagneticPolarizationCriteria( numSegmentsTotal, ftvalOrg, remainingSegmentsForThisFrequency );
			}

			int numOfRemainingSegments(0);
			for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
				if(remainingSegmentsForThisFrequency[iSeg]){
					++numOfRemainingSegments;
				}
			}
			if( numOfRemainingSegments <= 0 ){
				ptrOutputFiles->writeErrorMessage("Number of remaining segments is too small: " + Util::toString(numOfRemainingSegments));
			}

#ifdef _DEBUG_WRITE
			for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
				for( int iChan = 0; iChan < ptrControl->getNumberOfChannels(); ++iChan ){
					std::cout << ftvalOrg[iChan][iSeg] << " ";
				}
				std::cout << std::endl;
			}
#endif
			// Make data excluding some segments
			std::complex<double>** ftval = new std::complex<double>*[numChannels];
			for( int iChan = 0; iChan < numChannels; ++iChan ){
				ftval[iChan] = new std::complex<double>[numOfRemainingSegments];
				int icount(0);
				for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
					if(remainingSegmentsForThisFrequency[iSeg]){
						ftval[iChan][icount] = ftvalOrg[iChan][iSeg];
						++icount;
					}
				}
				assert( icount == numOfRemainingSegments );
			}
#ifdef _DEBUG_WRITE
			for( int iSeg = 0; iSeg < numOfRemainingSegments; ++iSeg ){
				for( int iChan = 0; iChan < ptrControl->getNumberOfChannels(); ++iChan ){
					std::cout << ftval[iChan][iSeg] << " ";
				}
				std::cout << std::endl;
			}
			for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
				std::cout << timesOrg[iSeg].first << " " << timesOrg[iSeg].second << " " << std::endl;
			}
#endif

			std::vector< std::pair<std::string, std::string> > times;
			times.reserve(numOfRemainingSegments);
			for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
				if(remainingSegmentsForThisFrequency[iSeg]){
					times.push_back(timesOrg[iSeg]);
				}
			}
#ifdef _DEBUG_WRITE
			for( int iSeg = 0; iSeg < numOfRemainingSegments; ++iSeg ){
				std::cout << times[iSeg].first << " " << times[iSeg].second << " " << std::endl;
			}
#endif
			assert( times.size() == numOfRemainingSegments );

			if( ptrControl->doesOutputFreqDomainDataToCsv() ){
				outputFrequencyDomainData( iSegLen, freqDegree, numOfRemainingSegments, ftval );
			}
			// Calculate response functions
			calculateResponseFunctions( iSegLen, freqDegree, timeLength, freq, numOfRemainingSegments, ftval, times, ofsResp, ofsRhoaPhs );

			// Delete arrays
			for( int iChan = 0; iChan < numChannels; ++iChan ){
				delete [] ftvalOrg[iChan];
			}
			delete [] ftvalOrg;
			for( int iChan = 0; iChan < numChannels; ++iChan ){
				delete [] ftval[iChan];
			}
			delete [] ftval;
		}

		// Delete arrays
		for( int iChan = 0; iChan < numChannels; ++iChan ){
			for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
				delete [] cdata[iChan][iSeg];
			}
			delete [] cdata[iChan];
		}
		delete [] cdata;
		delete [] meanSquares;
	}

}

// Set object of M-estimator with Huber weight
void Analysis::setRobustWeightHuber( const int index, const RobustWeightHuber* const ptrRobustWeight ){

	if( m_robustWeight[index] != NULL ){
		delete m_robustWeight[index];
	}

	m_robustWeight[index] = new RobustWeightHuber(*ptrRobustWeight);

}

// Set object of M-estimator with Thomson weight
void Analysis::setRobustWeightThomson( const int index, const RobustWeightThomson* const ptrRobustWeight ){

	if( m_robustWeight[index] != NULL ){
		delete m_robustWeight[index];
	}

	m_robustWeight[index] = new RobustWeightThomson(*ptrRobustWeight);

}

// Set object of M-estimator with Tukey's biweights
void Analysis::setRobustWeightTukeysBiweights( const int index, const RobustWeightTukeysBiweights* const ptrRobustWeight ){

	if( m_robustWeight[index] != NULL ){
		delete m_robustWeight[index];
	}

	m_robustWeight[index] = new RobustWeightTukeysBiweights(*ptrRobustWeight);

}

// Test function
void Analysis::test() const{

	const std::string fileName = "KrafftPoint.txt";
	std::ifstream ifs(fileName.c_str(), std::ios::in);
	int numOfData(0);
	ifs >> numOfData;
	double* x = new double[numOfData];
	double* y = new double[numOfData];
	for (int iData = 0; iData < numOfData; ++iData) {
		ifs >> y[iData] >> x[iData];
	}
	ifs.close();
	std::vector< std::pair<int, int> > candidatesOfPairsInit;
	const int numOfCandidatesOfPairsAll = (numOfData * numOfData - numOfData) / 2;
	for (int i1 = 0; i1 < numOfData; ++i1) {
		for (int i2 = i1 + 1; i2 < numOfData; ++i2) {
			candidatesOfPairsInit.push_back(std::make_pair(i1, i2));
		}
	}
	double slope(0.0);
	double intercept(0.0);
	double* weights = new double[numOfData];
	for (int iData = 0; iData < numOfData; ++iData) {
		weights[iData] = 1.0;
	}
	const int numOfCandidates = candidatesOfPairsInit.size() + 1;
	double* slopeCandidates = new double[numOfCandidates];
	double* interceptCandidates = new double[numOfCandidates];
	int icount(0);
	// Least square solution
	Util::calculateRegressionCoefficientsByWLSWithIntercept(numOfData, y, x, weights, slope, intercept);
	slopeCandidates[icount] = slope;
	interceptCandidates[icount] = intercept;
	++icount;
	std::cout << "OLS " << slope << " " << intercept << std::endl;
	for (std::vector< std::pair<int, int> >::const_iterator itr = candidatesOfPairsInit.begin(); itr != candidatesOfPairsInit.end(); ++itr, ++icount) {
		const int i1 = itr->first;
		const int i2 = itr->second;
		slope = (y[i2] - y[i1]) / (x[i2] - x[i1]);
		intercept = y[i1] - slope * x[i1];
		slopeCandidates[icount] = slope;
		interceptCandidates[icount] = intercept;
		std::cout << icount << " " << i1 << " " << i2 << " " << slope << " " << intercept << std::endl;
	}
	double scale(0.0);
	UtilRobust::computeSEstimatorForUnivariateLinearRegressionWithIntercept(numOfData, numOfCandidates,
		y, x, slopeCandidates, interceptCandidates, scale, slope, intercept);
	std::cout << "S-estimate " << slope << " " << intercept << " " <<scale << std::endl;
	delete[] x;
	delete[] y;
	delete[] slopeCandidates;
	delete[] interceptCandidates;
	return;

}

// Test function 2
void Analysis::test2( std::vector<CommonParameters::DataFileSet>& dataFileSets ){

	//const Control* const ptrControl = Control::getInstance();

	//for( int i = 1; i < 10000; ++i ){
	//	const double dsec = static_cast<double>(i) * 0.5;
	//	std::cout << std::setw(10) << dsec << " " << ptrControl->getTimeFromStartTimeOfEachSection(0, dsec) << std::endl;
	//}
	//for( int i = 1; i < 10000; ++i ){
	//	const double dsec = static_cast<double>(i) * 0.5;
	//	std::cout << std::setw(10) << dsec << " " << ptrControl->getTimeFromStartTimeOfEachSection(1, dsec) << std::endl;
	//}

	readTimeSeriesData(dataFileSets);

	//const int numChannels = m_numOutputVariables + 2 + 2 * m_numRemoteReferenceStations;
	const int numChannels = 5;
	for( std::vector<CommonParameters::DataFileSet>::const_iterator itr = dataFileSets.begin(); itr != dataFileSets.end(); ++itr ){
		const int numData = itr->numDataPoints;
		for( int iChan = 0; iChan < numChannels; ++iChan ){
			std::cout << "Var = " << iChan << std::endl;
			std::ostringstream fileName;
			fileName << "var_" << iChan <<".csv";
			std::ofstream ofs;
			ofs.open( fileName.str().c_str(), std::ios::out );
			if( ofs.fail() ){
				std::cerr << "File open error !! : " << fileName.str() << std::endl;
				exit(1);
			}
			std::complex<double>* cdata = new std::complex<double>[numData];
			for( int i = 0; i < numData; ++i ){
				cdata[i] = std::complex<double>(itr->dataFile[iChan].data[i], 0.0);
			}
			Util::fourierTransform(numData, cdata);
			const double samplingFrequency = Control::getInstance()->getSamplingFrequency();
			const double dt = 1.0 / samplingFrequency;
			const double T = static_cast<double>(numData) * dt;
			ofs << "i" << "," << "freq"
				<< "," << "amp" << "," << "amp*T/2"
				<< "," << "phs"
				<< "," << "power" << "," << "power*T"
				<< std::endl;
			for( int i = 0; i < numData / 2 + 1; ++i ){
				const double freq = static_cast<double>(i) / static_cast<double>(numData) * samplingFrequency;
				const double A =   2.0 * cdata[i].real();
				const double B = - 2.0 * cdata[i].imag();
				const double amp = hypot(A, B);
				const double phs = atan2(-B, A) * CommonParameters::RAD2DEG;
				double power = 2.0 * std::norm(cdata[i]);
				if( i == 0 || i == 8 ){
					power *= 0.5;
				}
				ofs << i << "," << freq
					<< "," << amp << "," << amp * 0.5 * T
					<< "," << phs
					<< "," << power << "," << power * T
					<< std::endl;
			}
			ofs.close();
			delete [] cdata;
		}
	}

}

// Test function 3
void Analysis::test3() const{

	std::ostringstream fileName;
	fileName << "temp.csv";
	std::ofstream ofs;
	ofs.open( fileName.str().c_str(), std::ios::out );
	if( ofs.fail() ){
		std::cerr << "File open error !! : " << fileName.str() << std::endl;
		exit(1);
	}

	const int numData = 4096;
	double data[4096];
	for( int iTime = 0; iTime < numData; ++iTime ){
		data[iTime] = 0.0;
	}
	const double samplingFreq = 32;
	const double dt = 1.0 / samplingFreq;
	for( int iFreq = 0; iFreq < 1000; ++iFreq ){
		const double freq = static_cast<double>(iFreq+1) / static_cast<double>(1000) * samplingFreq;
		double omega = 2.0 * CommonParameters::PI * freq;
		const double phase = static_cast<double>(iFreq*rand()) / static_cast<double>(RAND_MAX) * CommonParameters::DEG2RAD;
		for( int iTime = 0; iTime < numData; ++iTime ){
			const double t = static_cast<double>(iTime) * dt;
			data[iTime] += 2.0 * cos(omega*t + phase);
		}
	}
	std::complex<double> cdata[numData];
	for( int i = 0; i < numData; ++i ){
		cdata[i] = std::complex<double>(data[i], 0.0);
	}
	for( int i = 0; i < numData; ++i ){
		const double t = static_cast<double>(i) * dt;
		ofs << std::setw(15) << t << "," << std::setw(15) << cdata[i].real() << "," << std::setw(15) << cdata[i].imag() << std::endl;
	}
	ofs << std::endl;

	Util::fourierTransform(numData, cdata);
	ofs << "i" << "," << "freq"
		<< "," << "amp" << "," << "amp*T/2"
		<< "," << "phs"
		<< "," << "power" << "," << "power*T"
		<< std::endl;
	const double T = numData * dt;
	for( int i = 0; i < numData / 2 + 1; ++i ){
		const double freq = static_cast<double>(i) / static_cast<double>(numData) * samplingFreq;
		const double A =   2.0 * cdata[i].real();
		const double B = - 2.0 * cdata[i].imag();
		const double amp = hypot(A, B);
		const double phs = atan2(-B, A) * CommonParameters::RAD2DEG;
		double power = 2.0 * std::norm(cdata[i]);
		if( i == 0 || i == 8 ){
			power *= 0.5;
		}
		ofs << i << "," << freq
			<< "," << amp << "," << amp * 0.5 * T
			<< "," << phs
			<< "," << power << "," << power * T
			<< std::endl;
	};
	ofs << std::endl;

	const int numh = 200;
	double h[201];
	Util::calculateFIRFilterCoeffsByLeastSquare( numh, true, samplingFreq, 0.4, 0.5, 100.0, 1.0, h );

	for( int i = numh; i < numData; ++i ){
		cdata[i] = 0.0;
		for( int k = 0; k <= numh; ++k ){
			cdata[i] += h[k] * std::complex<double>(data[i-k], 0.0);
		}
	}
	for( int i = 0; i < numData; ++i ){
		ofs << i << "," << std::setw(15) << cdata[i].real() << "," << std::setw(15) << cdata[i].imag() << std::endl;
	}
	ofs << std::endl;

	Util::fourierTransform(numData, cdata);
	ofs << "i" << "," << "freq"
		<< "," << "amp" << "," << "amp*T/2"
		<< "," << "phs"
		<< "," << "power" << "," << "power*T"
		<< std::endl;
	for( int i = 0; i < numData / 2 + 1; ++i ){
		const double freq = static_cast<double>(i) / static_cast<double>(numData) * samplingFreq;
		const double A =   2.0 * cdata[i].real();
		const double B = - 2.0 * cdata[i].imag();
		const double amp = hypot(A, B);
		const double phs = atan2(-B, A) * CommonParameters::RAD2DEG;
		double power = 2.0 * std::norm(cdata[i]);
		if( i == 0 || i == 8 ){
			power *= 0.5;
		}
		ofs << i << "," << freq
			<< "," << amp << "," << amp * 0.5 * T
			<< "," << phs
			<< "," << power << "," << power * T
			<< std::endl;
	};
	ofs << std::endl;

	ofs.close();

	return;
}

// Decimation
void Analysis::decimation( std::vector<CommonParameters::DataFileSet>& dataFileSets ) const{

	Control* ptrControl = Control::getInstance();
	const Control::ParamsForDecimation params = ptrControl->getParamsForDecimation();
	if(!params.applyDecimation){
		return;
	}

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeLogMessage("Decimation");
	
	double* coeff = new double[params.filterLength + 1];
	const double samplingFreq = ptrControl->getSamplingFrequency();
	const double freqStop = 0.5 * samplingFreq / static_cast<double>(params.decimationInterval);
	const double freqPass = freqStop / pow(10.0, params.transitionBandWidthInLogarithmicScale);
	Util::calculateFIRFilterCoeffsByLeastSquare( params.filterLength, true, samplingFreq, freqPass, freqStop, 1.0, 1.0, coeff );

	std::string fileName = "lpf_for_decimation.csv";
	std::ofstream ofs;
	ofs.open( fileName.c_str(), std::ios::out );
	if( ofs.fail() ){
		ptrOutputFiles->writeErrorMessage("File open error : " + fileName );
	}
	const double freqNyquist = 0.5 * samplingFreq;
	const int numFreq = 100000;
	for( int i = 0; i < numFreq; ++i ){
		const double freq = freqNyquist * static_cast<double>(i+1) / static_cast<double>(numFreq);
		const std::complex<double> Hfreq = Util::calculateFrequencyCharacteristicsOfFIRFilter( 
			params.filterLength, samplingFreq, freq, coeff );
		const double angle = atan2( Hfreq.imag(), Hfreq.real() ) * CommonParameters::RAD2DEG;
		ofs << std::setw(20) << std::setprecision(12) << freq << ","
			<< std::setw(20) << std::setprecision(12) << std::abs(Hfreq) << "," 
			<< std::setw(20) << std::setprecision(12) << angle << std::endl;
	}
	ofs.close();

	const int numChannels = ptrControl->getNumberOfChannels();
	for( std::vector<CommonParameters::DataFileSet>::iterator itr = dataFileSets.begin(); itr != dataFileSets.end(); ++itr ){
		const int numDataOrg = itr->numDataPoints;
		for( int iChan = 0; iChan < numChannels; ++iChan ){
			double* dataOrg = new double[numDataOrg];
			for( int i = 0; i < numDataOrg; ++i ){
				dataOrg[i] = itr->dataFile[iChan].data[i];
			}
			const int numDataMod = numDataOrg - params.filterLength;
			double* dataMod = new double[numDataMod];
			for( int i = params.filterLength; i < numDataOrg; ++i ){
				double dbuf(0.0);
				for( int k = 0; k <= params.filterLength; ++k ){
					dbuf += coeff[k] * dataOrg[i-k];
				}
				dataMod[i - params.filterLength] = dbuf;
			}
			delete [] dataOrg;
			int numDataAfterDecimation(0);
			for( int i = 0; i < numDataMod; i += params.decimationInterval, ++numDataAfterDecimation ){
				itr->dataFile[iChan].data[numDataAfterDecimation] = dataMod[i];
			}
			for( int i = numDataAfterDecimation; i < numDataOrg; ++i ){
				itr->dataFile[iChan].data[i] = 0.0;
			}
			itr->numDataPoints = numDataAfterDecimation;
			delete [] dataMod;
		}
	}

	delete [] coeff;

	ptrControl->setSamplingFrequency(samplingFreq / static_cast<double>(params.decimationInterval));

}

// Read time-series data
void Analysis::readTimeSeriesData( std::vector<CommonParameters::DataFileSet>& dataFileSets ){

	const Control* const ptrControl = Control::getInstance();
	const int numTimeSeriesSections = ptrControl->getNumTimeSeriesSections();
	const int numChannels = ptrControl->getNumberOfChannels();

	for( int iSect = 0; iSect < numTimeSeriesSections; ++iSect ){
		const int numDataPoints = dataFileSets[iSect].numDataPoints;
		std::vector<CommonParameters::DataFile>& dataFileList = dataFileSets[iSect].dataFile;
		int iChan = 0;
		if (ptrControl->doesReadElogDualBinary()) {
			// Read ELOG binary data as channel 0 and 1
			ElogDual* ptrElogDual = ElogDual::getInstance();
			const std::string fileName = dataFileList[0].fileName;
			const int numSkipData = dataFileList[0].numSkipData;
			iChan = 2;
			for (int ch = 0; ch < iChan; ++ch) {
				dataFileList[ch].data = new double[numDataPoints];
			}
			if (fileName.find("*.dat") != std::string::npos) {
				const int directoryNameLength = fileName.find_last_of("\\") + 1;
				const std::string directoryName = fileName.substr(0, directoryNameLength);
				ptrElogDual->readElogBinaryFilesUnderADirectory(directoryName, numSkipData, numDataPoints, dataFileList[0].data, dataFileList[1].data);
			}
			else{
				int counter(0);
				ptrElogDual->readElogBinaryFile(fileName, numSkipData, numDataPoints, counter, dataFileList[0].data, dataFileList[1].data);
			}
		}

		if( ptrControl->doesReadElogMTBinary() ){
			// Read ELOG binary data as channel 0 and 1
			ElogMT* ptrElogMT = ElogMT::getInstance();
			const std::string fileName = dataFileList[0].fileName;
			const int numSkipData = dataFileList[0].numSkipData;
			const int elogMTReadingOption = ptrControl->getElogMTReadingOption();
			if (elogMTReadingOption == Control::READ_EX_EY_HX_HY_HZ_FROM_ELOGMT_DATA ||
				elogMTReadingOption == Control::READ_EX_EY_HX_HY_HZ_HRX_HRY_FROM_ELOGMT_DATA) {
				{
					iChan = 5;
					for (int ch = 0; ch < iChan; ++ch) {
						dataFileList[ch].data = new double[numDataPoints];
					}
					if (fileName.find("*.dat") != std::string::npos) {
						const int directoryNameLength = fileName.find_last_of("\\") + 1;
						const std::string directoryName = fileName.substr(0, directoryNameLength);
						ptrElogMT->readElogBinaryFilesUnderADirectory(directoryName, numSkipData, numDataPoints, dataFileList[0].data, dataFileList[1].data, dataFileList[2].data, dataFileList[3].data, dataFileList[4].data);
					}
					else {
						int counter(0);
						ptrElogMT->readElogBinaryFile(fileName, numSkipData, numDataPoints, counter, dataFileList[0].data, dataFileList[1].data, dataFileList[2].data, dataFileList[3].data, dataFileList[4].data);
					}
				}
				if (elogMTReadingOption == Control::READ_EX_EY_HX_HY_HZ_HRX_HRY_FROM_ELOGMT_DATA) {
					const std::string fileName = dataFileList[iChan].fileName;
					const int numSkipDataR= dataFileList[iChan].numSkipData;
					for (int ch = iChan; ch < iChan + 2; ++ch) {
						dataFileList[ch].data = new double[numDataPoints];
					}
					if (fileName.find("*.dat") != std::string::npos) {
						const int directoryNameLength = fileName.find_last_of("\\") + 1;
						const std::string directoryName = fileName.substr(0, directoryNameLength);
						ptrElogMT->readElogBinaryFilesUnderADirectoryHxHyOnly(directoryName, numSkipData, numDataPoints, dataFileList[iChan].data, dataFileList[iChan+1].data);
					}
					else {
						int counter(0);
						ptrElogMT->readElogBinaryFileHxHyOnly(fileName, numSkipData, numDataPoints, counter, dataFileList[iChan].data, dataFileList[iChan+1].data);
					}
					iChan += 2;
				}
			}
			else if (elogMTReadingOption == Control::READ_EX_EY_HX_HY_FROM_ELOGMT_DATA ||
					elogMTReadingOption == Control::READ_EX_EY_HX_HY_HRX_HRY_FROM_ELOGMT_DATA) {
				{
					iChan = 4;
					for (int ch = 0; ch < iChan; ++ch) {
						dataFileList[ch].data = new double[numDataPoints];
					}
					if (fileName.find("*.dat") != std::string::npos) {
						const int directoryNameLength = fileName.find_last_of("\\") + 1;
						const std::string directoryName = fileName.substr(0, directoryNameLength);
						ptrElogMT->readElogBinaryFilesUnderADirectoryExEyHxHyOnly(directoryName, numSkipData, numDataPoints, dataFileList[0].data, dataFileList[1].data, dataFileList[2].data, dataFileList[3].data);
					}
					else {
						int counter(0);
						ptrElogMT->readElogBinaryFileExEyHxHyOnly(fileName, numSkipData, numDataPoints, counter, dataFileList[0].data, dataFileList[1].data, dataFileList[2].data, dataFileList[3].data);
					}
				}
				if (elogMTReadingOption == Control::READ_EX_EY_HX_HY_HRX_HRY_FROM_ELOGMT_DATA) {
					const std::string fileName = dataFileList[iChan].fileName;
					const int numSkipData = dataFileList[iChan].numSkipData;
					for (int ch = iChan; ch < iChan + 2; ++ch) {
						dataFileList[ch].data = new double[numDataPoints];
					}
					if (fileName.find("*.dat") != std::string::npos) {
						const int directoryNameLength = fileName.find_last_of("\\") + 1;
						const std::string directoryName = fileName.substr(0, directoryNameLength);
						ptrElogMT->readElogBinaryFilesUnderADirectoryHxHyOnly(directoryName, numSkipData, numDataPoints, dataFileList[iChan].data, dataFileList[iChan + 1].data);
					}
					else {
						int counter(0);
						ptrElogMT->readElogBinaryFileHxHyOnly(fileName, numSkipData, numDataPoints, counter, dataFileList[iChan].data, dataFileList[iChan + 1].data);
					}
					iChan += 2;
				}
			}
			else if (elogMTReadingOption == Control::READ_HZ_HX_HY_FROM_ELOGMT_DATA ||
					elogMTReadingOption == Control::READ_HZ_HX_HY_HRX_HRY_FROM_ELOGMT_DATA) {
				{
					iChan = 3;
					for (int ch = 0; ch < iChan; ++ch) {
						dataFileList[ch].data = new double[numDataPoints];
					}
					if (fileName.find("*.dat") != std::string::npos) {
						const int directoryNameLength = fileName.find_last_of("\\") + 1;
						const std::string directoryName = fileName.substr(0, directoryNameLength);
						ptrElogMT->readElogBinaryFilesUnderADirectoryHzHxHyOnly(directoryName, numSkipData, numDataPoints, dataFileList[0].data, dataFileList[1].data, dataFileList[2].data);
					}
					else {
						int counter(0);
						ptrElogMT->readElogBinaryFileHzHxHyOnly(fileName, numSkipData, numDataPoints, counter, dataFileList[0].data, dataFileList[1].data, dataFileList[2].data);
					}
				}
				if (elogMTReadingOption == Control::READ_HZ_HX_HY_HRX_HRY_FROM_ELOGMT_DATA) {
					const std::string fileName = dataFileList[iChan].fileName;
					const int numSkipData = dataFileList[iChan].numSkipData;
					for (int ch = iChan; ch < iChan + 2; ++ch) {
						dataFileList[ch].data = new double[numDataPoints];
					}
					if (fileName.find("*.dat") != std::string::npos) {
						const int directoryNameLength = fileName.find_last_of("\\") + 1;
						const std::string directoryName = fileName.substr(0, directoryNameLength);
						ptrElogMT->readElogBinaryFilesUnderADirectoryHxHyOnly(directoryName, numSkipData, numDataPoints, dataFileList[iChan].data, dataFileList[iChan + 1].data);
					}
					else {
						int counter(0);
						ptrElogMT->readElogBinaryFileHxHyOnly(fileName, numSkipData, numDataPoints, counter, dataFileList[iChan].data, dataFileList[iChan + 1].data);
					}
					iChan += 2;
				}
			}
			else if (elogMTReadingOption == Control::READ_EX_EY_FROM_ELOGMT_DATA) {
				iChan = 2;
				for (int ch = 0; ch < iChan; ++ch) {
					dataFileList[ch].data = new double[numDataPoints];
				}
				if (fileName.find("*.dat") != std::string::npos) {
					const int directoryNameLength = fileName.find_last_of("\\") + 1;
					const std::string directoryName = fileName.substr(0, directoryNameLength);
					ptrElogMT->readElogBinaryFilesUnderADirectoryExEyOnly(directoryName, numSkipData, numDataPoints, dataFileList[0].data, dataFileList[1].data);
				}
				else {
					int counter(0);
					ptrElogMT->readElogBinaryFileExEyOnly(fileName, numSkipData, numDataPoints, counter, dataFileList[0].data, dataFileList[1].data);
				}
			}
			else if (elogMTReadingOption == Control::READ_HX_HY_HRX_HRY_FROM_ELOGMT_DATA ||
					elogMTReadingOption == Control::READ_HX_HY_FROM_ELOGMT_DATA) {
				{// Hx, Hy
					const std::string fileName = dataFileList[iChan].fileName;
					const int numSkipData = dataFileList[iChan].numSkipData;
					for (int ch = iChan; ch < iChan + 2; ++ch) {
						dataFileList[ch].data = new double[numDataPoints];
					}
					if (fileName.find("*.dat") != std::string::npos) {
						const int directoryNameLength = fileName.find_last_of("\\") + 1;
						const std::string directoryName = fileName.substr(0, directoryNameLength);
						ptrElogMT->readElogBinaryFilesUnderADirectoryHxHyOnly(directoryName, numSkipData, numDataPoints, dataFileList[iChan].data, dataFileList[iChan + 1].data);
					}
					else {
						int counter(0);
						ptrElogMT->readElogBinaryFileHxHyOnly(fileName, numSkipData, numDataPoints, counter, dataFileList[iChan].data, dataFileList[iChan + 1].data);
					}
					iChan += 2;
				}
				if(elogMTReadingOption == Control::READ_HX_HY_HRX_HRY_FROM_ELOGMT_DATA){
					// Hrx, Hry
					const std::string fileName = dataFileList[iChan].fileName;
					const int numSkipData = dataFileList[iChan].numSkipData;
					for (int ch = iChan; ch < iChan + 2; ++ch) {
						dataFileList[ch].data = new double[numDataPoints];
					}
					if (fileName.find("*.dat") != std::string::npos) {
						const int directoryNameLength = fileName.find_last_of("\\") + 1;
						const std::string directoryName = fileName.substr(0, directoryNameLength);
						ptrElogMT->readElogBinaryFilesUnderADirectoryHxHyOnly(directoryName, numSkipData, numDataPoints, dataFileList[iChan].data, dataFileList[iChan + 1].data);
					}
					else {
						int counter(0);
						ptrElogMT->readElogBinaryFileHxHyOnly(fileName, numSkipData, numDataPoints, counter, dataFileList[iChan].data, dataFileList[iChan + 1].data);
					}
					iChan += 2;
				}
			}
			else{
				OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
				ptrOutputFiles->writeErrorMessage("Unsupported type of ELOG-MT reading option: " + Util::toString(ptrControl->getElogMTReadingOption()));
			}
		}
		for( ; iChan < numChannels; ++iChan ){
			const std::string fileName = dataFileList[iChan].fileName;
			const int numSkipData = dataFileList[iChan].numSkipData;
			dataFileList[iChan].data = new double[numDataPoints];
			if (ptrControl->doesReadAtsBinary() && Util::extractExtensionOfFileName(fileName).find("ats") != std::string::npos) {
				Ats* ptrAts = Ats::getInstance();
				ptrAts->readAtsFile( fileName, numSkipData, numDataPoints, dataFileList[iChan].data );
			}else{
				readOneTimeSeriesData( fileName, numSkipData, numDataPoints, dataFileList[iChan].data );
			}
		}
	}

}

// Read calibration files
void Analysis::readCalibrationFiles( const std::vector<double>& freq ){

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	const Control* const ptrControl = Control::getInstance();

	if (ptrControl->doesMakeCalibrationFileForMFS()) {
		const int numFile = ptrControl->getNumCalibrationFilesForMFS();
		if (numFile > 0) {
			assert(numFile == ptrControl->getNumberOfChannels());
		}
		for (int iFile = 0; iFile < numFile; ++iFile) {
			const std::string fileName = ptrControl->getCalibrationFileNameForMFS(iFile);
			const Ats* ptrAts = Ats::getInstance();
			ptrAts->makeCalibrationFile(fileName, iFile, freq);
		}
	}

	if( ptrControl->doesMakeCalibrationFileForElogDual() ){
		const Control::ParamsForElogCalibration params = ptrControl->getParamsForElogDualCalibration();
		ElogDual* ptrElogDual = ElogDual::getInstance();
		if (ptrControl->doesMakeCalibrationFileForMFS()) {
			double dbuf(0.0);
			std::istringstream issX(ptrControl->getCalibrationFileNameForMFS(0));
			issX >> dbuf;
			const double dipoleLengthX = dbuf;
			std::istringstream issY(ptrControl->getCalibrationFileNameForMFS(1));
			issY >> dbuf;
			const double dipoleLengthY = dbuf;
			ptrElogDual->makeCalibrationFile(params.fileName, params.unitGroupDelay, 0, 1, dipoleLengthX, dipoleLengthY, freq);
		}
	}

	if (ptrControl->doesMakeCalibrationFileForElogMT()) {
		const Control::ParamsForElogCalibration params = ptrControl->getParamsForElogMTCalibration();
		ElogMT* ptrElogMT = ElogMT::getInstance();
		if (ptrControl->doesMakeCalibrationFileForMFS()) {
			double dbuf(0.0);
			std::istringstream issX(ptrControl->getCalibrationFileNameForMFS(0));
			issX >> dbuf;
			const double dipoleLengthX = dbuf;
			std::istringstream issY(ptrControl->getCalibrationFileNameForMFS(1));
			issY >> dbuf;
			const double dipoleLengthY = dbuf;
			std::vector<int> channelIndexes;
			channelIndexes.push_back(0);
			channelIndexes.push_back(1);
			ptrElogMT->makeCalibrationFile(params.fileName, params.unitGroupDelay, channelIndexes, dipoleLengthX, dipoleLengthY, freq);
		}
	}

	const int numCalFile = ptrControl->getNumCalibrationFiles();
	if( numCalFile > 0 ){
		assert( numCalFile == ptrControl->getNumberOfChannels() );
	}
	if( m_calibrationFunctions != NULL ){
		delete [] m_calibrationFunctions;
	}
	m_calibrationFunctions = new CalibrationFunction[numCalFile];
	for( int iFile = 0; iFile < numCalFile; ++iFile ){
		m_calibrationFunctions[iFile].readCalibrationFunction(ptrControl->getCalibrationFileName(iFile));
	}

}

// Read one time-series data
void Analysis::readOneTimeSeriesData( const std::string& fileName, const int numSkipData, const int numDataPoints, double* data ){

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeLogMessage("Read data from "+fileName);

	std::ifstream ifs( fileName.c_str(), std::ios::in );
	if( ifs.fail() ){
		ptrOutputFiles->writeErrorMessage("File open error : " + fileName);
	}

	for( int i = 0; i < numSkipData; ++i ){
		double dbuf(0.0);
		ifs >> dbuf;
	}

	for( int i = 0; i < numDataPoints; ++i ){
		if( ifs.eof() ){
			ptrOutputFiles->writeErrorMessage("End of file appears before reading all data points");
		}
		ifs >> data[i];
	}

	ifs.close();

}

// Perform preprocessing
void Analysis::preprocessing( std::vector<CommonParameters::DataFileSet>& dataFileSets ){

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeLogMessage("Perform preprocessing");

	const Control* const ptrControl = Control::getInstance();
	const double samplingFrequency = ptrControl->getSamplingFrequency();

	const int numChannels = ptrControl->getNumberOfChannels();
	int iSection(0);
	for( std::vector<CommonParameters::DataFileSet>::iterator itr = dataFileSets.begin(); itr != dataFileSets.end(); ++itr, ++iSection ){
		const int numData = itr->numDataPoints;
		for( int iChan = 0; iChan < numChannels; ++iChan ){
			ptrOutputFiles->writeLogMessage("Secton " + Util::toString(iSection) + ", Channel " + Util::toString(iChan) );
			// Hampel filter
			const Control::ParamsForHampelFilter paramsHampelFilter = ptrControl->getParamsForHampelFilter();
			if( paramsHampelFilter.applyHampelFilter ){
				const int numModifiedPoints = Util::hampelFilter( numData, paramsHampelFilter.numNeighborsOnEitherSide, paramsHampelFilter.nsigma, itr->dataFile[iChan].data );
				ptrOutputFiles->writeLogMessage("Number of the points modified by Hampel filter: " + Util::toString(numModifiedPoints) );
			}

			// Make new data records by subtracting mean
			const double mean = Util::calculateMeanValue(numData, itr->dataFile[iChan].data);
			ptrOutputFiles->writeLogMessage("Subtract mean (" + Util::toString(mean) + ")");
			for (int i = 0; i < numData; ++i) {
				itr->dataFile[iChan].data[i] -= mean;
			}
#ifdef _DEBUG_WRITE
			std::ostringstream oss;
			oss << "ts_sect_" << iSection << "_chan_" << iChan << "_mean.csv";
			std::ofstream ofs;
			ofs.open(oss.str().c_str(), std::ios::out);
			if (ofs.fail()) {
				ptrOutputFiles->writeLogMessage("File open error !! : " + oss.str());
			}
			for (int i = 0; i < numData; ++i) {
				const double elapsedTime = static_cast<double>(i) / samplingFrequency;
				ofs << ptrControl->getTimeFromStartTimeOfEachSection(iSection, elapsedTime) << ","
					<< std::setprecision(12) << std::scientific << itr->dataFile[iChan].data[i] << std::endl;
			}
			ofs.close();
#endif

			// Apply high-pass filter
			if( ptrControl->doesApplyIIRHighPassFilter() ){
				ptrOutputFiles->writeLogMessage("Apply high-pass filter");
				Util::applyIIRHighPassFilter(samplingFrequency, ptrControl->getCutoffFrequencyForIIRHighPassFilter(),
					numData, itr->dataFile[iChan].data);
				std::string fileName = "highpass.csv";
				std::ofstream ofs;
				ofs.open( fileName.c_str(), std::ios::out );
				if( ofs.fail() ){
					ptrOutputFiles->writeErrorMessage("File open error : " + fileName );
				}
				const double freqNyquist = 0.5 * samplingFrequency;
				const int numFreq = 100000;
				for( int i = 0; i < numFreq; ++i ){
					const double freq = freqNyquist * static_cast<double>(i+1) / static_cast<double>(numFreq);
					const std::complex<double> Hfreq = Util::calculateFrequencyCharacteristicsOfIIRHighPassFilter( freq,
							samplingFrequency, ptrControl->getCutoffFrequencyForIIRHighPassFilter() );
					const double angle = atan2( Hfreq.imag(), Hfreq.real() ) * CommonParameters::RAD2DEG;
					ofs << std::setw(20) << std::setprecision(12) << freq << ","
						<< std::setw(20) << std::setprecision(12) << std::abs(Hfreq) << "," 
						<< std::setw(20) << std::setprecision(12) << angle << std::endl;
				}
				ofs.close();		
			}

			// Apply low-pass filter
			if( ptrControl->doesApplyIIRLowPassFilter() ){
				ptrOutputFiles->writeLogMessage("Apply low-pass filter");
				Util::applyIIRLowPassFilter( samplingFrequency, ptrControl->getCutoffFrequencyForIIRLowPassFilter(),
					numData, itr->dataFile[iChan].data );
				std::string fileName = "lowpass.csv";
				std::ofstream ofs;
				ofs.open( fileName.c_str(), std::ios::out );
				if( ofs.fail() ){
					ptrOutputFiles->writeErrorMessage("File open error : " + fileName );
				}
				const double freqNyquist = 0.5 * samplingFrequency;
				const int numFreq = 100000;
				for( int i = 0; i < numFreq; ++i ){
					const double freq = freqNyquist * static_cast<double>(i+1) / static_cast<double>(numFreq);
					const std::complex<double> Hfreq = Util::calculateFrequencyCharacteristicsOfIIRLowPassFilter( freq,
						samplingFrequency, ptrControl->getCutoffFrequencyForIIRLowPassFilter() );
					const double angle = atan2( Hfreq.imag(), Hfreq.real() ) * CommonParameters::RAD2DEG;
					ofs << std::setw(20) << std::setprecision(12) << freq << ","
						<< std::setw(20) << std::setprecision(12) << std::abs(Hfreq) << "," 
						<< std::setw(20) << std::setprecision(12) << angle << std::endl;
				}
				ofs.close();		
			}

			// Apply notch filter
			if( ptrControl->getNumberOfCutoffFrequenciesForNotchFilter() > 0 ){
				ptrOutputFiles->writeLogMessage("Apply notch filter");
				const std::vector<double> cutoffFreqsForNotchFilter = ptrControl->getCutoffFrequenciesForNotchFilter();
				for( std::vector<double>::const_iterator itrCutoffFreq = cutoffFreqsForNotchFilter.begin();
					itrCutoffFreq != cutoffFreqsForNotchFilter.end(); ++itrCutoffFreq ){
					Util::applyNotchFilter(ptrControl->getParameterQForNotchFilter(), samplingFrequency, *itrCutoffFreq, 
						numData, itr->dataFile[iChan].data);
				}
				std::string fileName = "notch.csv";
				std::ofstream ofs;
				ofs.open( fileName.c_str(), std::ios::out );
				if( ofs.fail() ){
					ptrOutputFiles->writeErrorMessage("File open error : " + fileName );
				}
				const double freqNyquist = 0.5 * samplingFrequency;
				const int numFreq = 100000;
				for( int i = 0; i < numFreq; ++i ){
					const double freq = freqNyquist * static_cast<double>(i+1) / static_cast<double>(numFreq);
					std::complex<double> Hfreq = 1.0;
					for( std::vector<double>::const_iterator itrCutoffFreq = cutoffFreqsForNotchFilter.begin();
						itrCutoffFreq != cutoffFreqsForNotchFilter.end(); ++itrCutoffFreq ){
						Hfreq *= Util::calculateFrequencyCharacteristicsOfNotchFilter( freq, ptrControl->getParameterQForNotchFilter(),
							samplingFrequency, *itrCutoffFreq );
					}
					const double angle = atan2( Hfreq.imag(), Hfreq.real() ) * CommonParameters::RAD2DEG;
					ofs << std::setw(20) << std::setprecision(12) << freq << ","
						<< std::setw(20) << std::setprecision(12) << std::abs(Hfreq) << "," 
						<< std::setw(20) << std::setprecision(12) << angle << std::endl;
				}
				ofs.close();
			}
		}
	}
	
	if( (ptrControl->getParamsForPrewhitening()).applyPrewhitening ){
		m_coefficientsOfARModel = new std::vector<double>[numChannels];
		const Control::ParamsForPrewhitening params = ptrControl->getParamsForPrewhitening();
		if( params.typeOfEstimator == Control::MANUAL_INPUTS_AR_COEFFICIENTS_FOR_PREWHITENING ){
			( RobustPrewhitening::getInstance() )->prewhiteningUsingUserDefinedARCoeffs(dataFileSets, m_coefficientsOfARModel);
		}else{
			( RobustPrewhitening::getInstance() )->robustPrewhitening(dataFileSets, m_coefficientsOfARModel);
		}
	}

}

// Convert time-series data to frequency-domain data
void Analysis::convertToFrequencyData( const int segmentLength, const std::vector<CommonParameters::DataFileSet>& dataFileSets,
	int& numSegmentsTotal, std::complex<double>*** cdata, std::vector< std::pair<std::string, std::string> >& times, std::vector<double>* meanSquares ){

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeLogMessage("Convert time-series data to frequency-domain");

	const Control* const ptrControl = Control::getInstance();

	const int numChannels = ptrControl->getNumberOfChannels();
	const double overlappingRatio = ptrControl->getOverlappingRatio();
	const int shiftLength = static_cast<int>( (1.0 - overlappingRatio) * segmentLength );

	// Count total number of segments
	numSegmentsTotal = 0;
	for( std::vector<CommonParameters::DataFileSet>::const_iterator itr = dataFileSets.begin(); itr != dataFileSets.end(); ++itr ){
		if( itr->numDataPoints < segmentLength ){
			numSegmentsTotal += 1;
		}else{
			numSegmentsTotal += ( itr->numDataPoints - segmentLength ) / shiftLength + 1;
		}
	}
	if( numSegmentsTotal < 1 ){
		ptrOutputFiles->writeErrorMessage("Total number of segments is less than 1");
	}
	ptrOutputFiles->writeLogMessage("Total number of segments : " + Util::toString(numSegmentsTotal));
	
	times.reserve(numSegmentsTotal);

	// Make array of data segmetns
	const double samplingFrequency = ptrControl->getSamplingFrequency();
	double*** dataSegments = new double**[numChannels];
	for( int i = 0; i < numChannels; ++i ){
		dataSegments[i] = new double*[numSegmentsTotal];
	}
	int section(0);
	int counterSegment(0);
	for( std::vector<CommonParameters::DataFileSet>::const_iterator itr = dataFileSets.begin(); itr != dataFileSets.end(); ++itr, ++section ){
		const int numSegments = itr->numDataPoints < segmentLength ? 1 : ( itr->numDataPoints - segmentLength ) / shiftLength + 1;
		for( int iSeg = 0; iSeg < numSegments; ++iSeg, ++counterSegment ){
			for( int iChan = 0; iChan < numChannels; ++iChan ){
				dataSegments[iChan][counterSegment] = new double[segmentLength];
				const int index = iSeg * shiftLength;
				if( itr->numDataPoints < segmentLength ){
					memcpy(dataSegments[iChan][counterSegment], &(itr->dataFile[iChan].data[index]), sizeof(double)*itr->numDataPoints);
					for( int i = itr->numDataPoints; i < segmentLength; ++i ){
						dataSegments[iChan][counterSegment][i] = 0.0;
					}
				}else{
					memcpy(dataSegments[iChan][counterSegment], &(itr->dataFile[iChan].data[index]), sizeof(double)*segmentLength);
				}
//#ifdef _DEBUG_WRITE
//				std::ostringstream oss;
//				oss << "sect_" << section << "_seg_" << counterSegment << "_chan_" << iChan << ".csv"; 
//				std::ofstream ofs;
//				ofs.open( oss.str().c_str(), std::ios::out );
//				if( ofs.fail() ){
//					ptrOutputFiles->writeLogMessage("File open error !! : " + oss.str());
//				}
//				for( int i = 0; i < segmentLength; ++i ){
//					ofs << std::setprecision(12) << std::scientific << dataSegments[iChan][counterSegment][i] << std::endl;
//				}
//				ofs.close();
//#endif
			}
			// Start time
			const int index1 = iSeg * shiftLength;
			const double elapsedTime1 = static_cast<double>(index1) / samplingFrequency;
			std::string time1 = ptrControl->getTimeFromStartTimeOfEachSection(section, elapsedTime1);
			// End time
			const int index2 = index1 + segmentLength;
			const double elapsedTime2 = static_cast<double>(index2) / samplingFrequency;
			std::string time2 = ptrControl->getTimeFromStartTimeOfEachSection(section, elapsedTime2);
			// Insert
			times.push_back( std::make_pair(time1, time2) );
		}
	}
	assert( times.size() == numSegmentsTotal );

	// Calculate mean squares
	for( int iChan = 0; iChan < numChannels; ++iChan ){
		meanSquares[iChan].reserve(numSegmentsTotal);
		for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
			const double  meanSquare = Util::calculateMeanSquareValue(segmentLength, dataSegments[iChan][iSeg]);
			meanSquares[iChan].push_back(meanSquare);
		}
	}

	// Apply hanning window
	for( int iChan = 0; iChan < numChannels; ++iChan ){
		for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
			Util::hanningWindow(segmentLength, dataSegments[iChan][iSeg]);
//#ifdef _DEBUG_WRITE
//			std::ostringstream oss;
//			oss << "seg_" << iSeg << "_chan_" << iChan << "_hanning.csv"; 
//			std::ofstream ofs;
//			ofs.open( oss.str().c_str(), std::ios::out );
//			if( ofs.fail() ){
//				ptrOutputFiles->writeLogMessage("File open error !! : " + oss.str());
//			}
//			for( int i = 0; i < segmentLength; ++i ){
//				ofs << std::setprecision(12) << std::scientific << dataSegments[iChan][iSeg][i] << std::endl;
//			}
//			ofs.close();
//#endif
		}
	}

	// Fourier transform
	for( int iChan = 0; iChan < numChannels; ++iChan ){
		cdata[iChan] = new std::complex<double>*[numSegmentsTotal];
		for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
			cdata[iChan][iSeg] = new std::complex<double>[segmentLength];
			for( int i = 0; i < segmentLength; ++i ){
				cdata[iChan][iSeg][i] = std::complex<double>(dataSegments[iChan][iSeg][i], 0.0);
			}
			Util::fourierTransform(segmentLength, cdata[iChan][iSeg]);
		}
	}

	// Delete arrays
	for( int iChan = 0; iChan < numChannels; ++iChan ){
		for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
			delete [] dataSegments[iChan][iSeg];
		}
		delete [] dataSegments[iChan];
	}
	delete [] dataSegments;

}

// Perform calibration correction for all channels
void Analysis::calibrationCorrectionAllChannels( const int numChannels, const int numSegmentsTotal, 
									 const double freq, std::complex<double>** ftval ) const{

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeLogMessage("Perform calibration correction");

	// Calibration correction
	for( int iChan = 0; iChan < numChannels; ++iChan ){
		calibrationCorrection( iChan, numSegmentsTotal, freq, ftval[iChan] );
	}

}

// Perform calibration correction
void Analysis::calibrationCorrection( const int iChan, const int numSegmentsTotal, const double freq, std::complex<double>* ftval, const bool afterPreprocessing ) const{

	const Control* const ptrControl = Control::getInstance();

	// Interpolate calibration function
	std::complex<double> calCorrFunc(1.0, 0.0);
	if( ptrControl->getNumCalibrationFiles() == 0 ){
		// Calibration function is not used
	}else{
		calCorrFunc *= m_calibrationFunctions[iChan].calculateCalibrationFunction(freq);
	}
	
	if(ptrControl->getParamsForDecimation().applyDecimation){
		const Control::ParamsForDecimation params = ptrControl->getParamsForDecimation();
		double* coeff = new double[params.filterLength + 1];
		const double samplingFreqAfterDecimation = ptrControl->getSamplingFrequency();
		const double samplingFreqBeforeDecimation = samplingFreqAfterDecimation * static_cast<double>(params.decimationInterval);
		const double freqStop = 0.5 * samplingFreqBeforeDecimation;
		const double freqPass = freqStop / pow(10.0, params.transitionBandWidthInLogarithmicScale);
		Util::calculateFIRFilterCoeffsByLeastSquare( params.filterLength, true, samplingFreqBeforeDecimation, freqPass, freqStop, 1.0, 1.0, coeff );
		const std::complex<double> Hfreq = Util::calculateFrequencyCharacteristicsOfFIRFilter( params.filterLength, samplingFreqBeforeDecimation, freq, coeff );
		calCorrFunc /= Hfreq;
		delete [] coeff;
	}

	// Adjust the scaling factor for the loss due to the Hanning tapering
	calCorrFunc *= sqrt(8.0/3.0);

	if(afterPreprocessing){
		// Adjust the influences of preprocessing 
		const double samplingFrequency = ptrControl->getSamplingFrequency();

		// Adjust the effect of high-pass filter
		if( ptrControl->doesApplyIIRHighPassFilter() ){
			std::complex<double> Hfreq = Util::calculateFrequencyCharacteristicsOfIIRHighPassFilter( freq, samplingFrequency, 
				ptrControl->getCutoffFrequencyForIIRHighPassFilter() );
			calCorrFunc /= Hfreq;
		}

		// Adjust the effect of low-pass filter
		if( ptrControl->doesApplyIIRLowPassFilter() ){
			const double Q = 1.0/sqrt(2.0);
			std::complex<double> Hfreq = Util::calculateFrequencyCharacteristicsOfIIRLowPassFilter( freq, samplingFrequency, 
				ptrControl->getCutoffFrequencyForIIRLowPassFilter() );
			calCorrFunc /= Hfreq;
		}

		// Adjust the effect of notch filter
		if( ptrControl->getNumberOfCutoffFrequenciesForNotchFilter() > 0 ){
			const std::vector<double> cutoffFreqsForNotchFilter = ptrControl->getCutoffFrequenciesForNotchFilter();
			std::complex<double> Hfreq = std::complex<double>(1.0, 0.0);
			for( std::vector<double>::const_iterator itrCutoffFreq = cutoffFreqsForNotchFilter.begin();
				itrCutoffFreq != cutoffFreqsForNotchFilter.end(); ++itrCutoffFreq ){
				Hfreq *= Util::calculateFrequencyCharacteristicsOfNotchFilter( freq, ptrControl->getParameterQForNotchFilter(),
					samplingFrequency, *itrCutoffFreq );
			}
			if( std::abs(Hfreq) > CommonParameters::EPS ){
				calCorrFunc /= Hfreq;
			}
		}

		// Adjust prewhitening
		if( ptrControl->getParamsForPrewhitening().applyPrewhitening && m_coefficientsOfARModel != NULL ){
			const std::vector<double>& coeffs = m_coefficientsOfARModel[iChan];
			const double omega = 2.0 * CommonParameters::PI * freq / samplingFrequency;
			std::complex<double> H = std::complex<double>(1.0, 0.0);// Initialize 1
			int iAR = 1;// From 1
			for( std::vector<double>::const_iterator itr = coeffs.begin(); itr != coeffs.end(); ++itr, ++iAR ){
				const double coeff = *itr;
				const double arg = - static_cast<double>(iAR) * omega;
				// 1.0 - c1 * exp(-i*omega*1) - c2 * exp(-i*omega*2) - ... - ck * exp(-i*omega*k) - ...
				H -= std::complex<double>( cos(arg), sin(arg) ) * coeff;
			}
			calCorrFunc /= H;
		}
	}

	for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
		ftval[iSeg] *= calCorrFunc;
	}

}

// Calculate apparent resistivity
double Analysis::calcApparentResistivity( const double freq, const std::complex<double>& Z ) const{

	return 0.2 * std::norm(Z) / freq;

}

// Calculate error of apparent resistivity
double Analysis::calcApparentResistivityError( const double freq, const std::complex<double>& Z, const double dZ ) const{

	return 0.4 * std::abs(Z) * dZ / freq;

}

// Calculate diagonal components of hat matrix
double Analysis::calculateDiagonalComponentsOfHatMatrix( const int numSegments, const int channelX, const int channelY, 
	std::complex<double>** data, const double* const weights, double* hatDiagonals ) const{

	double GxxMean(0.0);
	double GyyMean(0.0);
	std::complex<double> GxyMean(0.0, 0.0);
	double sumWeights(0.0);
	for( int i = 0; i < numSegments; ++i ){
		GxxMean += std::norm(data[channelX][i]) * weights[i];
		GyyMean += std::norm(data[channelY][i]) * weights[i];
		GxyMean += data[channelX][i] * std::conj(data[channelY][i]) * weights[i];
		sumWeights += weights[i];
	}
	GxxMean /= sumWeights;
	GyyMean /= sumWeights;
	GxyMean /= sumWeights;
	const double coh = std::norm(GxyMean)/GxxMean/GyyMean;
	const double coh2 = coh * coh;

	double maxHatDiagonal(0.0);
	for( int i = 0; i < numSegments; ++i ){
		const double Gxx = std::norm(data[channelX][i]) * weights[i];
		const double Gyy = std::norm(data[channelY][i]) * weights[i];
		const std::complex<double> Gxy = data[channelX][i] * std::conj(data[channelY][i]) * weights[i];
		const double term0 = Gxx / GxxMean + Gyy / GyyMean;
		const double term1 = 2.0 * coh2 * (Gxy/GxyMean).real();
		hatDiagonals[i] = ( term0 - term1 ) / sumWeights / (1.0 - coh2);
		if( hatDiagonals[i] > maxHatDiagonal ){
			maxHatDiagonal =  hatDiagonals[i];
		}
#ifdef _DEBUG_WRITE
		std::cout << "Gxx Gyy Gxy term0 term1 hatDiagonals maxHatDiagonal: " << 
			Gxx << " " << Gyy << " " << Gxy << " " << term0 << " " << term1 << " " << hatDiagonals[i] << " " << maxHatDiagonal << std::endl;
#endif
	}

	return maxHatDiagonal;

}

// Calculate phase
double Analysis::calcPhase( const std::complex<double>& Z ) const{

	return std::arg(Z) * CommonParameters::RAD2DEG;

}

// Calculate error of phase
double Analysis::calcPhaseError( const std::complex<double>& Z, const double dZ ) const{

	const double ratio = dZ / std::abs(Z);
	if( ratio < 1.0 ){
		return asin(ratio) * CommonParameters::RAD2DEG;
	}else{
		return 360.0;
	}

}

// Calculate response functions by iteratively reweighted least squares
void Analysis::calculateResponseFunctionsByIRWLS( const int iRobustWeight, const std::complex<double>* const out, const std::complex<double>* const in0,
	const std::complex<double>* const in1, const int numSegments, const bool fixScale, double& scale, const double* const weightsPrior,
	double* weights, std::complex<double>* residuals, std::complex<double>& resp0, std::complex<double>& resp1, double& coherence, 
	std::vector<std::string>& titles, std::vector<double>* outputValues, const bool priorityOnFirst ) const{

	RobustWeight* const robustWeightCur = getPointerToRobustWeight(iRobustWeight);
	if (robustWeightCur == NULL) {
		return;
	}
	const std::string nameOfRobustWeight = robustWeightCur->getNameOfRobustWeight();

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeCvgAndLogMessage("Calculate response functions by iteratively reweighted least squares using " + nameOfRobustWeight + " weight");

	const int numIterationMax = robustWeightCur->getNumIterationMax();
	const double convergenceCriteria = robustWeightCur->getConvergenceCriteria();
	double GresPre(0.0);
	bool converge(false);
	for( int iter = 0; iter < numIterationMax; ++iter ){
		ptrOutputFiles->writeCvgMessage("Iteration number = " + Util::toString(iter));
		if (!fixScale) {
			// Calcuate scale
			scale = robustWeightCur->calculateScale(numSegments, residuals, scale);
		}
		// Calcuate weights
		const double sumOfWeights = robustWeightCur->calculateWeights(numSegments, residuals, scale, weightsPrior, weights);
		if( sumOfWeights < CommonParameters::EPS ){
			ptrOutputFiles->writeCvgMessage("	Sum of weights (" +  Util::toString(sumOfWeights) + ") is too small");
			ptrOutputFiles->writeWarningMessage("Sum of weights (" +  Util::toString(sumOfWeights) + ") is too small");
			break;
		}
		if( outputValues != NULL ){
			titles.push_back(nameOfRobustWeight + "_" + Util::toString(iter));
			for( int iSeg = 0 ; iSeg < numSegments; ++iSeg){
				outputValues[iSeg].push_back(residuals[iSeg].real());
				outputValues[iSeg].push_back(residuals[iSeg].imag());
				outputValues[iSeg].push_back(weights[iSeg]);
			}
		}
		const double Gres = calculateResponseFunctionByWLS( out, in0, in1, numSegments, weights, residuals, resp0, resp1, coherence, priorityOnFirst );
		ptrOutputFiles->writeCvgMessage("	Weighted residual power: " + Util::toString(Gres));
		if( Gres < CommonParameters::EPS ){
			ptrOutputFiles->writeCvgMessage("	Weighted residual power is too small (" +  Util::toString(Gres) + ")");
			ptrOutputFiles->writeWarningMessage("Weighted residual power is too small (" +  Util::toString(Gres) + ")");
			break;
		}
		if( fabs(Gres-GresPre) / fabs(GresPre) < convergenceCriteria ){
			// Converged
			ptrOutputFiles->writeCvgMessage("Iteration using " + nameOfRobustWeight + " weight converged");
			converge = true;
			break;
		}
		GresPre = Gres;
	}

	if(!converge){
		ptrOutputFiles->writeCvgMessage("	Iteration using " + nameOfRobustWeight + " weight does not converge");
	}

}

// Calculate response functions by iteratively reweighted remote reference
void Analysis::calculateResponseFunctionsByIRWLSRemoteReference(const int iRobustWeight, const std::complex<double>* const out,
	const std::complex<double>* const in0, const std::complex<double>* const in1,
	const std::complex<double>* const rr0, const std::complex<double>* const rr1,
	const int numSegments, const bool fixScale, double& scale, const double* const weightsPrior,
	double* weights, std::complex<double>* residuals, std::complex<double>& resp0, std::complex<double>& resp1, double& coherence,
	std::vector<std::string>& titles, std::vector<double>* outputValues) const {

	RobustWeight* const robustWeightCur = getPointerToRobustWeight(iRobustWeight);
	if (robustWeightCur == NULL) {
		return;
	}
	const std::string nameOfRobustWeight = robustWeightCur->getNameOfRobustWeight();

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeCvgAndLogMessage("Calculate response functions by iteratively reweighted remote reference using " + nameOfRobustWeight + " weight");

	const int numIterationMax = robustWeightCur->getNumIterationMax();
	const double convergenceCriteria = robustWeightCur->getConvergenceCriteria();
	double GresPre(0.0);
	bool converge(false);
	for (int iter = 0; iter < numIterationMax; ++iter) {
		ptrOutputFiles->writeCvgMessage("Iteration number = " + Util::toString(iter));
		if (!fixScale) {
			// Calcuate scale
			scale = robustWeightCur->calculateScale(numSegments, residuals, scale);
		}
		// Calcuate weights
		const double sumOfWeights = robustWeightCur->calculateWeights(numSegments, residuals, scale, weightsPrior, weights);
		if (sumOfWeights < CommonParameters::EPS) {
			ptrOutputFiles->writeCvgMessage("	Sum of weights (" + Util::toString(sumOfWeights) + ") is too small");
			ptrOutputFiles->writeWarningMessage("Sum of weights (" + Util::toString(sumOfWeights) + ") is too small");
			break;
		}
		if (outputValues != NULL) {
			titles.push_back(nameOfRobustWeight + "_" + Util::toString(iter));
			for (int iSeg = 0; iSeg < numSegments; ++iSeg) {
				outputValues[iSeg].push_back(residuals[iSeg].real());
				outputValues[iSeg].push_back(residuals[iSeg].imag());
				outputValues[iSeg].push_back(weights[iSeg]);
			}
		}
		const double Gres = calculateResponseFunctionByWLSRemoteReference(out, in0, in1, rr0, rr1, numSegments, weights, residuals, resp0, resp1, coherence);
		ptrOutputFiles->writeCvgMessage("	Weighted residual power: " + Util::toString(Gres));
		if (Gres < CommonParameters::EPS) {
			ptrOutputFiles->writeCvgMessage("	Weighted residual power is too small (" + Util::toString(Gres) + ")");
			ptrOutputFiles->writeWarningMessage("Weighted residual power is too small (" + Util::toString(Gres) + ")");
			break;
		}
		if (fabs(Gres - GresPre) / fabs(GresPre) < convergenceCriteria) {
			// Converged
			ptrOutputFiles->writeCvgMessage("Iteration using " + nameOfRobustWeight + " weight converged");
			converge = true;
			break;
		}
		GresPre = Gres;
	}

	if (!converge) {
		ptrOutputFiles->writeCvgMessage("	Iteration using " + nameOfRobustWeight + " weight does not converge");
	}
}

// Calculate response function by the ordinary remote reference method
double Analysis::calculateResponseFunctionByOrdinaryRemoteReference( const int numSegments, const std::complex<double>* const out,
	const std::complex<double>* const in1, const std::complex<double>* const in2, 
	const std::complex<double>* const rr1, const std::complex<double>* const rr2, 
	std::complex<double>& resp1, std::complex<double>& resp2 ) const{

	double* weights = new double[numSegments];
	for (int iSeg = 0; iSeg < numSegments; ++iSeg) {
		weights[iSeg] = 1.0;
	}
	calculateResponseFunctionByWLSRemoteReferenceAux(numSegments, out, in1, in2, rr1, rr2, weights, resp1, resp2);
	delete[] weights;

	double Gvar(0.0);
	double Gsyn(0.0);
	double Gres(0.0);
	for( int i = 0; i < numSegments; ++i ){
		Gvar += std::norm(out[i]);
		const std::complex<double> outSyn = resp1 * in1[i] + resp2 * in2[i];
		Gsyn += std::norm(outSyn);
		const std::complex<double> residual = out[i] - outSyn;
		Gres += std::norm(residual);
	}

	double coherence = 1.0 - Gres / Gvar;
	if( coherence > 1.0 ){
		coherence = 1.0;
	}else if( coherence < 0.0 ){
		coherence = 0.0;
	}
	return coherence;

}

// Calculate response function by the ordinary weighted square method
double Analysis::calculateResponseFunctionByWLS( const std::complex<double>* const out, const std::complex<double>* const in0,
	const std::complex<double>* const in1, const int numSegments, const double* const weights, std::complex<double>* residuals,
	std::complex<double>& resp0, std::complex<double>& resp1, double& coherence, const bool priorityOnFirst ) const{

	const Control* const ptrControl = Control::getInstance();
	calculateResponseFunctionByWLSAux( numSegments, out, in0, in1, weights, resp0, resp1, priorityOnFirst );
	double Gvar(0.0);
	double Gsyn(0.0);
	double Gres(0.0);
	double sumWeights(0.0);

#ifdef _USE_OMP
	int numThreads(1);
	#pragma omp parallel
	{
		numThreads = omp_get_num_threads();
	}
	double* GvarThread = new double[numThreads];
	double* GsynThread = new double[numThreads];
	double* GresThread = new double[numThreads];
	double* sumWeightsThread = new double[numThreads];
	for( int iThread = 0; iThread < numThreads; ++iThread ){
		// Zero clear
		GvarThread[iThread] = 0.0;
		GsynThread[iThread] = 0.0;
		GresThread[iThread] = 0.0;
		sumWeightsThread[iThread] = 0.0;
	}
	int i(0);
	int iThread(0);
	std::complex<double> outSyn(0.0, 0.0);
	#pragma omp parallel default(shared) private(i, iThread, outSyn)
	{
		iThread = omp_get_thread_num();
		#pragma omp for
		for( i = 0; i < numSegments; ++i ){
			GvarThread[iThread] += std::norm(out[i]) * weights[i];
			outSyn = resp0 * in0[i] + resp1 * in1[i];
			GsynThread[iThread] += std::norm(outSyn) * weights[i];
			residuals[i] = out[i] - outSyn;
			GresThread[iThread] += std::norm(residuals[i]) * weights[i];
			sumWeightsThread[iThread] += weights[i];
		}
	}
	for( int iThread = 0; iThread < numThreads; ++iThread ){
		Gvar += GvarThread[iThread];
		Gsyn += GsynThread[iThread];
		Gres += GresThread[iThread];
		sumWeights += sumWeightsThread[iThread];
	}
	delete [] GvarThread;
	delete [] GsynThread;
	delete [] GresThread;
	delete [] sumWeightsThread;
#ifdef _DEBUG_WRITE
	for( int i = 0; i < numSegments; ++i ){
		std::cout << "iSeg weights residuals: "
			<< i << " " << weights[i] << " " << residuals[i] << std::endl;
	}
#endif
#else
	for( int i = 0; i < numSegments; ++i ){
		Gvar += std::norm(out[i]) * weights[i];
		const std::complex<double> outSyn = resp0 * in0[i] + resp1 * in1[i];
		Gsyn += std::norm(outSyn) * weights[i];
		residuals[i] = out[i] - outSyn;
		Gres += std::norm(residuals[i]) * weights[i];
#ifdef _DEBUG_WRITE
		std::cout << "iSeg outSyn weights Gvar residuals Gres: "
			<< i << " " << outSyn << " " << weights[i] << " " <<Gvar << " " << residuals[i] << " " << Gres << std::endl;
#endif
		sumWeights += weights[i];
	}
#endif

	if( sumWeights < CommonParameters::EPS ){
		Gvar = 0.0;
		Gsyn = 0.0;
		Gres = 0.0;
		coherence = 0.0;
	}else{
		Gvar /= sumWeights;
		Gsyn /= sumWeights;
		Gres /= sumWeights;
		coherence = 1.0 - Gres / Gvar;
#ifdef _DEBUG_WRITE
		std::cout << "coherence : " << coherence << " " << 1.0 - Gres / Gvar << std::endl;
#endif
	}

	if (coherence > 1.0) {
		coherence = 1.0;
	}
	else if (coherence < 0.0) {
		coherence = 0.0;
	}

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeCvgMessage("	Sum of weights: " + Util::toString(sumWeights) );
	ptrOutputFiles->writeCvgMessage("	Squared coherence: " + Util::toString(coherence) );
	std::ostringstream msg1;
	msg1 << "	Estimated response function: ";
	msg1 << "(" << std::setw(12) << std::setprecision(4) << std::scientific << resp0.real() << "," 
				<< std::setw(12) << std::setprecision(4) << std::scientific << resp0.imag() << ")";
	msg1 << ", ("  << std::setw(12) << std::setprecision(4) << std::scientific << resp1.real() << "," 
				<< std::setw(12) << std::setprecision(4) << std::scientific << resp1.imag() << ")";
	ptrOutputFiles->writeCvgMessage(msg1.str());
	std::ostringstream msg2;
	msg2 << "	Amplitude of the estimated response function:";
	msg2 << std::setw(12) << std::setprecision(4) << std::scientific << std::abs(resp0) << "," 
		<< std::setw(12) << std::setprecision(4) << std::scientific << std::abs(resp1);
	ptrOutputFiles->writeCvgMessage(msg2.str());
	std::ostringstream msg3;
	msg3 << "	Phase(deg.) of the estimated response function:";
	msg3 << std::setw(7) << std::setprecision(1) << std::fixed << std::arg(resp0)*CommonParameters::RAD2DEG << "," 
		<< std::setw(7) << std::setprecision(1) << std::fixed << std::arg(resp1)*CommonParameters::RAD2DEG;
	ptrOutputFiles->writeCvgMessage(msg3.str());

	return Gres;

}

// Auxiliary function for calculating response function by the weighted leaset square method
void Analysis::calculateResponseFunctionByWLSAux( const int numSegments, const std::complex<double>* const out,
	const std::complex<double>* const in1, const std::complex<double>* const in2, const double* const weights, 
	std::complex<double>& resp1, std::complex<double>& resp2, const bool priorityOnFirst ) const{

	std::complex<double> out_in1(0.0, 0.0);
	std::complex<double> out_in2(0.0, 0.0);
	std::complex<double> in2_in1(0.0, 0.0);
	std::complex<double> in1_in1(0.0, 0.0);
	std::complex<double> in2_in2(0.0, 0.0);
#ifdef _USE_OMP
	int numThreads(1);
	#pragma omp parallel
	{
		numThreads = omp_get_num_threads();
	}
	std::complex<double>* out_in1Thread = new std::complex<double>[numThreads];
	std::complex<double>* out_in2Thread = new std::complex<double>[numThreads];
	std::complex<double>* in2_in1Thread = new std::complex<double>[numThreads];
	std::complex<double>* in1_in1Thread = new std::complex<double>[numThreads];
	std::complex<double>* in2_in2Thread = new std::complex<double>[numThreads];
	const std::complex<double> czero = std::complex<double>(0.0,0.0);
	for( int iThread = 0; iThread < numThreads; ++iThread ){
		// Zero clear
		out_in1Thread[iThread] = czero;
		out_in2Thread[iThread] = czero;
		in2_in1Thread[iThread] = czero;
		in1_in1Thread[iThread] = czero;
		in2_in2Thread[iThread] = czero;
	}
	int i(0);
	int iThread(0);
	#pragma omp parallel default(shared) private(i, iThread)
	{
		iThread = omp_get_thread_num();
		#pragma omp for
		for( i = 0; i < numSegments; ++i ){
			out_in1Thread[iThread] += out[i] * std::conj(in1[i]) * weights[i];
			out_in2Thread[iThread] += out[i] * std::conj(in2[i]) * weights[i];
			in2_in1Thread[iThread] += in2[i] * std::conj(in1[i]) * weights[i];
			in1_in1Thread[iThread] += in1[i] * std::conj(in1[i]) * weights[i];
			in2_in2Thread[iThread] += in2[i] * std::conj(in2[i]) * weights[i];
		}
	}
	for( int iThread = 0; iThread < numThreads; ++iThread ){
		out_in1 += out_in1Thread[iThread];
		out_in2 += out_in2Thread[iThread];
		in2_in1 += in2_in1Thread[iThread];
		in1_in1 += in1_in1Thread[iThread];
		in2_in2 += in2_in2Thread[iThread];
	}
	delete [] out_in1Thread;
	delete [] out_in2Thread;
	delete [] in2_in1Thread;
	delete [] in1_in1Thread;
	delete [] in2_in2Thread;
#else
	for( int i = 0; i < numSegments; ++i ){
		out_in1 += out[i] * std::conj(in1[i]) * weights[i];
		out_in2 += out[i] * std::conj(in2[i]) * weights[i];
		in2_in1 += in2[i] * std::conj(in1[i]) * weights[i];
		in1_in1 += in1[i] * std::conj(in1[i]) * weights[i];
		in2_in2 += in2[i] * std::conj(in2[i]) * weights[i];
	}
#endif

	const double coherence12 = std::norm(in2_in1)/Util::calculateAbsoluteValue(in1_in1)/Util::calculateAbsoluteValue(in2_in2);
	if( coherence12 > 1.0 - CommonParameters::EPS ){
		if( priorityOnFirst ){
			resp1 = out_in1 / in1_in1;
			resp2 = std::complex<double>(0.0,0.0);
		}else{
			resp1 = std::complex<double>(0.0,0.0);
			resp2 = out_in2 / in2_in2;
		}
		return;
	}

	const std::complex<double> in1_in2 = std::conj(in2_in1);
	const std::complex<double> det = in1_in1 * in2_in2 - in1_in2 * in2_in1;

	resp1 = ( out_in1 * in2_in2 - out_in2 * in2_in1 ) / det;
	resp2 = ( out_in2 * in1_in1 - out_in1 * in1_in2 ) / det;

}

// Calculate response function by the weighted leaset square method for bootstrap
void Analysis::calculateResponseFunctionByWLSForBootstrap( const int numSegments, const int* const segmentIndexes,
	const std::complex<double>* const out, const std::complex<double>* const in1, const std::complex<double>* const in2, 
	const double* const weights, std::complex<double>& resp1, std::complex<double>& resp2, const bool priorityOnFirst ) const{

	std::complex<double> out_in1(0.0, 0.0);
	std::complex<double> out_in2(0.0, 0.0);
	std::complex<double> in2_in1(0.0, 0.0);
	std::complex<double> in1_in1(0.0, 0.0);
	std::complex<double> in2_in2(0.0, 0.0);

#ifdef _USE_OMP
	int numThreads(1);
	#pragma omp parallel
	{
		numThreads = omp_get_num_threads();
	}
	std::complex<double>* out_in1Thread = new std::complex<double>[numThreads];
	std::complex<double>* out_in2Thread = new std::complex<double>[numThreads];
	std::complex<double>* in2_in1Thread = new std::complex<double>[numThreads];
	std::complex<double>* in1_in1Thread = new std::complex<double>[numThreads];
	std::complex<double>* in2_in2Thread = new std::complex<double>[numThreads];
	const std::complex<double> czero = std::complex<double>(0.0,0.0);
	for( int iThread = 0; iThread < numThreads; ++iThread ){
		// Zero clear
		out_in1Thread[iThread] = czero;
		out_in2Thread[iThread] = czero;
		in2_in1Thread[iThread] = czero;
		in1_in1Thread[iThread] = czero;
		in2_in2Thread[iThread] = czero;
	}
	int icount(0);
	int iThread(0);
	int iSeg(0);
	#pragma omp parallel default(shared) private(icount, iThread, iSeg)
	{
		iThread = omp_get_thread_num();
		#pragma omp for
		for( icount = 0; icount < numSegments; ++icount ){
			iSeg = segmentIndexes[icount];
			out_in1Thread[iThread] += out[iSeg] * std::conj(in1[iSeg]) * weights[iSeg];
			out_in2Thread[iThread] += out[iSeg] * std::conj(in2[iSeg]) * weights[iSeg];
			in2_in1Thread[iThread] += in2[iSeg] * std::conj(in1[iSeg]) * weights[iSeg];
			in1_in1Thread[iThread] += in1[iSeg] * std::conj(in1[iSeg]) * weights[iSeg];
			in2_in2Thread[iThread] += in2[iSeg] * std::conj(in2[iSeg]) * weights[iSeg];
		}
	}
	for( int iThread = 0; iThread < numThreads; ++iThread ){
		out_in1 += out_in1Thread[iThread];
		out_in2 += out_in2Thread[iThread];
		in2_in1 += in2_in1Thread[iThread];
		in1_in1 += in1_in1Thread[iThread];
		in2_in2 += in2_in2Thread[iThread];
	}
	delete [] out_in1Thread;
	delete [] out_in2Thread;
	delete [] in2_in1Thread;
	delete [] in1_in1Thread;
	delete [] in2_in2Thread;
#else
	for( int icount = 0; icount < numSegments; ++icount ){
		const int iSeg = segmentIndexes[icount];
		out_in1 += out[iSeg] * std::conj(in1[iSeg]) * weights[iSeg];
		out_in2 += out[iSeg] * std::conj(in2[iSeg]) * weights[iSeg];
		in2_in1 += in2[iSeg] * std::conj(in1[iSeg]) * weights[iSeg];
		in1_in1 += in1[iSeg] * std::conj(in1[iSeg]) * weights[iSeg];
		in2_in2 += in2[iSeg] * std::conj(in2[iSeg]) * weights[iSeg];
	}
#endif

	const double coherence12 = std::norm(in2_in1)/Util::calculateAbsoluteValue(in1_in1)/Util::calculateAbsoluteValue(in2_in2);
	if( coherence12 > 1.0 - CommonParameters::EPS ){
		if( priorityOnFirst ){
			resp1 = out_in1 / in1_in1;
			resp2 = std::complex<double>(0.0,0.0);
		}else{
			resp1 = std::complex<double>(0.0,0.0);
			resp2 = out_in2 / in2_in2;
		}
		return;
	}

	const std::complex<double> in1_in2 = std::conj(in2_in1);
	const std::complex<double> det = in1_in1 * in2_in2 - in1_in2 * in2_in1;

	resp1 = ( out_in1 * in2_in2 - out_in2 * in2_in1 ) / det;
	resp2 = ( out_in2 * in1_in1 - out_in1 * in1_in2 ) / det;

}

// Calculate response function by the weighted remote reference method
double Analysis::calculateResponseFunctionByWLSRemoteReference(const std::complex<double>* const out,
	const std::complex<double>* const in1, const std::complex<double>* const in2,
	const std::complex<double>* const rr1, const std::complex<double>* const rr2, 
	const int numSegments, const double* const weights, std::complex<double>* residuals,
	std::complex<double>& resp1, std::complex<double>& resp2, double& coherence) const {

	const Control* const ptrControl = Control::getInstance();
	calculateResponseFunctionByWLSRemoteReferenceAux(numSegments, out, in1, in2, rr1, rr2, weights, resp1, resp2);
	double Gvar(0.0);
	double Gsyn(0.0);
	double Gres(0.0);
	double sumWeights(0.0);

#ifdef _USE_OMP
	int numThreads(1);
#pragma omp parallel
	{
		numThreads = omp_get_num_threads();
	}
	double* GvarThread = new double[numThreads];
	double* GsynThread = new double[numThreads];
	double* GresThread = new double[numThreads];
	double* sumWeightsThread = new double[numThreads];
	for (int iThread = 0; iThread < numThreads; ++iThread) {
		// Zero clear
		GvarThread[iThread] = 0.0;
		GsynThread[iThread] = 0.0;
		GresThread[iThread] = 0.0;
		sumWeightsThread[iThread] = 0.0;
	}
	int i(0);
	int iThread(0);
	std::complex<double> outSyn(0.0, 0.0);
#pragma omp parallel default(shared) private(i, iThread, outSyn)
	{
		iThread = omp_get_thread_num();
#pragma omp for
		for (i = 0; i < numSegments; ++i) {
			GvarThread[iThread] += std::norm(out[i]) * weights[i];
			outSyn = resp1 * in1[i] + resp2 * in2[i];
			GsynThread[iThread] += std::norm(outSyn) * weights[i];
			residuals[i] = out[i] - outSyn;
			GresThread[iThread] += std::norm(residuals[i]) * weights[i];
			sumWeightsThread[iThread] += weights[i];
		}
	}
	for (int iThread = 0; iThread < numThreads; ++iThread) {
		Gvar += GvarThread[iThread];
		Gsyn += GsynThread[iThread];
		Gres += GresThread[iThread];
		sumWeights += sumWeightsThread[iThread];
	}
	delete[] GvarThread;
	delete[] GsynThread;
	delete[] GresThread;
	delete[] sumWeightsThread;
#ifdef _DEBUG_WRITE
	for (int i = 0; i < numSegments; ++i) {
		std::cout << "iSeg weights residuals: "
			<< i << " " << weights[i] << " " << residuals[i] << std::endl;
	}
#endif
#else
	for (int i = 0; i < numSegments; ++i) {
		Gvar += std::norm(out[i]) * weights[i];
		const std::complex<double> outSyn = resp1 * in1[i] + resp2 * in2[i];
		Gsyn += std::norm(outSyn) * weights[i];
		residuals[i] = out[i] - outSyn;
		Gres += std::norm(residuals[i]) * weights[i];
#ifdef _DEBUG_WRITE
		std::cout << "iSeg outSyn weights Gvar residuals Gres: "
			<< i << " " << outSyn << " " << weights[i] << " " << Gvar << " " << residuals[i] << " " << Gres << std::endl;
#endif
		sumWeights += weights[i];
	}
#endif

	if (sumWeights < CommonParameters::EPS) {
		Gvar = 0.0;
		Gsyn = 0.0;
		Gres = 0.0;
		coherence = 0.0;
	}
	else {
		Gvar /= sumWeights;
		Gsyn /= sumWeights;
		Gres /= sumWeights;
		coherence = 1.0 - Gres / Gvar;
#ifdef _DEBUG_WRITE
		std::cout << "coherence : " << coherence << " " << 1.0 - Gres / Gvar << std::endl;
#endif
	}

	if (coherence > 1.0) {
		coherence = 1.0;
	}
	else if (coherence < 0.0) {
		coherence = 0.0;
	}

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeCvgMessage("	Sum of weights: " + Util::toString(sumWeights));
	ptrOutputFiles->writeCvgMessage("	Squared coherence: " + Util::toString(coherence));
	std::ostringstream msg1;
	msg1 << "	Estimated response function: ";
	msg1 << "(" << std::setw(12) << std::setprecision(4) << std::scientific << resp1.real() << ","
		<< std::setw(12) << std::setprecision(4) << std::scientific << resp1.imag() << ")";
	msg1 << ", (" << std::setw(12) << std::setprecision(4) << std::scientific << resp2.real() << ","
		<< std::setw(12) << std::setprecision(4) << std::scientific << resp2.imag() << ")";
	ptrOutputFiles->writeCvgMessage(msg1.str());
	std::ostringstream msg2;
	msg2 << "	Amplitude of the estimated response function:";
	msg2 << std::setw(12) << std::setprecision(4) << std::scientific << std::abs(resp1) << ","
		<< std::setw(12) << std::setprecision(4) << std::scientific << std::abs(resp2);
	ptrOutputFiles->writeCvgMessage(msg2.str());
	std::ostringstream msg3;
	msg3 << "	Phase(deg.) of the estimated response function:";
	msg3 << std::setw(7) << std::setprecision(1) << std::fixed << std::arg(resp1) * CommonParameters::RAD2DEG << ","
		<< std::setw(7) << std::setprecision(1) << std::fixed << std::arg(resp2) * CommonParameters::RAD2DEG;
	ptrOutputFiles->writeCvgMessage(msg3.str());

	return Gres;
}

// Auxiliary function for calculating response function by the weighted remote reference method
void Analysis::calculateResponseFunctionByWLSRemoteReferenceAux(const int numSegments, const std::complex<double>* const out,
	const std::complex<double>* const in1, const std::complex<double>* const in2,
	const std::complex<double>* const rr1, const std::complex<double>* const rr2,
	const double* const weights, std::complex<double>& resp1, std::complex<double>& resp2) const {

	std::complex<double> out_rr1(0.0, 0.0);
	std::complex<double> out_rr2(0.0, 0.0);
	std::complex<double> in2_rr1(0.0, 0.0);
	std::complex<double> in1_rr1(0.0, 0.0);
	std::complex<double> in2_rr2(0.0, 0.0);
	std::complex<double> in1_rr2(0.0, 0.0);

#ifdef _USE_OMP
	int numThreads(1);
#pragma omp parallel
	{
		numThreads = omp_get_num_threads();
	}
	std::complex<double>* out_rr1Thread = new std::complex<double>[numThreads];
	std::complex<double>* out_rr2Thread = new std::complex<double>[numThreads];
	std::complex<double>* in2_rr1Thread = new std::complex<double>[numThreads];
	std::complex<double>* in1_rr1Thread = new std::complex<double>[numThreads];
	std::complex<double>* in2_rr2Thread = new std::complex<double>[numThreads];
	std::complex<double>* in1_rr2Thread = new std::complex<double>[numThreads];

	const std::complex<double> czero = std::complex<double>(0.0, 0.0);
	for (int iThread = 0; iThread < numThreads; ++iThread) {
		// Zero clear
		out_rr1Thread[iThread] = czero;
		out_rr2Thread[iThread] = czero;
		in2_rr1Thread[iThread] = czero;
		in1_rr1Thread[iThread] = czero;
		in2_rr2Thread[iThread] = czero;
		in1_rr2Thread[iThread] = czero;
	}
	int i(0);
	int iThread(0);
#pragma omp parallel default(shared) private(i, iThread)
	{
		iThread = omp_get_thread_num();
#pragma omp for
		for (i = 0; i < numSegments; ++i) {
			out_rr1Thread[iThread] += out[i] * std::conj(rr1[i]) * weights[i];
			out_rr2Thread[iThread] += out[i] * std::conj(rr2[i]) * weights[i];
			in2_rr1Thread[iThread] += in2[i] * std::conj(rr1[i]) * weights[i];
			in1_rr1Thread[iThread] += in1[i] * std::conj(rr1[i]) * weights[i];
			in2_rr2Thread[iThread] += in2[i] * std::conj(rr2[i]) * weights[i];
			in1_rr2Thread[iThread] += in1[i] * std::conj(rr2[i]) * weights[i];
		}
	}
	for (int iThread = 0; iThread < numThreads; ++iThread) {
		out_rr1 += out_rr1Thread[iThread];
		out_rr2 += out_rr2Thread[iThread];
		in2_rr1 += in2_rr1Thread[iThread];
		in1_rr1 += in1_rr1Thread[iThread];
		in2_rr2 += in2_rr2Thread[iThread];
		in1_rr2 += in1_rr2Thread[iThread];
	}
	delete[] out_rr1Thread;
	delete[] out_rr2Thread;
	delete[] in2_rr1Thread;
	delete[] in1_rr1Thread;
	delete[] in2_rr2Thread;
	delete[] in1_rr2Thread;
#else
	for (int i = 0; i < numSegments; ++i) {
		out_rr1 += out[i] * std::conj(rr1[i]) * weights[i];
		out_rr2 += out[i] * std::conj(rr2[i]) * weights[i];
		in2_rr1 += in2[i] * std::conj(rr1[i]) * weights[i];
		in1_rr1 += in1[i] * std::conj(rr1[i]) * weights[i];
		in2_rr2 += in2[i] * std::conj(rr2[i]) * weights[i];
		in1_rr2 += in1[i] * std::conj(rr2[i]) * weights[i];
	}
#endif

	const std::complex<double> det = in1_rr1 * in2_rr2 - in1_rr2 * in2_rr1;

	resp1 = (out_rr1 * in2_rr2 - out_rr2 * in2_rr1) / det;
	resp2 = (out_rr2 * in1_rr1 - out_rr1 * in1_rr2) / det;

	//double Gvar(0.0);
	//double Gsyn(0.0);
	//for (int i = 0; i < numSegments; ++i) {
	//	Gvar += std::norm(out[i]) * weights[i];
	//	const std::complex<double> outSyn = resp1 * in1[i] + resp2 * in2[i];
	//	Gsyn += std::norm(outSyn) * weights[i];
	//}

	//double coherence = Gsyn / Gvar;
	//if (coherence > 1.0) {
	//	coherence = 1.0;
	//}
	//else if (coherence < 0.0) {
	//	coherence = 0.0;
	//}
	//return coherence;

}

// Calculate response function by the weighted remote reference method for bootstrap
void Analysis::calculateResponseFunctionByWLSRemoteReferenceForBootstrap(const int numSegments, const int* const segmentIndexes,
	const std::complex<double>* const out,
	const std::complex<double>* const in1, const std::complex<double>* const in2,
	const std::complex<double>* const rr1, const std::complex<double>* const rr2,
	const double* const weights, std::complex<double>& resp1, std::complex<double>& resp2) const {


	std::complex<double> out_rr1(0.0, 0.0);
	std::complex<double> out_rr2(0.0, 0.0);
	std::complex<double> in2_rr1(0.0, 0.0);
	std::complex<double> in1_rr1(0.0, 0.0);
	std::complex<double> in2_rr2(0.0, 0.0);
	std::complex<double> in1_rr2(0.0, 0.0);

#ifdef _USE_OMP
	int numThreads(1);
#pragma omp parallel
	{
		numThreads = omp_get_num_threads();
	}
	std::complex<double>* out_rr1Thread = new std::complex<double>[numThreads];
	std::complex<double>* out_rr2Thread = new std::complex<double>[numThreads];
	std::complex<double>* in2_rr1Thread = new std::complex<double>[numThreads];
	std::complex<double>* in1_rr1Thread = new std::complex<double>[numThreads];
	std::complex<double>* in2_rr2Thread = new std::complex<double>[numThreads];
	std::complex<double>* in1_rr2Thread = new std::complex<double>[numThreads];

	const std::complex<double> czero = std::complex<double>(0.0, 0.0);
	for (int iThread = 0; iThread < numThreads; ++iThread) {
		// Zero clear
		out_rr1Thread[iThread] = czero;
		out_rr2Thread[iThread] = czero;
		in2_rr1Thread[iThread] = czero;
		in1_rr1Thread[iThread] = czero;
		in2_rr2Thread[iThread] = czero;
		in1_rr2Thread[iThread] = czero;
	}
	int icount(0);
	int iThread(0);
	int iSeg(0);
#pragma omp parallel default(shared) private(icount, iThread, iSeg)
	{
		iThread = omp_get_thread_num();
#pragma omp for
		for (icount = 0; icount < numSegments; ++icount) {
			iSeg = segmentIndexes[icount];
			out_rr1Thread[iThread] += out[iSeg] * std::conj(rr1[iSeg]) * weights[iSeg];
			out_rr2Thread[iThread] += out[iSeg] * std::conj(rr2[iSeg]) * weights[iSeg];
			in2_rr1Thread[iThread] += in2[iSeg] * std::conj(rr1[iSeg]) * weights[iSeg];
			in1_rr1Thread[iThread] += in1[iSeg] * std::conj(rr1[iSeg]) * weights[iSeg];
			in2_rr2Thread[iThread] += in2[iSeg] * std::conj(rr2[iSeg]) * weights[iSeg];
			in1_rr2Thread[iThread] += in1[iSeg] * std::conj(rr2[iSeg]) * weights[iSeg];
		}
	}
	for (int iThread = 0; iThread < numThreads; ++iThread) {
		out_rr1 += out_rr1Thread[iThread];
		out_rr2 += out_rr2Thread[iThread];
		in2_rr1 += in2_rr1Thread[iThread];
		in1_rr1 += in1_rr1Thread[iThread];
		in2_rr2 += in2_rr2Thread[iThread];
		in1_rr2 += in1_rr2Thread[iThread];
	}
	delete[] out_rr1Thread;
	delete[] out_rr2Thread;
	delete[] in2_rr1Thread;
	delete[] in1_rr1Thread;
	delete[] in2_rr2Thread;
	delete[] in1_rr2Thread;
#else
	for (int icount = 0; icount < numSegments; ++icount) {
		const int iSeg = segmentIndexes[icount];
		out_rr1 += out[iSeg] * std::conj(rr1[iSeg]) * weights[iSeg];
		out_rr2 += out[iSeg] * std::conj(rr2[iSeg]) * weights[iSeg];
		in2_rr1 += in2[iSeg] * std::conj(rr1[iSeg]) * weights[iSeg];
		in1_rr1 += in1[iSeg] * std::conj(rr1[iSeg]) * weights[iSeg];
		in2_rr2 += in2[iSeg] * std::conj(rr2[iSeg]) * weights[iSeg];
		in1_rr2 += in1[iSeg] * std::conj(rr2[iSeg]) * weights[iSeg];
	}
#endif

	const std::complex<double> det = in1_rr1 * in2_rr2 - in1_rr2 * in2_rr1;

	resp1 = (out_rr1 * in2_rr2 - out_rr2 * in2_rr1) / det;
	resp2 = (out_rr2 * in1_rr1 - out_rr1 * in1_rr2) / det;

}

// Get pointer to M-estimators
RobustWeight* Analysis::getPointerToRobustWeight(const int iRobustWeight) const {
		return m_robustWeight[iRobustWeight];
}

// Output spectral density functions to cvg file
void Analysis::outputSpectralDensityFunctionsToCvgFile( const int numSegments, const double timeLength, const std::complex<double>* const out,
	const std::complex<double>* const in1, const std::complex<double>* const in2, const double* const weights ) const{

	std::complex<double> out_in1(0.0, 0.0);
	std::complex<double> out_in2(0.0, 0.0);
	std::complex<double> in2_in1(0.0, 0.0);
	std::complex<double> in1_in1(0.0, 0.0);
	std::complex<double> in2_in2(0.0, 0.0);
	std::complex<double> out_out(0.0, 0.0);
	double sumOfWeights(0.0);
	for( int i = 0; i < numSegments; ++i ){
		out_in1 += out[i] * std::conj(in1[i]) * weights[i];
		out_in2 += out[i] * std::conj(in2[i]) * weights[i];
		in2_in1 += in2[i] * std::conj(in1[i]) * weights[i];
		in1_in1 += in1[i] * std::conj(in1[i]) * weights[i];
		in2_in2 += in2[i] * std::conj(in2[i]) * weights[i];
		out_out += out[i] * std::conj(out[i]) * weights[i];
		sumOfWeights += weights[i];
	}

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	if( sumOfWeights < CommonParameters::EPS ){
		ptrOutputFiles->writeCvgMessage("	Sum of weights (" +  Util::toString(sumOfWeights) + ") is too small");
		ptrOutputFiles->writeWarningMessage("Sum of weights (" +  Util::toString(sumOfWeights) + ") is too small");
	}

	const double factor = 2.0 / timeLength / sumOfWeights;
	out_in1 *= factor;
	out_in2 *= factor;
	in2_in1 *= factor;
	in1_in1 *= factor;
	in2_in2 *= factor;
	out_out *= factor;

	ptrOutputFiles->writeCvgMessage("Spectral density function <Out *Inp1>: " + Util::toString(out_in1));
	ptrOutputFiles->writeCvgMessage("Spectral density function <Out *Inp2>: " + Util::toString(out_in2));
	ptrOutputFiles->writeCvgMessage("Spectral density function <Inp2*Inp1>: " + Util::toString(in2_in1));
	ptrOutputFiles->writeCvgMessage("Spectral density function <Inp1*Inp1>: " + Util::toString(in1_in1));
	ptrOutputFiles->writeCvgMessage("Spectral density function <Inp2*Inp2>: " + Util::toString(in2_in2));
	ptrOutputFiles->writeCvgMessage("Spectral density function <Out *Out >: " + Util::toString(out_out));

}

// Output spectral density functions to cvg file
void Analysis::outputSpectralDensityFunctionsToCvgFile( const int numSegments, const double timeLength, const std::complex<double>* const out,
	const std::complex<double>* const in1, const std::complex<double>* const in2, 
	const std::complex<double>* const rr1, const std::complex<double>* const rr2 ) const{

	std::complex<double> out_rr1(0.0, 0.0);
	std::complex<double> out_rr2(0.0, 0.0);
	std::complex<double> in1_rr1(0.0, 0.0);
	std::complex<double> in1_rr2(0.0, 0.0);
	std::complex<double> in2_rr1(0.0, 0.0);
	std::complex<double> in2_rr2(0.0, 0.0);
	double sumOfWeights(0.0);
	for( int i = 0; i < numSegments; ++i ){
		out_rr1 += out[i] * std::conj(rr1[i]);
		out_rr2 += out[i] * std::conj(rr2[i]);
		in1_rr1 += in1[i] * std::conj(rr1[i]);
		in1_rr2 += in1[i] * std::conj(rr2[i]);
		in2_rr1 += in2[i] * std::conj(rr1[i]);
		in2_rr2 += in2[i] * std::conj(rr2[i]);
		sumOfWeights += 1.0;
	}

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	if( sumOfWeights < CommonParameters::EPS ){
		ptrOutputFiles->writeCvgMessage("	Sum of weights (" +  Util::toString(sumOfWeights) + ") is too small");
		ptrOutputFiles->writeWarningMessage("Sum of weights (" +  Util::toString(sumOfWeights) + ") is too small");
	}

	const double factor = 2.0 / timeLength / sumOfWeights;
	out_rr1 *= factor;
	out_rr2 *= factor;
	in1_rr1 *= factor;
	in1_rr2 *= factor;
	in2_rr1 *= factor;
	in2_rr2 *= factor;

	ptrOutputFiles->writeCvgMessage("Spectral density function <Out *RR1>: " + Util::toString(out_rr1));
	ptrOutputFiles->writeCvgMessage("Spectral density function <Out *RR2>: " + Util::toString(out_rr2));
	ptrOutputFiles->writeCvgMessage("Spectral density function <Inp1*RR1>: " + Util::toString(in1_rr1));
	ptrOutputFiles->writeCvgMessage("Spectral density function <Inp1*RR2>: " + Util::toString(in1_rr2));
	ptrOutputFiles->writeCvgMessage("Spectral density function <Inp2*RR1>: " + Util::toString(in2_rr1));
	ptrOutputFiles->writeCvgMessage("Spectral density function <Inp2*RR2>: " + Util::toString(in2_rr2));
}

// Merge sections
void Analysis::mergeSections( std::vector<CommonParameters::DataFileSet>& dataFileSets ) const{

	const Control* const ptrControl = Control::getInstance();
	std::vector<CommonParameters::DataFileSet> dataFileSetsAfterMerge;
	const int numSectionsAfterMerge = ptrControl->getNumRangeOfSectionsForMerge();
	if( numSectionsAfterMerge < 1 ){
		return;
	}
	dataFileSetsAfterMerge.reserve(numSectionsAfterMerge);

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeLogMessage("Merge sections");

	for( int iSectionAfterMerge = 0; iSectionAfterMerge < numSectionsAfterMerge; ++iSectionAfterMerge ){
		const std::pair<int, int> rangeOfSections = ptrControl->getRangeOfSectionsForMerge(iSectionAfterMerge);
		int numDataAfterMerge(0);
		int iSection(0);
		for( std::vector<CommonParameters::DataFileSet>::const_iterator itrOrg = dataFileSets.begin(); itrOrg != dataFileSets.end(); ++itrOrg, ++iSection ){
			if( iSection >= rangeOfSections.first && iSection <= rangeOfSections.second ){
				const int numDataOrg = itrOrg->numDataPoints;
				numDataAfterMerge += numDataOrg;
			}
		}
		std::ostringstream oss;
		oss << "Merged section " << iSectionAfterMerge << " : number of data after merging is " << numDataAfterMerge;
		ptrOutputFiles->writeLogMessage(oss.str());
		CommonParameters::DataFileSet dataFileSetNew;
		dataFileSetNew.numDataPoints = numDataAfterMerge;
		const int numChannels = ptrControl->getNumberOfChannels();
		dataFileSetNew.dataFile.reserve(numChannels);
		for( int iChan = 0; iChan < numChannels; ++iChan ){
			CommonParameters::DataFile dataFileNew;
			dataFileNew.fileName = "merged_section_" + Util::toString(iSectionAfterMerge) + "_channel_" + Util::toString(iChan);
			dataFileNew.numSkipData = 0;
			dataFileNew.data = new double [numDataAfterMerge];
			const std::pair<int, int> rangeOfSections = ptrControl->getRangeOfSectionsForMerge(iSectionAfterMerge);
			int iSection(0);
			int dataCounter(0);
			for( std::vector<CommonParameters::DataFileSet>::const_iterator itrOrg = dataFileSets.begin(); itrOrg != dataFileSets.end(); ++itrOrg, ++iSection ){
				if( iSection >= rangeOfSections.first && iSection <= rangeOfSections.second ){
					const int numDataOrg = itrOrg->numDataPoints;
					for( int i = 0; i < numDataOrg; ++i, ++dataCounter ){
						dataFileNew.data[dataCounter] = itrOrg->dataFile[iChan].data[i];
					}
				}
			}
			assert( dataCounter == numDataAfterMerge );
			dataFileSetNew.dataFile.push_back(dataFileNew);
		}
		dataFileSetsAfterMerge.push_back(dataFileSetNew);
	}

	// Delete allocated memory for time-series data
	for( std::vector<CommonParameters::DataFileSet>::const_iterator itrDataFileSet = dataFileSets.begin(); itrDataFileSet != dataFileSets.end(); ++itrDataFileSet ){
		for( std::vector<CommonParameters::DataFile>::const_iterator itrDataFile = itrDataFileSet->dataFile.begin(); itrDataFile != itrDataFileSet->dataFile.end(); ++itrDataFile ){
			delete [] itrDataFile->data;
		}
	}

	dataFileSets.swap(dataFileSetsAfterMerge);

}
// Calculate rotated fields
void Analysis::calculateRotatedFields( const int numSegmentsTotal, std::complex<double>** ftval ) const{

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeLogMessage("Calculate rotated fields for segment length");

	const Control* const ptrControl = Control::getInstance();
	const double rotationAngle = ptrControl->getRotationAngle();

	// Output variables
	const int numOutputVariables = ptrControl->getNumOutputVariables();		
	for( int i = 0; i < numOutputVariables / 2; ++i ){
		const int iChan0 = ptrControl->getChannelIndex(CommonParameters::OUTPUT, i);
		const int iChan1 = ptrControl->getChannelIndex(CommonParameters::OUTPUT, i+1);
		const double azimuth0 = ptrControl->getAzimuth(iChan0);
		const double azimuth1 = ptrControl->getAzimuth(iChan1);
		if( fabs(azimuth0 - azimuth1) < CommonParameters::EPS ){
			ptrOutputFiles->writeErrorMessage("Measured directions of channels " + Util::toString(iChan0) +
				" and " + Util::toString(iChan1) + " are the same (" +  Util::toString(azimuth0) + ")");
		}
		for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
			const std::complex<double> v0 = ftval[iChan0][iSeg];
			const std::complex<double> v1 = ftval[iChan1][iSeg];
			ftval[iChan0][iSeg] = Util::calculateRotatedField( azimuth0, azimuth1, rotationAngle, v0, v1 );
			ftval[iChan1][iSeg] = Util::calculateRotatedField( azimuth0, azimuth1, rotationAngle + 90, v0, v1 );
		}
	}
	// Input variables
	const int numInputVariables = ptrControl->getNumInputVariables();
	for( int i = 0; i < numInputVariables / 2; ++i ){
		const int iChan0 = ptrControl->getChannelIndex(CommonParameters::INPUT, i);
		const int iChan1 = ptrControl->getChannelIndex(CommonParameters::INPUT, i+1);
		const double azimuth0 = ptrControl->getAzimuth(iChan0);
		const double azimuth1 = ptrControl->getAzimuth(iChan1);
		if( fabs(azimuth0 - azimuth1) < CommonParameters::EPS ){
			ptrOutputFiles->writeErrorMessage("Measured directions of channels " + Util::toString(iChan0) +
				" and " + Util::toString(iChan1) + " are the same (" +  Util::toString(azimuth0) + ")");
		}
		for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
			const std::complex<double> v0 = ftval[iChan0][iSeg];
			const std::complex<double> v1 = ftval[iChan1][iSeg];
			ftval[iChan0][iSeg] = Util::calculateRotatedField( azimuth0, azimuth1, rotationAngle, v0, v1 );
			ftval[iChan1][iSeg] = Util::calculateRotatedField( azimuth0, azimuth1, rotationAngle + 90, v0, v1 );
		}
	}
	// Remote reference variables
	const int numRemoteReferenceVariables = ptrControl->getNumRemoteReferenceVariables();
	for( int i = 0; i < numRemoteReferenceVariables / 2; ++i ){
		const int iChan0 = ptrControl->getChannelIndex(CommonParameters::REMOTE_REFERENCE, i);
		const int iChan1 = ptrControl->getChannelIndex(CommonParameters::REMOTE_REFERENCE, i+1);
		const double azimuth0 = ptrControl->getAzimuth(iChan0);
		const double azimuth1 = ptrControl->getAzimuth(iChan1);
		if( fabs(azimuth0 - azimuth1) < CommonParameters::EPS ){
			ptrOutputFiles->writeErrorMessage("Measured directions of channels " + Util::toString(iChan0) +
				" and " + Util::toString(iChan1) + " are the same (" +  Util::toString(azimuth0) + ")");
		}
		for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
			const std::complex<double> v0 = ftval[iChan0][iSeg];
			const std::complex<double> v1 = ftval[iChan1][iSeg];
			ftval[iChan0][iSeg] = Util::calculateRotatedField( azimuth0, azimuth1, rotationAngle, v0, v1 );
			ftval[iChan1][iSeg] = Util::calculateRotatedField( azimuth0, azimuth1, rotationAngle + 90, v0, v1 );
		}
	}

}

// Output average spectrum
void Analysis::outputAverageSpectrum( std::complex<double>** cdata, const int numData, const int section, 
	const bool afterPreprocessing, const bool afterCalibration ) const{

	const Control* const ptrControl = Control::getInstance();
	const int numChannels = ptrControl->getNumberOfChannels();

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();

	std::ostringstream fileNameForAverage;
	if(afterPreprocessing){
		if(afterCalibration){
			fileNameForAverage << "stat_freq_average_after_preprocessing_sect" << section <<".csv";
		}else{
			fileNameForAverage << "stat_freq_average_after_preprocessing_before_calibration_sect" << section <<".csv";
		}
	}else{
		if(afterCalibration){
			fileNameForAverage << "stat_freq_average_sect" << section <<".csv";
		}else{
			fileNameForAverage << "stat_freq_average_before_calibration_sect" << section <<".csv";
		}
	}
	std::ofstream ofsAvg;
	ofsAvg.open( fileNameForAverage.str().c_str(), std::ios::out );
	if( ofsAvg.fail() ){
		ptrOutputFiles->writeErrorMessage( "File open error : " + fileNameForAverage.str() );
	}
	// Write header
	ofsAvg << "frequency";
	for( int iChan = 0; iChan < numChannels; ++iChan ){
		ofsAvg << ",auto_power_" << iChan;
	}
	ofsAvg << std::endl;

	const double samplingFrequency = ptrControl->getSamplingFrequency();
	const double T =  static_cast<double>(numData) / samplingFrequency;
	const Control::ParamsForFreqDomainEvaluation params = ptrControl->getParamsForFreqDomainEvaluation();
	const int numSegmentLengths = ptrControl->getNumSegmentLengths();

	const double logFreqStartAvg = log10(params.startOutputFrequency);
	const double logFreqEndAvg = log10(params.endOutputFrequency);
	const double logFreqIncAvg = ( logFreqEndAvg - logFreqStartAvg ) / static_cast<double>(params.numOfOutputFrequencyForAverage);
	int outputFreqIndexPre = -1;
	int counterOfFreqs(0);
	double logFreqAverage = 0.0;// Zero clear
	std::vector<double> logPowerAverage;
	logPowerAverage.reserve(numChannels);
	for( int iChan = 0; iChan < numChannels; ++iChan ){
		logPowerAverage.push_back(0.0);// Zero clear
	}
	for( int i = 1; i < numData / 2 + 1; ++i ){
		const double freq = static_cast<double>(i) / static_cast<double>(numData) * samplingFrequency;
		if( freq < params.startOutputFrequency || freq > params.endOutputFrequency ){
			continue;
		}
		const double logFreq = log10(freq);
		logFreqAverage += logFreq;
		const double factor = ( i == 0 || i == numData / 2 ) ? 1.0 : 2.0;
		for( int iChan = 0; iChan < numChannels; ++iChan ){
			const double power = factor * std::norm(cdata[iChan][i]) * T;
			logPowerAverage[iChan] += log10(power);
		}
		++counterOfFreqs;
		const int outputFreqIndexCur = static_cast<int>( std::floor( ( logFreq - logFreqStartAvg ) / logFreqIncAvg ) );
		if( outputFreqIndexCur <= outputFreqIndexPre ){
			continue;
		}
		logFreqAverage /= static_cast<double>(counterOfFreqs);
		ofsAvg << std::setw(10)  << std::scientific << pow(10,logFreqAverage);
		logFreqAverage = 0.0;
		for( int iChan = 0; iChan < numChannels; ++iChan ){
			logPowerAverage[iChan] /= static_cast<double>(counterOfFreqs);
			ofsAvg << ", ";
			ofsAvg << std::setw(10)  << std::scientific << pow(10,logPowerAverage[iChan]);
			logPowerAverage[iChan] = 0.0;
		}
		ofsAvg << std::endl;
		counterOfFreqs = 0;
		outputFreqIndexPre = outputFreqIndexCur;
	}
	ofsAvg.close();

}

// Output average spectrum for the frequencies where response functions are estimated
// @note This function has not yet been debugged
void Analysis::outputAverageSpectrum2( std::complex<double>** cdata, const int numData, const int section, 
	const bool afterPreprocessing, const bool afterCalibration ) const{

	const Control* const ptrControl = Control::getInstance();
	const int numChannels = ptrControl->getNumberOfChannels();

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();

	std::ostringstream fileNameForAverage;
	if(afterPreprocessing){
		if(afterCalibration){
			fileNameForAverage << "stat_freq_average_after_preprocessing_sect" << section <<".csv";
		}else{
			fileNameForAverage << "stat_freq_average_after_preprocessing_before_calibration_sect" << section <<".csv";
		}
	}else{
		if(afterCalibration){
			fileNameForAverage << "stat_freq_average_sect" << section <<".csv";
		}else{
			fileNameForAverage << "stat_freq_average_before_calibration_sect" << section <<".csv";
		}
	}
	std::ofstream ofsAvg;
	ofsAvg.open( fileNameForAverage.str().c_str(), std::ios::out );
	if( ofsAvg.fail() ){
		ptrOutputFiles->writeErrorMessage( "File open error : " + fileNameForAverage.str() );
	}
	// Write header
	ofsAvg << "frequency";
	for( int iChan = 0; iChan < numChannels; ++iChan ){
		ofsAvg << ",auto_power_" << iChan;
	}
	ofsAvg << std::endl;

	const double samplingFrequency = ptrControl->getSamplingFrequency();
	const double T =  static_cast<double>(numData) / samplingFrequency;
	const Control::ParamsForFreqDomainEvaluation params = ptrControl->getParamsForFreqDomainEvaluation();
	const int numSegmentLengths = ptrControl->getNumSegmentLengths();
	std::vector<double> targetFreqs;
	for( int iSegLen = 0; iSegLen < numSegmentLengths; ++iSegLen ){
		const int segmentLength = ptrControl->getSegmentLength(iSegLen);
		const int numTargetFreqsInSegment = ptrControl->getNumTargetFrequencyInSegment(iSegLen);
		for( int iFreq = 0; iFreq < numTargetFreqsInSegment; ++iFreq ){
			ptrOutputFiles->writeLogMessage("-------------------------------------------------------------------------------",false);
			const int freqDegree = ptrControl->getTargetFrequencyDegreesInSegment(iSegLen, iFreq);
			const double freqTarget = samplingFrequency * static_cast<double>(freqDegree) / static_cast<double>(segmentLength);
			targetFreqs.push_back(freqTarget);
		}
	}

	const int numOfTargetFreqs = static_cast<int>(targetFreqs.size());
	int* ids = new int[numData / 2];
	double* diffFreqs = new double[numData / 2];
	for( int iTargetFreqs = 0; iTargetFreqs < numOfTargetFreqs; ++iTargetFreqs ){
		int counter(0);
		for( int i = 1; i < numData / 2 + 1; ++i ){
			const double freq = static_cast<double>(i) / static_cast<double>(numData) * samplingFrequency;
			//if( iTargetFreqs + 1 < numOfTargetFreqs && freq >= targetFreqs[iTargetFreqs + 1] ){
			//	// Exclude the frequencies larger than or equal to the next target frequency
			//	break;
			//}
			//if( iTargetFreqs > 0  && freq <= targetFreqs[iTargetFreqs - 1] ){
			//	// Exclude the frequencies lest than or equal to the previous target frequency
			//	break;
			//}
			ids[i] = i;
			diffFreqs[i] = fabs(freq - targetFreqs[iTargetFreqs]);
			++counter;
		}
		Util::quickSort(counter, ids, diffFreqs);
		int numForAverage = params.numOfOutputFrequencyForAverage;
		if( counter < numForAverage ){
			numForAverage = counter;
		}
		std::vector<double> logPowerAverage;
		logPowerAverage.reserve(numChannels);
		for( int iChan = 0; iChan < numChannels; ++iChan ){
			logPowerAverage.push_back(0.0);// Zero clear
		}
		for( int iSample = 0; iSample < numForAverage; ++iSample ){
			for( int iChan = 0; iChan < numChannels; ++iChan ){
				const int index = ids[iSample];
				const double power = 2.0 * std::norm(cdata[iChan][index]) * T;
				logPowerAverage[iChan] += log10(power);
			}
		}
		if( numForAverage > 1 ){
			ofsAvg << targetFreqs[iTargetFreqs] << std::endl;
			for( int iChan = 0; iChan < numChannels; ++iChan ){
				logPowerAverage[iChan] /= static_cast<double>(numForAverage);
				ofsAvg << ", ";
				ofsAvg << std::setw(10)  << std::scientific << pow(10,logPowerAverage[iChan]);
			}
			ofsAvg << std::endl;
		}
	}
	delete [] ids;
	delete [] diffFreqs;

}

// Output spectrum
void Analysis::outputSpectrum( std::complex<double>** cdata, const int numData, const int section, 
	const bool afterPreprocessing, const bool afterCalibration ) const{

	const Control* const ptrControl = Control::getInstance();
	const int numChannels = ptrControl->getNumberOfChannels();
	const int numOutputVariables = ptrControl->getNumOutputVariables();		
	const int numInputVariables = ptrControl->getNumInputVariables();		
	const int numRemoteReferenceVariables = ptrControl->getNumRemoteReferenceVariables();
	const double samplingFrequency = ptrControl->getSamplingFrequency();
	const double T = static_cast<double>(numData) / samplingFrequency;
	const Control::ParamsForFreqDomainEvaluation params = ptrControl->getParamsForFreqDomainEvaluation();

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();

	std::ostringstream fileName;
	if(afterPreprocessing){
		if(afterCalibration){
			fileName << "stat_freq_after_preprocessing_sect" << section <<".csv";
		}else{
			fileName << "stat_freq_after_preprocessing_before_calibration_sect" << section <<".csv";
		}
	}else{
		if(afterCalibration){
			fileName << "stat_freq_sect" << section <<".csv";
		}else{
			fileName << "stat_freq_before_calibration_sect" << section <<".csv";
		}
	}
	std::ofstream ofs;
	ofs.open( fileName.str().c_str(), std::ios::out );
	if( ofs.fail() ){
		ptrOutputFiles->writeErrorMessage( "File open error : " + fileName.str() );
	}
	// Write header
	ofs << "index,frequency";
	for( int iChan = 0; iChan < numChannels; ++iChan ){
		ofs << ",amplitude_" << iChan;
		ofs << ",phase_" << iChan;
		ofs << ",power_" << iChan;
	}
	for( int iOut = 0; iOut < numOutputVariables; ++iOut ){
		const int iChan0 = ptrControl->getChannelIndex(CommonParameters::OUTPUT, iOut);
		for( int iRR = 0; iRR < numRemoteReferenceVariables; ++iRR ){
			const int iChan1 = ptrControl->getChannelIndex(CommonParameters::REMOTE_REFERENCE, iRR);
			ofs << ",cross_power_gain_" << iChan0 << "_" << iChan1;
			ofs << ",cross_power_phase_" << iChan0 << "_" << iChan1;
		}
	}
	for( int iInp = 0; iInp < numInputVariables; ++iInp ){
		const int iChan0 = ptrControl->getChannelIndex(CommonParameters::INPUT, iInp);
		for( int iRR = 0; iRR < numRemoteReferenceVariables; ++iRR ){
			const int iChan1 = ptrControl->getChannelIndex(CommonParameters::REMOTE_REFERENCE, iRR);
			ofs << ",cross_power_gain_" << iChan0 << "_" << iChan1;
			ofs << ",cross_power_phase_" << iChan0 << "_" << iChan1;
		}
	}
	for( int iOut = 0; iOut < numOutputVariables; ++iOut ){
		const int iChan0 = ptrControl->getChannelIndex(CommonParameters::OUTPUT, iOut);
		for( int iInp = 0; iInp < numInputVariables; ++iInp ){
			const int iChan1 = ptrControl->getChannelIndex(CommonParameters::INPUT, iInp);
			ofs << ",cross_power_gain_" << iChan0 << "_" << iChan1;
			ofs << ",cross_power_phase_" << iChan0 << "_" << iChan1;
		}
	}
	ofs << std::endl;

	const double logFreqStart = log10(params.startOutputFrequency);
	const double logFreqEnd = log10(params.endOutputFrequency);
	const double logFreqInc = ( logFreqEnd - logFreqStart ) / static_cast<double>(params.numOfOutputFrequency);
	int outputFreqIndexPre = -1;
	for( int i = 1; i < numData / 2 + 1; ++i ){
		const double freq = static_cast<double>(i) / static_cast<double>(numData) * samplingFrequency;
		if( freq < params.startOutputFrequency || freq > params.endOutputFrequency ){
			continue;
		}
		const double logFreq = log10(freq);
		const int outputFreqIndexCur = static_cast<int>( std::floor( ( logFreq - logFreqStart ) / logFreqInc ) );
		if( outputFreqIndexCur <= outputFreqIndexPre ){
			continue;
		}
		outputFreqIndexPre = outputFreqIndexCur;
		const double factor = ( i == 0 || i == numData / 2 ) ? 1.0 : 2.0;
		ofs << i << "," << freq;
		for( int iChan = 0; iChan < numChannels; ++iChan ){
			const double amp = std::abs(cdata[iChan][i]) * 0.5 * T;
			const double phs = std::arg(cdata[iChan][i]) * CommonParameters::RAD2DEG;
			const double power = factor * std::norm(cdata[iChan][i]) * T;
			ofs << "," << std::setprecision(10) << std::scientific << amp;
			ofs << "," << std::setprecision(10) << std::scientific << phs;
			ofs << "," << std::setprecision(10) << std::scientific << power;
		}
		for( int iOut = 0; iOut < numOutputVariables; ++iOut ){
			const int iChan0 = ptrControl->getChannelIndex(CommonParameters::OUTPUT, iOut);
			for( int iRR = 0; iRR < numRemoteReferenceVariables; ++iRR ){
				const int iChan1 = ptrControl->getChannelIndex(CommonParameters::REMOTE_REFERENCE, iRR);
				std::complex<double> cross = factor / T * cdata[iChan0][i] * std::conj(cdata[iChan1][i]);
				ofs << "," << std::setprecision(10) << std::scientific << std::abs(cross);
				ofs << "," << std::setprecision(10) << std::scientific << std::arg(cross) * CommonParameters::RAD2DEG;
			}
		}
		for( int iInp = 0; iInp < numInputVariables; ++iInp ){
			const int iChan0 = ptrControl->getChannelIndex(CommonParameters::INPUT, iInp);
			for( int iRR = 0; iRR < numRemoteReferenceVariables; ++iRR ){
				const int iChan1 = ptrControl->getChannelIndex(CommonParameters::REMOTE_REFERENCE, iRR);
				std::complex<double> cross = factor / T * cdata[iChan0][i] * std::conj(cdata[iChan1][i]);
				ofs << "," << std::setprecision(10) << std::scientific << std::abs(cross);
				ofs << "," << std::setprecision(10) << std::scientific << std::arg(cross) * CommonParameters::RAD2DEG;
			}
		}
		for( int iOut = 0; iOut < numOutputVariables; ++iOut ){
			const int iChan0 = ptrControl->getChannelIndex(CommonParameters::OUTPUT, iOut);
			for( int iInp = 0; iInp < numInputVariables; ++iInp ){
				const int iChan1 = ptrControl->getChannelIndex(CommonParameters::INPUT, iInp);
				std::complex<double> cross = factor / T * cdata[iChan0][i] * std::conj(cdata[iChan1][i]);
				ofs << "," << std::setprecision(10) << std::scientific << std::abs(cross);
				ofs << "," << std::setprecision(10) << std::scientific << std::arg(cross) * CommonParameters::RAD2DEG;
			}
		}
		ofs << std::endl;
	}
	ofs.close();

}

// Output frequency-domain data
void Analysis::outputFrequencyDomainData( const int iSegLen, const int freqDegree, const int numSegmentsTotal, std::complex<double>** ftval ) const{

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeLogMessage("Output frequency-domain data");

	const Control* const ptrControl = Control::getInstance();
	for( int iChan = 0; iChan < ptrControl->getNumberOfChannels(); ++iChan ){
		std::ostringstream fileName;
		fileName << "freq_domain_data_segm_" << iSegLen << "_index_" << freqDegree << "_chan_" << iChan << ".csv";
		std::ofstream ofs;
		ofs.open( fileName.str().c_str(), std::ios::out );
		if( ofs.fail() ){
			ptrOutputFiles->writeLogMessage("File open error !! : " + fileName.str());
		}
		ofs << "real,imaginary" << std::endl;
		for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
			ofs << std::setprecision(12) << std::scientific << ftval[iChan][iSeg].real() << ",";
			ofs << std::setprecision(12) << std::scientific << ftval[iChan][iSeg].imag() << std::endl;
		}
		ofs.close();
	}

}

// Output time-series data
void Analysis::outputTimeSeriesData( const std::vector<CommonParameters::DataFileSet>& dataFileSets, const bool afterPreprocessing ) const{

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	if(afterPreprocessing){
		ptrOutputFiles->writeLogMessage("Output time-series data after preprocessing");
	}else{
		ptrOutputFiles->writeLogMessage("Output time-series data");
	}

	const Control* const ptrControl = Control::getInstance();
	const int numChannels = ptrControl->getNumberOfChannels();
	int iSection(0);
	for( std::vector<CommonParameters::DataFileSet>::const_iterator itr = dataFileSets.begin(); itr != dataFileSets.end(); ++itr, ++iSection ){
		const int numData = itr->numDataPoints;
		for( int iChan = 0; iChan < numChannels; ++iChan ){
			// Output time serieas
			if( ptrControl->doesOutputTimeSeriesToCsv() ){
				std::ostringstream oss;
				if(afterPreprocessing){
					oss << "time_series_after_preprocessing_sect_" << iSection << "_chan_" << iChan << ".csv"; 
				}else{
					oss << "time_series_sect_" << iSection << "_chan_" << iChan << ".csv"; 
				}
				std::ofstream ofs;
				ofs.open( oss.str().c_str(), std::ios::out );
				if( ofs.fail() ){
					ptrOutputFiles->writeLogMessage("File open error !! : " + oss.str());
				}
				for( int i = 0; i < numData; ++i ){
					ofs << std::setprecision(12) << std::scientific << itr->dataFile[iChan].data[i] << std::endl;
				}
				ofs.close();
			}
		}
	}

}

// Evaluate characteristics of time-series data prior to the estimation of the response functions
void Analysis::priorEvaluationOfTimeSeriesData( const int interval, const std::vector<CommonParameters::DataFileSet>& dataFileSets,
	const bool afterPreprocessing ) const{

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	if(afterPreprocessing){
		ptrOutputFiles->writeLogMessage("Perform prior evaluation of time-domain characteristics after preprocessing");
	}else{
		ptrOutputFiles->writeLogMessage("Perform prior evaluation of time-domain characteristics");
	}

	const Control* const ptrControl = Control::getInstance();
	const double samplingFrequency = ptrControl->getSamplingFrequency();

	const int numChannels = ptrControl->getNumberOfChannels();
	std::vector<double>* means = new std::vector<double>[numChannels];
	std::vector<double>* meanSquares = new std::vector<double>[numChannels];
	std::vector<double>* variances = new std::vector<double>[numChannels];
	std::vector<int> sections;
	std::vector<std::string> times;

	int section(0);
	for( std::vector<CommonParameters::DataFileSet>::const_iterator itr = dataFileSets.begin(); itr != dataFileSets.end(); ++itr, ++section ){
		const int numData = itr->numDataPoints;
		for( int iChan = 0; iChan < numChannels; ++iChan ){
			int icount(0);
			while(icount + interval < numData){
				const double mean = Util::calculateMeanValue( interval, &(itr->dataFile[iChan].data[icount]) );
				const double mean2 = Util::calculateMeanSquareValue( interval, &(itr->dataFile[iChan].data[icount]) );
				const double var = Util::calculateVariance( interval, mean , &(itr->dataFile[iChan].data[icount]) );
				means[iChan].push_back(mean);
				meanSquares[iChan].push_back(mean2);
				variances[iChan].push_back(var);
				icount += interval;
			}
		}
		int icount(0);
		while(icount + interval < numData){
			sections.push_back(section);
			const double elapsedTime = static_cast<double>(icount) / samplingFrequency;
			times.push_back( ptrControl->getTimeFromStartTimeOfEachSection(section, elapsedTime) );
			icount += interval;
		}
	}

	const std::string fileName = afterPreprocessing ? "stat_time_after_preprocessing.csv" : "stat_time.csv";
	static std::ofstream ofs;
	ofs.open( fileName.c_str(), std::ios::out );
	if( ofs.fail() ){
		ptrOutputFiles->writeErrorMessage( "File open error : " + fileName );
	}

	ofs << "index,section#";
	ofs << ",time";
	for( int iChan = 0; iChan < numChannels; ++iChan ){
		ofs << ",mean_channel" << iChan;
		ofs << ",meanSquare_channel" << iChan;
		ofs << ",variance_channel" << iChan;
	}
	ofs << std::endl;

	const int numSegment = static_cast<int>( means[0].size() );
	for( int iSeg = 0; iSeg < numSegment; ++iSeg ){
		ofs << iSeg << "," << sections[iSeg];
		ofs <<  "," << times[iSeg];
		for( int iChan = 0; iChan < numChannels; ++iChan ){
			ofs << "," << std::setprecision(10) << std::scientific << means[iChan][iSeg];
			ofs	<< "," << std::setprecision(10) << std::scientific << meanSquares[iChan][iSeg];
			ofs	<< "," << std::setprecision(10) << std::scientific << variances[iChan][iSeg];
		}
		ofs << std::endl;
	}

	ofs.close();

}

// Evaluate characteristics of frequency data from all data prior to the estimation of the response functions
void Analysis::priorEvaluationOfFrequencyDataFromAllData( const std::vector<CommonParameters::DataFileSet>& dataFileSets,
	const bool afterPreprocessing ) const{

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	if(afterPreprocessing){
		ptrOutputFiles->writeLogMessage("Perform prior evaluation of frequency-domain characteristics from all data after preprocessing");
	}else{
		ptrOutputFiles->writeLogMessage("Perform prior evaluation of frequency-domain characteristics from all data");
	}

	const Control* const ptrControl = Control::getInstance();
	const int numChannels = ptrControl->getNumberOfChannels();
	const int numOutputVariables = ptrControl->getNumOutputVariables();		
	const int numInputVariables = ptrControl->getNumInputVariables();		
	const int numRemoteReferenceVariables = ptrControl->getNumRemoteReferenceVariables();
	const double samplingFrequency = ptrControl->getSamplingFrequency();
	const Control::ParamsForFreqDomainEvaluation params = ptrControl->getParamsForFreqDomainEvaluation();

	int section(0);
	for( std::vector<CommonParameters::DataFileSet>::const_iterator itr = dataFileSets.begin(); itr != dataFileSets.end(); ++itr, ++section ){
		const int numData = itr->numDataPoints;
		int numDataForFFT(-1);	
		for( int i = 1; i < 1000; ++i ){
			const int vpow2 = static_cast<int>( pow(2,i) );
			if( vpow2 >= numData ){
				numDataForFFT = vpow2;
				break;
			}
		}

		// FFT
		std::complex<double>** cdata = new std::complex<double>*[numChannels];
		for( int iChan = 0; iChan < numChannels; ++iChan ){
			cdata[iChan] = new std::complex<double>[numDataForFFT];
			for( int i = 0; i < numData; ++i ){
				cdata[iChan][i] = std::complex<double>(itr->dataFile[iChan].data[i], 0.0);
			}
			for( int i = numData; i < numDataForFFT; ++i ){
				cdata[iChan][i] = std::complex<double>(0.0, 0.0);
			}
			Util::fourierTransform(numDataForFFT, cdata[iChan]);
		}
		outputSpectrum( cdata, numDataForFFT, section, afterPreprocessing, false );
		outputAverageSpectrum( cdata, numDataForFFT, section, afterPreprocessing, false );

		// Calibration
		for( int iChan = 0; iChan < numChannels; ++iChan ){
			// Stop to write warning otherwise huge number of warnings can be outputed because the frequency
			// ranges of calibration files are usually too narrow to cover the whole frequency range.
			ptrOutputFiles->stopToWriteWarningMessage();
			for( int i = 1; i < numDataForFFT / 2 + 1; ++i ){
				// Exclude frequency = 0 (Hz) because log(frequency) is caculated for calibration correcton
				const double freq = static_cast<double>(i) / static_cast<double>(numDataForFFT) * samplingFrequency;
				if( freq < params.startOutputFrequency || freq > params.endOutputFrequency ){
					continue;
				}
				calibrationCorrection( iChan, 1, freq, &cdata[iChan][i], afterPreprocessing );
			}
			ptrOutputFiles->restartToWriteWarningMessage();
		}
		outputSpectrum( cdata, numDataForFFT, section, afterPreprocessing, true );
		outputAverageSpectrum( cdata, numDataForFFT, section, afterPreprocessing, true );

		// Release memory 
		for( int iChan = 0; iChan < numChannels; ++iChan ){
			delete [] cdata[iChan];
		}
		delete [] cdata;	
	}

}

// Evaluate characteristics of data segments
void Analysis::priorEvaluationOfDataSegments( const int numSegmentsTotal, std::complex<double>** ftval,
	const std::vector< std::pair<std::string, std::string> >& times, const double timeLength, const std::string& fileName ) const{
	
	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeLogMessage("Prior evaluation of data segments");

	const Control* const ptrControl = Control::getInstance();
	const Control::ParamsForDataSegmentEvaluation params = ptrControl->getParamsForDataSegmentEvaluation();
	const int numSegments = params.numDataSegments;
	const int numChannels = ptrControl->getNumberOfChannels();

	std::ofstream ofs;
	ofs.open( fileName.c_str(), std::ios::out );
	if( ofs.fail() ){
		ptrOutputFiles->writeErrorMessage( "File open error : " + fileName );
	}
	ofs << "start_index,end_index,start_time,end_time";
	const int numOutputVariables = ptrControl->getNumOutputVariables();
	const int numInputVariables = ptrControl->getNumInputVariables();
	const int numRemoteReferenceVariables = ptrControl->getNumRemoteReferenceVariables();
	if( numSegments >= 2 ){
		for( int iVar = 0; iVar < numOutputVariables; ++iVar ){
			const int var = ptrControl->getChannelIndex( CommonParameters::OUTPUT, iVar );
			ofs << ",coherence_" << var;
		}
		for( int iVar = 0; iVar < numOutputVariables; ++iVar ){
			const int out = ptrControl->getChannelIndex( CommonParameters::OUTPUT, iVar );
			const int in0 = ptrControl->getChannelIndex( CommonParameters::INPUT, 0 );
			const int in1 = ptrControl->getChannelIndex( CommonParameters::INPUT, 1 );
			ofs << ",resp_real_" << out << "_" << in0;
			ofs << ",resp_imag_" << out << "_" << in0;
			ofs << ",resp_real_" << out << "_" << in1;
			ofs << ",resp_imag_" << out << "_" << in1;
		}
	}
	for( int iVar = 0; iVar < numOutputVariables / 2; ++iVar ){
		const int var0 = ptrControl->getChannelIndex( CommonParameters::OUTPUT, iVar );
		const int var1 = ptrControl->getChannelIndex( CommonParameters::OUTPUT, iVar+1 );
		ofs << ",degree_of_polarization_" << var0 << "_" << var1;
		ofs << ",polarization_direction_" << var0 << "_" << var1;
		ofs << ",ellipticity_" << var0 << "_" << var1;
	}
	for( int iVar = 0; iVar < numInputVariables / 2; ++iVar ){
		const int var0 = ptrControl->getChannelIndex( CommonParameters::INPUT, iVar );
		const int var1 = ptrControl->getChannelIndex( CommonParameters::INPUT, iVar+1 );
		ofs << ",degree_of_polarization_" << var0 << "_" << var1;
		ofs << ",polarization_direction_" << var0 << "_" << var1;
		ofs << ",ellipticity_" << var0 << "_" << var1;
	}
	for( int iVar = 0; iVar < numRemoteReferenceVariables / 2; ++iVar ){
		const int var0 = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, iVar );
		const int var1 = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, iVar+1 );
		ofs << ",degree_of_polarization_" << var0 << "_" << var1;
		ofs << ",polarization_direction_" << var0 << "_" << var1;
		ofs << ",ellipticity_" << var0 << "_" << var1;
	}
	for( int iVar = 0; iVar < numOutputVariables; ++iVar ){
		const int var = ptrControl->getChannelIndex( CommonParameters::OUTPUT, iVar );
		ofs << ",auto_power_" << var;
	}
	for( int iVar = 0; iVar < numInputVariables; ++iVar ){
		const int var = ptrControl->getChannelIndex( CommonParameters::INPUT, iVar );
		ofs << ",auto_power_" << var;
	}
	for( int iVar = 0; iVar < numRemoteReferenceVariables; ++iVar ){
		const int var = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, iVar );
		ofs << ",auto_power_" << var;
	}
	if(numOutputVariables == 3){
		ofs << ",Poynting_vector_x,Poynting_vector_y,Poynting_vector_z";
	}
	ofs << std::endl;

	std::vector<int> segmentIndexes;
	segmentIndexes.reserve(numSegments);
	for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
		segmentIndexes.push_back(iSeg);
		if( static_cast<int>(segmentIndexes.size()) == numSegments ){
			priorEvaluationOfDataSegmentsAux( segmentIndexes, ftval, times, timeLength, numSegments, ofs );
			segmentIndexes.clear();
		}
	}

	ofs.close();

}

// Auxiliary function for evaluating characteristics of data segments
void Analysis::priorEvaluationOfDataSegmentsAux( const std::vector<int>& segmentIndexes, std::complex<double>** ftval,
	const std::vector< std::pair<std::string, std::string> >& times, const double timeLength, const int numSegments,
	std::ofstream& ofs ) const{

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	const Control* const ptrControl = Control::getInstance();
	const int numChannels = ptrControl->getNumberOfChannels();

	// Copy data
	std::complex<double>** data = new std::complex<double>*[numChannels];
	for( int iChan = 0; iChan < numChannels; ++iChan ){
		data[iChan] = new std::complex<double>[numSegments];
		int index(0);
		for( std::vector<int>::const_iterator itr = segmentIndexes.begin(); itr != segmentIndexes.end(); ++itr, ++index ){
			data[iChan][index] = ftval[iChan][*itr];
		}
	}
	// Main calculation
	const int index0 = segmentIndexes.front();
	const int index1 = segmentIndexes.back();
	const std::string timeStart = times[index0].first;
	const std::string timeEnd = times[index1].second;
	ofs << index0 << "," << index1 << "," << timeStart << "," << timeEnd;

	const int numOutputVariables = ptrControl->getNumOutputVariables();
	const int numInputVariables = ptrControl->getNumInputVariables();
	const int numRemoteReferenceVariables = ptrControl->getNumRemoteReferenceVariables();
	std::complex<double>* resp0 = NULL;
	std::complex<double>* resp1 = NULL;
	if( numSegments >= 2 ){
		resp0 = new std::complex<double>[numOutputVariables];
		resp1 = new std::complex<double>[numOutputVariables];
		// Squared coherence
		for( int iVar = 0; iVar < numOutputVariables; ++iVar ){
			const int var = ptrControl->getChannelIndex( CommonParameters::OUTPUT, iVar );
			const int in0 = ptrControl->getChannelIndex( CommonParameters::INPUT, 0 );
			const int in1 = ptrControl->getChannelIndex( CommonParameters::INPUT, 1 );
			assert( numRemoteReferenceVariables >= 2 );
			const int rr0 = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, 0 );
			const int rr1 = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, 1 );
			const double coherence = calculateResponseFunctionByOrdinaryRemoteReference( numSegments, 
				data[var], data[in0], data[in1], data[rr0], data[rr1], resp0[iVar], resp1[iVar] );
			ofs << "," << std::setprecision(10) << std::scientific << coherence;
		}
		// Response functions
		for( int iVar = 0; iVar < numOutputVariables; ++iVar ){
			ofs << "," << std::setprecision(10) << std::scientific << resp0[iVar].real();
			ofs << "," << std::setprecision(10) << std::scientific << resp0[iVar].imag();
			ofs << "," << std::setprecision(10) << std::scientific << resp1[iVar].real();
			ofs << "," << std::setprecision(10) << std::scientific << resp1[iVar].imag();
		}
	}

	// Polarization direction, ellipticity, and degree of polarization
	for( int iVar = 0; iVar < numOutputVariables / 2; ++iVar ){
		const int var0 = ptrControl->getChannelIndex( CommonParameters::OUTPUT, iVar );
		const int var1 = ptrControl->getChannelIndex( CommonParameters::OUTPUT, iVar+1 );
		double Gvar0(0.0);	
		double Gvar1(0.0);	
		std::complex<double> Gvar01(0.0, 0.0);
		for( int i = 0; i < numSegments; ++i ){
			Gvar0 += std::norm(data[var0][i]);
			Gvar1 += std::norm(data[var1][i]);
			Gvar01 += data[var0][i] * std::conj(data[var1][i]);
		}
#if 0
		const double polarizationDirection = 0.5 * atan(2.0 * Gvar01.real() / (Gvar0 - Gvar1)) * CommonParameters::RAD2DEG;
		ofs << "," << std::setprecision(10) << std::scientific << polarizationDirection;
		const double chi = 0.5 * asin(2.0 * Gvar01.imag() / (Gvar0 + Gvar1));
		const double ellipticity = fabs(tan(chi));
		ofs << "," << std::setprecision(10) << std::scientific << ellipticity;
#endif
		const double det = Gvar0 * Gvar1 - std::norm(Gvar01);
		const double tr = Gvar0 + Gvar1;
		const double degreeOfPolarization = sqrt(1.0 - 4.0 * det / (tr * tr));
		ofs << "," << std::setprecision(10) << std::scientific << degreeOfPolarization;
		const double D = 0.5 * tr - 0.5 * sqrt(tr * tr - 4.0 * det);
		const double P0 = Gvar0 - D;
		const double P1 = Gvar1 - D;
		const std::complex<double> P01 = Gvar01;
		const double polarizationDirection2 = 0.5 * atan(2.0 * P01.real() / (P0 - P1)) * CommonParameters::RAD2DEG;
		ofs << "," << std::setprecision(10) << std::scientific << polarizationDirection2;
		const double chi2 = 0.5 * asin(2.0 * P01.imag() / (P0 + P1));
		const double ellipticity2 = fabs(tan(chi2));
		ofs << "," << std::setprecision(10) << std::scientific << ellipticity2;
	}
	for( int iVar = 0; iVar < numInputVariables / 2; ++iVar ){
		const int var0 = ptrControl->getChannelIndex( CommonParameters::INPUT, iVar );
		const int var1 = ptrControl->getChannelIndex( CommonParameters::INPUT, iVar+1 );
		double Gvar0(0.0);	
		double Gvar1(0.0);	
		std::complex<double> Gvar01(0.0, 0.0);
		for( int i = 0; i < numSegments; ++i ){
			Gvar0 += std::norm(data[var0][i]);
			Gvar1 += std::norm(data[var1][i]);
			Gvar01 += data[var0][i] * std::conj(data[var1][i]);
		}
#if 0
		const double polarizationDirection = 0.5 * atan(2.0 * Gvar01.real() / (Gvar0 - Gvar1)) * CommonParameters::RAD2DEG;
		ofs << "," << std::setprecision(10) << std::scientific << polarizationDirection;
		const double chi = 0.5 * asin(2.0 * Gvar01.imag() / (Gvar0 + Gvar1));
		const double ellipticity = fabs(tan(chi));
		ofs << "," << std::setprecision(10) << std::scientific << ellipticity;
#endif
		const double det = Gvar0 * Gvar1 - std::norm(Gvar01);
		const double tr = Gvar0 + Gvar1;
		const double degreeOfPolarization = sqrt(1.0 - 4.0 * det / (tr * tr));
		ofs << "," << std::setprecision(10) << std::scientific << degreeOfPolarization;
		const double D = 0.5 * tr - 0.5 * sqrt(tr * tr - 4.0 * det);
		const double P0 = Gvar0 - D;
		const double P1 = Gvar1 - D;
		const std::complex<double> P01 = Gvar01;
		const double polarizationDirection2 = 0.5 * atan(2.0 * P01.real() / (P0 - P1)) * CommonParameters::RAD2DEG;
		ofs << "," << std::setprecision(10) << std::scientific << polarizationDirection2;
		const double chi2 = 0.5 * asin(2.0 * P01.imag() / (P0 + P1));
		const double ellipticity2 = fabs(tan(chi2));
		ofs << "," << std::setprecision(10) << std::scientific << ellipticity2;
	}
	for( int iVar = 0; iVar < numRemoteReferenceVariables / 2; ++iVar ){
		const int var0 = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, iVar );
		const int var1 = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, iVar+1 );
		double Gvar0(0.0);
		double Gvar1(0.0);
		std::complex<double> Gvar01(0.0, 0.0);
		for( int i = 0; i < numSegments; ++i ){
			Gvar0 += std::norm(data[var0][i]);
			Gvar1 += std::norm(data[var1][i]);
			Gvar01 += data[var0][i] * std::conj(data[var1][i]);
		}
#if 0
		const double polarizationDirection = 0.5 * atan(2.0 * Gvar01.real() / (Gvar0 - Gvar1)) * CommonParameters::RAD2DEG;
		ofs << "," << std::setprecision(10) << std::scientific << polarizationDirection;
		const double chi = 0.5 * asin(2.0 * Gvar01.imag() / (Gvar0 + Gvar1));
		const double ellipticity = fabs(tan(chi));
		ofs << "," << std::setprecision(10) << std::scientific << ellipticity;
#endif
		const double det = Gvar0 * Gvar1 - std::norm(Gvar01);
		const double tr = Gvar0 + Gvar1;
		const double degreeOfPolarization = sqrt(1.0 - 4.0 * det / (tr * tr));
		ofs << "," << std::setprecision(10) << std::scientific << degreeOfPolarization;
		const double D = 0.5 * tr - 0.5 * sqrt(tr * tr - 4.0 * det);
		const double P0 = Gvar0 - D;
		const double P1 = Gvar1 - D;
		const std::complex<double> P01 = Gvar01;
		const double polarizationDirection2 = 0.5 * atan(2.0 * P01.real() / (P0 - P1)) * CommonParameters::RAD2DEG;
		ofs << "," << std::setprecision(10) << std::scientific << polarizationDirection2;
		const double chi2 = 0.5 * asin(2.0 * P01.imag() / (P0 + P1));
		const double ellipticity2 = fabs(tan(chi2));
		ofs << "," << std::setprecision(10) << std::scientific << ellipticity2;
	}

	// Auto power
	const double factor = 2.0 / timeLength / static_cast<double>(numSegments);
	for( int iVar = 0; iVar < numOutputVariables; ++iVar ){
		double Gvar(0.0);
		const int var = ptrControl->getChannelIndex( CommonParameters::OUTPUT, iVar );
		for( int i = 0; i < numSegments; ++i ){
			Gvar += std::norm(data[var][i]);
		}
		ofs << "," << std::setprecision(10) << std::scientific << Gvar * factor;
	}
	for( int iVar = 0; iVar < numInputVariables; ++iVar ){
		double Gvar(0.0);
		const int var = ptrControl->getChannelIndex( CommonParameters::INPUT, iVar );
		for( int i = 0; i < numSegments; ++i ){
			Gvar += std::norm(data[var][i]);
		}
		ofs << "," << std::setprecision(10) << std::scientific << Gvar * factor;
	}
	for( int iVar = 0; iVar < numRemoteReferenceVariables; ++iVar ){
		double Gvar(0.0);
		const int var = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, iVar );
		for( int i = 0; i < numSegments; ++i ){
			Gvar += std::norm(data[var][i]);
		}
		ofs << "," << std::setprecision(10) << std::scientific << Gvar * factor;
	}

	// Poynting vector
	if( numOutputVariables == 3 ){
		const int varEx = ptrControl->getChannelIndex( CommonParameters::OUTPUT, 0 );
		const int varEy = ptrControl->getChannelIndex( CommonParameters::OUTPUT, 1 );
		const int varBz = ptrControl->getChannelIndex( CommonParameters::OUTPUT, 2 );
		const int varBx = ptrControl->getChannelIndex( CommonParameters::INPUT, 0 );
		const int varBy = ptrControl->getChannelIndex( CommonParameters::INPUT, 1 );
		std::complex<double> compX(0.0, 0.0);
		std::complex<double> compY(0.0, 0.0);
		std::complex<double> compZ(0.0, 0.0);
		for( int i = 0; i < numSegments; ++i ){
			compX += data[varEy][i] * std::conj(data[varBz][i]);// assuming Ez = 0
			compY -= data[varEx][i] * std::conj(data[varBz][i]);// assuming Ez = 0
			compZ += data[varEx][i] * std::conj(data[varBy][i]) -  data[varEy][i] * std::conj(data[varBx][i]);
		}
		compX /= static_cast<double>(numSegments);
		compY /= static_cast<double>(numSegments);
		compZ /= static_cast<double>(numSegments);
		ofs << "," << std::setprecision(10) << std::scientific << 0.5 * compX.real();
		ofs << "," << std::setprecision(10) << std::scientific << 0.5 * compY.real();
		ofs << "," << std::setprecision(10) << std::scientific << 0.5 * compZ.real();
	}

	ofs << std::endl;

	if( resp0 != NULL ){
		delete [] resp0;
	}
	if( resp1 != NULL ){
		delete [] resp1;
	}
	for( int iChan = 0; iChan < numChannels; ++iChan ){
		delete [] data[iChan];
	}
	delete [] data;

}

// Select segments to be excluded by degree of magnetic polarization criteria
void Analysis::selectSegmentsToBeExcludedByDegreeOfMagneticPolarizationCriteria( const int numSegmentsTotal, std::complex<double>** ftval, std::vector<bool>& remainingSegments ) const{

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeLogMessage("Select segments to be excluded by degree of magnetic polarization criteria");

	const Control* const ptrControl = Control::getInstance();
	const int in0 = ptrControl->getChannelIndex( CommonParameters::INPUT, 0 );
	const int in1 = ptrControl->getChannelIndex( CommonParameters::INPUT, 1 );

	// Make copy
	std::vector<bool> remainingSegmentsOrg = remainingSegments;

	const Control::ParamsForDegreeOfMagneticPolarizationCriteria params = ptrControl->getParamsForDegreeOfMagneticPolarizationCriteria();
	const int numSegments = params.numSegments;
	const double threshold = params.threshold;
	std::vector<int> segmentIndexes;
	segmentIndexes.reserve(numSegments);
	for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
		segmentIndexes.push_back(iSeg);
		if( static_cast<int>(segmentIndexes.size()) == numSegments ){
			double Jxx(0.0);
			double Jyy(0.0);
			std::complex<double> Jxy(0.0, 0.0);
			for( std::vector<int>::const_iterator itr = segmentIndexes.begin(); itr != segmentIndexes.end(); ++itr ){
				Jxx += std::norm(ftval[in0][*itr]);
				Jyy += std::norm(ftval[in1][*itr]);
				Jxy += ftval[in0][*itr] * std::conj(ftval[in1][*itr]);
			}
			const double det = Jxx * Jyy - std::norm(Jxy);
			const double tr = Jxx + Jyy;
			const double degreeOfPolarization = sqrt(1.0 - 4.0 * det / (tr * tr));
			if( degreeOfPolarization > threshold ){
				for( std::vector<int>::const_iterator itr = segmentIndexes.begin(); itr != segmentIndexes.end(); ++itr ){
					remainingSegments[*itr] = false;
				}
			}
			segmentIndexes.clear();
		}
	}

	int numOfRemainingSegments(0);
	for( int iSeg = 0 ; iSeg < numSegmentsTotal; ++iSeg ){
		if(remainingSegments[iSeg]){
			++numOfRemainingSegments;
		}
	}
	const double percent = static_cast<double>(numOfRemainingSegments) / static_cast<double>(numSegmentsTotal) * 100.0;
	ptrOutputFiles->writeLogMessage("Number of remaining segments: " + Util::toString(numOfRemainingSegments) + " (" + Util::toString(percent) + " %)");

}

// Select segments to be excluded by mean square criteria
void Analysis::selectSegmentsToBeExcludedByMeanSquareCriteria( const int numSegmentsTotal, const std::vector<double>* const meanSquares, 
	const std::vector< std::pair<std::string, std::string> >& times, std::vector<bool>& remainingSegments ) const{

	const Control* const ptrControl = Control::getInstance();
	const Control::ParamsForMeanSquareCriteria params = ptrControl->getParamsForMeanSquareCriteria();
	if( !params.applyMeanSquareCriteria ){
		return;
	}
	const double nsigma = params.nsigma;

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeLogMessage("Select segments to be excluded by mean square criteria");

	const int numChannels = ptrControl->getNumberOfChannels();
	for( int iChan = 0; iChan < numChannels; ++iChan ){
		const double median = Util::calculateMedian(meanSquares[iChan]);
		const double MADN = Util::calculateMADN(meanSquares[iChan]);
		int numRemoved(0);
		for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
			if( fabs( meanSquares[iChan][iSeg] - median ) > nsigma * MADN ){
#ifdef _DEBUG_WRITE
				const std::string timeStart = times[iSeg].first;
				const std::string timeEnd = times[iSeg].second;
				const double absDif = fabs( meanSquares[iChan][iSeg] - median );
				const double ratio = absDif / MADN;
#endif
				remainingSegments[iSeg] = false;
				++numRemoved;
			}
		}
		ptrOutputFiles->writeLogMessage( "At channel " + Util::toString(iChan) +
			", mean squares of " + Util::toString(numRemoved) + " segments are greater than the criteria.");
	}

	int numOfRemainingSegments(0);
	for( int iSeg = 0 ; iSeg < numSegmentsTotal; ++iSeg ){
		if(remainingSegments[iSeg]){
			++numOfRemainingSegments;
		}
	}
	const double percent = static_cast<double>(numOfRemainingSegments) / static_cast<double>(numSegmentsTotal) * 100.0;
	ptrOutputFiles->writeLogMessage("Number of remaining segments: " + Util::toString(numOfRemainingSegments) + " (" + Util::toString(percent) + " %)");

}

// Select segments to be excluded by polarizatiton direction criteria
void Analysis::selectSegmentsToBeExcludedByMagneticPolarizatitonDirectionCriteria( const int numSegmentsTotal, std::complex<double>** ftval, std::vector<bool>& remainingSegments ) const{

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeLogMessage("Select segments to be excluded by magnetic polarizatiton direction criteria");

	const Control* const ptrControl = Control::getInstance();
	const int in0 = ptrControl->getChannelIndex( CommonParameters::INPUT, 0 );
	const int in1 = ptrControl->getChannelIndex( CommonParameters::INPUT, 1 );

	double* numMPDs = new double[90];
	for( int iDeg = 0; iDeg < 90; ++iDeg ){
		numMPDs[iDeg] = 0.0;
	}
	std::vector<int>* segmentIndexes = new std::vector<int>[90];

	for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
		const double Jxx = std::norm(ftval[in0][iSeg]);
		const double Jyy = std::norm(ftval[in1][iSeg]);
		const std::complex<double> Jxy = ftval[in0][iSeg] * std::conj(ftval[in1][iSeg]);
		//const double polarizationDirection = 0.5 * atan(2.0 * Jxy.real() / (Jxx - Jyy)) * CommonParameters::RAD2DEG;
		const double det = Jxx * Jyy - std::norm(Jxy);
		const double tr = Jxx + Jyy;
		const double D = 0.5 * tr - 0.5 * sqrt(tr * tr - 4.0 * det);
		const double Pxx = Jxx - D;
		const double Pyy = Jyy - D;
		const std::complex<double> Pxy = Jxy;
		const double polarizationDirection = 0.5 * atan(2.0 * Pxy.real() / (Pxx - Pyy)) * CommonParameters::RAD2DEG;
		int iDeg = static_cast<int>(polarizationDirection + 45.0);
		if( iDeg < 0 ){
			iDeg = 0;
		}else if( iDeg > 89 ){
			iDeg = 89;
		}
		++numMPDs[iDeg];
		segmentIndexes[iDeg].push_back(iSeg);
	}

	const double threshold = ptrControl->getParamsForMagneticPolarizatitonDirectionCriteria().threshold;
	const double expectedValue = static_cast<double>(numSegmentsTotal) / 90.0;// assuming uniform distribution
	for( int iDeg = 0; iDeg < 90; ++iDeg ){
#ifdef _DEBUG_WRITE
		std::cout << numMPDs[iDeg] << std::endl;
#endif
		if( numMPDs[iDeg] > threshold * expectedValue ){
			for( std::vector<int>::const_iterator itr = segmentIndexes[iDeg].begin(); itr != segmentIndexes[iDeg].end(); ++itr ){
				remainingSegments[*itr] = false;
			}
		}
	}
	delete [] numMPDs;
	delete [] segmentIndexes;

	int numOfRemainingSegments(0);
	for( int iSeg = 0 ; iSeg < numSegmentsTotal; ++iSeg ){
		if(remainingSegments[iSeg]){
			++numOfRemainingSegments;
		}
	}
	const double percent = static_cast<double>(numOfRemainingSegments) / static_cast<double>(numSegmentsTotal) * 100.0;
	ptrOutputFiles->writeLogMessage("Number of remaining segments: " + Util::toString(numOfRemainingSegments) + " (" + Util::toString(percent) + " %)");

}

// Select segments to be excluded by square coherence criteria
void Analysis::selectSegmentsToBeExcludedBySquareCoherenceCriteria( int numSegmentsTotal, std::complex<double>** ftval, 
	std::vector<bool>& remainingSegments ) const{

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeLogMessage("Select segments to be excluded by square coherence criteria");

	const Control* const ptrControl = Control::getInstance();
	const int numChannels = ptrControl->getNumberOfChannels();
	const Control::ParamsForSquareCoherenceCriteria params = ptrControl->getParamsForSquareCoherenceCriteria();
	
	std::vector<int> segmentIndexes;
	segmentIndexes.reserve(params.numSegments);
	std::complex<double>** data = new std::complex<double>*[numChannels];
	for( int iChan = 0; iChan < numChannels; ++iChan ){
		data[iChan] = new std::complex<double>[params.numSegments];
	}
	for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
		segmentIndexes.push_back(iSeg);
		if( static_cast<int>(segmentIndexes.size()) == params.numSegments ){
			for( int iChan = 0; iChan < numChannels; ++iChan ){
				int index(0);
				for( std::vector<int>::const_iterator itr = segmentIndexes.begin(); itr != segmentIndexes.end(); ++itr, ++index ){
					data[iChan][index] = ftval[iChan][*itr];
				}
			}
			const double coherence = selectSegmentsToBeExcludedBySquareCoherenceCriteriaAux(params.numSegments, data);
			if( coherence < params.threshold ){
				for( std::vector<int>::const_iterator itr = segmentIndexes.begin(); itr != segmentIndexes.end(); ++itr ){
					remainingSegments[*itr] = false;
				}
			}
			segmentIndexes.clear();
		}
	}
	for( int iChan = 0; iChan < numChannels; ++iChan ){
		delete [] data[iChan];
	}
	delete [] data;

	int numOfRemainingSegments(0);
	for( int iSeg = 0 ; iSeg < numSegmentsTotal; ++iSeg ){
		if(remainingSegments[iSeg]){
			++numOfRemainingSegments;
		}
	}
	const double percent = static_cast<double>(numOfRemainingSegments) / static_cast<double>(numSegmentsTotal) * 100.0;
	ptrOutputFiles->writeLogMessage("Number of remaining segments: " + Util::toString(numOfRemainingSegments) + " (" + Util::toString(percent) + " %)");

}

// Auxiliary function for selecting segments to be excluded by square coherence criteria
double Analysis::selectSegmentsToBeExcludedBySquareCoherenceCriteriaAux( const int numSegments, std::complex<double>** data ) const{

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	const Control* const ptrControl = Control::getInstance();
	const int numChannels = ptrControl->getNumberOfChannels();

	double coherence(0.0);

	const int numOutputVariables = ptrControl->getNumOutputVariables();
	const int numRemoteReferenceVariables = ptrControl->getNumRemoteReferenceVariables();
	std::complex<double>* resp0 = new std::complex<double>[numOutputVariables];
	std::complex<double>* resp1 = new std::complex<double>[numOutputVariables];
	for( int iVar = 0; iVar < numOutputVariables; ++iVar ){
		const int var = ptrControl->getChannelIndex( CommonParameters::OUTPUT, iVar );
		const int in0 = ptrControl->getChannelIndex( CommonParameters::INPUT, 0 );
		const int in1 = ptrControl->getChannelIndex( CommonParameters::INPUT, 1 );
		assert( numRemoteReferenceVariables >= 2 );
		const int rr0 = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, 0 );
		const int rr1 = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, 1 );
		coherence = calculateResponseFunctionByOrdinaryRemoteReference( numSegments,
			data[var], data[in0], data[in1], data[rr0], data[rr1], resp0[iVar], resp1[iVar] );
	}
	
	delete [] resp0;
	delete [] resp1;

	return coherence;

}

// Select segments to be excluded by square coherence criteria with random sampling
void Analysis::selectSegmentsToBeExcludedBySquareCoherenceCriteriaWithRandomSampling( const int numSegmentsTotal, 
	std::complex<double>** ftval, std::vector<bool>& remainingSegments ) const{

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeLogMessage("Select segments to be excluded by square coherence criteria with random sampling");

	const Control* const ptrControl = Control::getInstance();
	const int numChannels = ptrControl->getNumberOfChannels();
	const Control::ParamsForSquareCoherenceCriteriaWithRandomSampling params = ptrControl->getParamsForSquareCoherenceCriteriaWithRandomSampling ();

	if( numSegmentsTotal < params.numSegments ){
		ptrOutputFiles->writeWarningMessage("Number of segments is less that number of segments used square coherence criteria");
		return;
	}

	std::set<int> segmentIndexes;

#ifdef _RAND
	srand(1234);
#else
#ifdef _MERSENNE_TWISTER_ORIGINAL
	init_genrand64(1234);
#else
	std::mt19937_64 gen(1234);
	std::uniform_int_distribution<int> uniformDistibution(0, numSegmentsTotal - 1);
#endif
#endif
	const int maxInterationNum = 10000;

	std::complex<double>** data = new std::complex<double>*[numChannels];
	for( int iChan = 0; iChan < numChannels; ++iChan ){
		data[iChan] = new std::complex<double>[params.numSegments];
	}
	std::vector<double>* squareCoherenceVectors = new std::vector<double>[numSegmentsTotal];
	for( int iSamples = 0; iSamples < params.numRandomSamples; ++iSamples ){
		// Make random samples
		for( int iter = 0; iter < maxInterationNum; ++iter ){
#ifdef _RAND
			const int iSeg = (rand() / RAND_MAX) * (numSegmentsTotal - 1);
#else
#ifdef _MERSENNE_TWISTER_ORIGINAL
			const int iSeg = static_cast<int>(genrand64_real1() * numSegmentsTotal);
#else
			const int iSeg = uniformDistibution(gen);
#endif
#endif
			segmentIndexes.insert(iSeg);
			if( segmentIndexes.size() == params.numSegments ){
				break;
			}
		}
		if( segmentIndexes.size() < params.numSegments ){
			ptrOutputFiles->writeWarningMessage("Intertion number reachs the upper limit in determing segments square coherence criteria");
			break;
		}
		for( int iChan = 0; iChan < numChannels; ++iChan ){
			int index(0);
			for( std::set<int>::const_iterator itr = segmentIndexes.begin(); itr != segmentIndexes.end(); ++itr, ++index ){
				data[iChan][index] = ftval[iChan][*itr];
			}
		}
		const double coherence = selectSegmentsToBeExcludedBySquareCoherenceCriteriaAux(params.numSegments, data);
		for( std::set<int>::const_iterator itr = segmentIndexes.begin(); itr != segmentIndexes.end(); ++itr ){
			squareCoherenceVectors[*itr].push_back(coherence);
		}
		segmentIndexes.clear();
		assert(segmentIndexes.empty());
	}
	for( int iChan = 0; iChan < numChannels; ++iChan ){
		delete [] data[iChan];
	}
	delete [] data;

	for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
		if( squareCoherenceVectors[iSeg].empty() ){
			continue;
		}
		const double medianCoherence = Util::calculateMedian(squareCoherenceVectors[iSeg]);
		if( medianCoherence < params.threshold ){
			remainingSegments[iSeg] = false;
		}
	}
	delete [] squareCoherenceVectors;

	int numOfRemainingSegmentsMod(0);
	for( int iSeg = 0 ; iSeg < numSegmentsTotal; ++iSeg ){
		if(remainingSegments[iSeg]){
			++numOfRemainingSegmentsMod;
		}
	}
	const double percent = static_cast<double>(numOfRemainingSegmentsMod) / static_cast<double>(numSegmentsTotal) * 100.0;
	ptrOutputFiles->writeLogMessage("Number of remaining segments: " + Util::toString(numOfRemainingSegmentsMod) + " (" + Util::toString(percent) + " %)");

}

// Write header to the output file for apparent resistivity and phase
void Analysis::writeHeaderToOutputFileForApparentResistivityAndPhase( std::ofstream& ofs ) const{

	const Control* const ptrControl = Control::getInstance();
	const int numOutputVariables = ptrControl->getNumOutputVariables();		
	const int numInputVariables = ptrControl->getNumRemoteReferenceVariables();

	ofs << "frequency,period";
	assert( numInputVariables == 2 );
	for( int iVar = 0; iVar < numOutputVariables; ++iVar ){
		const int out = ptrControl->getChannelIndex( CommonParameters::OUTPUT, iVar );
		const int in0 = ptrControl->getChannelIndex( CommonParameters::INPUT, 0 );
		const int in1 = ptrControl->getChannelIndex( CommonParameters::INPUT, 1 );
		ofs << ",app_res_" << out << "_" << in0;
		ofs << ",phase_" << out << "_" << in0;
		ofs << ",app_res_" << out << "_" << in1;
		ofs << ",phase_" << out << "_" << in1;
		ofs << ",coherence_" << out << "_" << in0 << "+" << in1;
	}
	for( int iVar = 0; iVar < numOutputVariables; ++iVar ){
		const int out = ptrControl->getChannelIndex( CommonParameters::OUTPUT, iVar );
		const int in0 = ptrControl->getChannelIndex( CommonParameters::INPUT, 0 );
		const int in1 = ptrControl->getChannelIndex( CommonParameters::INPUT, 1 );
		ofs << ",dapp_res_" << out << "_" << in0;
		ofs << ",dphase_" << out << "_" << in0;
		ofs << ",dapp_res_" << out << "_" << in1;
		ofs << ",dphase_" << out << "_" << in1;
	}
	ofs << std::endl;

}

// Write header to the output file of response function
void Analysis::writeHeaderToOutputFileForResponseFunctions( std::ofstream& ofs ) const{

	const Control* const ptrControl = Control::getInstance();
	const int numOutputVariables = ptrControl->getNumOutputVariables();		
	const int numInputVariables = ptrControl->getNumRemoteReferenceVariables();

	ofs << "frequency,period";
	assert( numInputVariables == 2 );
	for( int iVar = 0; iVar < numOutputVariables; ++iVar ){
		const int out = ptrControl->getChannelIndex( CommonParameters::OUTPUT, iVar );
		const int in0 = ptrControl->getChannelIndex( CommonParameters::INPUT, 0 );
		const int in1 = ptrControl->getChannelIndex( CommonParameters::INPUT, 1 );
		ofs << ",resp_real_" << out << "_" << in0;
		ofs << ",resp_imag_" << out << "_" << in0;
		ofs << ",resp_real_" << out << "_" << in1;
		ofs << ",resp_imag_" << out << "_" << in1;
		ofs << ",coherence_" << out << "_" << in0 << "+" << in1;
	}
	for( int iVar = 0; iVar < numOutputVariables; ++iVar ){
		const int out = ptrControl->getChannelIndex( CommonParameters::OUTPUT, iVar );
		const int in0 = ptrControl->getChannelIndex( CommonParameters::INPUT, 0 );
		const int in1 = ptrControl->getChannelIndex( CommonParameters::INPUT, 1 );
		ofs << ",dresp_" << out << "_" << in0;
		ofs << ",dresp_" << out << "_" << in1;
	}
	ofs << std::endl;

}