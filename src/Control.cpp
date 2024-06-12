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
#include "Control.h"
#include "OutputFiles.h"
#include "Util.h"
#include "TableOfTukeysBiweightParameters.h"
#include "AnalysisMultivariateRegression.h"
#include "AnalysisTwoStage.h"
#include "AnalysisRepeatedMedian.h"
#include "AnalysisOrdinaryRemoteReference.h"
#include "AnalysisTest.h"
#include "RobustWeightHuber.h"
#include "RobustWeightThomson.h"
#include "RobustWeightTukeysBiweights.h"
#include "Ats.h"

#include <cmath>
#include <stdlib.h>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>
#include <assert.h>
#include <iomanip>

#ifdef _USE_OMP
#include <omp.h>
#endif

// Return the instance of the class
Control* Control::getInstance(){
   	static Control instance;// The only instance
  	return &instance;
}

// Default constructer
Control::Control() :
	m_cutoffFrequencyForIIRHighPassFilter(-1.0),
	m_cutoffFrequencyForIIRLowPassFilter(-1.0),
	m_doesApplyIIRHighPassFilter(false),
	m_doesApplyIIRLowPassFilter(false),
	m_elogMTReadingOption(Control::READ_EX_EY_HX_HY_HZ_HRX_HRY_FROM_ELOGMT_DATA),
	m_errorEstimationMethod(Control::FIXED_WEIGHTS_BOOTSTRAP),
	m_numOutputVariables(-1),
	m_numInputVariables(2),
	m_numRemoteReferenceVariables(2),
	m_samplingFrequency(-1),
	m_numThreads(1),
	m_numTimeSeriesSections(-1),
	m_numRepetitionsOfBootstrap(1000),
	m_overlappingRatio(0.5),
	m_calibForMFS(false),
	m_calibForElogDual(false),
	m_calibForElogMT(false),
	m_rotationAngle(0.0),
	m_analysis(NULL),
	m_outputApparentResistivityAndPhase(false),
	m_outputLevel(0),
	m_outputFreqDomainDataToCsv(false),
	m_outputTimeSeriesToCsv(false),
	m_parameterQForNotchFilter(10.0),
	m_percentageOfOmmitedDataSubsetDeletionJackknife(5.0),
	m_procedureType(Control::ORDINARY_REMOTE_REFERENCE),
	m_readAtsBinary(false),
	m_readElogDualBinary(false),
	m_readElogMTBinary(false),
	m_typeOfElogDual(Control::ELOG1K),
	m_typeOfElogMT(Control::ELOGMT_ADU_MODE)
{
	// Parameters for ELOG calibration
	m_paramsForElogDualCalibration.unitGroupDelay = 5.17;
	m_paramsForElogMTCalibration.unitGroupDelay = 5.17;

	// Parameters for the data-segment evaluation prior to the calculation of the response functions
	m_paramsForDataSegmentEvaluation.doEvaluation = false;
	m_paramsForDataSegmentEvaluation.numDataSegments = 10;

	// Parameters for decimation
	m_paramsForDecimation.applyDecimation = false;
	m_paramsForDecimation.decimationInterval = 1;
	m_paramsForDecimation.filterLength = 100;
	m_paramsForDecimation.transitionBandWidthInLogarithmicScale = 0.5;

	// Parameters for the frequency-domain evaluation prior to the calculation of the response functions
	m_paramsForFreqDomainEvaluation.doEvaluation = false;
	m_paramsForFreqDomainEvaluation.startOutputFrequency = 1.0e-10;
	m_paramsForFreqDomainEvaluation.endOutputFrequency = 1.0e+10;
	m_paramsForFreqDomainEvaluation.numOfOutputFrequency = 100000;
	m_paramsForFreqDomainEvaluation.numOfOutputFrequencyForAverage = 100;

	// Parameters for the time-domain evaluation prior to the calculation of the response functions
	m_paramsForTimeDomainEvaluation.doEvaluation = false;
	m_paramsForTimeDomainEvaluation.timeSeriesInterval = 10;

	// Parameters for Hampel filter
	m_paramsForParamsForHampelFilter.applyHampelFilter = false;
	m_paramsForParamsForHampelFilter.numNeighborsOnEitherSide = 50;
	m_paramsForParamsForHampelFilter.nsigma = 100.0;

	// Parameters for degree of polarizatiton direction criteria
	m_paramsForDegreeOfMagneticPolarizationCriteria.applyDegreeOfMagneticPolarizationCriteria = false;
	m_paramsForDegreeOfMagneticPolarizationCriteria.numSegments = 10;
	m_paramsForDegreeOfMagneticPolarizationCriteria.threshold = 0.9;

	// Parameters for magnetic polarizatiton direction criteria
	m_paramsForMagneticPolarizatitonDirectionCriteria.applyMagneticPolarizatitonDirectionCriteria = false;
	m_paramsForMagneticPolarizatitonDirectionCriteria.threshold = 5.0;

	// Parameters for mean square criteria
	m_paramsForMeanSquareCriteria.applyMeanSquareCriteria = false;
	m_paramsForMeanSquareCriteria.nsigma = 5.0;

	// Parameters for square coherence criteria
	m_paramsForSquareCoherenceCriteria.applySquareCoherenceCriteria = false;
	m_paramsForSquareCoherenceCriteria.numSegments = 10;
	m_paramsForSquareCoherenceCriteria.threshold = 0.3;

	// Parameters for square coherence criteria with random sampling
	m_paramsForSquareCoherenceCriteriaWithRandomSampling.applySquareCoherenceCriteria = false;
	m_paramsForSquareCoherenceCriteriaWithRandomSampling.numSegments = 10;
	m_paramsForSquareCoherenceCriteriaWithRandomSampling.numRandomSamples = 1000;
	m_paramsForSquareCoherenceCriteriaWithRandomSampling.threshold = 0.3;

	// Parameters for robust calculation of autocovariance matrix
	m_paramsForRobustAutoCovariance.percentageOfSmallNumberAddedToDiagonals = 0.0;
	m_paramsForRobustAutoCovariance.maxNumOfIterations = 10;
	m_paramsForRobustAutoCovariance.convergenceCriteria = 0.05;

	// Parameters for robust prewhitening
	m_paramsForPrewhitening.applyPrewhitening = false;
	m_paramsForPrewhitening.typeOfEstimator = Control::USE_LEAST_SQUARE_ESTIMATOR_FOR_PREWHITENING;
	m_paramsForPrewhitening.maxDegreeOfARModel = 10;
	m_paramsForPrewhitening.numCandidatesOfPartialAutocorrelationFunction = 5;

	// Parameters for robust filter
	m_paramsForRobustFilter.applyRobustFilter = false;
	m_paramsForRobustFilter.replaceTimeSeriesWithFilteredValues = false;

	// Parameters for robust multivariate regression
	m_paramsForRobustMultivariateRegression.selectInitialCandidatesByRandomSamplingAtEachFrequency = true;
	m_paramsForRobustMultivariateRegression.numOfMaxInitialCandidates = 100;
	m_paramsForRobustMultivariateRegression.numOfMaxIterationsOfFirstIstep = 3;
	m_paramsForRobustMultivariateRegression.convergenceCriteriaOfFirstIstep = 0.05;
	m_paramsForRobustMultivariateRegression.numOfMaxCandidatesOfSecondIstep = 10;
	m_paramsForRobustMultivariateRegression.numOfMaxIterationsOfSecondIstep = 16;
	m_paramsForRobustMultivariateRegression.convergenceCriteriaOfSecondIstep = 0.01;
	m_paramsForRobustMultivariateRegression.startOfTimeRange = "00:00:00";
	m_paramsForRobustMultivariateRegression.endOfTimeRange = "24:00:00";
	
	// Parameters for deciding candicates for subsequent frequencies
	m_paramsForDecidingCandicatesForSubsequentFrequencies.useResponseFunctionsOfPreviousFrequency = true;

	// Parameters for the treatment of hat matrix
	m_paramsForTreatmentOfHatMatrix.applyLeverageWeights = false;
	m_paramsForTreatmentOfHatMatrix.threshold = 3.0;
	m_paramsForTreatmentOfHatMatrix.maxNumberOfOuterIteration = 5;

}

// Destructer
Control::~Control(){

	if( m_analysis != NULL ){
		delete m_analysis;
	}

	// Delete allocated memory for time-series data
	for( std::vector<CommonParameters::DataFileSet>::const_iterator itrDataFileSet = m_dataFileSets.begin();
		itrDataFileSet != m_dataFileSets.end(); ++itrDataFileSet ){
		for( std::vector<CommonParameters::DataFile>::const_iterator itrDataFile = itrDataFileSet->dataFile.begin();
			itrDataFile != itrDataFileSet->dataFile.end(); ++itrDataFile ){
			delete [] itrDataFile->data;
		}
	}

}

// Run analysis
void Control::run(const bool outputToConsole){

	////----- For test >>>>>
	//m_analysis = new AnalysisTwoStage();
	//m_analysis->test();
	//return;
	////----- For test <<<<<

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->setOutputToConsole(outputToConsole);
	readParameterFile();

	////----- For test >>>>>
	//m_analysis->test2(m_dataFileSets);
	//return;
	////----- For test <<<<<

	m_analysis->run(m_dataFileSets);

	ptrOutputFiles->writeLogMessage("End " + Util::toString(CommonParameters::programName) );

}

// Get flag specifing whether IIR high-pass filter is applied
bool Control::doesApplyIIRHighPassFilter() const{
	return m_doesApplyIIRHighPassFilter;
}

// Get flag specifing whether IIR low-pass filter is applied
bool Control::doesApplyIIRLowPassFilter() const{
	return m_doesApplyIIRLowPassFilter;
}

// Get flag specifing whether calibration for MFS is performed
bool Control::doesMakeCalibrationFileForMFS() const {
	return m_calibForMFS;
}

// Get flag specifing whether calibration for ELOG-Dual file is performed
bool Control::doesMakeCalibrationFileForElogDual() const {
	return m_calibForElogDual;
}

// Get flag specifing whether calibration for ELOG-MT file is performed
bool Control::doesMakeCalibrationFileForElogMT() const{
	return m_calibForElogMT;
}

// Get flag specifing whether apparent resistivity and phase are outputed in a seperate file
bool Control::doesOutputApparentResistivityAndPhase() const{
	return m_outputApparentResistivityAndPhase;
}

// Get flag specifing wheter output frequency domain data as csv file
bool Control::doesOutputFreqDomainDataToCsv() const{
	return m_outputFreqDomainDataToCsv;
}

// Get flag specifing wheter output time series data as csv file
bool Control::doesOutputTimeSeriesToCsv() const{
	return m_outputTimeSeriesToCsv;
}

// Get flag specifing whether input file is ATS binary file
bool Control::doesReadAtsBinary() const{
	return m_readAtsBinary;
}

// Get flag specifing whether input file is ELOG-Dual binary file
bool Control::doesReadElogDualBinary() const {
	return m_readElogDualBinary;
}

// Get flag specifing whether input file is ELOG-MT binary file
bool Control::doesReadElogMTBinary() const{
	return m_readElogMTBinary;
}

// Get azimuth
double Control::getAzimuth( const int iChan ) const{
	return m_azimuths[iChan];
}

// Get channel index
int Control::getChannelIndex( const int dataType, const int index ) const{

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();

	switch (dataType)
	{
	case CommonParameters::OUTPUT:
		if( index < 0 || index > getNumOutputVariables() ){
			ptrOutputFiles->writeErrorMessage("Index of output variable is wrong: " + Util::toString(index));
		}
		return index;
	case CommonParameters::INPUT:
		if( index < 0 || index > getNumInputVariables() ){
			ptrOutputFiles->writeErrorMessage("Index of input variable is wrong: " + Util::toString(index));
		}
		return index + getNumOutputVariables();
	case CommonParameters::REMOTE_REFERENCE:
		if( index < 0 || index > getNumRemoteReferenceVariables() ){
			ptrOutputFiles->writeErrorMessage("Index of remote reference variable is wrong: " + Util::toString(index));
		}
		return index + getNumOutputVariables() + getNumInputVariables();
	default:
		ptrOutputFiles->writeErrorMessage("Wrong type of data: " + Util::toString(dataType));
		break;
	}

	return -1;

}

// Get option of ELOG-Dual binary data reading
int Control::getElogDualReadingOption() const {
	return m_elogDualReadingOption;
}

// Get option of ELOG-MT binary data reading
int Control::getElogMTReadingOption() const {
	return m_elogMTReadingOption;
}

// Get error estimation method
int Control::getErrorEstimationMethod() const{
	return m_errorEstimationMethod;
}

// Get number of output variables
int Control::getNumOutputVariables() const{
	return m_numOutputVariables;
}

// Get number of input variables
int Control::getNumInputVariables() const{
	return m_numInputVariables;
}

// Get number of remote reference variables
int Control::getNumRemoteReferenceVariables() const{
	return m_numRemoteReferenceVariables;
}

// Get sampling frequency
double Control::getSamplingFrequency() const{
	return m_samplingFrequency;
}

// Get number of threads
int Control::getNumThreads() const{
	return m_numThreads;
}

// Get number of time-series sections
int Control::getNumTimeSeriesSections() const{
	return m_numTimeSeriesSections;
}

// Get number of segment lengths
int Control::getNumSegmentLengths() const{
	return static_cast<int>( m_segments.size() );
}

// Get number of target frequency in an input segment
int Control::getNumTargetFrequencyInSegment( const int iSeg ) const{
	return static_cast<int>( m_segments[iSeg].degrees.size() );
}

// Get number of the ranges of sections for merge
int Control::getNumRangeOfSectionsForMerge() const{
	return static_cast<int>( m_rangeOfSectionsForMerge.size() );
}

// Get number of repetitions of bootstrap method
int Control::getNumRepetitionsOfBootstrap() const{
	return m_numRepetitionsOfBootstrap;
}

// Get target frequency degrees in an input segment
int Control::getTargetFrequencyDegreesInSegment( const int iSeg, const int index ) const{
	return m_segments[iSeg].degrees[index];
}

// Get length of each segment
int Control::getSegmentLength( const int iSeg ) const{
	return m_segments[iSeg].length;
}

// Get ratio of overlapping part to whole segment length
double Control::getOverlappingRatio() const{
	return m_overlappingRatio;
}

// Get rotation angle
double Control::getRotationAngle() const{
	return m_rotationAngle;
}

// Get numebur of calibration files for MFS
int Control::getNumCalibrationFilesForMFS() const {
	return static_cast<int>(m_calibrationFilesForMFS.size());
}

// Get name of calibration file for MFS
std::string Control::getCalibrationFileNameForMFS(const int iFile) const {
	return m_calibrationFilesForMFS[iFile];
}
// Get numebur of calibration files
int Control::getNumCalibrationFiles() const{
	return static_cast<int>( m_calibrationFiles.size() );
}

// Get name of calibration file
std::string Control::getCalibrationFileName( const int iFile ) const{
	return m_calibrationFiles[iFile];
}

// Get cutoff frequency for IIR high-pass filter
double Control::getCutoffFrequencyForIIRHighPassFilter() const{
	return m_cutoffFrequencyForIIRHighPassFilter;
}

// Get cutoff frequency for IIR low-pass filter
double Control::getCutoffFrequencyForIIRLowPassFilter() const{
	return m_cutoffFrequencyForIIRLowPassFilter;
}

// Get cutoff frequencies for notch filter
std::vector<double> Control::getCutoffFrequenciesForNotchFilter() const{
	return m_cutoffFrequenciesForNotchFilter;
}

// Get directory of logger calibration files
std::string Control::getDirectoryOfLoggerCalibrationFiles() const {
	return m_directoryOfLoggerCalibrationFiles;
}

// Get number of channels
int Control::getNumberOfChannels() const{
	return getNumOutputVariables() + getNumInputVariables() + getNumRemoteReferenceVariables();
}

// Get number of frequencies
int Control::getNumberOfFrequencies() const{
	return static_cast<int>( m_frequencies.size() );
}

// Get number of cutoff frequencies for notch filter
int Control::getNumberOfCutoffFrequenciesForNotchFilter() const{
	return static_cast<int>( m_cutoffFrequenciesForNotchFilter.size() );
}

// Get frequency
double Control::getFrequency( const int iFreq ) const{
	return m_frequencies[iFreq];
}

// Get all frequencies without duplications
std::vector<double> Control::getFrequenciesAllWithoutDuplications() const{

	std::vector<double> temp = m_frequencies;
	std::sort(temp.begin(), temp.end());

	std::vector<double> freqOut;
	double freqPre(-1.0);
	for( std::vector<double>::const_iterator itr = temp.begin(); itr != temp.end(); ++itr ){
		if( fabs(*itr - freqPre) > 1.0e-6 ){
			freqOut.push_back(*itr);
			freqPre = *itr;
		}
	}

	return freqOut;

}

// Get output level
int Control::getOutputLevel() const{
	return m_outputLevel;
}

// Get parameters for ELOG-Dual calibration
Control::ParamsForElogCalibration Control::getParamsForElogDualCalibration() const {
	return m_paramsForElogDualCalibration;
}

// Get parameters for ELOG-MT calibration
Control::ParamsForElogCalibration Control::getParamsForElogMTCalibration() const{
	return m_paramsForElogMTCalibration;
}

// Get parameters for the data-segment evaluation prior to the calculation of the response functions
Control::ParamsForDataSegmentEvaluation Control::getParamsForDataSegmentEvaluation() const{
	return m_paramsForDataSegmentEvaluation;
}

// Get parameters for decimation
Control::ParamsForDecimation Control::getParamsForDecimation() const{
	return m_paramsForDecimation;
}

// Get parameters for the frequency-domain evaluation prior to the calculation of the response functions
Control::ParamsForFreqDomainEvaluation Control::getParamsForFreqDomainEvaluation() const{
	return m_paramsForFreqDomainEvaluation;
}

// Get parameters for the time-domain evaluation prior to the calculation of the response functions
Control::ParamsForTimeDomainEvaluation Control::getParamsForTimeDomainEvaluation() const{
	return m_paramsForTimeDomainEvaluation;
}

// Get parameters for Hampel filter
Control::ParamsForHampelFilter Control::getParamsForHampelFilter() const{
	return m_paramsForParamsForHampelFilter;
}

// Get parameters for degree of magnetic polarizatiton criteria
Control::ParamsForDegreeOfMagneticPolarizationCriteria Control::getParamsForDegreeOfMagneticPolarizationCriteria() const{
	return m_paramsForDegreeOfMagneticPolarizationCriteria;
}

// Get parameters for magnetic polarizatiton direction criteria
Control::ParamsForMagneticPolarizatitonDirectionCriteria Control::getParamsForMagneticPolarizatitonDirectionCriteria() const{
	return m_paramsForMagneticPolarizatitonDirectionCriteria;
}

// Get parameters for mean square criteria
Control::ParamsForMeanSquareCriteria Control::getParamsForMeanSquareCriteria() const{
	return m_paramsForMeanSquareCriteria;
}

// Get parameters for square coherence thresholding
Control::ParamsForSquareCoherenceCriteria Control::getParamsForSquareCoherenceCriteria() const{
	return m_paramsForSquareCoherenceCriteria;
}

// Get parameters for square coherence criteria with random sampling
Control::ParamsForSquareCoherenceCriteriaWithRandomSampling Control::getParamsForSquareCoherenceCriteriaWithRandomSampling() const{
	return m_paramsForSquareCoherenceCriteriaWithRandomSampling;
}

// Get parameters for robust calculation of autocovariance matrix
Control::ParamsForRobustAutoCovariance Control::getParamsForRobustAutoCovariance() const{
	return m_paramsForRobustAutoCovariance;
}

// Get parameters for robust prewhitening
Control::ParamsForPrewhitening Control::getParamsForPrewhitening() const{
	return m_paramsForPrewhitening;
}

// Get parameters for robust filter
Control::ParamsForRobustFilter Control::getParamsForRobustFilter() const{
	return m_paramsForRobustFilter;
}

// Get parameters for the treatment of hat matrix
Control::ParamsForTreatmentOfHatMatrix Control::getParamsForTreatmentOfHatMatrix() const{
	return m_paramsForTreatmentOfHatMatrix;
}

// Get parameters for robust multivariate regression
Control::ParamsForRobustMultivariateRegression Control::getParamsForRobustMultivariateRegression() const{
	return m_paramsForRobustMultivariateRegression;
}

// Get parameters for deciding candicates for subsequent frequencies
Control::ParamsForDecidingCandicatesForSubsequentFrequencies Control::getParamsForDecidingCandicatesForSubsequentFrequencies() const{
	return m_paramsForDecidingCandicatesForSubsequentFrequencies;
}

// Get parameter Q for notch filter
double Control::getParameterQForNotchFilter() const{
	return m_parameterQForNotchFilter;
}

// Get percentage of ommited data in subset deletion jackknife
double Control::getPercentageOfOmmitedDataSubsetDeletionJackknife() const{
	return m_percentageOfOmmitedDataSubsetDeletionJackknife;
}

// Get ranges of sections for merge
std::pair<int, int> Control::getRangeOfSectionsForMerge( const int iSectionAfterMerge ) const{
	return m_rangeOfSectionsForMerge[iSectionAfterMerge];
}

// Get time from start time and elapsed time (seconds)
std::string Control::getTimeFromStartTimeOfEachSection( const int sectionIndex, const double elapsedTime, const bool forAfterMergingSections ) const{

	assert( sectionIndex >=0 && sectionIndex < m_numTimeSeriesSections );

	int sectionIndexMod(sectionIndex);
	if( forAfterMergingSections && getNumRangeOfSectionsForMerge() > 0 ){
		sectionIndexMod = getRangeOfSectionsForMerge(sectionIndex).first;
	}
	const int hour0 = Util::stringToInt( m_startTimeOfEachSection[sectionIndexMod].substr(0,2) );
	const int min0 = Util::stringToInt( m_startTimeOfEachSection[sectionIndexMod].substr(3,2) );
	const double sec0 = Util::stringToDouble( m_startTimeOfEachSection[sectionIndexMod].substr(6) );

	double temp = sec0 + elapsedTime;
	if( m_paramsForDecimation.applyDecimation ){
		const double samplingFreq = getSamplingFrequency();
		const double samplingFreqBeforeDecimation = samplingFreq * static_cast<double>(m_paramsForDecimation.decimationInterval);
		const int halfOfDimension = m_paramsForDecimation.filterLength / 2;
		temp += static_cast<double>(halfOfDimension) / samplingFreqBeforeDecimation;
	}
	int dhour = static_cast<int>(floor(temp)) / 3600;
	int dmin = ( static_cast<int>(floor(temp)) - dhour * 3600 ) / 60;
	const double sec1 = temp - static_cast<double>(dhour) * 3600  - static_cast<double>(dmin) * 60;
	if( min0 + dmin >= 60 ){
		const int ddhour = ( min0 + dmin ) / 60;
		assert( ddhour == 1 );
		dhour += ddhour;
		dmin -= ddhour * 60;
	}

	const int hour1 = ( hour0 + dhour ) % 24;
	const int min1 = min0 + dmin;

	std::ostringstream(oss);
	if( hour1 < 10 ){
		oss << "0";
	}
	oss << hour1 << ":";
	if( min1 < 10 ){
		oss << "0";
	}
	oss << min1 << ":";
	if( sec1 < 10 ){
		oss << "0";
	}
	oss << static_cast<int>(std::floor(sec1));

	return oss.str();

}

// Get number of start times of sections
int Control::getNumStartTimesSections() const{

	return static_cast<int>( m_startTimeOfEachSection.size() );

}

// Get procedure type
int Control::getProcedureType () const{

	return m_procedureType;

}

// Get type of ELOG-Dual
int Control::getTypeOfElogDual() const {
	return m_typeOfElogDual;
}

// Get type of ELOG-MT
int Control::getTypeOfElogMT() const {
	return m_typeOfElogMT;
}

// Set sampling frequency
void Control::setSamplingFrequency( const double samplingFrequency ){
	m_samplingFrequency = samplingFrequency;
}

// Read control parameters from input file
void Control::readParameterFile(){

	const std::string fileName = "param.dat";

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeLogMessage("Read parameters.");

	std::ifstream ifs( fileName.c_str(), std::ios::in );
	if( ifs.fail() ){
		ptrOutputFiles->writeErrorMessage("File open error : " + fileName);
	}

	RobustWeightHuber objRobustWeightHuber;
	RobustWeightThomson objRobustWeightThomson;
	RobustWeightTukeysBiweights objRobustWeightTukeysBiweights;
	int typesOfRobustWeight[2] = { HUBER, TUKEYS_BIWEIGHTS };

	std::string line;
	while( getline( ifs, line ) ){
#ifdef _DEBUG_WRITE
		std::cout << "line : " << line << std::endl;
#endif
		if( line.substr(0,1).compare("#") == 0 ){
			continue;
		}else if( line.find("NUM_OUT") != std::string::npos ){
			ifs >> m_numOutputVariables;
		}
		else if( line.find("NUM_RR") != std::string::npos ){
			ifs >> m_numRemoteReferenceVariables;
		}
		else if( line.find("SAMPLING_FREQ") != std::string::npos ){
			ifs >> m_samplingFrequency;
		}
		else if( line.find("NUM_SECTION") != std::string::npos ){
			ifs >> m_numTimeSeriesSections;
		}
		else if( line.substr(0,7).compare("SEGMENT") == 0 ){
			int numSegments(0);
			ifs >> numSegments;
			for( int iSeg = 0 ; iSeg < numSegments; ++iSeg ){
				SegmentInfo segBuf;
				int ibuf(0);
				ifs >> ibuf;
				if(!Util::isPow2(ibuf)){
					ptrOutputFiles->writeErrorMessage("Length of each segment must be power of two.");
				}
				segBuf.length = ibuf;
				int num(0);
				ifs >> num;
				for( int i = 0 ; i < num; ++i ){
					ifs >> ibuf;
					if( ibuf < 0 ){
						std::ostringstream msg;
						msg << "The frequency index (" << ibuf << ") of the " << iSeg << "-th segment (" << segBuf.length << ") is too small.";
						ptrOutputFiles->writeErrorMessage(msg.str());
					}
					if( ibuf > segBuf.length / 2 ){
						std::ostringstream msg;
						msg << "The " << ibuf << "-th frequency of the " << iSeg << "-th segment (" << segBuf.length << ") is greater than the Nyquists frequency.";
						ptrOutputFiles->writeErrorMessage(msg.str());
					}
					segBuf.degrees.push_back(ibuf);
				}
				m_segments.push_back(segBuf);
			}
		}
		else if( line.find("OVERLAP") != std::string::npos ){
			ifs >> m_overlappingRatio;
		}
		else if( line.find("DATA_FILES") != std::string::npos ){
			if( m_numOutputVariables < 0 ){
				ptrOutputFiles->writeErrorMessage("You must define NOUT data before DATA_FILES data.");
			}
			if( m_numTimeSeriesSections < 0 ){
				ptrOutputFiles->writeErrorMessage("You must define NSECT data before DATA_FILES data.");
			}
			m_dataFileSets.reserve(m_numTimeSeriesSections);
			for( int iSect = 0; iSect < m_numTimeSeriesSections; ++iSect ){
				CommonParameters::DataFileSet dataFileSetTemp;
				ifs >> dataFileSetTemp.numDataPoints;
				const int numChans = getNumberOfChannels();
				dataFileSetTemp.dataFile.reserve(numChans);
				for( int iChan = 0; iChan < numChans; ++iChan ){
					CommonParameters::DataFile dataFileTemp;
					ifs >> dataFileTemp.fileName;
					ifs >> dataFileTemp.numSkipData;
					dataFileTemp.data = NULL;
					dataFileSetTemp.dataFile.push_back(dataFileTemp);
				}
				m_dataFileSets.push_back(dataFileSetTemp);
			}
		}
		else if( line.find("AZIMUTH") != std::string::npos ){
			const int numChannels = getNumberOfChannels();
			m_azimuths.reserve(numChannels);
			for( int iChan = 0; iChan < numChannels; ++iChan ){
				double dbuf(0.0);
				ifs >> dbuf;
				m_azimuths.push_back(dbuf);
			}
		}
		else if( line.find("ROTATION") != std::string::npos ){
			double dbuf(0.0);
			ifs >> dbuf;
			m_rotationAngle = dbuf;
		}
		else if (line.find("MFS_CAL") != std::string::npos) {
			m_calibForMFS = true;
			const int numChannels = getNumberOfChannels();
			m_calibrationFilesForMFS.clear();
			m_calibrationFiles.clear();
			m_calibrationFilesForMFS.reserve(numChannels);
			m_calibrationFiles.reserve(numChannels);
			const Ats* ptrAts = Ats::getInstance();
			for (int iChan = 0; iChan < numChannels; ++iChan) {
				std::string sbuf;
				ifs >> sbuf;
				m_calibrationFilesForMFS.push_back(sbuf);
				m_calibrationFiles.push_back(ptrAts->getCalibrationFileName(iChan));
			}
		}
		else if (line.find("ATS_CAL") != std::string::npos) {
			m_calibForMFS = true;
			const int numChannels = getNumberOfChannels();
			m_calibrationFilesForMFS.clear();
			m_calibrationFiles.clear();
			m_calibrationFilesForMFS.reserve(numChannels);
			m_calibrationFiles.reserve(numChannels);
			Ats* ptrAts = Ats::getInstance();
			ptrAts->setIsCalADUResp(true);
			for (int iChan = 0; iChan < numChannels; ++iChan) {
				std::string sbuf;
				ifs >> sbuf;
				m_calibrationFilesForMFS.push_back(sbuf);
				m_calibrationFiles.push_back(ptrAts->getCalibrationFileName(iChan));
			}
		}
		else if( line.find("ATS_BINARY") != std::string::npos ){
			m_readAtsBinary = true;
		}
		else if (line.find("ELOGDUAL_CAL") != std::string::npos) {
			m_calibForElogDual = true;
			int ibuf(0);
			ifs >> ibuf;
			m_typeOfElogDual = ibuf;
			std::string sbuf;
			ifs >> sbuf;
			m_paramsForElogDualCalibration.fileName = sbuf;
			double dbuf(0.0);
			ifs >> dbuf;
			m_paramsForElogDualCalibration.unitGroupDelay = dbuf;
			}
		else if (line.find("ELOGDUAL_BINARY") != std::string::npos) {
			m_readElogDualBinary = true;
		}
		else if( line.find("ELOGMT_CAL") != std::string::npos ){
			m_calibForElogMT = true;
			int ibuf(0);
			ifs >> ibuf;
			m_typeOfElogMT = ibuf;
			std::string sbuf;
			ifs >> sbuf;
			m_paramsForElogMTCalibration.fileName = sbuf;
			double dbuf(0.0);
			ifs >> dbuf;
			m_paramsForElogMTCalibration.unitGroupDelay = dbuf;
		}
		else if( line.find("ELOGMT_BINARY") != std::string::npos ){
			m_readElogMTBinary = true;
		}
		else if (line.find("ELOGMT_READ_OPTION") != std::string::npos) {
			int ibuf(0);
			ifs >> ibuf;
			if (ibuf < 0 && ibuf > Control::NUM_TYPE_OF_ELOGMT_READING_OPTION) {
				ptrOutputFiles->writeErrorMessage("Unsupported type of ELOG-MT reading option: " + Util::toString(ibuf));
			}
			m_elogMTReadingOption = ibuf;
		}
		else if (line.find("LOGGER_CAL_DIRECTORY") != std::string::npos) {
			std::string sbuf;
			ifs >> sbuf;
			m_directoryOfLoggerCalibrationFiles = sbuf;
		}
		else if( line.find("CAL_FILES") != std::string::npos ){
			const int numChannels = getNumberOfChannels();
			m_calibrationFiles.clear();
			m_calibrationFiles.reserve(numChannels);
			for( int iChan = 0; iChan < numChannels; ++iChan ){
				std::string sbuf;
				ifs >> sbuf;
				m_calibrationFiles.push_back(sbuf);
			}
		}
		else if( line.find("OUTPUT_LEVEL") != std::string::npos ){
			int ibuf(0);
			ifs >> ibuf;
			m_outputLevel = ibuf;
		}
		else if( line.find("DATA_SEGMENT_EVAL") != std::string::npos ){
			m_paramsForDataSegmentEvaluation.doEvaluation = true;
			int ibuf(0);
			ifs >> ibuf;
			m_paramsForDataSegmentEvaluation.numDataSegments = ibuf;
		}
		else if( line.find("FREQ_DOMAIN_EVAL") != std::string::npos ){
			m_paramsForFreqDomainEvaluation.doEvaluation = true;
			double dbuf(0);
			ifs >> dbuf;
			m_paramsForFreqDomainEvaluation.startOutputFrequency = dbuf;
			ifs >> dbuf;
			m_paramsForFreqDomainEvaluation.endOutputFrequency = dbuf;
			int ibuf(0);
			ifs >> ibuf;
			m_paramsForFreqDomainEvaluation.numOfOutputFrequency = ibuf;
			ifs >> ibuf;
			m_paramsForFreqDomainEvaluation.numOfOutputFrequencyForAverage = ibuf;
		}
		else if( line.find("TIME_DOMAIN_EVAL") != std::string::npos ){
			m_paramsForTimeDomainEvaluation.doEvaluation = true;
			int ibuf(0);
			ifs >> ibuf;
			m_paramsForTimeDomainEvaluation.timeSeriesInterval = ibuf;
		}
		else if( line.find("DECIMATION") != std::string::npos ){
			m_paramsForDecimation.applyDecimation = true;
			int ibuf(0);
			ifs >> ibuf;
			m_paramsForDecimation.decimationInterval = ibuf;
			ifs >> ibuf;
			m_paramsForDecimation.filterLength = ibuf;
			double dbuf(0);
			ifs >> dbuf;
			m_paramsForDecimation.transitionBandWidthInLogarithmicScale = dbuf;
		}
		else if( line.find("START_TIMES") != std::string::npos ){
			for( int iSect = 0 ; iSect < m_numTimeSeriesSections; ++iSect ){
				std::string sbuf;
				ifs >> sbuf;
				m_startTimeOfEachSection.push_back(sbuf);
			}
		}
		else if( line.find("COHERENCE_CRITERIA") != std::string::npos && line.length() < 25){
			m_paramsForSquareCoherenceCriteria.applySquareCoherenceCriteria = true;
			int ibuf(0);
			ifs >> ibuf;
			if( ibuf < 2 ){
				ptrOutputFiles->writeErrorMessage("Segment size for square coherence criteria is too small: " + Util::toString(ibuf));
			}
			m_paramsForSquareCoherenceCriteria.numSegments = ibuf;
			double dbuf(0);
			ifs >> dbuf;
			if( dbuf < 0.0 ){
				ptrOutputFiles->writeErrorMessage("Square coherence criteria is less than 0.0: " + Util::toString(dbuf));
			}
			if( dbuf > 1.0 ){
				ptrOutputFiles->writeErrorMessage("Square coherence criteria is larger than 1.0: " + Util::toString(dbuf));
			}
			m_paramsForSquareCoherenceCriteria.threshold = dbuf;
		}
		else if( line.find("COHERENCE_CRITERIA_RANDOM") != std::string::npos ){
			m_paramsForSquareCoherenceCriteriaWithRandomSampling.applySquareCoherenceCriteria = true;
			int ibuf(0);
			ifs >> ibuf;
			if( ibuf < 2 ){
				ptrOutputFiles->writeErrorMessage("Segment size for square coherence criteria is too small: " + Util::toString(ibuf));
			}
			m_paramsForSquareCoherenceCriteriaWithRandomSampling.numSegments = ibuf;
			ifs >> ibuf;
			m_paramsForSquareCoherenceCriteriaWithRandomSampling.numRandomSamples = ibuf;
			double dbuf(0);
			ifs >> dbuf;
			if( dbuf < 0.0 ){
				ptrOutputFiles->writeErrorMessage("Square coherence criteria is less than 0.0: " + Util::toString(dbuf));
			}
			if( dbuf > 1.0 ){
				ptrOutputFiles->writeErrorMessage("Square coherence criteria is larger than 1.0: " + Util::toString(dbuf));
			}
			m_paramsForSquareCoherenceCriteriaWithRandomSampling.threshold = dbuf;
		}
		else if( line.find("HAMPEL_FILTER") != std::string::npos ){
			m_paramsForParamsForHampelFilter.applyHampelFilter = true;
			int ibuf(0);
			ifs >> ibuf;// Number of neighbor points on either side for Hampel filter
			if( ibuf < 1 ){
				ptrOutputFiles->writeErrorMessage("Number of neighbor points on either side for Hampel filter is less than 1: " + Util::toString(ibuf));
			}
			m_paramsForParamsForHampelFilter.numNeighborsOnEitherSide = ibuf;
			double dbuf(0);
			ifs >> dbuf;
			if( dbuf <= 0 ){
				ptrOutputFiles->writeErrorMessage("Criteria for Hampel filter must be positive: " + Util::toString(dbuf));
			}
			m_paramsForParamsForHampelFilter.nsigma = dbuf;
		}
		else if( line.find("MS_CRITERIA") != std::string::npos ){
			m_paramsForMeanSquareCriteria.applyMeanSquareCriteria = true;
			double dbuf(0);
			ifs >> dbuf;
			if( dbuf <= 0 ){
				ptrOutputFiles->writeErrorMessage("Mean square criteria must be positive: " + Util::toString(dbuf));
			}
			m_paramsForMeanSquareCriteria.nsigma = dbuf;
		}
		else if( line.find("PROCEDURE") != std::string::npos ){
			int ibuf(0);
			ifs >> ibuf;
			if( ibuf < 0 || ibuf >= Control::NUM_TYPE_OF_PROCEDURE_TYPE ){
				ptrOutputFiles->writeErrorMessage("Unsupported type of procedure : " + Util::toString(ibuf));
			}
			m_procedureType = ibuf;
		}
		else if( line.find("MESTIMATORS") != std::string::npos ){
			int ibuf(0);
			ifs >> ibuf;
			if( ibuf < Control::NOT_USED || ibuf >= Control::NUM_TYPE_OF_MESTIMATORS ){
				ptrOutputFiles->writeErrorMessage("Unsupported type of M-estimator : " + Util::toString(ibuf));
			}
			typesOfRobustWeight[0] = ibuf;
			ifs >> ibuf;
			if( ibuf < Control::NOT_USED || ibuf >= Control::NUM_TYPE_OF_MESTIMATORS ){
				ptrOutputFiles->writeErrorMessage("Unsupported type of M-estimator : " + Util::toString(ibuf));
			}
			typesOfRobustWeight[1] = ibuf;
		}
		else if( line.find("HUBER") != std::string::npos ){
			double dbuf(0);
			ifs >> dbuf;
			objRobustWeightHuber.setThreshould(dbuf);
			int ibuf(0);
			ifs >> ibuf;
			objRobustWeightHuber.setNumIterationMax(ibuf);
			ifs >> dbuf;
			objRobustWeightHuber.setConvergenceCriteria(dbuf);
		}
		else if( line.find("THOMSON") != std::string::npos ){
			double dbuf(0);
			ifs >> dbuf;
			objRobustWeightThomson.setParameterForDetermingProbability(dbuf);
			int ibuf(0);
			ifs >> ibuf;
			objRobustWeightThomson.setNumIterationMax(ibuf);
			ifs >> dbuf;
			objRobustWeightThomson.setConvergenceCriteria(dbuf);
		}
		else if( line.find("TUKEYS_BIWEIGHTS") != std::string::npos ){
			int ibuf(0);
			ifs >> ibuf;
			objRobustWeightTukeysBiweights.setNumIterationMax(ibuf);
			double dbuf(0);
			ifs >> dbuf;
			objRobustWeightTukeysBiweights.setConvergenceCriteria(dbuf);
		}
		else if( line.find("HAT_MATRIX") != std::string::npos ){
			m_paramsForTreatmentOfHatMatrix.applyLeverageWeights = true;
			double dbuf(0);
			ifs >> dbuf;
			m_paramsForTreatmentOfHatMatrix.threshold = dbuf;
			int ibuf(0);
			ifs >> ibuf;
			if( ibuf < 1 ){
				ptrOutputFiles->writeErrorMessage("Maximum number of the outer iteration is less than 1: " + Util::toString(ibuf));
			}
			m_paramsForTreatmentOfHatMatrix.maxNumberOfOuterIteration = ibuf;
		}
		else if( line.find("OUTPUT_TIME_SERIES") != std::string::npos ){
			m_outputTimeSeriesToCsv = true;
		}
		else if( line.find("JACKKNIFE") != std::string::npos ){
			double dbuf(0);
			ifs >> dbuf;
			m_percentageOfOmmitedDataSubsetDeletionJackknife = dbuf;
		}
		else if( line.find("HIGH_PASS") != std::string::npos ){
			m_doesApplyIIRHighPassFilter = true;
			double dbuf(0.0);
			ifs >> dbuf;
			if( dbuf <= 0.0 ){
				ptrOutputFiles->writeErrorMessage("Cutoff frequency must be positive");
			}
			m_cutoffFrequencyForIIRHighPassFilter = dbuf;
		}
		else if( line.find("LOW_PASS") != std::string::npos ){
			m_doesApplyIIRLowPassFilter = true;
			double dbuf(0.0);
			ifs >> dbuf;
			if( dbuf <= 0.0 ){
				ptrOutputFiles->writeErrorMessage("Cutoff frequency must be positive");
			}
			m_cutoffFrequencyForIIRLowPassFilter = dbuf;
		}
		else if( line.find("ROBUST_AUTO_COVARIANCE") != std::string::npos ){
			// Parameters for robust calculation of autocovariance matrix
			double dbuf(0.0);
			ifs >> dbuf;
			if( dbuf < 0 ){
				ptrOutputFiles->writeErrorMessage("Percentage of the small number added to diagonals of auto-covariance matrix is negative");
			}
			m_paramsForRobustAutoCovariance.percentageOfSmallNumberAddedToDiagonals = dbuf;
			int ibuf(0);
			ifs >> ibuf;
			if( ibuf < 0 ){
				ptrOutputFiles->writeErrorMessage("Maximum number of the iterations of robust auto-covariance calculation is negative");
			}
			m_paramsForRobustAutoCovariance.maxNumOfIterations = ibuf;
			ifs >> dbuf;
			if( dbuf < 0 ){
				ptrOutputFiles->writeErrorMessage("Convergence criteria of the iterations of robust auto-covariance calculation");
			}
			m_paramsForRobustAutoCovariance.convergenceCriteria = dbuf;
		}
		else if( line.find("PREWHITENING") != std::string::npos ){
			// Parameters for robust prewhitening
			m_paramsForPrewhitening.applyPrewhitening = true;
			int ibuf(0);
			ifs >> ibuf;
			if( ibuf >= Control::NUM_TYPE_OF_MESTIMATORS_FOR_PREWHITENING ){
				ptrOutputFiles->writeErrorMessage("Unsupported type of pre-whitening : " + Util::toString(ibuf));
			}
			m_paramsForPrewhitening.typeOfEstimator = ibuf;
			ifs >> ibuf;
			if( ibuf < 0 ){
				ptrOutputFiles->writeErrorMessage("Degree of AR model is negative");
			}
			if( ibuf > TableOfTukeysBiweightParameters::numDimensions ){
				ptrOutputFiles->writeErrorMessage("Degree of AR model is too large (" + Util::toString(ibuf) + ")");
			}
			m_paramsForPrewhitening.maxDegreeOfARModel = ibuf;
			ifs >> ibuf;
			if( ibuf % 2 == 0 ){
				ptrOutputFiles->writeErrorMessage("Numbef of candidates of partial autocorrelation function should be an odd number");
			}
			m_paramsForPrewhitening.numCandidatesOfPartialAutocorrelationFunction = ibuf;
		}
		else if( line.find("ROBUST_FILTER") != std::string::npos ){
			// Parameters for robust filter
			m_paramsForRobustFilter.applyRobustFilter = true;
			int ibuf(0);
			ifs >> ibuf;
			if( ibuf != 0 ){
				m_paramsForRobustFilter.replaceTimeSeriesWithFilteredValues = true;
			}else{
				m_paramsForRobustFilter.replaceTimeSeriesWithFilteredValues = false;
			}
			const int numChannels = getNumberOfChannels();
			m_paramsForRobustFilter.thresholds.clear();
			m_paramsForRobustFilter.thresholds.reserve(numChannels);
			for( int iChan = 0; iChan < numChannels; ++iChan ){
				ThresholdsForRobustFilter thresholds;
				ifs >> thresholds.firstThresholdFactor;
				if( thresholds.firstThresholdFactor < 0 ){
					ptrOutputFiles->writeErrorMessage("Threshold factor must be positive");
				}
				ifs >> thresholds.secondThresholdFactor;
				if( thresholds.secondThresholdFactor < 0 ){
					ptrOutputFiles->writeErrorMessage("Threshold factor must be positive");
				}
				if( thresholds.firstThresholdFactor >= thresholds.secondThresholdFactor ){
					ptrOutputFiles->writeErrorMessage("The second threshold factor must be larger than the first threshould factor");
				}
				ifs >> thresholds.maxNumOfConsecutiveReplacements;
				m_paramsForRobustFilter.thresholds.push_back(thresholds);
			}
		}
		else if( line.find("ROBUST_MULTIVARIATE_REGRESSION") != std::string::npos || line.find("RRMS") != std::string::npos){
			// Parameters for robust multivariate regression
			int ibuf(0);
			ifs >> ibuf;
			if( ibuf == 1 ){
				m_paramsForRobustMultivariateRegression.selectInitialCandidatesByRandomSamplingAtEachFrequency = true;
			}
			ifs >> ibuf;
			if( ibuf < 0 ){
				ptrOutputFiles->writeErrorMessage("Maximum number of initial candidates is negative");
			}
			m_paramsForRobustMultivariateRegression.numOfMaxInitialCandidates = ibuf;
			ifs >> ibuf;
			if( ibuf < 0 ){
				ptrOutputFiles->writeErrorMessage("Maximum number of iteration of the first imporvements is negative");
			}
			m_paramsForRobustMultivariateRegression.numOfMaxIterationsOfFirstIstep = ibuf;
			double dbuf(0.0);
			ifs >> dbuf;
			if( dbuf < 0 ){
				ptrOutputFiles->writeErrorMessage("Convegence criteria of the first imporvements is negative");
			}
			m_paramsForRobustMultivariateRegression.convergenceCriteriaOfFirstIstep = dbuf;
			ifs >> ibuf;
			if( ibuf < 0 ){
				ptrOutputFiles->writeErrorMessage("Maximum number of the candidates of the second imporvements is negative");
			}
			m_paramsForRobustMultivariateRegression.numOfMaxCandidatesOfSecondIstep = ibuf;
			ifs >> ibuf;
			if( ibuf < 0 ){
				ptrOutputFiles->writeErrorMessage("Maximum number of iteration of the second imporvements is negative");
			}
			m_paramsForRobustMultivariateRegression.numOfMaxIterationsOfSecondIstep = ibuf;
			ifs >> dbuf;
			if( dbuf < 0 ){
				ptrOutputFiles->writeErrorMessage("Convegence criteria of the second imporvements is negative");
			}
			m_paramsForRobustMultivariateRegression.convergenceCriteriaOfSecondIstep = dbuf;
		}
		else if( line.find("DIFFERENCE_THRESHOLD") != std::string::npos ){	
			m_paramsForDecidingCandicatesForSubsequentFrequencies.useResponseFunctionsOfPreviousFrequency = false;
			const int numChannels = getNumberOfChannels() - getNumRemoteReferenceVariables();
			m_paramsForDecidingCandicatesForSubsequentFrequencies.thresholdOfDifferences.clear();
			m_paramsForDecidingCandicatesForSubsequentFrequencies.thresholdOfDifferences.reserve(numChannels);
			for( int iChan = 0; iChan < numChannels; ++iChan ){
				ParamsOfThresholdOfDifferenceResponseFuctions params;
				int ibuf(0);
				ifs >> ibuf;
				if( ibuf < 0 && ibuf > Control::NUM_OF_THRESHOLD_TYPE_OF_DIFFERENCE_OF_RESPONSE_FUCTIONS ){
					ptrOutputFiles->writeErrorMessage("Unsupported type of difference threshold of response fuctions: " + Util::toString(ibuf));
				}
				params.type = ibuf;
				double dbuf(0.0);
				ifs >> dbuf;
				if( dbuf < 0 ){
					ptrOutputFiles->writeErrorMessage("Difference threshold is negative");
				}
				params.threshold = dbuf;
				m_paramsForDecidingCandicatesForSubsequentFrequencies.thresholdOfDifferences.push_back(params);
			}
		}
		else if( line.find("TIME_RANGE_OF_CANDIDATES") != std::string::npos ){	
			std::string sbuf;
			ifs >> sbuf;
			m_paramsForRobustMultivariateRegression.startOfTimeRange = sbuf;
			ifs >> sbuf;
			m_paramsForRobustMultivariateRegression.endOfTimeRange = sbuf;
		}
		else if( line.find("BOOTSTRAP") != std::string::npos ){	
			int ibuf(0);
			ifs >> ibuf;
			if( ibuf < 2 ){
				ptrOutputFiles->writeErrorMessage("Number of repetitions of bootstrap is too small: " + Util::toString(ibuf));
			}
			m_numRepetitionsOfBootstrap = ibuf;
		}
		else if( line.find("ERROR_ESTIMATION") != std::string::npos ){	
			int ibuf(0);
			ifs >> ibuf;
			m_errorEstimationMethod = ibuf;
		}
		else if( line.find("OUTPUT_FREQ_DOMAIN_DATA") != std::string::npos ){
			m_outputFreqDomainDataToCsv = true;
		}
		else if( line.find("NUM_THREADS") != std::string::npos ){
			int ibuf(0);
			ifs >> ibuf;
			m_numThreads = ibuf;
		}
		else if( line.find("NOTCH") != std::string::npos && line.length() < 12 ){
			int ibuf(0);
			ifs >> ibuf;
			const int numFreqs = ibuf;
			m_cutoffFrequenciesForNotchFilter.reserve(numFreqs);
			for( int iFreq = 0; iFreq < numFreqs; ++iFreq ){
				double dbuf(0.0);
				ifs >> dbuf;
				if( dbuf <= 0.0 ){
					ptrOutputFiles->writeErrorMessage("Cutoff frequency must be positive");
				}
				m_cutoffFrequenciesForNotchFilter.push_back(dbuf);
			}
		}
		else if( line.find("NOTCH_PARAM_Q") != std::string::npos ){
			double dbuf(0);
			ifs >> dbuf;
			if( dbuf <= 0.0 ){
				ptrOutputFiles->writeErrorMessage("Parameter Q must be positive");
			}
			m_parameterQForNotchFilter = dbuf;
		}
		else if( line.find("DMP_CRITERIA") != std::string::npos ){
			m_paramsForDegreeOfMagneticPolarizationCriteria.applyDegreeOfMagneticPolarizationCriteria = true;
			int ibuf(0);
			ifs >> ibuf;
			if( ibuf < 2 ){
				ptrOutputFiles->writeErrorMessage("Segment size for degree of magnetic polarization criteria is too small: " + Util::toString(ibuf));
			}
			m_paramsForDegreeOfMagneticPolarizationCriteria.numSegments = ibuf;
			double dbuf(0);
			ifs >> dbuf;
			if( dbuf <= 0.0 ){
				ptrOutputFiles->writeErrorMessage("Threshold must be positive");
			}
			m_paramsForDegreeOfMagneticPolarizationCriteria.threshold = dbuf;
		}
		else if( line.find("MPD_CRITERIA") != std::string::npos ){
			m_paramsForMagneticPolarizatitonDirectionCriteria.applyMagneticPolarizatitonDirectionCriteria = true;
			double dbuf(0);
			ifs >> dbuf;
			if( dbuf <= 0.0 ){
				ptrOutputFiles->writeErrorMessage("Threshold must be positive");
			}
			m_paramsForMagneticPolarizatitonDirectionCriteria.threshold = dbuf;
		}
		else if( line.find("OUTPUT_RHOA_PHS") != std::string::npos ){
			m_outputApparentResistivityAndPhase = true;
		}
		else if( line.find("MERGE_SECTIONS") != std::string::npos ){
			int ibuf(0);
			ifs >> ibuf;
			if( ibuf < 1 ){
				ptrOutputFiles->writeErrorMessage("Number of sections after merging is too small: " + Util::toString(ibuf));
			}
			const int numSectionsAfterMerge(ibuf);
			for( int iSection = 0; iSection < numSectionsAfterMerge; ++iSection ){
				int startSectionIndex(-1);
				int endSectionIndex(-1);
				ifs >> startSectionIndex >> endSectionIndex;
				m_rangeOfSectionsForMerge.push_back ( std::make_pair(startSectionIndex, endSectionIndex ) );
			}
		}
		else if( line.find("END") != std::string::npos ){
			break;
		}
	}

	ifs.close();

	if (getNumCalibrationFilesForMFS() > 0) {
		assert(getNumCalibrationFilesForMFS() == getNumberOfChannels());
	}
	if( getNumCalibrationFiles() > 0 ){
		assert( getNumCalibrationFiles() == getNumberOfChannels() );
	}

	const double samplingFrequencyAfterDecimation = m_samplingFrequency / static_cast<double>(m_paramsForDecimation.decimationInterval);

	ptrOutputFiles->writeLogMessage("================================================================================",false);
	ptrOutputFiles->writeLogMessage("Summary of control parameters",false);
	ptrOutputFiles->writeLogMessage("================================================================================",false);
	ptrOutputFiles->writeLogMessage("Number of threads : " + Util::toString(getNumThreads()),false);
#ifdef _USE_OMP
	omp_set_num_threads( getNumThreads() );
#endif
	switch (getProcedureType()){
		case Control::TWO_STAGE:
			ptrOutputFiles->writeLogMessage("Procedure type : two stage procedure", false);
			break;
		case Control::MULTIVARIATE_REGRESSION:
			ptrOutputFiles->writeLogMessage("Procedure type : multivariate regression (RRMS)", false);
			break;
		case Control::ORDINARY_REMOTE_REFERENCE:
			ptrOutputFiles->writeLogMessage("Procedure type : ordinary remote reference", false);
			break;
		case Control::REPEATED_MEDIAN:
			ptrOutputFiles->writeLogMessage("Procedure type : repeated median estimator", false);
			break;
		case Control::TEST:
			ptrOutputFiles->writeLogMessage("Procedure type : test", false);
			break;
		default:
			ptrOutputFiles->writeErrorMessage("Unsupported procedure type : " + Util::toString(getProcedureType()));
			break;
	}
	// Make instance for analysis class
	switch (getProcedureType())
	{
		case Control::TWO_STAGE:
			m_analysis = new AnalysisTwoStage();
			break;
		case Control::MULTIVARIATE_REGRESSION:
			m_analysis = new AnalysisMultivariateRegression();
			break;
		case Control::ORDINARY_REMOTE_REFERENCE:
			m_analysis = new AnalysisOrdinaryRemoteReference();
			break;
		case Control::REPEATED_MEDIAN:
			m_analysis = new AnalysisRepeatedMedian();
			break;
		case Control::TEST:
			m_analysis = new AnalysisTest();
			break;
		default:
			ptrOutputFiles->writeErrorMessage("Unsupported procedure type : " + Util::toString(getProcedureType()));
			break;
	}
	// Output parameters of procedures
	if (getProcedureType() == Control::TWO_STAGE || getProcedureType() == Control::ORDINARY_REMOTE_REFERENCE || getProcedureType() == Control::TEST) {
		for( int index = 0; index < 2; ++index ){
			std::string sbuf;
			if( index == 0 ){
				sbuf = "The first M-estimator : ";
			}else{
				sbuf = "The second M-estimator : ";
			}
			switch (typesOfRobustWeight[index])
			{
			case Control::HUBER:
				m_analysis->setRobustWeightHuber( index, &objRobustWeightHuber );
				sbuf += "Hubur";
				ptrOutputFiles->writeLogMessage(sbuf, false);
				ptrOutputFiles->writeLogMessage("	Threshould value for downweighting : " + Util::toString(objRobustWeightHuber.getThreshould()), false);
				ptrOutputFiles->writeLogMessage("	Maximum number of iteration : " + Util::toString(objRobustWeightHuber.getNumIterationMax()), false);
				ptrOutputFiles->writeLogMessage("	Convergence criteria : " + Util::toString(objRobustWeightHuber.getConvergenceCriteria()), false);
				break;
			case Control::THOMSON:
				m_analysis->setRobustWeightThomson( index, &objRobustWeightThomson );
				sbuf += "Thomson";
				ptrOutputFiles->writeLogMessage(sbuf, false);
				ptrOutputFiles->writeLogMessage("	Parameter for determing the probability used for outlier rejection : " + Util::toString(objRobustWeightThomson.getParameterForDetermingProbability()), false);
				ptrOutputFiles->writeLogMessage("	Maximum number of iteration : " + Util::toString(objRobustWeightThomson.getNumIterationMax()), false);
				ptrOutputFiles->writeLogMessage("	Convergence criteria : " + Util::toString(objRobustWeightThomson.getConvergenceCriteria()), false);
				break;
			case Control::TUKEYS_BIWEIGHTS:
				m_analysis->setRobustWeightTukeysBiweights( index, &objRobustWeightTukeysBiweights );
				sbuf += "Tukey's biweights";
				ptrOutputFiles->writeLogMessage(sbuf, false);
				ptrOutputFiles->writeLogMessage("	Maximum number of iteration : " + Util::toString(objRobustWeightTukeysBiweights.getNumIterationMax()), false);
				ptrOutputFiles->writeLogMessage("	Convergence criteria : " + Util::toString(objRobustWeightTukeysBiweights.getConvergenceCriteria()), false);
				break;
			case Control::NOT_USED:
				// The pointer remains to be NULL
				sbuf += "not to be used";
				ptrOutputFiles->writeLogMessage(sbuf, false);
				break;
			default:
				ptrOutputFiles->writeErrorMessage("Unsupported type of M-estimator : " + Util::toString(typesOfRobustWeight[index]));
				break;
			}
		}
		if (getProcedureType() == Control::TWO_STAGE && m_paramsForTreatmentOfHatMatrix.applyLeverageWeights) {
			ptrOutputFiles->writeLogMessage("Parameters for leverage weights based on hat matrix : ", false);
			ptrOutputFiles->writeLogMessage("     Threshold value for determing leverage points: " + Util::toString(m_paramsForTreatmentOfHatMatrix.threshold), false);
			ptrOutputFiles->writeLogMessage("     Maximum number of the outer iteration: " + Util::toString(m_paramsForTreatmentOfHatMatrix.maxNumberOfOuterIteration), false);
		}
	}
	else if( getProcedureType() == Control::MULTIVARIATE_REGRESSION ){
		if( m_paramsForRobustMultivariateRegression.selectInitialCandidatesByRandomSamplingAtEachFrequency ){
			ptrOutputFiles->writeLogMessage("At each frequency, initial candidates are selected by random sampling", false);
		}
		ptrOutputFiles->writeLogMessage("Maximum number of initial candidates: "
			+ Util::toString(m_paramsForRobustMultivariateRegression.numOfMaxInitialCandidates), false);
		ptrOutputFiles->writeLogMessage("Maximum number of iteration of the first imporvements: "
			+ Util::toString(m_paramsForRobustMultivariateRegression.numOfMaxIterationsOfFirstIstep), false);
		ptrOutputFiles->writeLogMessage("Convegence criteria of the first imporvements: "
			+ Util::toString(m_paramsForRobustMultivariateRegression.convergenceCriteriaOfFirstIstep), false);
		ptrOutputFiles->writeLogMessage("Maximum number of the candidates of the second imporvements: "
			+ Util::toString(m_paramsForRobustMultivariateRegression.numOfMaxCandidatesOfSecondIstep), false);
		ptrOutputFiles->writeLogMessage("Maximum number of iteration of the second imporvements: "
			+ Util::toString(m_paramsForRobustMultivariateRegression.numOfMaxIterationsOfSecondIstep), false);
		ptrOutputFiles->writeLogMessage("Convegence criteria of the second imporvements: "
			+ Util::toString(m_paramsForRobustMultivariateRegression.convergenceCriteriaOfSecondIstep), false);
		ptrOutputFiles->writeLogMessage("Time range for selecting intial candidates: " + 
			m_paramsForRobustMultivariateRegression.startOfTimeRange + " - " +
			m_paramsForRobustMultivariateRegression.endOfTimeRange, false);
		if( !m_paramsForDecidingCandicatesForSubsequentFrequencies.useResponseFunctionsOfPreviousFrequency ){
			const int numChannels = getNumberOfChannels() - getNumRemoteReferenceVariables();
			ptrOutputFiles->writeLogMessage("Threshold type of differences of response functions:", false);
			ptrOutputFiles->writeLogMessage("  Channel#      Type Threshold", false);
			for( int iChan = 0; iChan < numChannels; ++iChan ){
				const ParamsOfThresholdOfDifferenceResponseFuctions& params = m_paramsForDecidingCandicatesForSubsequentFrequencies.thresholdOfDifferences[iChan];
				std::ostringstream oss;
				oss << std::setw(10) << iChan;
				oss << std::setw(10) << params.type;
				oss << std::setw(10) << params.threshold;
				ptrOutputFiles->writeLogMessage(oss.str(), false);
			}
		}
	}
	switch (getErrorEstimationMethod())
	{
	case SUBSET_DELETION_JACKKNIFE:
		ptrOutputFiles->writeLogMessage("Error estimation method : subset deletion jackknife", false);
		ptrOutputFiles->writeLogMessage("Percentage of ommited data in subset deletion jackknife : " + Util::toString(m_percentageOfOmmitedDataSubsetDeletionJackknife), false);
		break;
	case FIXED_WEIGHTS_JACKKNIFE:
		ptrOutputFiles->writeLogMessage("Error estimation method : fixed-weights jackknife", false);
		break;
	case FIXED_WEIGHTS_BOOTSTRAP:
		ptrOutputFiles->writeLogMessage("Error estimation method : fixed-weights bootstrap", false);
		ptrOutputFiles->writeLogMessage("Number or repetitions in bootstrap : " + Util::toString(getNumRepetitionsOfBootstrap()), false);
		break;
	case STRICT_BOOTSTRAP:
		ptrOutputFiles->writeLogMessage("Error estimation method : strict bootstrap", false);
		ptrOutputFiles->writeLogMessage("Number or repetitions in bootstrap : " + Util::toString(getNumRepetitionsOfBootstrap()), false);
		break;
	case PARAMETRIC:
		ptrOutputFiles->writeLogMessage("Error estimation method : parametric", false);
		break;
	default:
		ptrOutputFiles->writeErrorMessage("Unsupported error estimation method : " + Util::toString(getErrorEstimationMethod()));
		break;
	}

	ptrOutputFiles->writeLogMessage("Number of output variables : " + Util::toString(getNumOutputVariables()), false);
	ptrOutputFiles->writeLogMessage("Number of input variables : " + Util::toString(getNumInputVariables()), false);
	ptrOutputFiles->writeLogMessage("Number of remote reference variables : " + Util::toString(getNumRemoteReferenceVariables()), false);
	ptrOutputFiles->writeLogMessage("Sampling frequency (Hz) : " + Util::toString(getSamplingFrequency()), false);
	ptrOutputFiles->writeLogMessage("Number of time-series sections : " + Util::toString(getNumTimeSeriesSections()), false);
	ptrOutputFiles->writeLogMessage("Ratio of overlapping part to whole segment length : " + Util::toString(getOverlappingRatio()), false);
	ptrOutputFiles->writeLogMessage("Output level : " + Util::toString(getOutputLevel()), false);
	if( doesOutputTimeSeriesToCsv() ){
		ptrOutputFiles->writeLogMessage("Output time-series data to csv files", false);
	}
	if( doesOutputFreqDomainDataToCsv() ){
		ptrOutputFiles->writeLogMessage("Output frequency-domain data to csv files", false);
	}
	if( doesOutputApparentResistivityAndPhase() ){
		ptrOutputFiles->writeLogMessage("Output apparent resistivity and phase to a seperate csv file", false);
	}
	ptrOutputFiles->writeLogMessage("Information about the segment lengths and frequencies : ", false);
	{
		ptrOutputFiles->writeLogMessage("  Segment#    Length     Index       Frequency(Hz)         Period(sec)", false);
		int iSeg(0);
		for( std::vector<SegmentInfo>::const_iterator itrSeg = m_segments.begin(); itrSeg != m_segments.end(); ++itrSeg, ++iSeg ){
			for( std::vector<int>::const_iterator itr = itrSeg->degrees.begin(); itr != itrSeg->degrees.end(); ++itr ){
				std::ostringstream oss;
				oss << std::setw(10) << iSeg;
				oss << std::setw(10) << itrSeg->length;
				oss << std::setw(10) << *itr;
				const double freq = samplingFrequencyAfterDecimation * static_cast<double>(*itr) / static_cast<double>(itrSeg->length);
				const double period = 1.0 / freq;
				oss << std::setw(20) << std::scientific << std::setprecision(9) << freq;
				oss << std::setw(20) << std::scientific << std::setprecision(9) << period;
				ptrOutputFiles->writeLogMessage(oss.str(), false);
				m_frequencies.push_back(freq);// Construct array for frequencies
			}
		}
	}
	ptrOutputFiles->writeLogMessage("Information about the time-series data : ", false);
	ptrOutputFiles->writeLogMessage("  Section#  Channel#      Type          NSkip          NData     File", false);
	for( int iSect = 0; iSect < m_numTimeSeriesSections; ++iSect ){
		const int numDataPoints = m_dataFileSets[iSect].numDataPoints;
		std::vector<CommonParameters::DataFile>& dataFileList = m_dataFileSets[iSect].dataFile;
		for( int i = 0; i < m_numOutputVariables; ++i ){
			const int iChan = getChannelIndex( CommonParameters::OUTPUT, i );
			const std::string fileName = dataFileList[iChan].fileName;
			const int numSkipData = dataFileList[iChan].numSkipData;
			std::ostringstream oss;
			oss << std::setw(10) << iSect;
			oss << std::setw(10) << iChan;
			oss << std::setw(10) << "Out"+Util::toString(i);
			oss << std::setw(15) << numSkipData;
			oss << std::setw(15) << numDataPoints;
			oss << std::setw(5) << "";
			oss << fileName;
			ptrOutputFiles->writeLogMessage(oss.str(), false);
		}
		for( int i = 0; i < m_numInputVariables; ++i ){
			const int iChan = getChannelIndex( CommonParameters::INPUT, i );
			const std::string fileName = dataFileList[iChan].fileName;
			const int numSkipData = dataFileList[iChan].numSkipData;
			std::ostringstream oss;
			oss << std::setw(10) << iSect;
			oss << std::setw(10) << iChan;
			oss << std::setw(10) << "Inp"+Util::toString(i);
			oss << std::setw(15) << numSkipData;
			oss << std::setw(15) << numDataPoints;
			oss << std::setw(5) << "";
			oss << fileName;
			ptrOutputFiles->writeLogMessage(oss.str(), false);
		}
		for( int i = 0; i < m_numRemoteReferenceVariables; ++i ){
			const int iChan = getChannelIndex( CommonParameters::REMOTE_REFERENCE, i );
			const std::string fileName = dataFileList[iChan].fileName;
			const int numSkipData = dataFileList[iChan].numSkipData;
			std::ostringstream oss;
			oss << std::setw(10) << iSect;
			oss << std::setw(10) << iChan;
			oss << std::setw(10) << "RR"+Util::toString(i);
			oss << std::setw(15) << numSkipData;
			oss << std::setw(15) << numDataPoints;
			oss << std::setw(5) << "";
			oss << fileName;
			ptrOutputFiles->writeLogMessage(oss.str(), false);
		}
	}
	if( getNumRangeOfSectionsForMerge() > 0 ){
		const int numRangeOfSectionsForMerge = getNumRangeOfSectionsForMerge();
		ptrOutputFiles->writeLogMessage("Number of sections after merging : " + Util::toString(numRangeOfSectionsForMerge), false);
		ptrOutputFiles->writeLogMessage("  Original section#  Merged section#", false);
		int endSectionIndexPre(-1);
		for( int iSectionAfterMerge = 0; iSectionAfterMerge < numRangeOfSectionsForMerge; ++iSectionAfterMerge ){
			const std::pair<int, int> range = getRangeOfSectionsForMerge(iSectionAfterMerge);
			std::ostringstream oss;
			oss << std::setw(19) << range.first << std::setw(17) << range.second;
			ptrOutputFiles->writeLogMessage(oss.str(), false);
			if( range.first < 0 ){
				ptrOutputFiles->writeErrorMessage("Start section index is too small: " + Util::toString(range.first));
			}
			if( range.second < 0 ){
				ptrOutputFiles->writeErrorMessage("End section index is too small: " + Util::toString(range.second));
			}
			if( range.first >= getNumTimeSeriesSections() ){
				ptrOutputFiles->writeErrorMessage("Start section index is too large: " + Util::toString(range.first));
			}
			if( range.second >= getNumTimeSeriesSections() ){
				ptrOutputFiles->writeErrorMessage("End section index is too large: " + Util::toString(range.second));
			}
			if( range.second < range.first ){
				ptrOutputFiles->writeErrorMessage("End section index should be larger than start section index!!");
			}
			if( range.first <= endSectionIndexPre ){
				ptrOutputFiles->writeErrorMessage("Start section index should be greater than the previous end section index.");
			}
			endSectionIndexPre = range.second;
		}
	}
	if( !m_startTimeOfEachSection.empty() ){
		if( getNumStartTimesSections() != getNumTimeSeriesSections() ){
			ptrOutputFiles->writeErrorMessage("Number of start times is not equal to the number of sections");
		}
		ptrOutputFiles->writeLogMessage("Start time of each section : ", false);
		ptrOutputFiles->writeLogMessage("  Section#     Time(hh:mm::ss)", false);
		for( int iSect = 0; iSect < m_numTimeSeriesSections; ++iSect ){
			std::ostringstream oss;
			oss << std::setw(10) << iSect;
			oss << std::setw(20) << m_startTimeOfEachSection[iSect];
			ptrOutputFiles->writeLogMessage(oss.str(), false);
		}
	}else{
		// If start times are not specified, all start times are forced to be 00:00:00
		for( int iSect = 0 ; iSect < m_numTimeSeriesSections; ++iSect ){
			const std::string sbuf = "00:00:00";
			m_startTimeOfEachSection.push_back(sbuf);
		}
	}
	const int numChannels = getNumberOfChannels();
	if( m_azimuths.empty() ){
		ptrOutputFiles->writeErrorMessage("Azimuth is not defined.");
	}
	ptrOutputFiles->writeLogMessage("Rotation angle (deg.) : " + Util::toString(getRotationAngle()), false);
	ptrOutputFiles->writeLogMessage("  Channel#      Type  Azimuth(deg.)", false);
	{
		for( int i = 0; i < m_numOutputVariables; ++i ){
			const int iChan = getChannelIndex( CommonParameters::OUTPUT, i );
			std::ostringstream oss;
			oss << std::setw(10) << iChan;
			oss << std::setw(10) << "Out"+Util::toString(i);
			oss << std::setw(15) << getAzimuth(iChan);
			ptrOutputFiles->writeLogMessage(oss.str(), false);
		}
		for( int i = 0; i < m_numInputVariables; ++i ){
			const int iChan = getChannelIndex( CommonParameters::INPUT, i );
			std::ostringstream oss;
			oss << std::setw(10) << iChan;
			oss << std::setw(10) << "Inp"+Util::toString(i);
			oss << std::setw(15) << getAzimuth(iChan);
			ptrOutputFiles->writeLogMessage(oss.str(), false);
		}
		for( int i = 0; i < m_numRemoteReferenceVariables; ++i ){
			const int iChan = getChannelIndex( CommonParameters::REMOTE_REFERENCE, i );
			std::ostringstream oss;
			oss << std::setw(10) << iChan;
			oss << std::setw(10) << "RR"+Util::toString(i);
			oss << std::setw(15) << getAzimuth(iChan);
			ptrOutputFiles->writeLogMessage(oss.str(), false);
		}
	}
	if( doesReadAtsBinary() ){
		ptrOutputFiles->writeLogMessage("Read binary ATS files", false);
	}
	if (doesMakeCalibrationFileForMFS()) {
		ptrOutputFiles->writeLogMessage("Information about the inputs for calibration files : ", false);
		ptrOutputFiles->writeLogMessage("      Type     Input", false);
		const Ats* ptrAts = Ats::getInstance();
		for (int i = 0; i < m_numOutputVariables; ++i) {
			const int iChan = getChannelIndex(CommonParameters::OUTPUT, i);
			std::ostringstream oss;
			oss << std::setw(10) << "Out" + Util::toString(i) << "     ";
			oss << getCalibrationFileNameForMFS(iChan);
			ptrOutputFiles->writeLogMessage(oss.str(), false);
		}
		for (int i = 0; i < m_numInputVariables; ++i) {
			const int iChan = getChannelIndex(CommonParameters::INPUT, i);
			std::ostringstream oss;
			oss << std::setw(10) << "Inp" + Util::toString(i) << "     ";
			oss << getCalibrationFileNameForMFS(iChan);
			ptrOutputFiles->writeLogMessage(oss.str(), false);
		}
		for (int i = 0; i < m_numRemoteReferenceVariables; ++i) {
			const int iChan = getChannelIndex(CommonParameters::REMOTE_REFERENCE, i);
			std::ostringstream oss;
			oss << std::setw(10) << "RR" + Util::toString(i) << "     ";
			oss << getCalibrationFileNameForMFS(iChan);
			ptrOutputFiles->writeLogMessage(oss.str(), false);
		}
	}
	if (doesReadElogDualBinary()){
		ptrOutputFiles->writeLogMessage("Read two electric field data from ELOG1K/ELOG-Dual binary files", false);
	}
	if (doesMakeCalibrationFileForElogDual()) {
		ptrOutputFiles->writeLogMessage("Information about the inputs for calibration of ELOG-Dual files : ", false);
		switch (getTypeOfElogDual())
		{
			case Control::ELOG1K:
				ptrOutputFiles->writeLogMessage("      Type : ELOG1K" , false);
				break;
			case Control::ELOGDUAL_ADU_MODE:
				ptrOutputFiles->writeLogMessage("      Type : ELOG-Dual (ADU mode)", false);
				break;
			case Control::ELOGDUAL_PHX_MODE:
				ptrOutputFiles->writeLogMessage("      Type : ELOG-Dual (PHX mode)", false);
				break;
			default:
				break;
		}
		ptrOutputFiles->writeLogMessage("      Calibration file : " + Util::toString(m_paramsForElogDualCalibration.fileName), false);
		ptrOutputFiles->writeLogMessage("      Unit group delay : " + Util::toString(m_paramsForElogDualCalibration.unitGroupDelay), false);
	}
	if( doesReadElogMTBinary() ){
		switch (getElogMTReadingOption()){
			case Control::READ_EX_EY_HX_HY_HZ_FROM_ELOGMT_DATA:
				ptrOutputFiles->writeLogMessage("Read Ex, Ey, Hx, Hy, and Hz from ELOG-MT binary files", false);
				break;
			case Control::READ_EX_EY_HX_HY_FROM_ELOGMT_DATA:
				ptrOutputFiles->writeLogMessage("Read Ex, Ey, Hx, and Hy from ELOG-MT binary files", false);
				break;
			case Control::READ_HZ_HX_HY_FROM_ELOGMT_DATA:
				ptrOutputFiles->writeLogMessage("Read Hx, Hy, and Hz from ELOG-MT binary files", false);
				break;
			case Control::READ_EX_EY_FROM_ELOGMT_DATA:
				ptrOutputFiles->writeLogMessage("Read Ex and Ey from ELOG-MT binary files", false);
				break;
			case Control::READ_EX_EY_HX_HY_HZ_HRX_HRY_FROM_ELOGMT_DATA:
				ptrOutputFiles->writeLogMessage("Read Ex, Ey, Hx, Hy, Hz, Hrx, and Hry from ELOG-MT binary files", false);
				break;
			case Control::READ_EX_EY_HX_HY_HRX_HRY_FROM_ELOGMT_DATA:
				ptrOutputFiles->writeLogMessage("Read Ex, Ey, Hx, Hy, Hrx, and Hry from ELOG-MT binary files", false);
				break;
			case Control::READ_HZ_HX_HY_HRX_HRY_FROM_ELOGMT_DATA:
				ptrOutputFiles->writeLogMessage("Read Hx, Hy, Hz, Hrx, and Hry from ELOG-MT binary files", false);
				break;
			case Control::READ_HX_HY_HRX_HRY_FROM_ELOGMT_DATA:
				ptrOutputFiles->writeLogMessage("Read Hx, Hy, Hrx, and Hry from ELOG-MT binary files", false);
				break;
			case Control::READ_HX_HY_FROM_ELOGMT_DATA:
				ptrOutputFiles->writeLogMessage("Read Hx and Hy from ELOG-MT binary files", false);
				break;
			default:
				ptrOutputFiles->writeErrorMessage("Unsupported type of ELOG-MT reading option: " + Util::toString(getElogMTReadingOption()));
				break;
		}
	}
	if( doesMakeCalibrationFileForElogMT() ){
		ptrOutputFiles->writeLogMessage("Information about the inputs for calibration of ELOG-MT files : ", false);
		switch (getTypeOfElogMT())
		{
		case Control::ELOGMT_ADU_MODE:
			ptrOutputFiles->writeLogMessage("      Type : ELOG-MT (ADU mode)", false);
			break;
		case Control::ELOGMT_PHX_MODE:
			ptrOutputFiles->writeLogMessage("      Type : ELOG-MT (PHX mode)", false);
			break;
		default:
			break;
		}
		ptrOutputFiles->writeLogMessage("      Calibration file : " + Util::toString(m_paramsForElogMTCalibration.fileName), false);
		ptrOutputFiles->writeLogMessage("      Unit group delay : " + Util::toString(m_paramsForElogMTCalibration.unitGroupDelay), false);
	}
	if( !m_calibrationFiles.empty() ){
		assert( m_calibrationFiles.size() == getNumberOfChannels() );
		ptrOutputFiles->writeLogMessage("Information about calibration files : ", false);
		ptrOutputFiles->writeLogMessage("      Type     Calibration file", false);
		{
			for( int i = 0; i < m_numOutputVariables; ++i ){
				const int iChan = getChannelIndex( CommonParameters::OUTPUT, i );
				std::ostringstream oss;
				oss << std::setw(10) << "Out"+Util::toString(i);
				oss << std::setw(5) << "";
				oss << getCalibrationFileName(iChan);
				ptrOutputFiles->writeLogMessage(oss.str(), false);
			}
			for( int i = 0; i < m_numInputVariables; ++i ){
				const int iChan = getChannelIndex( CommonParameters::INPUT, i );
				std::ostringstream oss;
				oss << std::setw(10) << "Inp"+Util::toString(i);
				oss << std::setw(5) << "";
				oss << getCalibrationFileName(iChan);
				ptrOutputFiles->writeLogMessage(oss.str(), false);
			}
			for( int i = 0; i < m_numRemoteReferenceVariables; ++i ){
				const int iChan = getChannelIndex( CommonParameters::REMOTE_REFERENCE, i );
				std::ostringstream oss;
				oss << std::setw(10) << "RR"+Util::toString(i);
				oss << std::setw(5) << "";
				oss << getCalibrationFileName(iChan);
				ptrOutputFiles->writeLogMessage(oss.str(), false);
			}
		}
	}
	if (!getDirectoryOfLoggerCalibrationFiles().empty()) {
		ptrOutputFiles->writeLogMessage("Directory of logger calibration files : " + getDirectoryOfLoggerCalibrationFiles(), false);
	}
	if( m_paramsForParamsForHampelFilter.applyHampelFilter ){
		ptrOutputFiles->writeLogMessage("Parameters for Hampel filter : ", false);
		ptrOutputFiles->writeLogMessage("     Number of neighbor points on either side for Hampel filter: " + Util::toString(m_paramsForParamsForHampelFilter.numNeighborsOnEitherSide), false);
		ptrOutputFiles->writeLogMessage("     Criteria for Hampel filter: " + Util::toString(m_paramsForParamsForHampelFilter.nsigma), false);
	}
	if( doesApplyIIRHighPassFilter() ){
		ptrOutputFiles->writeLogMessage("Parameters for high-pass filter : ", false);
		ptrOutputFiles->writeLogMessage("     Cutoff frequency (Hz): " + Util::toString(getCutoffFrequencyForIIRHighPassFilter()), false);
	}
	if( doesApplyIIRLowPassFilter() ){
		ptrOutputFiles->writeLogMessage("Parameters for low-pass filter : ", false);
		ptrOutputFiles->writeLogMessage("     Cutoff frequency (Hz): " + Util::toString(getCutoffFrequencyForIIRLowPassFilter()), false);
	}
	if( !m_cutoffFrequenciesForNotchFilter.empty() ){
		ptrOutputFiles->writeLogMessage("Parameters for notch filter : ", false);
		ptrOutputFiles->writeLogMessage("     Q: " + Util::toString(m_parameterQForNotchFilter), false);
		std::ostringstream oss;
		oss << "     Cutoff frequencies (Hz): ";
		for( std::vector<double>::const_iterator itr = m_cutoffFrequenciesForNotchFilter.begin(); 
			itr != m_cutoffFrequenciesForNotchFilter.end(); ++itr ){
			oss << *itr << " ";
		}
		ptrOutputFiles->writeLogMessage(oss.str(), false);
	}
	if( m_paramsForDecimation.applyDecimation ){
		ptrOutputFiles->writeLogMessage("Parameters for decimation : ", false);
		ptrOutputFiles->writeLogMessage("     Decimation interval: " + Util::toString(m_paramsForDecimation.decimationInterval), false);
		ptrOutputFiles->writeLogMessage("     Filter length: " + Util::toString(m_paramsForDecimation.filterLength), false);
		ptrOutputFiles->writeLogMessage("     Transition band width in logarithmic scale: " + Util::toString(m_paramsForDecimation.transitionBandWidthInLogarithmicScale), false);
	}
	if( m_paramsForPrewhitening.applyPrewhitening ){
		ptrOutputFiles->writeLogMessage("Parameters for robust prewhitening : ", false);
		switch (m_paramsForPrewhitening.typeOfEstimator)
		{
			case Control::USE_LEAST_SQUARE_ESTIMATOR_FOR_PREWHITENING:
				ptrOutputFiles->writeLogMessage("     Least square estimator is used", false);
				ptrOutputFiles->writeLogMessage("     Maximum degree of AR model: " + Util::toString(m_paramsForPrewhitening.maxDegreeOfARModel), false);
				break;
			case Control::USE_S_ESTIMATOR_FOR_PREWHITENING:
				ptrOutputFiles->writeLogMessage("     S-estimator is used", false);
				ptrOutputFiles->writeLogMessage("     Maximum degree of AR model: " + Util::toString(m_paramsForPrewhitening.maxDegreeOfARModel), false);
				ptrOutputFiles->writeLogMessage("     Numbef of candidates of partial autocorrelation function: " + Util::toString(m_paramsForPrewhitening.numCandidatesOfPartialAutocorrelationFunction), false);
				break;
			case Control::MANUAL_INPUTS_AR_COEFFICIENTS_FOR_PREWHITENING:
				ptrOutputFiles->writeLogMessage("     AR coefficients are specified by an input file", false);
				break;
			default:
				ptrOutputFiles->writeErrorMessage("Unsupported type of estimator for prewhitening: " + Util::toString(m_paramsForPrewhitening.typeOfEstimator));
				break;
		}
	}
	if( m_paramsForRobustFilter.applyRobustFilter ){
		ptrOutputFiles->writeLogMessage("Parameters for robust filter : ", false);
		if( m_paramsForRobustFilter.replaceTimeSeriesWithFilteredValues ){
			ptrOutputFiles->writeLogMessage("     Time-series data are repladed with filtered values", false);
		}
		ptrOutputFiles->writeLogMessage("       Channel#   The 1st threshold   The 2nd threshold  Max consecutive replacements", false);
		if( m_paramsForRobustFilter.thresholds.size() != numChannels ){
			ptrOutputFiles->writeErrorMessage("Number of thresholds for robust filter (" + Util::toString(m_paramsForRobustFilter.thresholds.size()) 
				+ ") is not equal to the number of channel (" + Util::toString(numChannels) + ")" );
		}
		for( int iChan = 0; iChan < numChannels; ++iChan ){
			const ThresholdsForRobustFilter& thresholds = m_paramsForRobustFilter.thresholds[iChan];
			std::ostringstream oss;
			oss << "     ";
			oss << std::setw(10) << iChan;
			oss << std::setw(20) << thresholds.firstThresholdFactor;
			oss << std::setw(20) << thresholds.secondThresholdFactor;
			oss << std::setw(30) << thresholds.maxNumOfConsecutiveReplacements;
			ptrOutputFiles->writeLogMessage(oss.str(), false);
		}
		ptrOutputFiles->writeLogMessage("Parameters for robust auto-covariance calculation : ", false);
		ptrOutputFiles->writeLogMessage("     Percentage of the small number added to diagonals: " + Util::toString(m_paramsForRobustAutoCovariance.percentageOfSmallNumberAddedToDiagonals), false);
		ptrOutputFiles->writeLogMessage("     Maximum number of the iterations: " + Util::toString(m_paramsForRobustAutoCovariance.maxNumOfIterations), false);
		ptrOutputFiles->writeLogMessage("     Convergence criteria of the iterations: " + Util::toString(m_paramsForRobustAutoCovariance.convergenceCriteria), false);
	}
	if( m_paramsForMeanSquareCriteria.applyMeanSquareCriteria ){
		ptrOutputFiles->writeLogMessage("Mean square criteria : " + Util::toString(m_paramsForMeanSquareCriteria.nsigma), false);
	}
	if( m_paramsForSquareCoherenceCriteria.applySquareCoherenceCriteria ){
		ptrOutputFiles->writeLogMessage("Parameters for square coherence criteria : ", false);
		ptrOutputFiles->writeLogMessage("     Number of segments: " + Util::toString(m_paramsForSquareCoherenceCriteria.numSegments), false);
		ptrOutputFiles->writeLogMessage("     Square coherence criteria: " + Util::toString(m_paramsForSquareCoherenceCriteria.threshold), false);
	}
	if( m_paramsForSquareCoherenceCriteriaWithRandomSampling.applySquareCoherenceCriteria ){
		ptrOutputFiles->writeLogMessage("Parameters for square coherence criteria with random sampling: ", false);
		ptrOutputFiles->writeLogMessage("     Number of segments: " + Util::toString(m_paramsForSquareCoherenceCriteriaWithRandomSampling.numSegments), false);
		ptrOutputFiles->writeLogMessage("     Number of random samples: " + Util::toString(m_paramsForSquareCoherenceCriteriaWithRandomSampling.numRandomSamples), false);
		ptrOutputFiles->writeLogMessage("     Square coherence criteria: " + Util::toString(m_paramsForSquareCoherenceCriteriaWithRandomSampling.threshold), false);
	}
	if( m_paramsForDegreeOfMagneticPolarizationCriteria.applyDegreeOfMagneticPolarizationCriteria ){
		ptrOutputFiles->writeLogMessage("Parameters for degree of magnetic polarizatiton criteria: ", false);
		ptrOutputFiles->writeLogMessage("     Number of segments: " + Util::toString(m_paramsForDegreeOfMagneticPolarizationCriteria.numSegments), false);
		ptrOutputFiles->writeLogMessage("     Threshould: " + Util::toString(m_paramsForDegreeOfMagneticPolarizationCriteria.threshold), false);
	}
	if( m_paramsForMagneticPolarizatitonDirectionCriteria.applyMagneticPolarizatitonDirectionCriteria ){
		ptrOutputFiles->writeLogMessage("Parameters for magnetic polarizatiton direction criteria: ", false);
		ptrOutputFiles->writeLogMessage("     Threshould: " + Util::toString(m_paramsForMagneticPolarizatitonDirectionCriteria.threshold), false);
	}
	if( m_paramsForDataSegmentEvaluation.doEvaluation ){
		ptrOutputFiles->writeLogMessage("Parameters for evaluations of data segments: ", false);
		ptrOutputFiles->writeLogMessage("     Number of data segments: " + Util::toString(m_paramsForDataSegmentEvaluation.numDataSegments), false);
	}
	if( m_paramsForFreqDomainEvaluation.doEvaluation ){
		if( m_paramsForFreqDomainEvaluation.startOutputFrequency < CommonParameters::EPS ){
			ptrOutputFiles->writeErrorMessage("Start of output frequency (" + Util::toString(m_paramsForFreqDomainEvaluation.startOutputFrequency)
				+ ") is too small" );
		}
		if( m_paramsForFreqDomainEvaluation.endOutputFrequency > 0.5 * samplingFrequencyAfterDecimation ){
			ptrOutputFiles->writeWarningMessage("End of output frequency (" + Util::toString(m_paramsForFreqDomainEvaluation.endOutputFrequency)
				+ ") is changed to be the Nyquist frequency ("+ Util::toString(0.5 * samplingFrequencyAfterDecimation) +")" );
			m_paramsForFreqDomainEvaluation.endOutputFrequency = 0.5 * samplingFrequencyAfterDecimation;
		}
		ptrOutputFiles->writeLogMessage("Parameters for prior evaluations of frequency-domain characteristics: ", false);
		ptrOutputFiles->writeLogMessage("     Start of output frequency: " + Util::toString(m_paramsForFreqDomainEvaluation.startOutputFrequency), false);
		ptrOutputFiles->writeLogMessage("     End of output frequency: " + Util::toString(m_paramsForFreqDomainEvaluation.endOutputFrequency), false);
		ptrOutputFiles->writeLogMessage("     Number of output frequencies: " + Util::toString(m_paramsForFreqDomainEvaluation.numOfOutputFrequency), false);
		ptrOutputFiles->writeLogMessage("     Number of output frequencies for average values: " + Util::toString(m_paramsForFreqDomainEvaluation.numOfOutputFrequencyForAverage), false);
	}
	if( m_paramsForTimeDomainEvaluation.doEvaluation ){
		ptrOutputFiles->writeLogMessage("Parameters for prior evaluations of time-domain characteristics: ", false);
		ptrOutputFiles->writeLogMessage("     Time-series interval: " + Util::toString(m_paramsForTimeDomainEvaluation.timeSeriesInterval), false);
	}
	ptrOutputFiles->writeLogMessage("================================================================================",false);

}
