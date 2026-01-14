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
#ifndef DBLDEF_CONTROL
#define DBLDEF_CONTROL

#include <vector>
#include "CommonParameters.h"
#include "Analysis.h"

// Class of the control parameters
class Control{

public:

	enum ElogDualType {
		ELOG1K = 0,
		ELOGDUAL_ADU_MODE,
		ELOGDUAL_PHX_MODE,
		NUM_TYPE_OF_ELOGDUAL_TYPE,
	};

	enum ElogMTType {
		ELOGMT_ADU_MODE,
		ELOGMT_PHX_MODE,
		NUM_TYPE_OF_ELOGMT_TYPE,
	};

	enum ElogMTReadingOption {
		NOT_SPECIFIED = -1,
		READ_EX_EY_HX_HY_HZ_HRX_HRY_FROM_ELOGMT_DATA = 0,
		READ_EX_EY_HX_HY_HRX_HRY_FROM_ELOGMT_DATA,
		READ_HZ_HX_HY_HRX_HRY_FROM_ELOGMT_DATA,
		READ_EX_EY_HX_HY_HZ_FROM_ELOGMT_DATA,
		READ_EX_EY_HX_HY_FROM_ELOGMT_DATA,
		READ_HZ_HX_HY_FROM_ELOGMT_DATA,
		READ_EX_EY_FROM_ELOGMT_DATA,
		READ_HX_HY_HRX_HRY_FROM_ELOGMT_DATA,
		READ_HX_HY_FROM_ELOGMT_DATA,
		NUM_TYPE_OF_ELOGMT_READING_OPTION,
	};

	enum ErrorEstimationMethod{
		PARAMETRIC = 0,
		FIXED_WEIGHTS_BOOTSTRAP,
		STRICT_BOOTSTRAP,
		ROBUST_BOOTSTRAP,
		FIXED_WEIGHTS_JACKKNIFE,
		SUBSET_DELETION_JACKKNIFE,
		NUM_TYPE_OF_ERROR_ESTIMATION_METHOD,
	};

	enum ProcedureType{
		ORDINARY_REMOTE_REFERENCE = 0,
		MULTIVARIATE_REGRESSION,
		MODIFIED_MULTIVARIATE_REGRESSION,
		TWO_STAGE,
		REPEATED_MEDIAN,
		TEST,
		NUM_TYPE_OF_PROCEDURE_TYPE,
	};

	enum TypeOfEstimatorForPrewhitening{
		MANUAL_INPUTS_AR_COEFFICIENTS_FOR_PREWHITENING = -1,
		USE_LEAST_SQUARE_ESTIMATOR_FOR_PREWHITENING = 0,
		USE_S_ESTIMATOR_FOR_PREWHITENING,
		NUM_TYPE_OF_MESTIMATORS_FOR_PREWHITENING,
	};

	enum TypeOfRobustWeight {
		NOT_USED = -1,
		HUBER = 0,
		TUKEYS_BIWEIGHTS,
		THOMSON,
		NUM_TYPE_OF_MESTIMATORS,
	};

	// Type of the threshold of the difference of response fuctions 
	enum ThresholdTypeOfDifferenceOfResponseFuctions{
		DIFFERENCE_OF_RESPONSES_DIVIDED_BY_SQUARE_ROOT_FREQUENCY,
		DIFFERENCE_OF_RAW_RESPONSES,
		NUM_OF_THRESHOLD_TYPE_OF_DIFFERENCE_OF_RESPONSE_FUCTIONS
	};

	// Segment length and the index of the frequency where response functions are estimated
	enum TimingEOFBasedDenoising{
		BEFORE_DECIMATION = 0 ,
		AFTER_DEGITAL_FILTERS,
	};

	// Segment length and the index of the frequency where response functions are estimated
	struct SegmentInfo{
		int length;
		std::vector<int> degrees;
	};

	// Parameters for ELOG calibration
	struct ParamsForElogCalibration{
		std::string fileName;
		double unitGroupDelay;
	};

	// Parameters for the data-segment evaluation prior to the calculation of the response functions
	struct ParamsForDataSegmentEvaluation{
		bool doEvaluation;
		int numDataSegments;
	};

	// Parameters for decimation
	struct ParamsForDecimation{
		bool applyDecimation;
		int decimationInterval;
		int filterLength;
		double transitionBandWidthInLogarithmicScale;
	};

	// Parameters for the frequency-domain evaluation prior to the calculation of the response functions
	struct ParamsForFreqDomainEvaluation{
		bool doEvaluation;
		double startOutputFrequency;
		double endOutputFrequency;
		int numOfOutputFrequency;
		int numOfOutputFrequencyForAverage;
	};

	// Parameters for the time-domain evaluation prior to the calculation of the response functions
	struct ParamsForTimeDomainEvaluation{
		bool doEvaluation;
		int timeSeriesInterval;
	};

	// Parameters for Hampel filter
	struct ParamsForHampelFilter{
		bool applyHampelFilter;
		int numNeighborsOnEitherSide;
		double nsigma;
	};

	// Parameters for mean square criteria
	struct ParamsForMeanSquareCriteria{
		bool applyMeanSquareCriteria;
		double nsigma;
	};

	// Parameters for square coherence criteria
	struct ParamsForSquareCoherenceCriteria{
		bool applySquareCoherenceCriteria;
		int numSegments;
		double threshold;
	};

	// Parameters for square coherence criteria with random sampling
	struct ParamsForSquareCoherenceCriteriaWithRandomSampling{
		bool applySquareCoherenceCriteria;
		int numSegments;
		int numRandomSamples;
		double threshold;
	};

	// Parameters for robust calculation of autocovariance matrix
	struct ParamsForRobustAutoCovariance{
		double percentageOfSmallNumberAddedToDiagonals;
		int maxNumOfIterations;
		double convergenceCriteria;
	};

	// Parameters of the threshold of the difference of response fuctions 
	struct ParamsOfThresholdOfDifferenceResponseFuctions{
		int type;
		double threshold;
	};

	// Parameters for prewhitening
	struct ParamsForPrewhitening{
		bool applyPrewhitening;
		int typeOfEstimator;
		int maxDegreeOfARModel;
		int numCandidatesOfPartialAutocorrelationFunction;
	};

	// Thresholds for robust filter
	struct ThresholdsForRobustFilter{
		double firstThresholdFactor;
		double secondThresholdFactor;
		int maxNumOfConsecutiveReplacements;
	};

	// Parameters for robust filter
	struct ParamsForRobustFilter{
		bool applyRobustFilter;
		bool replaceTimeSeriesWithFilteredValues;
		std::vector<ThresholdsForRobustFilter> thresholds;
	};

	// Parameters for robust multivariate regression
	struct ParamsForRobustMultivariateRegression{
		bool selectInitialCandidatesByRandomSamplingAtEachFrequency;
		int numOfMaxInitialCandidates;
		int numOfMaxIterationsOfFirstIstep;
		double convergenceCriteriaOfFirstIstep;
		int numOfMaxCandidatesOfSecondIstep;
		int numOfMaxIterationsOfSecondIstep;
		double convergenceCriteriaOfSecondIstep;
		std::string startOfTimeRange;
		std::string endOfTimeRange;
	};

	// Parameters for deciding candicates for subsequent frequencies
	struct ParamsForDecidingCandicatesForSubsequentFrequencies{
		bool useResponseFunctionsOfPreviousFrequency;
		std::vector<ParamsOfThresholdOfDifferenceResponseFuctions> thresholdOfDifferences;
	};

	// Parameters for degree of polarizatiton direction criteria
	struct ParamsForDegreeOfMagneticPolarizationCriteria{
		bool applyDegreeOfMagneticPolarizationCriteria;
		int numSegments;
		double threshold;
	};

	// Parameters for magnetic polarizatiton direction criteria
	struct ParamsForMagneticPolarizatitonDirectionCriteria{
		bool applyMagneticPolarizatitonDirectionCriteria;
		double threshold;
	};

	// Parameters for the treatment of hat matrix
	struct ParamsForTreatmentOfHatMatrix{
		bool applyLeverageWeights;
		double threshold;
		int maxNumberOfOuterIteration;
	};

	// Return the the instance of the class
    static Control* getInstance();

	// Run analysis
	void run(const bool outputToConsole);

	// Get flag specifing whether IIR high-pass filter is applied
	bool doesApplyIIRHighPassFilter() const; 

	// Get flag specifing whether IIR low-pass filter is applied
	bool doesApplyIIRLowPassFilter() const; 

	// Get flag specifing whether calibration for MFS is performed
	bool doesMakeCalibrationFileForMFS() const;

	// Get flag specifing whether calibration for ELOG-Dual file is performed
	bool doesMakeCalibrationFileForElogDual() const;

	// Get flag specifing whether calibration for ELOG-MT file is performed
	bool doesMakeCalibrationFileForElogMT() const; 

	// Get flag specifing whether apparent resistivity and phase are outputed in a seperate file
	bool doesOutputApparentResistivityAndPhase() const;

	// Get flag specifing wheter output frequency domain data as csv file
	bool doesOutputFreqDomainDataToCsv() const;

	// Get flag specifing wheter output time series data as csv file
	bool doesOutputTimeSeriesToCsv() const;

	// Get flag specifing whether input file is ATS binary file
	bool doesOutputCalibratedTimeSeriesToCsv() const;

	// Get flag specifing whether input file is ATS binary file
	bool doesReadAtsBinary() const;

	// Get flag specifing whether input file is ELOG-Dual binary file
	bool doesReadMTH5() const;

	// Get flag specifing whether input file is ELOG-Dual binary file
	bool doesReadElogDualBinary() const;

	// Get flag specifing whether input file is ELOG-MT binary file
	bool doesReadElogMTBinary() const;

	// Get azimuth
	bool doesPeformEOFBasedDenoising() const;

	// Get azimuth
	double getAzimuth( const int iChan ) const;

	// Get channel index
	int getChannelIndex( const int dataType, const int index ) const;

	// Get option of ELOG-Dual binary data reading
	int getElogDualReadingOption() const;

	// Get option of ELOG-MT binary data reading
	int getElogMTReadingOption() const;

	// Get error estimation method
	int getErrorEstimationMethod() const;

	// Get number of output variables
	int getNumOutputVariables() const;

	// Get number of input variables
	int getNumInputVariables() const;

	// Get number of remote reference variables
	int getNumRemoteReferenceVariables() const;

	// Get sampling frequency
	double getSamplingFrequency() const;

	// Get number of threads
	double getSamplingFrequencyOrg() const;

	// Get number of threads
	int getNumThreads() const;

	// Get number of time-series sections
	int getNumTimeSeriesSections() const;

	// Get number of segment lengths
	int getNumSegmentLengths() const;

	// Get number of target frequency degrees in an input segment
	int getNumTargetFrequencyInSegment( const int iSeg ) const;

	// Get number of the ranges of sections for merge
	int getNumRangeOfSectionsForMerge() const;

	// Get number of repetitions of bootstrap method
	int getNumRepetitionsOfBootstrap() const;

	// Get target frequency degrees in an input segment
	int getTargetFrequencyDegreesInSegment( const int iSeg, const int index ) const;

	// Get length of each segment
	int getSegmentLength( const int iSeg ) const;

	// Get ratio of overlapping part to whole segment length
	double getOverlappingRatio() const;

	// Get rotation angle
	double getRotationAngle() const;

	// Get numebur of calibration files for MFS
	int getNumCalibrationFilesForMFS() const;

	// Get name of calibration file for MFS
	std::string getCalibrationFileNameForMFS(const int iFile) const;

	// Get numebur of calibration files
	int getNumCalibrationFiles() const;

	// Get name of calibration file
	std::string getCalibrationFileName( const int iFile ) const;

	// Get cutoff frequency for IIR high-pass filter
	double getCutoffFrequencyForIIRHighPassFilter() const;

	// Get cutoff frequency for IIR low-pass filter
	double getCutoffFrequencyForIIRLowPassFilter() const;

	// Get cutoff frequencies for notch filter
	std::vector<double> getCutoffFrequenciesForNotchFilter() const;

	// Get directory of logger calibration files
	std::string getDirectoryOfLoggerCalibrationFiles() const;

	// Get number of ichannels
	int getNumberOfChannels() const;

	// Get number of cutoff frequencies for notch filter
	int getNumberOfCutoffFrequenciesForNotchFilter() const;

	// Get number of frequencies
	int getNumberOfFrequencies() const;

	// Get frequency
	double getFrequency( const int iFreq ) const;

	// Get all frequencies without duplications
	std::vector<double> getFrequenciesAllWithoutDuplications() const;

	// Get output level
	int getOutputLevel() const;

	// Get parameter Q for notch filter
	double getParameterQForNotchFilter() const;

	// Get parameters for ELOG-Dual calibration
	ParamsForElogCalibration getParamsForElogDualCalibration() const;

	// Get parameters for ELOG-MT calibration
	ParamsForElogCalibration getParamsForElogMTCalibration() const;

	// Get parameters for the data-segment evaluation prior to the calculation of the response functions
	ParamsForDataSegmentEvaluation getParamsForDataSegmentEvaluation() const;

	// Get parameters for decimation
	ParamsForDecimation getParamsForDecimation() const;

	// Get parameters for the frequency-domain evaluation prior to the calculation of the response functions
	ParamsForFreqDomainEvaluation getParamsForFreqDomainEvaluation() const;

	// Get parameters for the time-domain evaluation prior to the calculation of the response functions
	ParamsForTimeDomainEvaluation getParamsForTimeDomainEvaluation() const;

	// Get parameters for Hampel filter
	ParamsForHampelFilter getParamsForHampelFilter() const;

	// Get parameters for mean square criteria
	ParamsForMeanSquareCriteria getParamsForMeanSquareCriteria() const;

	// Get parameters for degree of polarizatiton direction criteria
	ParamsForDegreeOfMagneticPolarizationCriteria getParamsForDegreeOfMagneticPolarizationCriteria() const;

	// Get parameters for magnetic polarizatiton direction criteria
	ParamsForMagneticPolarizatitonDirectionCriteria getParamsForMagneticPolarizatitonDirectionCriteria() const;

	// Get parameters for square coherence criteria
	ParamsForSquareCoherenceCriteria getParamsForSquareCoherenceCriteria() const;

	// Get parameters for square coherence criteria with random sampling
	ParamsForSquareCoherenceCriteriaWithRandomSampling getParamsForSquareCoherenceCriteriaWithRandomSampling() const;

	// Get parameters for robust calculation of autocovariance matrix
	ParamsForRobustAutoCovariance getParamsForRobustAutoCovariance() const;

	// Get parameters for robust prewhitening
	ParamsForPrewhitening getParamsForPrewhitening() const;

	// Get parameters for robust filter
	ParamsForRobustFilter getParamsForRobustFilter() const;

	// Get parameters for robust multivariate regression
	ParamsForRobustMultivariateRegression getParamsForRobustMultivariateRegression() const;

	// Get parameters for deciding candicates for subsequent frequencies
	ParamsForDecidingCandicatesForSubsequentFrequencies getParamsForDecidingCandicatesForSubsequentFrequencies() const;

	// Get parameters for the treatment of hat matrix
	ParamsForTreatmentOfHatMatrix getParamsForTreatmentOfHatMatrix() const;

	// Get percentage of ommited data in subset deletion jackknife
	double getPercentageOfOmmitedDataSubsetDeletionJackknife() const;

	// Get ranges of sections for merge
	std::pair<int, int> getRangeOfSectionsForMerge( const int iSectionAfterMerge ) const;

	// Get time from start time and elapsed time (seconds)
	std::string getTimeFromStartTimeOfEachSection( const int sectionIndex, const double elapsedTime, const bool forAfterMergingSections = true ) const;

	// Get number of start times of sections
	int getNumStartTimesSections() const;

	// Get procedure type
	int getProcedureType() const;

	// Get timing of denoising based on EOF
	int getTimingEOFBasedDenoising() const;

	// Get type of ELOG-Dual
	int getTypeOfElogDual() const;

	// Get type of ELOG-MT
	int getTypeOfElogMT() const;

	// Set sampling frequency
	void setSamplingFrequency( const double samplingFrequency );

private:

	// Constructer
	Control();

	// Destructer
	~Control();

	// Copy constructer
	Control(const Control& rhs);

	// Assignment operator
	Control& operator=(const Control& rhs);

	// Azimuths
	std::vector<double> m_azimuths;

	// Cutoff frequency for IIR high-pass filter
	double m_cutoffFrequencyForIIRHighPassFilter;

	// Cutoff frequency for IIR low-pass filter
	double m_cutoffFrequencyForIIRLowPassFilter;

	// Cutoff frequencies for notch filter
	std::vector<double> m_cutoffFrequenciesForNotchFilter;

	// Directory of logger calibration files
	std::string m_directoryOfLoggerCalibrationFiles;

	// Flag specifing whether IIR high-pass filter is applied
	bool m_doesApplyIIRHighPassFilter;

	// Flag specifing whether IIR low-pass filter is applied
	bool m_doesApplyIIRLowPassFilter;

	// Option of ELOG-Dual binary data reading
	bool m_doesPeformEOFBasedDenoising;

	// Option of ELOG-Dual binary data reading
	int m_elogDualReadingOption;

	// Option of ELOG-MT binary data reading
	int m_elogMTReadingOption;

	// Error estimation method
	int m_errorEstimationMethod;

	// Number of output variables
	int m_numOutputVariables;

	// Number of input variables
	int m_numInputVariables;

	// Number of remote reference variables
	int m_numRemoteReferenceVariables;

	// Sampling frequency
	double m_samplingFrequencyOrg;

	// Sampling frequency
	double m_samplingFrequency;

	// Number of threads
	int m_numThreads;

	// Number of time-series sections
	int m_numTimeSeriesSections;

	// Number of repetitions of bootstrap method
	int m_numRepetitionsOfBootstrap;

	// Array for segment length and the index of the frequency where response functions are estimated
	std::vector<SegmentInfo> m_segments;

	// Ratio of overlapping part to whole segment length
	double m_overlappingRatio;

	// Frequencies where response functions are estimated
	std::vector<double> m_frequencies;

	// Flag specifing whether input calibration for MFS is performed
	bool m_calibForMFS;

	// Flag specifing whether input calibration for ELOG-Dual file is performed
	bool m_calibForElogDual;

	// Flag specifing whether input calibration for ELOG-MT file is performed
	bool m_calibForElogMT;

	// Pointer to Analysis class
	Analysis* m_analysis;

	// Data file sets
	std::vector<CommonParameters::DataFileSet> m_dataFileSets;

	// Rotation angles
	double m_rotationAngle;

	// Calibration files
	std::vector<std::string> m_calibrationFiles;

	// Calibration file for MFS
	std::vector<std::string> m_calibrationFilesForMFS;

	// Percentage of ommited data in subset deletion jackknife
	double m_percentageOfOmmitedDataSubsetDeletionJackknife;

	// Start time of each section
	std::vector<std::string> m_startTimeOfEachSection;

	// Flag specifing whether apparent resistivity and phase are outputed in a seperate file
	bool m_outputApparentResistivityAndPhase;

	// Output level
	int m_outputLevel;

	// Flag specifing wheter output frequency domain data as csv file
	bool m_outputFreqDomainDataToCsv;

	// Flag specifing wheter output time series data as csv file
	bool m_outputTimeSeriesToCsv;

	// Flag specifing wheter output calibrated time series data as csv file
	bool m_outputCalibratedTimeSeriesToCsv;

	// Parameter Q for notch filter
	double m_parameterQForNotchFilter;

	// Parameters for ELOGDual calibration
	ParamsForElogCalibration m_paramsForElogDualCalibration;

	// Parameters for ELOG-MT calibration
	ParamsForElogCalibration m_paramsForElogMTCalibration;

	// Parameters for the data-segment evaluation prior to the calculation of the response functions
	ParamsForDataSegmentEvaluation m_paramsForDataSegmentEvaluation;

	// Parameters for decimation
	ParamsForDecimation m_paramsForDecimation;

	// Parameters for the frequency-domain evaluation prior to the calculation of the response functions
	ParamsForFreqDomainEvaluation m_paramsForFreqDomainEvaluation;

	// Parameters for the time-domain evaluation prior to the calculation of the response functions
	ParamsForTimeDomainEvaluation m_paramsForTimeDomainEvaluation;

	// Parameters for square coherence criteria
	ParamsForSquareCoherenceCriteria m_paramsForSquareCoherenceCriteria;

	// Parameters for square coherence criteria with random sampling
	ParamsForSquareCoherenceCriteriaWithRandomSampling m_paramsForSquareCoherenceCriteriaWithRandomSampling;

	// Parameters for Hampel filter
	ParamsForHampelFilter m_paramsForParamsForHampelFilter;

	// Parameters for mean square criteria
	ParamsForMeanSquareCriteria m_paramsForMeanSquareCriteria;

	// Parameters for degree of polarizatiton direction criteria
	ParamsForDegreeOfMagneticPolarizationCriteria m_paramsForDegreeOfMagneticPolarizationCriteria;

	// Parameters for magnetic polarizatiton direction criteria
	ParamsForMagneticPolarizatitonDirectionCriteria m_paramsForMagneticPolarizatitonDirectionCriteria;

	// Parameters for robust calculation of autocovariance matrix
	ParamsForRobustAutoCovariance m_paramsForRobustAutoCovariance;

	// Parameters for robust prewhitening
	ParamsForPrewhitening m_paramsForPrewhitening;

	// Parameters for robust filter
	ParamsForRobustFilter m_paramsForRobustFilter;

	// Parameters for robust multivariate regression
	ParamsForRobustMultivariateRegression m_paramsForRobustMultivariateRegression;

	// Parameters for deciding candicates for subsequent frequencies
	ParamsForDecidingCandicatesForSubsequentFrequencies m_paramsForDecidingCandicatesForSubsequentFrequencies;

	// Parameters for the treatment of hat matrix
	ParamsForTreatmentOfHatMatrix m_paramsForTreatmentOfHatMatrix;

	// Procedure type
	int m_procedureType;

	// Ranges of sections for merge
	std::vector< std::pair<int, int> > m_rangeOfSectionsForMerge;

	// Flag specifing whether input file is ATS binary file
	bool m_readAtsBinary;

	// Flag specifing whether input file is MTH5 file
	bool m_readMTH5;

	// Flag specifing whether ELOG-Dual binary is read
	bool m_readElogDualBinary;

	// Flag specifing whether ELOG-MT binary is read
	bool m_readElogMTBinary;

	// Timing of denoising based on EOF
	int m_timingEOFBasedDenoising;

	// Type of ELOG-Dual
	int m_typeOfElogDual;

	// Type of ELOG-MT
	int m_typeOfElogMT;

	// Read control parameters from input file
	void readParameterFile();

};

#endif
