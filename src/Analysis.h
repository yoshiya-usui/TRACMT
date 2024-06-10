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
#ifndef DBLDEF_ANALYSIS
#define DBLDEF_ANALYSIS

#include "CommonParameters.h"
#include "CalibrationFunction.h"
#include "RobustWeightHuber.h"
#include "RobustWeightThomson.h"
#include "RobustWeightTukeysBiweights.h"
#include <vector>

// Class of analysis
class Analysis{

public:

	// Default constructer
	Analysis();

	// Destructer
	virtual ~Analysis();

	// Exexute analysis
	void run( std::vector<CommonParameters::DataFileSet>& dataFileSets );

	// Set object of M-estimator with Huber weight
	void setRobustWeightHuber( const int index, const RobustWeightHuber* const ptrRobustWeight );

	// Set object of M-estimator with Thomson weight
	void setRobustWeightThomson( const int index, const RobustWeightThomson* const ptrRobustWeight );

	// Set object of M-estimator with Tukey's biweights
	void setRobustWeightTukeysBiweights( const int index, const RobustWeightTukeysBiweights* const ptrRobustWeight );

	// Test function
	void test() const;

	// Test function 2
	void test2( std::vector<CommonParameters::DataFileSet>& dataFileSets );

	// Test function 3
	void test3() const;

protected:

	// Calculate apparent resistivity
	double calcApparentResistivity( const double freq, const std::complex<double>& Z ) const;

	// Calculate error of apparent resistivity
	double calcApparentResistivityError( const double freq, const std::complex<double>& Z, const double dZ ) const;

	// Calculate diagonal components of hat matrix
	double calculateDiagonalComponentsOfHatMatrix( const int numSegments, const int channelX, const int channelY,
		std::complex<double>** data, const double* const weights, double* hatDiagonals ) const;

	// Calculate phase
	double calcPhase( const std::complex<double>& Z ) const;

	// Calculate error of phase
	double calcPhaseError( const std::complex<double>& Z, const double dZ ) const;

	// Calculate response functions by iteratively reweighted least squares
	void calculateResponseFunctionsByIRWLS( const int iRobustWeight, const std::complex<double>* const out, const std::complex<double>* const in0,
		const std::complex<double>* const in1, const int numSegments, const bool fixScale, double& scale, const double* const weightsPrior,
		double* weights, std::complex<double>* residuals, std::complex<double>& resp0, std::complex<double>& resp1, double& coherence,
		std::vector<std::string>& titles, std::vector<double>* outputValues, const bool priorityOnFirst = true ) const;

	// Calculate response functions by iteratively reweighted remote reference
	void calculateResponseFunctionsByIRWLSRemoteReference(const int iRobustWeight, const std::complex<double>* const out, 
		const std::complex<double>* const in0, const std::complex<double>* const in1, 
		const std::complex<double>* const rr1, const std::complex<double>* const rr2,
		const int numSegments, const bool fixScale, double& scale, const double* const weightsPrior,
		double* weights, std::complex<double>* residuals, std::complex<double>& resp0, std::complex<double>& resp1, double& coherence,
		std::vector<std::string>& titles, std::vector<double>* outputValues) const;

	// Calculate response function by the orinary remote reference method
	double calculateResponseFunctionByOrdinaryRemoteReference( const int numSegments, const std::complex<double>* const out,
		const std::complex<double>* const in1, const std::complex<double>* const in2, 
		const std::complex<double>* const rr1, const std::complex<double>* const rr2, 
		std::complex<double>& resp1, std::complex<double>& resp2 ) const;

	// Calculate response function by the ordinary weighted square method
	double calculateResponseFunctionByWLS( const std::complex<double>* const out, const std::complex<double>* const in0,
		const std::complex<double>* const in1, const int numSegments, const double* const weights, std::complex<double>* residuals,
		std::complex<double>& resp0, std::complex<double>& resp1, double& coherence, const bool priorityOnFirst = true ) const;

	// Auxiliary function for calculating response function by the weighted leaset square method
	void calculateResponseFunctionByWLSAux( const int numSegments, const std::complex<double>* const out,
		const std::complex<double>* const in1, const std::complex<double>* const in2, const double* const weights, 
		std::complex<double>& resp1, std::complex<double>& resp2, const bool priorityOnFirst = true ) const;

	// Calculate response function by the weighted leaset square method for bootstrap
	void calculateResponseFunctionByWLSForBootstrap( const int numSegments, const int* const segmentIndexes,
		const std::complex<double>* const out, const std::complex<double>* const in1, const std::complex<double>* const in2,
		const double* const weights, std::complex<double>& resp1, std::complex<double>& resp2, const bool priorityOnFirst = true ) const;

	// Calculate response function by the weighted remote reference method
	double calculateResponseFunctionByWLSRemoteReference(const std::complex<double>* const out,
		const std::complex<double>* const in1, const std::complex<double>* const in2,
		const std::complex<double>* const rr1, const std::complex<double>* const rr2,
		const int numSegments, const double* const weights, std::complex<double>* residuals, 
		std::complex<double>& resp1, std::complex<double>& resp2, double& coherence ) const;

	// Auxiliary function for calculating response function by the weighted remote reference method
	void calculateResponseFunctionByWLSRemoteReferenceAux(const int numSegments, const std::complex<double>* const out,
		const std::complex<double>* const in1, const std::complex<double>* const in2,
		const std::complex<double>* const rr1, const std::complex<double>* const rr2,
		const double* const weights, std::complex<double>& resp1, std::complex<double>& resp2) const;

	// Calculate response function by the weighted remote reference method for bootstrap
	void calculateResponseFunctionByWLSRemoteReferenceForBootstrap(const int numSegments, const int* const segmentIndexes,
		const std::complex<double>* const out,
		const std::complex<double>* const in1, const std::complex<double>* const in2,
		const std::complex<double>* const rr1, const std::complex<double>* const rr2,
		const double* const weights, std::complex<double>& resp1, std::complex<double>& resp2) const;

	// Get pointer to M-estimators
	RobustWeight* getPointerToRobustWeight(const int iRobustWeight) const;

	// Output spectral density functions to cvg file
	void outputSpectralDensityFunctionsToCvgFile( const int numSegments, const double timeLength, const std::complex<double>* const out,
		const std::complex<double>* const in1, const std::complex<double>* const in2, const double* const weights ) const;

	// Output spectral density functions to cvg file
	void outputSpectralDensityFunctionsToCvgFile( const int numSegments, const double timeLength, const std::complex<double>* const out,
		const std::complex<double>* const in1, const std::complex<double>* const in2, const std::complex<double>* const rr1, const std::complex<double>* const rr2 ) const;

	// Merge sections
	void mergeSections( std::vector<CommonParameters::DataFileSet>& dataFileSets ) const;

private:
	
	// Copy constructer
	Analysis(const Analysis& rhs);

	// Assignment operator
	Analysis& operator=(const Analysis& rhs);

	// Decimation
	void decimation( std::vector<CommonParameters::DataFileSet>& dataFileSets ) const;

	// Read time-series data
	void readTimeSeriesData( std::vector<CommonParameters::DataFileSet>& dataFileSets );

	// Read calibration files
	void readCalibrationFiles( const std::vector<double>& freq );

	// Read one time-series data
	void readOneTimeSeriesData( const std::string& fileName, const int numSkipData, const int numDataPoints, double* data );
	
	// Perform preprocessing
	void preprocessing( std::vector<CommonParameters::DataFileSet>& dataFileSets );

	// Convert time-series data to frequency-domain data
	void convertToFrequencyData( const int segmentLength, const std::vector<CommonParameters::DataFileSet>& dataFileSets,
		int& numSegmentsTotal, std::complex<double>*** cdata, std::vector< std::pair<std::string, std::string> >& times, std::vector<double>* meanSquares );

	// Perform calibration correction for all channels
	void calibrationCorrectionAllChannels( const int numChannels, const int numSegmentsTotal, const double freq, std::complex<double>** ftval ) const;

	// Perform calibration correction for main analyss
	void calibrationCorrection( const int iChan, const int numSegmentsTotal, const double freq, std::complex<double>* ftval, const bool afterPreprocessing = true ) const;

	// Calculate rotated fields
	void calculateRotatedFields( const int numSegmentsTotal, std::complex<double>** ftval ) const;

	virtual void calculateResponseFunctions( const int iSegLen, const int freqDegree, const double timeLength, const double freq,
		const int numSegmentsTotal, std::complex<double>** ftval, const std::vector< std::pair<std::string, std::string> >& times,
		std::ofstream& ofsResp, std::ofstream& ofsRhoaPhs ) = 0;

	// Output average spectrum
	void outputAverageSpectrum( std::complex<double>** cdata, const int numData, const int section, 
		const bool afterPreprocessing, const bool afterCalibration ) const;

	// Output average spectrum for the frequencies where response functions are estimated
	void outputAverageSpectrum2( std::complex<double>** cdata, const int numData, const int section, 
		const bool afterPreprocessing, const bool afterCalibration ) const;

	// Output spectrum
	void outputSpectrum( std::complex<double>** cdata, const int numData, const int section, 
		const bool afterPreprocessing, const bool afterCalibration ) const;

	// Output frequency-domain data
	void outputFrequencyDomainData( const int iSegLen, const int freqDegree, const int numSegmentsTotal, std::complex<double>** ftval ) const;

	// Output time-series data
	void outputTimeSeriesData( const std::vector<CommonParameters::DataFileSet>& dataFileSets, const bool afterPreprocessing ) const;

	// Evaluate characteristics of time-series data prior to the estimation of the response functions
	void priorEvaluationOfTimeSeriesData( const int interval, const std::vector<CommonParameters::DataFileSet>& dataFileSets,
		const bool afterPreprocessing ) const;

	// Evaluate characteristics of frequency data from all data prior to the estimation of the response functions
	void priorEvaluationOfFrequencyDataFromAllData( const std::vector<CommonParameters::DataFileSet>& dataFileSets,
		const bool afterPreprocessing ) const;

	// Evaluate characteristics of data segments
	void priorEvaluationOfDataSegments( const int numSegmentsTotal, std::complex<double>** ftval,
		const std::vector< std::pair<std::string, std::string> >& times, const double timeLength, const std::string& fileName ) const;

	// Auxiliary function for evaluating characteristics of data segments
	void priorEvaluationOfDataSegmentsAux( const std::vector<int>& segmentIndexes, std::complex<double>** ftval,
		const std::vector< std::pair<std::string, std::string> >& times, const double timeLength, const int numSegments,
		std::ofstream& ofs ) const;

	// Select segments to be excluded by degree of magnetic polarization criteria
	void selectSegmentsToBeExcludedByDegreeOfMagneticPolarizationCriteria( const int numSegmentsTotal, std::complex<double>** ftval, std::vector<bool>& remainingSegments ) const;

	// Select segments to be excluded by mean square criteria
	void selectSegmentsToBeExcludedByMeanSquareCriteria( const int numSegmentsTotal, const std::vector<double>* const meanSquares, 
		const std::vector< std::pair<std::string, std::string> >& times, std::vector<bool>& remainingSegments ) const;

	// Select segments to be excluded by magnetic polarizatiton direction criteria
	void selectSegmentsToBeExcludedByMagneticPolarizatitonDirectionCriteria( const int numSegmentsTotal, std::complex<double>** ftval, std::vector<bool>& remainingSegments ) const;

	// Select segments to be excluded by square coherence criteria
	void selectSegmentsToBeExcludedBySquareCoherenceCriteria( const int numSegmentsTotal, std::complex<double>** ftval,
		std::vector<bool>& remainingSegments ) const;

	// Auxiliary function for selecting segments to be excluded by square coherence criteria
	double selectSegmentsToBeExcludedBySquareCoherenceCriteriaAux( const int numSegments, std::complex<double>** data ) const;

	// Select segments to be excluded by square coherence criteria with random sampling
	void selectSegmentsToBeExcludedBySquareCoherenceCriteriaWithRandomSampling( const int numSegmentsTotal, std::complex<double>** ftval,
		std::vector<bool>& remainingSegments ) const;

	// Auxiliary function for selecting segments to be excluded by square coherence criteria with random sampling
	bool selectSegmentsToBeExcludedBySquareCoherenceCriteriaWithRandomSamplingAux( const std::vector<int>& segmentIndexes, std::complex<double>** ftval,
		const std::vector< std::pair<std::string, std::string> >& times, const double timeLength ) const;

	// Write header to the output file for apparent resistivity and phase
	void writeHeaderToOutputFileForApparentResistivityAndPhase( std::ofstream& ofs ) const;

	// Write header to the output file of response function
	void writeHeaderToOutputFileForResponseFunctions( std::ofstream& ofs ) const;

	// Arrays of calibration function
	CalibrationFunction* m_calibrationFunctions;
	
	// Coefficients of AR model
	std::vector<double>* m_coefficientsOfARModel;

	// M-estimators
	RobustWeight* m_robustWeight[2];

};

#endif
