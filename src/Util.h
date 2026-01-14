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
#ifndef DBLDEF_UTIL
#define DBLDEF_UTIL

#include <complex>
#include <vector>
#include <sstream>

namespace Util
{

// Apply IIR high-pass filter
void applyIIRHighPassFilter( const double samplingFreq, const double cutoffFreq, const int numData, double* data );

// Apply IIR low-pass filter
void applyIIRLowPassFilter( const double samplingFreq, const double cutoffFreq, const int numData, double* data );

// Apply notch filter
void applyNotchFilter( const double Q, const double samplingFreq, const double cutoffFreq, const int numData, double* data );

// Convert string "HH:MM:SS" to seconds
int convertHHMMSSToSeconds( const std::string& hhmmss );

// Debug write real matrix
// Matrix must be stored by column major order
void debugWriteRealMatrix( const int numRows, const int numColumns, const double* const matrix );

// Sort elements by its key value with quick sort
void quickSort( const int numOfIDs, int* ids, const double* const values );

// Sort elements by its key value with quick sort and replace original arrays
void quickSort( const int numOfIDs, double* values );

// Calculate Absolute value fo a complex value
double calculateAbsoluteValue( const std::complex<double> value );

// Calculate calibration function for FIR filter (Type1)
std::complex<double> calculateCalibrationForFIRFilterType1( const std::string& fileName, const double samplingFreq, const double freq, const bool groupDelay );

// Calculate calibration function for FIR filter (Type2)
std::complex<double> calculateCalibrationForFIRFilterType2( const std::string& fileName, const double samplingFreq, const double freq, const bool isELOG, 
	int nfstart, int nfend );

// Calculate determinant of real square matrix
double calculateDeterminantOfMatrix( const int dimension, const double* const matrix );

// Calcualte all eigenvalues and eigenvectors of a real symmetric
void calculateEigenValuesAndVectorsOfRealSymmetricMatrix( const int dimension, const double* const matrix, 
	double* eigenValues, double* eigenVectors );

// Calculate IQR
double calculateIQR( const int num, const double* const data );

// Calculate MADN
double calculateMADN( const int num, const double* const data );

// Calculate MADN
double calculateMADN( const std::vector<double>& data );

// calculate mean value
double calculateMeanValue( const int num, const double* const data );

// calculate mean square value
double calculateMeanSquareValue( const int num, const double* const data );

// Calculate median
double calculateMedian( const int num, const double* const data );

// Calculate median
double calculateMedian( const std::vector<double>& data );

// Calculate median absolute deviation
double calculateMAD(const int num, const double median, const double* const data);

// Calculate median absolute deviation
double calculateMAD( const int num, const double* const data );

// Calculate median absolute deviation
double calculateMAD( const std::vector<double>& data );

// Calculate regression coefficients by weighted least square method
double calculateRegressionCoefficientsByWLS( const int numOfData, const double* const y, const double* const x,
	const double* const weights );

// Calculate regression coefficients by weighted least square method with intercept
void calculateRegressionCoefficientsByWLSWithIntercept(const int numOfData, const double* const y, const double* const x,
	const double* const weights, double& slope, double& intercept);

// calculate variance
double calculateVariance( const int num, const double mean, const double* const data );

// calculate linear trend factor by Hino (1985)
double calculateLinearTrendFactorByHino1985( const int num, const double* const data );

// calculate linear trend factor by the least square
void calculateLinearTrendFactorByLeastSquare( const int num, const double* const data, double& b0, double& b1 );

// Apply Hampel filter
int hampelFilter( const int num, const int numNeighborsOnEitherSide, const double nsigma, double* data );

// Apply Hanning window
void hanningWindow( const int num, double* data );

// Perform FFT
void fft( const int num, std::complex<double>* data, const int isign );

//  Fourier transform
void fourierTransform( const int num, std::complex<double>* data );

// Inverse fourier transform
void inverseFourierTransform( const int num, std::complex<double>* data );

// Calculate field at an angle from those of two different directions
std::complex<double> calculateRotatedField( const double d1, const double d2, const double rot, const std::complex<double>& v1, const std::complex<double>& v2 );

// Check wether the input value is power of two
bool isPow2( const int val );

// Interpolation by the algorithm of Akima (1970)
double interpolationAkima(  const int num, const double* const dataX, const double* const dataY, const double x );

// Calculate slope t for the interpolation by the algorithm of Akima (1970)
double calculateSlopeForInterpolationAkima( const double* const dataX, const double* const dataY );

// Linear interpolation
double interpolationLinear( const double x1, const double x2, const double y1, const double y2, const double x );

// Calculate FIR filter coefficients by the least square method
void calculateFIRFilterCoeffsByLeastSquare( const int dimension, const bool isLowPass, const double samplingFreq,
	const double freq1, const double freq2, const double weight1, const double weight2, double* coeff );

// Calculate frequency characteristics of FIR filter
std::complex<double> calculateFrequencyCharacteristicsOfFIRFilter( const int dimension, const double samplingFreq, const double freq, const double* coeff );

// Calculate frequency characteristics of IIR high-pass filter
std::complex<double> calculateFrequencyCharacteristicsOfIIRHighPassFilter( const double freq, const double samplingFrequency, const double cutoffFreq );

// Calculate frequency characteristics of IIR low-pass filter
std::complex<double> calculateFrequencyCharacteristicsOfIIRLowPassFilter( const double freq, const double samplingFrequency, const double cutoffFreq );

// Calculate frequency characteristics of notch filter
std::complex<double> calculateFrequencyCharacteristicsOfNotchFilter( const double freq, const double Q, const double samplingFrequency, const double cutoffFreq );

// Factrize and solve a linear equation with real coefficents matrix
void factorizeAndSolveLinearEquationRealMatrix( const int dimension, const int nRhs, const double* const matrix, const double* const rhsVectors, double* result );

// Factrize and solve a linear equation with real symmetric matrix
void factorizeAndSolveLinearEquationRealSymmetricMatrix( const int dimension, const int nRhs, const double* const matrix, const double* const rhsVectors, double* result );

// Factrize and solve a linear equation with real symmetric positive definite matrix
void factorizeAndSolveLinearEquationRealSymmetricPositiveDefiniteMatrix( const int dimension, const int nRhs, const double* const matrix, const double* const rhsVectors, double* result );

// Factorize a real square matrix
void factorizeRealSquareMatrix( const int dimension, const double* const matrix, double* factorizedMatrix, int* ipiv );

// Factorize a real symmetric matrix
void factorizeRealSymmetricMatrix( const int dimension, const double* const matrix, double* factorizedMatrix, int* ipiv );

// Factorize a real symmetric positive definite matrix
void factorizeRealSymmetricPositiveDefiniteMatrix( const int dimension, const double* const matrix, double* factorizedMatrix );

// Solve a linear equation with real square matrix
void solveLinearEquationRealSquareMatrix( const int dimension, const int nRhs, const int* const ipivInt, double* factorizedMatrix, const double* const rhsVectors, double* result );

// Solve a linear equation with real symmetric matrix
void solveLinearEquationRealSymmetricMatrix( const int dimension, const int nRhs, const int* const ipivInt, double* factorizedMatrix, const double* const rhsVectors, double* result );

// Solve a linear equation with real symmetric positive definite matrix
void solveLinearEquationRealSymmetricPositiveDefiniteMatrix( const int dimension, const int nRhs, double* factorizedMatrix, const double* const rhsVectors, double* result );

// Extract filename without extension
std::string extractFileNameWithoutExtension( const std::string& fileNameFull );

// Extract extension of filename
std::string extractExtensionOfFileName(const std::string& fileNameFull);

// The second degree Lagrange interpolation
double interpolation2ndOrderLagrange( const double x1, const double x2, const double x3, const double y1, const double y2, const double y3, const double x );

// Convert a decimal number to a binary-coded form
int decimalToBinary(int dec );

// Convert to string
template < typename T >
std::string toString( const T x ){

	std::ostringstream oss;
	oss << x;
	return oss.str();

}

// Convert string to integer
int stringToInt( const std::string& sbuf );

// Convert string to double
double stringToDouble( const std::string& sbuf );

}

#endif
