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
#include "Util.h"
#include "OutputFiles.h"
#include "CommonParameters.h"
#include "LapackInterface.h"
#include <algorithm>
#include <complex>
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <vector>
#include <assert.h>

// Apply IIR high-pass filter
void Util::applyIIRHighPassFilter( const double samplingFreq, const double cutoffFreq, const int numData, double* data ){

	double* dataOrg = new double[numData];
	const double lamda = tan(CommonParameters::PI * cutoffFreq / samplingFreq);
	const double beta  = 1.0 / (lamda + 1.0);
	const double alpha = (lamda - 1.0) / (lamda + 1.0);
	for( int i = 0; i < numData; ++i ){
		dataOrg[i] = data[i];
	}
	data[0] = dataOrg[0];
	double yPre1 = dataOrg[0];
	for( int i = 1; i < numData; ++i ){
		data[i] = beta * dataOrg[i] - beta * dataOrg[i-1] - alpha * yPre1;
		yPre1 = data[i];
	}
	delete [] dataOrg;

}

// Apply IIR low-pass filter
void Util::applyIIRLowPassFilter( const double samplingFreq, const double cutoffFreq, const int numData, double* data ){

	double* dataOrg = new double[numData];
	const double Q = 1.0/sqrt(2.0);
	const double lamda = tan(CommonParameters::PI * cutoffFreq / samplingFreq);
	const double delta = lamda * lamda + lamda / Q + 1.0;
	const double beta = lamda * lamda / delta;
	const double alpha1 = 2.0 * (lamda * lamda - 1.0) / delta;
	const double alpha2 = 1.0 - 2.0 * lamda / Q / delta;
	for( int i = 0; i < numData; ++i ){
		dataOrg[i] = data[i];
	}
	data[0] = dataOrg[0];
	data[1] = dataOrg[1];
	double yPre2 = dataOrg[0];// 2-2=0
	double yPre1 = dataOrg[1];// 2-1=1
	for( int i = 2; i < numData; ++i ){
		data[i] = beta * dataOrg[i] + 2.0 * beta * dataOrg[i-1] + beta * dataOrg[i-2] - alpha1 * yPre1 - alpha2 * yPre2;
		yPre2 = yPre1;
		yPre1 = data[i];
	}
	delete [] dataOrg;

}

// Apply notch filter
void Util::applyNotchFilter( const double Q, const double samplingFreq, const double cutoffFreq, const int numData, double* data ){

	double* dataOrg = new double[numData];
	const double lamda = tan(CommonParameters::PI * cutoffFreq / samplingFreq);
	const double delta = lamda * lamda + lamda / Q + 1.0;
	const double beta0 =       (lamda * lamda + 1.0) / delta;
	const double beta1 = 2.0 * (lamda * lamda - 1.0) / delta;
	const double beta2 = beta0;
	const double alpha1 = beta1;
	const double alpha2 = 1.0 - 2.0 * lamda / Q / delta;
	for( int i = 0; i < numData; ++i ){
		dataOrg[i] = data[i];
	}
	data[0] = dataOrg[0];
	data[1] = dataOrg[1];
	double yPre2 = dataOrg[0];// 2-2=0
	double yPre1 = dataOrg[1];// 2-1=1
	for( int i = 2; i < numData; ++i ){
		data[i] = beta0 * dataOrg[i] + beta1 * dataOrg[i-1] + beta2 * dataOrg[i-2] - alpha1 * yPre1 - alpha2 * yPre2;
		yPre2 = yPre1;
		yPre1 = data[i];
	}
	delete [] dataOrg;

}

// Convert string "HH:MM:SS" to seconds
int Util::convertHHMMSSToSeconds( const std::string& hhmmss ){

	//const int hour = stoi( hhmmss.substr(0,2) );
	//const int min = stoi( hhmmss.substr(3,2) );
	//const int sec = stoi( hhmmss.substr(6,2) );

	std::istringstream issHour(hhmmss.substr(0,2));
	int hour(0);
	issHour >> hour;

	std::istringstream issMin(hhmmss.substr(3,2));
	int min(0);
	issMin >> min;

	std::istringstream issSec(hhmmss.substr(6,2));
	int sec(0);
	issSec >> sec;

	return hour * 3600 + min * 60 + sec;

}

// Debug write real matrix
// Matrix must be stored by column major order
void Util::debugWriteRealMatrix( const int numRows, const int numColumns, const double* const matrix ){

	std::cout << "[";
	for( int row = 0; row < numRows; ++row ){
		for( int col = 0; col < numColumns; ++col ){		
			// Column major
			const int index = col * numRows + row;
			std::cout << matrix[index] << " ";
		}
		if( row + 1 < numRows ){
			std::cout << ";";
		}
	}
	std::cout << "]" << std::endl;

}

// Sort elements by its key value with quick sort
// cf) Numerical Recipes in C++ Second Edition, p336-p339.
// [Input]:
//   1) numOfIDs: Number of elements to be sorted
//   2)   values: Values of each element. Elements are sorted by this values.
// [Input/Output]:
//   1)      ids: Array of elements to be sorted
void Util::quickSort( const int numOfIDs, int* ids, const double* const values ){

	const int minimumSizeOfQuickSort = 7;
	std::vector< std::pair<int, int> > stack;

	if( ids == NULL ){
		OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
		ptrOutputFiles->writeErrorMessage("ids is NULL in quickSort");
	}
	if( values == NULL ){
		OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
		ptrOutputFiles->writeErrorMessage("values is NULL quickSort");
	}
	if( numOfIDs <= 0 ){
		OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
		ptrOutputFiles->writeErrorMessage("numOfIDs is equal to or less than zero");
	}
	int iLeft(0);
	int iRight(numOfIDs - 1);

	while(1){

		if( iRight - iLeft < minimumSizeOfQuickSort ){
			//----- Insertion sort >>>>>
			for( int j = iLeft + 1; j <= iRight; ++j ){
				const int moved = ids[j];
				const double valueOfElemnentMoved = values[moved];
				int i = j - 1;
				for( ; i >= iLeft; --i ){
					if( values[ids[i]] <= valueOfElemnentMoved ){
						break;
					}
					ids[i+1] = ids[i];
				}
				ids[i+1] = moved;
			}

			if( stack.empty() ){
				break;// Finish while loop
			}
			iLeft  = stack.back().first;
			iRight = stack.back().second;
			stack.pop_back();
			//----- Insertion sort <<<<<
		}else{
			//----- Quick sort >>>>>
			std::swap( ids[(iLeft+iRight)/2] , ids[iLeft+1] );// Exchange element
		    if( values[ids[iLeft]] > values[ids[iRight]] ){
				std::swap( ids[iLeft], ids[iRight] );// Exchange element
			}
			if( values[ids[iLeft+1]] > values[ids[iRight]] ){
				std::swap( ids[iLeft+1], ids[iRight] );// Exchange element
			}
			if( values[ids[iLeft]] > values[ids[iLeft+1]] ){
				std::swap( ids[iLeft], ids[iLeft+1] );// Exchange element
			}
			int i = iLeft + 1;
			int j = iRight;
			const int partitioningElement = ids[iLeft+1];// Partitioning element
			const double valueOfPartitioningElement = values[ids[iLeft+1]];// values of partitioning element

			while(1){
				++i;
				while( values[ids[i]] < valueOfPartitioningElement ){
					++i;// Scan up to find element which >= partitioning element
				}
				--j;
				while( values[ids[j]] > valueOfPartitioningElement ){
					--j;// Scan up to find element which <= partitioning element
				}
				if( j < i ){// i and j are crossed, so that partitioning complete
					break;
				}
				std::swap( ids[i], ids[j] );// Exchange element
			}

			 // Exchange IDs of j and iLeft+1 ( partitioning element )
			ids[iLeft+1] = ids[j];
			ids[j] = partitioningElement;

			if( iLeft - i + 1 >= j - iLeft ){// Elements in the right side >= Elements in the left side
				stack.push_back( std::pair<int, int>( i, iRight ) );
				iRight = j - 1;
			}else{// Elements in the right side < Elements in the left side
				stack.push_back( std::pair<int, int>( iLeft, j - 1) );
				iLeft = i;
			}
			//----- Quick sort <<<<<
		}
		
	}

	return;

}

// Sort elements by its key value with quick sort and replace input values
void Util::quickSort( const int numOfIDs, double* values ){

	int* ids = new int[numOfIDs];
	double* valuesOrg = new double[numOfIDs];
	for( int i = 0; i < numOfIDs; ++i ){
		ids[i] = i;
		valuesOrg[i] = values[i];
	}
	Util::quickSort(numOfIDs, ids, valuesOrg);
	for( int i = 0; i < numOfIDs; ++i ){
		values[i] = valuesOrg[ids[i]];
	}
	delete [] ids;
	delete [] valuesOrg;

}

// Calculate Absolute value fo a complex value
double Util::calculateAbsoluteValue( const std::complex<double> value ){

	return hypot(value.real(), value.imag());

}


// Calculate calibration function for FIR filter (Type1)
std::complex<double> Util::calculateCalibrationForFIRFilterType1( const std::string& fileName, const double samplingFreq, const double freq, const bool groupDelay ){

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();

	std::ifstream ifs( fileName.c_str(), std::ios::in );
	if( ifs.fail() ){
		ptrOutputFiles->writeErrorMessage("File open error : " + fileName);
	}
	//ptrOutputFiles->writeLogMessage("Read FIR filter coefficients from " + fileName );

	std::string sbuf;
	//for( int i = 0; i < nskip; ++i){
	//	getline(ifs,sbuf);
	//}

	std::vector<double> h;
	double sumh(0.0); 
	while(getline(ifs,sbuf)){
		if (sbuf.substr(0, 1).compare("#") == 0) {
			// Skip comment lines
			continue;
		}
		std::istringstream iss(sbuf);
		double dbuf;
		iss >> dbuf;
		h.push_back(dbuf);
		sumh += dbuf;
	}
	
	const int numh = static_cast<int>( h.size() );
	//ptrOutputFiles->writeLogMessage("FIR filter length: " + Util::toString(numh) );

	for( int i = 0; i < numh; ++i ){
		h[i] /= sumh;
	}

	const double omega = 2.0 * CommonParameters::PI * freq / samplingFreq;
	std::complex<double> H = std::complex<double>(0.0, 0.0);
	if(groupDelay){
		// Timer is delayed according to group delay
		for( int i = 1; i <= numh; ++i ){
			const double arg = static_cast<double>(i - numh) * omega;
			H += h[i-1] * std::complex<double>( cos(arg), sin(arg) );
		}
	}else{
		// Timer group delay in convolution is adjusted
		const double adjust = (1.0 - static_cast<double>(numh)) / 2.0;
		for( int i = 1; i <= numh; ++i ){
			const double arg = ( static_cast<double>(i - numh) - adjust ) * omega;
			H += h[i-1] * std::complex<double>( cos(arg), sin(arg) );
		}
	}

//#ifdef _DEBUG_WRITE
//	std::cout << std::setw(20) << std::setprecision(12) << 1.0/freq;
//	std::cout << std::setw(20) << std::setprecision(12) << hypot(H.real(), H.imag());
//	std::cout << std::setw(20) << std::setprecision(12) << atan2(H.imag(), H.real()) * CommonParameters::RAD2DEG;
//	std::cout << std::endl;
//#endif

	return H;

}

// Calculate calibration function for FIR filter (Type2)
std::complex<double> Util::calculateCalibrationForFIRFilterType2( const std::string& fileName, const double samplingFreq, const double freq, const bool isELOG,
	int nfstart, int nfend ){

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();

	std::ifstream ifs( fileName.c_str(), std::ios::in );
	if( ifs.fail() ){
		ptrOutputFiles->writeErrorMessage("File open error : " + fileName);
	}

	std::string sbuf;
	//for( int i = 0; i < nskip; ++i){
	//	getline(ifs,sbuf);
	//}

	std::vector<double> h;
	double sumh(0.0);
	while(getline(ifs,sbuf)){
		if (sbuf.substr(0, 1).compare("#") == 0) {
			// Skip comment lines
			continue;
		}
		std::istringstream iss(sbuf);
		double dbuf;
		iss >> dbuf;
		h.push_back(dbuf);
		sumh += dbuf;
	}
	
	const int numh = static_cast<int>( h.size() );

	if (isELOG) {
		nfend = 0;
		nfstart = nfend - numh + 1;
	}
	else {
		if( numh != nfend- nfstart + 1 ){
			ptrOutputFiles->writeErrorMessage("nfend-nfstart+1 is not equal to filter length: " + Util::toString(nfend-nfstart+1) );
		}
	}
	
	const double omega = 2.0 * CommonParameters::PI * freq / samplingFreq;
	std::complex<double> H = std::complex<double>(0.0, 0.0);
	int icount(0);
	for( int i = nfstart; i <= nfend; ++i, ++icount ){
		const double arg = static_cast<double>(i) * omega;
		H += h[icount] * std::complex<double>( cos(arg), sin(arg) );
	}

	return H;

}

// Calculate determinant of real square matrix
double Util::calculateDeterminantOfMatrix( const int dimension, const double* const matrixOrg ){

	const int numComps = dimension * dimension;
	double* matrixWork = new double[numComps];
	// Copy matrix
	for( int i = 0; i < numComps; ++i ){
		matrixWork[i] = matrixOrg[i];
	}

//#ifdef _DEBUG_WRITE
//	Util::debugWriteRealMatrix(dimension, dimension, matrixWork);
//#endif

	int numPivoting(0);
	// Change to triangular matrix
	for( int col = 0; col < dimension - 1; ++col ){
		// Search maximum component
		double maxAbsDiag(0.0);
		int rowMax(-1);
		for( int row = col; row < dimension; ++row ){
			const int index = col * dimension + row;
			if( fabs(matrixWork[index]) > maxAbsDiag ){
				maxAbsDiag = fabs(matrixWork[index]);
				rowMax = row;
			}
		}
		if( maxAbsDiag < 1.0e-20 ){
			delete [] matrixWork;
			return 0.0;
		}
		// Pivoting
		if( rowMax != col ){
			for( int i = col; i < dimension; ++i ){
				const int indexCur = i * dimension + col;
				const int indexMax = i * dimension + rowMax;
				const double valCur = matrixWork[indexCur];
				const double valMax = matrixWork[indexMax];
				matrixWork[indexCur] = valMax;
				matrixWork[indexMax] = valCur;
			}
			++numPivoting;
//#ifdef _DEBUG_WRITE
//			Util::debugWriteRealMatrix(dimension, dimension, matrixWork);
//#endif
		}
		const int indexDiag = col * dimension + col;
		const double diagComp = matrixWork[indexDiag];
		for( int row = col + 1; row < dimension; ++row ){
			// Column major
			const int index = col * dimension + row;
			const double factor = matrixWork[index] / diagComp;
			for( int i = 0; i < dimension; ++i ){
				// Column major
				const int index1 = i * dimension + row;
				const int index2 = i * dimension + col;
				matrixWork[index1] -= matrixWork[index2] * factor; 
			}
		}
//#ifdef _DEBUG_WRITE
//		Util::debugWriteRealMatrix(dimension, dimension, matrixWork);
//#endif
	}
	
	double determinant(1.0);
	for( int i = 0; i < dimension; ++i ){
		// Column major
		const int index = i * dimension + i;
		determinant *= matrixWork[index]; 
	}
	delete [] matrixWork;
	determinant *= pow(-1.0, numPivoting);

	return determinant;

}

// Calcualte all eigenvalues and eigenvectors of a real symmetric
void Util::calculateEigenValuesAndVectorsOfRealSymmetricMatrix(const int dimension, const double* const matrix,
	double* eigenValues, double* eigenVectors) {

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	if (dimension < 1) {
		ptrOutputFiles->writeErrorMessage("Dimension of linear equation is less than 1");
	}

	int icount(0);
	// Column major
	for (int col = 0; col < dimension; ++col) {
		for (int row = 0; row < col; ++row) {
			const int index = col * dimension + row;
			eigenVectors[index] = 0.0;
		}
		for (int row = col; row < dimension; ++row) {
			const int index = col * dimension + row;
			eigenVectors[index] = matrix[icount];
			++icount;
		}
	}

	LapackInterface::calculateEigenValuesAndVectorsOfRealSymmetricMatrix(dimension, eigenVectors, eigenValues);

}

// Calculate IQR
double Util::calculateIQR( const int num, const double* const data ){

	int* ids = new int[num];
	for( int i = 0; i < num; ++i ){
		ids[i] = i;
	}
	Util::quickSort(num, ids, data);
	double* dataSort = new double[num];
	for( int i = 0; i < num; ++i ){
		dataSort[i] = data[ids[i]];
	}
	delete [] ids;
	
	const int quartile1 = num / 4;
	const int quartile3 = num / 4 * 3;
	const double result = dataSort[quartile3] - dataSort[quartile1];
	delete [] dataSort;
	return result;

}

// Calculate MADN
double Util::calculateMADN( const int num, const double* const data ){
	return Util::calculateMAD( num, data ) / 0.675;
}

// Calculate MADN
double Util::calculateMADN( const std::vector<double>& data ){
	return Util::calculateMAD( data ) / 0.675;
}

// Calculate mean value
double Util::calculateMeanValue( const int num, const double* const data ){

	double sum(0.0);
	for( int i = 0; i < num; ++i ){
		sum += data[i];
	}

	return sum / static_cast<double>(num);

}

// Calculate mean square value
double Util::calculateMeanSquareValue( const int num, const double* const data ){

	double sum(0.0);
	for( int i = 0; i < num; ++i ){
		sum += data[i] * data[i];
	}

	return sum / static_cast<double>(num);

}

// Calculate median
double Util::calculateMedian( const int num, const double* const data ){

	int* ids = new int[num];
	for( int i = 0; i < num; ++i ){
		ids[i] = i;
	}
	Util::quickSort(num, ids, data);
//#ifdef _DEBUG_WRITE
//	for( int i = 0; i < num; ++i ){
//		std::cout << "ids[" << i << "]=" << ids[i] << std::endl;
//	}
//#endif
	double* dataSort = new double[num];
	for( int i = 0; i < num; ++i ){
		dataSort[i] = data[ids[i]];
	}
	delete [] ids;
	
//#ifdef _DEBUG_WRITE
//	for( int i = 0; i < num; ++i ){
//		std::cout << "dataSort[" << i << "]=" << dataSort[i] << std::endl;
//	}
//#endif

	if( num % 2 == 1 ){
		// Add number
		const int center = num / 2;
		const double result = dataSort[center];
		delete [] dataSort;
		return result;
	}else{
		// Even number
		const int center0 = num / 2 - 1;
		const int center1 = num / 2;
		const double result = 0.5 * ( dataSort[center0] + dataSort[center1] );
		delete [] dataSort;
		return result;
	}

}

// Calculate median
double Util::calculateMedian( const std::vector<double>& data ){

	std::vector<double> vec = data;
	std::sort(vec.begin(), vec.end());

	const int num = static_cast<int>(data.size());

	if( num % 2 == 1 ){
		// Add number
		const int center = num / 2;
		const double result = vec[center];
		return result;
	}else{
		// Even number
		const int center0 = num / 2 - 1;
		const int center1 = num / 2;
		const double result = 0.5 * ( vec[center0] + vec[center1] );
		return result;
	}

}

// Calculate median absolute deviation
double Util::calculateMAD(const int num, const double median, const double* const data) {

	double* diffAbs = new double[num];
	for (int i = 0; i < num; ++i) {
		diffAbs[i] = fabs(data[i] - median);
	}
	const double result = Util::calculateMedian(num, diffAbs);
	delete[] diffAbs;

	return result;

}

// Calculate median absolute deviation
double Util::calculateMAD(const int num, const double* const data) {

	const double median = Util::calculateMedian(num, data);
	return calculateMAD(num, median, data);

}

// Calculate median absolute deviation
double Util::calculateMAD( const std::vector<double>& data ){

	const double median = Util::calculateMedian( data );

	const int num = static_cast<int>(data.size());
	double* diffAbs = new double[num];
	for( int i = 0; i < num; ++i ){
		diffAbs[i] = fabs( data[i] - median );
	}

	const double result = Util::calculateMedian( num, diffAbs );
	delete [] diffAbs;

	return result;

}

// Calculate regression coefficients by weighted least square method
double Util::calculateRegressionCoefficientsByWLS( const int numOfData, const double* const y, const double* const x,
	const double* const weights ){

	double xx = 0.0;
	double xy = 0.0;
	for( int iData = 0; iData < numOfData; ++iData ){
		xx += x[iData] * x[iData] * weights[iData]; 
		xy += x[iData] * y[iData] * weights[iData]; 
	}

	if( fabs(xx) < CommonParameters::EPS ){
		OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
		ptrOutputFiles->writeWarningMessage("Denominator is too small ("+ Util::toString(xx) + ") in the weighted least square method");
		return xy / CommonParameters::EPS;
	}

	return xy / xx;

}

// Calculate regression coefficients by weighted least square method with intercept
void Util::calculateRegressionCoefficientsByWLSWithIntercept(const int numOfData, const double* const y, const double* const x,
	const double* const weights, double& slope, double& intercept) {

	double wSum = 0.0;
	double wx = 0.0;
	double wy = 0.0;
	double wxx = 0.0;
	double wxy = 0.0;
	for (int iData = 0; iData < numOfData; ++iData) {
		wSum += weights[iData];
		wx += x[iData] * weights[iData];
		wy += y[iData] * weights[iData];
		wxx += x[iData] * x[iData] * weights[iData];
		wxy += x[iData] * y[iData] * weights[iData];
	}

	const double det = wSum * wxx - wx * wx;
	if (fabs(det) < CommonParameters::EPS) {
		OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
		ptrOutputFiles->writeWarningMessage("Denominator is too small (" + Util::toString(det) + ") in the weighted least square method");
	}

	intercept = ( wxx * wy - wx * wxy ) / det;
	slope = ( wSum * wxy - wx * wy ) / det;

}

// Calculate variance
double Util::calculateVariance( const int num, const double mean, const double* const data ){

	double sum(0.0);
	for( int i = 0; i < num; ++i ){
		sum += ( data[i] - mean ) * ( data[i] - mean );
	}

	return sum / static_cast<double>(num-1);

}

// Calculate linear trend factor by Hino (1985)
// [References] 
// Mikio Hino, Spectral analysis, Asakura publishing, 1985 (Japanese)
double Util::calculateLinearTrendFactorByHino1985( const int num, const double* const data ){

	const int numDiv3 = num / 3; 
	const double mean0 = Util::calculateMeanValue( numDiv3, data );
	const double mean1 = Util::calculateMeanValue( numDiv3, &data[num-numDiv3] );
	const double alpha = ( mean1 - mean0 ) / ( num - numDiv3 );
	return alpha;

}

// Calculate linear trend factor by the least square
void Util::calculateLinearTrendFactorByLeastSquare( const int num, const double* const data, double& b0, double& b1 ){

	const double mean = Util::calculateMeanValue( num, data );
	const double sum1 = static_cast<double>( num * (num + 1) ) / 2.0;
	const double sum2 = static_cast<double>( num * (num + 1) * (2*num + 1) ) / 6.0;
	double sum3(0.0);
	for( int i = 0; i < num; ++i ){
		sum3 += data[i] * static_cast<double>(i+1);
	}
	const double dnum = static_cast<double>(num);
	const double deno = dnum * sum2 - sum1 * sum1;
	b0 = ( dnum * mean * sum2 - sum1 * sum3 ) / deno;
	b1 = ( dnum * sum3 - dnum * mean * sum1 ) / deno;

}

// Apply Hampel filter
int Util::hampelFilter( const int num, const int numNeighborsOnEitherSide, const double nsigma, double* data ){

	double* dataMod = new double[num];
	for( int i = 0; i < num; ++i ){
		dataMod[i] = data[i];
	}

	int numModifiedPoints(0);

	const int windowLength = 2 * numNeighborsOnEitherSide + 1;
	for( int iData = numNeighborsOnEitherSide; iData < num - numNeighborsOnEitherSide; ++iData ){
		const double median = Util::calculateMedian( windowLength, &data[iData-numNeighborsOnEitherSide] );
		const double MADN = Util::calculateMADN( windowLength, &data[iData-numNeighborsOnEitherSide] );
		if( fabs( data[iData] - median ) > nsigma * MADN ){
#ifdef _DEBUG_WRITE
			const double absErr = fabs( data[iData] - median );
			const double ratio = absErr / MADN;
#endif
			dataMod[iData] = median;
			++numModifiedPoints;
		}
	}
	// Near the start point
	for( int iData = 0; iData < numNeighborsOnEitherSide; ++iData ){
		const double median = Util::calculateMedian( numNeighborsOnEitherSide, data );
		const double MADN = Util::calculateMADN( numNeighborsOnEitherSide, data );
		if( fabs( data[iData] - median ) > nsigma * MADN ){
#ifdef _DEBUG_WRITE
			const double absErr = fabs( data[iData] - median );
			const double ratio = absErr / MADN;
#endif
			dataMod[iData] = median;
			++numModifiedPoints;
		}
	}
	// Near the end point
	for( int iData = num - numNeighborsOnEitherSide; iData < num; ++iData ){
		const double median = Util::calculateMedian( numNeighborsOnEitherSide, &data[num-numNeighborsOnEitherSide] );
		const double MADN = Util::calculateMADN( numNeighborsOnEitherSide, &data[num-numNeighborsOnEitherSide] );
		if( fabs( data[iData] - median ) > nsigma * MADN ){
#ifdef _DEBUG_WRITE
			const double absErr = fabs( data[iData] - median );
			const double ratio = absErr / MADN;
#endif
			dataMod[iData] = median;
			++numModifiedPoints;
		}
	}

	for( int i = 0; i < num; ++i ){
		data[i] = dataMod[i];
	}
	
	delete [] dataMod;

	return numModifiedPoints;

}

// Apply Hanning window
void Util::hanningWindow( const int num, double* data ){
	
	for( int i = 0; i < num; ++i ){
		const double ratio = static_cast<double>(i) / static_cast<double>(num-1);
		data[i] *= 0.5 * ( 1.0 - cos(2.0 * CommonParameters::PI * ratio) );
	}
	
}

// Perform FFT
// [References] 
// 1) Numerical Recipes in C++ Second Edition, p504-p510.
// 2) Yorihiko Osaki, A new guide to spectral analysis of earthquake motions, Kajima publishing, 1994 (Japanese)
void Util::fft( const int num, std::complex<double>* data, const int isign ){

	if(!Util::isPow2(num)){
		OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
		ptrOutputFiles->writeErrorMessage("Data length for FFT routine must be power of two.");
	}

	// Bit-reversal section
//	int j = 1;
//	for( int i = 1; i <= num; ++i ){
//		if( j > i ){
//#ifdef _DEBUG_WRITE
//			const int bini = decimalToBinary(i-1);
//			const int binj = decimalToBinary(j-1);
//			std::cout << i-1 << " " << j-1 << " " << bini << " " << binj << std::endl;
//#endif
//			const std::complex<double> temp = data[j-1];
//			data[j-1] = data[i-1];
//			data[i-1] = temp;
//		}
//		int m = num / 2;
//		while( m >= 2 && j > m ){
//			j -= m;
//			m /= 2;
//		}
//		j += m;
//	}
	int j = 0;
	for( int i = 0; i < num; ++i ){
		if( j > i ){
//#ifdef _DEBUG_WRITE
//			const int bini = Util::decimalToBinary(i);
//			const int binj = Util::decimalToBinary(j);
//			std::cout << i << " " << j << " " << bini << " " << binj << std::endl;
//#endif
			const std::complex<double> temp = data[j];
			data[j] = data[i];
			data[i] = temp;
		}
		int m = num / 2;
		while( j >= m ){
			j -= m;
			m /= 2;
			if( m < 2 ){
				break;
			}
		}
		j += m;
	}

	// Danielson-Lanczos section
	int kmax = 1;
	while( num > kmax ){
		const int istep = kmax * 2;
		for( int k = 0; k < kmax; ++k ){
			const double ratio = static_cast<double>(isign*k) / static_cast<double>(kmax);
			const std::complex<double> theta = std::complex<double>(0.0, CommonParameters::PI*ratio);
			for( int i = k; i < num; i += istep ){
				const int j = i + kmax;
				const std::complex<double> temp = data[j] * exp(theta);
				data[j] = data[i] - temp; 
				data[i] = data[i] + temp;
			}
		}
		kmax = istep;
	}

}

// Fourier transform
void Util::fourierTransform( const int num, std::complex<double>* data ){

	Util::fft( num, data, -1 );

	const double factor = 1.0 / static_cast<double>(num);
	for( int i = 0; i < num; ++i ){
		data[i] *= factor;
	}	

}

// Inverse fourier transform
void Util::inverseFourierTransform( const int num, std::complex<double>* data ){

	Util::fft( num, data, 1 );

}

// Calculate field at an angle from those of two different directions
std::complex<double> Util::calculateRotatedField( const double direction1, const double direction2, const double rotation, 
							 const std::complex<double>& v1, const std::complex<double>& v2 ){

	const double d1 = direction1 * CommonParameters::DEG2RAD;
	const double d2 = direction2 * CommonParameters::DEG2RAD;
	const double rot = rotation * CommonParameters::DEG2RAD;
	const double factor = 1.0 / sin(d2 - d1);
	const double factor1 = - sin(rot - d2) * factor;
	const double factor2 =   sin(rot - d1) * factor;
	
	return factor1 * v1 + factor2 * v2;

}

// Check wether the input value is power of two
bool Util::isPow2( const int val ){

	for( int i = 1; i < 1000; ++i ){
		const int vpow2 = static_cast<int>( pow(2,i) );
		if( val == vpow2 ){
			return true;
		}
		if( val < vpow2 ){
			return false;
		}
	}

	return false;
}

// Interpolation by the algorithm of Akima (1970)
// [References] 
// Akima, Hiroshi, A new method of interpolation and smooth curve: fitting based on local procedures. J. ACM 17, 4 (Oct. 1970), 589-602.
double Util::interpolationAkima( const int num, const double* const dataX, const double* const dataY, const double x ){

#ifdef _DEBUG_WRITE
	for( int i = 0; i < num; ++i ){
		std::cout << "dataX[" << i << "]=" << dataX[i] << ", dataY[" << i << "]=" << dataY[i] << std::endl;
	}
#endif

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();

	if( num < 1 ){
		ptrOutputFiles->writeErrorMessage("Number of data is less than one in the interpolation by the algorithm of Akima (1970)");
		return -99999;
	}else if( num == 1 ){
		return dataY[0];
	}else if( num == 2 ){
		const double xMin = dataX[0] < dataX[1] ? dataX[0] : dataX[1];
		const double xMax = dataX[0] > dataX[1] ? dataX[0] : dataX[1];
		if( x < xMin ){
			std::ostringstream msg;
			msg << "Target x (" << x << ") is less than the minimum value of data (" << xMin << ") in interpolation. The end point value (" << xMin << ") is used as the result of the interpolation.";
			ptrOutputFiles->writeWarningMessage(msg.str());
			if( dataX[0] < dataX[1] ){
				return dataY[0];
			}else{
				return dataY[1];
			}
		}else if( x > xMax ){
			std::ostringstream msg;
			msg << "Target x (" << x << ") is larger than the maximum value of data (" << xMax << ") in interpolation. The end point value (" << xMax << ") is used as the result of the interpolation.";
			ptrOutputFiles->writeWarningMessage(msg.str());
			if( dataX[0] > dataX[1] ){
				return dataY[0];
			}else{
				return dataY[1];
			}
		}
		return Util::interpolationLinear( dataX[0], dataX[1], dataY[0], dataY[1], x );
	}

	// Arrange the data in ascending order
	int* ids = new int[num];
	for( int i = 0; i < num; ++i ){
		ids[i] = i;
	}
	quickSort(num, ids, dataX);
#ifdef _DEBUG_WRITE
	for( int i = 0; i < num; ++i ){
		std::cout << "ids[" << i << "]=" << ids[i] << std::endl;
	}
#endif
	double* dataXAdd = new double[num+4];
	double* dataYAdd = new double[num+4];
	for( int i = 0; i < num; ++i ){
		dataXAdd[i+2] = dataX[ids[i]];
		dataYAdd[i+2] = dataY[ids[i]];
	}
	delete [] ids;

	const int orgIndexMin = 2;
	const int orgIndexMax = num+1;
	if( x < dataXAdd[orgIndexMin] ){
		std::ostringstream msg;
		msg << "Target x (" << x << ") is less than the minimum value of data (" << dataXAdd[orgIndexMin] << ") in interpolation. The end point value (" << dataYAdd[orgIndexMin] << ") is used as the result of the interpolation.";
		ptrOutputFiles->writeWarningMessage(msg.str());
		return dataYAdd[orgIndexMin];
	}else if( x > dataXAdd[orgIndexMax] ){
		std::ostringstream msg;
		msg << "Target x (" << x << ") is larger than the maximum value of data (" << dataXAdd[orgIndexMax] << ") in interpolation. The end point value (" << dataYAdd[orgIndexMax] << ") is used as the result of the interpolation.";
		ptrOutputFiles->writeWarningMessage(msg.str());
		return dataYAdd[orgIndexMax];
	}
	for( int i = orgIndexMin; i < orgIndexMax; ++i ){
		const double x1 = dataXAdd[i];
		const double x2 = dataXAdd[i+1];
		if( fabs(x1 - x2) < CommonParameters::EPS ){
			ptrOutputFiles->writeErrorMessage("Adjacent two x values are nearly identical in the interpolation");
		}
	}

	// Add two points to left and right ends of the data
	const double dx2Left = dataXAdd[orgIndexMin+2] - dataXAdd[orgIndexMin];
	const double dx2Right = dataXAdd[orgIndexMax] - dataXAdd[orgIndexMax-2];
	dataXAdd[orgIndexMin-2] = dataXAdd[orgIndexMin  ] - dx2Left;
	dataXAdd[orgIndexMin-1] = dataXAdd[orgIndexMin+1] - dx2Left;
	dataXAdd[orgIndexMax+1] = dataXAdd[orgIndexMax-1] + dx2Right;
	dataXAdd[orgIndexMax+2] = dataXAdd[orgIndexMax  ] + dx2Right;
	dataYAdd[orgIndexMin-2] = Util::interpolation2ndOrderLagrange( dataXAdd[orgIndexMin], dataXAdd[orgIndexMin+1], dataXAdd[orgIndexMin+2], 
		dataYAdd[orgIndexMin], dataYAdd[orgIndexMin+1], dataYAdd[orgIndexMin+2], dataXAdd[orgIndexMin-2] );
	dataYAdd[orgIndexMin-1] = Util::interpolation2ndOrderLagrange( dataXAdd[orgIndexMin], dataXAdd[orgIndexMin+1], dataXAdd[orgIndexMin+2], 
		dataYAdd[orgIndexMin], dataYAdd[orgIndexMin+1], dataYAdd[orgIndexMin+2], dataXAdd[orgIndexMin-1] );
	dataYAdd[orgIndexMax+1] = Util::interpolation2ndOrderLagrange( dataXAdd[orgIndexMax], dataXAdd[orgIndexMax-1], dataXAdd[orgIndexMax-2],
		dataYAdd[orgIndexMax], dataYAdd[orgIndexMax-1], dataYAdd[orgIndexMax-2], dataXAdd[orgIndexMax+1] );
	dataYAdd[orgIndexMax+2] = Util::interpolation2ndOrderLagrange( dataXAdd[orgIndexMax], dataXAdd[orgIndexMax-1], dataXAdd[orgIndexMax-2],
		dataYAdd[orgIndexMax], dataYAdd[orgIndexMax-1], dataYAdd[orgIndexMax-2], dataXAdd[orgIndexMax+2] );

#ifdef _DEBUG_WRITE
	for( int i = 0; i < num+4; ++i ){
		std::cout << "dataXAdd[" << i << "]=" << dataXAdd[i] << ", dataYAdd[" << i << "]=" << dataYAdd[i] << std::endl;
	}
#endif

	for( int i = orgIndexMin; i < orgIndexMax; ++i ){
		const double x1 = dataXAdd[i];
		const double x2 = dataXAdd[i+1];
		if( x >= x1 - CommonParameters::EPS && x <= x2 + CommonParameters::EPS ){
			const double y1 = dataYAdd[i];
			const double y2 = dataYAdd[i+1];
			const double t1 = Util::calculateSlopeForInterpolationAkima( &dataXAdd[i-2], &dataYAdd[i-2] );
			const double t2 = Util::calculateSlopeForInterpolationAkima( &dataXAdd[i-1], &dataYAdd[i-1] );
			const double p0 = y1;
			const double p1 = t1;
			const double p2 = ( 3.0 * (y2 - y1) / (x2 - x1) - 2.0 * t1 - t2 ) / ( x2 - x1 );
			const double p3 = ( t1 + t2 - 2.0 * (y2 - y1) / (x2 - x1) ) / ( x2 - x1 ) / ( x2 - x1 );
			const double dx = x-x1;
			const double y = p0 + p1 * dx + p2 * dx * dx + p3 * dx * dx * dx; 
			return y;
		}
	}

	delete [] dataXAdd;
	delete [] dataYAdd;

	std::ostringstream msg;
	msg << "Interpolation was not performed for the target x (" << x << ").";
	ptrOutputFiles->writeErrorMessage(msg.str());
	return -99999;

}

// Calculate slope t for the interpolation by the algorithm of Akima (1970)
// [References] 
// Akima, Hiroshi, A new method of interpolation and smooth curve: fitting based on local procedures. J. ACM 17, 4 (Oct. 1970), 589-602.
double Util::calculateSlopeForInterpolationAkima( const double* const dataX, const double* const dataY ){

	const double m1 = ( dataY[1] - dataY[0] ) / ( dataX[1] - dataX[0] ); 
	const double m2 = ( dataY[2] - dataY[1] ) / ( dataX[2] - dataX[1] ); 
	const double m3 = ( dataY[3] - dataY[2] ) / ( dataX[3] - dataX[2] ); 
	const double m4 = ( dataY[4] - dataY[3] ) / ( dataX[4] - dataX[3] ); 

	if( fabs(m1-m2) < CommonParameters::EPS && fabs(m3-m4) < CommonParameters::EPS ){
		return 0.5 * (m1 + m4);
	}else{
		const double absm43 = fabs( m4 - m3 );
		const double absm21 = fabs( m2 - m1 );
		return ( absm43 * m2 + absm21 * m3 ) / ( absm43 + absm21 );
	}

}

// Linear interpolation
double Util::interpolationLinear( const double x1, const double x2, const double y1, const double y2, const double x ){

	if( fabs(x1 - x2) < CommonParameters::EPS ){
		OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
		ptrOutputFiles->writeErrorMessage("x1 and x2 are nearly identical in the linear interpolation");
	}

	const double term1 = y1 * (x - x2) / (x1 - x2);
	const double term2 = y2 * (x - x1) / (x2 - x1);

	return term1 + term2;

}

// Calculate FIR filter coefficients by the least square method
void Util::calculateFIRFilterCoeffsByLeastSquare( const int dimension, const bool isLowPass, const double samplingFreq,
	const double freq1, const double freq2, const double weight1, const double weight2, double* coeff ){

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	if( dimension <= 0 ){
		ptrOutputFiles->writeErrorMessage("Dimension of FIR filter must be positive");
	}
	if( dimension % 2 != 0 ){
		ptrOutputFiles->writeErrorMessage("Dimension must be an even number when FIR filter is calculated by the least square method.");
	}
	if( samplingFreq <= 0.0 ){
		ptrOutputFiles->writeErrorMessage("Sampling frequency must be positive");
	}
	if( freq1 <= 0.0 || freq2 <= 0.0 ){
		ptrOutputFiles->writeErrorMessage("Frequency must be positive");
	}
	if( freq1 >= freq2 ){
		ptrOutputFiles->writeErrorMessage("The first frequency must be smaller than the second frequency");
	}

	const double normOmega1 = 2.0 * CommonParameters::PI * freq1 / samplingFreq;
	const double normOmega2 = 2.0 * CommonParameters::PI * freq2 / samplingFreq;

	const int num = dimension / 2 + 1;
	double* Amat = new double [ num * num ];
	double* bvec = new double [ num ];

	if( isLowPass ){
		// Low pass
		bvec[0] = weight1 * normOmega1;
		for( int iRow = 1; iRow < num; ++iRow ){
			const double m = static_cast<double>(iRow);
			bvec[iRow] = weight1 * sin(m * normOmega1)/m;
		}
	}else{
		// High pass
		bvec[0] = weight2 * ( CommonParameters::PI - normOmega2 );
		for( int iRow = 1; iRow < num; ++iRow ){
			const double m = static_cast<double>(iRow);
			bvec[iRow] = - weight2 * sin(m * normOmega2)/m;
		}
	}

	int icount(0);
	for( int iCol = 0; iCol < num; ++iCol ){
		for( int iRow = 0; iRow < num; ++iRow ){
			if( iRow == 0 && iCol == 0 ){
				Amat[icount] = weight1 * normOmega1 + weight2 * ( CommonParameters::PI - normOmega2 );
			}else if( iRow == iCol ){
				const double nn = 2.0 * static_cast<double>(iRow);
				const double term1 = 0.5 * weight1 * ( sin(nn * normOmega1) / nn + normOmega1 );
				const double term2 = 0.5 * weight2 * sin(nn * normOmega2) / nn;
				const double term3 = 0.5 * weight2 * ( CommonParameters::PI - normOmega2 );
				Amat[icount] = term1 - term2 + term3;
			}else{
				const double nm1 = static_cast<double>(iRow + iCol);
				const double nm2 = static_cast<double>(iRow - iCol);
				const double term1 = 0.5 * weight1 * sin(nm1 * normOmega1) / nm1;
				const double term2 = 0.5 * weight1 * sin(nm2 * normOmega1) / nm2;
				const double term3 = 0.5 * weight2 * sin(nm1 * normOmega2) / nm1;
				const double term4 = 0.5 * weight2 * sin(nm2 * normOmega2) / nm2;
				Amat[icount] = term1 + term2 - term3 - term4;
			}
			++icount;
		}
	}

	double* result = new double [ num ];
	Util::factorizeAndSolveLinearEquationRealMatrix( num, 1, Amat, bvec, result );
	delete [] Amat;
	delete [] bvec;
	coeff[dimension/2] = result[0];
	for( int iRow = 1; iRow < num; ++iRow ){
		const int index = dimension / 2 - iRow;
		const double h = 0.5 * result[iRow];
		coeff[index] = h;
		coeff[dimension - index] = h;
	}
#ifdef _DEBUG_WRITE
	for( int iRow = 0; iRow < num; ++iRow ){
		std::cout << "a[" << iRow << "]=" << std::setw(20) << std::setprecision(12) << result[iRow] << std::endl;
	}
#endif
	delete [] result;

}

// Calculate frequency characteristics of FIR filter
std::complex<double> Util::calculateFrequencyCharacteristicsOfFIRFilter( const int dimension, const double samplingFreq, const double freq, const double* coeff ){

	const double omega = 2.0 * CommonParameters::PI * freq / samplingFreq;

	std::complex<double> H = std::complex<double>(0.0, 0.0);

	const int halfOfDimension = dimension / 2;
	for( int i = 0; i <= dimension; ++i ){
		const double arg = - static_cast<double>(i - halfOfDimension) * omega;// Group delay is adjusted
		H += std::complex<double>( cos(arg), sin(arg) ) * coeff[i];
	}

	return H;

}

// Calculate frequency characteristics of IIR high-pass filter
std::complex<double> Util::calculateFrequencyCharacteristicsOfIIRHighPassFilter( const double freq, const double samplingFrequency, const double cutoffFreq ){

	const double omega = 2.0 * CommonParameters::PI * freq / samplingFrequency;
	const double lamda = tan(CommonParameters::PI * cutoffFreq / samplingFrequency);
	const double beta  = 1.0 / (lamda + 1.0);
	const double alpha = (lamda - 1.0) / (lamda + 1.0);

	const double arg = - omega;

	std::complex<double> numerator = std::complex<double>(beta, 0.0);
	numerator -= std::complex<double>( cos(arg), sin(arg) ) * beta;

	std::complex<double> denominator = std::complex<double>(1.0, 0.0);
	denominator += std::complex<double>( cos(arg), sin(arg) ) * alpha;

	return numerator / denominator;

}

// Calculate frequency characteristics of IIR low-pass filter
std::complex<double> Util::calculateFrequencyCharacteristicsOfIIRLowPassFilter( const double freq, const double samplingFrequency, const double cutoffFreq ){

	const double Q = 1.0/sqrt(2.0);
	const double omega = 2.0 * CommonParameters::PI * freq / samplingFrequency;
	const double lamda = tan(CommonParameters::PI * cutoffFreq / samplingFrequency);
	const double delta = lamda * lamda + lamda / Q + 1.0;
	const double beta0 = lamda * lamda / delta;
	const double alpha1 = 2.0 * (lamda * lamda - 1.0) / delta;
	const double alpha2 = 1.0 - 2.0 * lamda / Q / delta;
	const double beta[2]  = { 2.0 * beta0,  beta0 };
	const double alpha[2] = { alpha1, alpha2 };

	std::complex<double> numerator = std::complex<double>(beta0, 0.0);
	for( int i = 0; i < 2; ++i ){
		const double arg = - static_cast<double>(i+1) * omega;
		numerator += std::complex<double>( cos(arg), sin(arg) ) * beta[i];
	}

	std::complex<double> denominator = std::complex<double>(1.0, 0.0);
	for( int i = 0; i < 2; ++i ){
		const double arg = - static_cast<double>(i+1) * omega;
		denominator += std::complex<double>( cos(arg), sin(arg) ) * alpha[i];
	}

	return numerator / denominator;

}

// Calculate frequency characteristics of notch filter
std::complex<double> Util::calculateFrequencyCharacteristicsOfNotchFilter( const double freq, const double Q, 	const double samplingFrequency, const double cutoffFreq ){

	const double omega = 2.0 * CommonParameters::PI * freq / samplingFrequency;
	const double lamda = tan(CommonParameters::PI * cutoffFreq / samplingFrequency);
	const double delta = lamda * lamda + lamda / Q + 1.0;
	const double beta0 = (lamda * lamda + 1.0) / delta;
	const double beta1 = 2.0 * (lamda * lamda - 1.0) / delta;
	const double beta2 = beta0;
	const double alpha1 = beta1;
	const double alpha2 = 1.0 - 2.0 * lamda / Q / delta;
	const double beta[2]  = { beta1,  beta2 };
	const double alpha[2] = { alpha1, alpha2 };

	std::complex<double> numerator = std::complex<double>(beta0, 0.0);
	for( int i = 0; i < 2; ++i ){
		const double arg = - static_cast<double>(i+1) * omega;
		numerator += std::complex<double>( cos(arg), sin(arg) ) * beta[i];
	}

	std::complex<double> denominator = std::complex<double>(1.0, 0.0);
	for( int i = 0; i < 2; ++i ){
		const double arg = - static_cast<double>(i+1) * omega;
		denominator += std::complex<double>( cos(arg), sin(arg) ) * alpha[i];
	}

	return numerator / denominator;

}

// Factrize and solve a linear equation with real coefficents and vectors
void Util::factorizeAndSolveLinearEquationRealMatrix( const int dimension, const int nRhs, const double* const matrix, const double* const rhsVectors, double* result ){

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	if( dimension < 1 ){
		ptrOutputFiles->writeErrorMessage("Dimension of linear equation is less than 1" );
	}
	if( nRhs < 1 ){
		ptrOutputFiles->writeErrorMessage("Number of the right-hand-side vectors is less than 1" );
	}

	double* matrixWork = new double[dimension * dimension];
	memcpy(matrixWork, matrix, sizeof(double)*dimension*dimension);
	memcpy(result, rhsVectors, sizeof(double)*nRhs*dimension);

	LapackInterface::factorizeAndSolveLinearEquationRealMatrix(dimension, nRhs, matrixWork, result);
	delete [] matrixWork;

}

// Factrize and solve a linear equation with real symmetric matrix
void Util::factorizeAndSolveLinearEquationRealSymmetricMatrix( const int dimension, const int nRhs, const double* const matrix, const double* const rhsVectors, double* result ){

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	if( dimension < 1 ){
		ptrOutputFiles->writeErrorMessage("Dimension of linear equation is less than 1" );
	}
	if( nRhs < 1 ){
		ptrOutputFiles->writeErrorMessage("Number of the right-hand-side vectors is less than 1" );
	}

	double* matrixWork = new double[dimension * dimension];

	int icount(0);
	// Column major
	for( int col = 0; col < dimension; ++col ){
		for( int row = 0; row < col; ++row ){
			const int index = col * dimension + row;
			matrixWork[index] = 0.0;
		}
		for( int row = col; row < dimension; ++row ){
			const int index = col * dimension + row;
			matrixWork[index] = matrix[icount];
			++icount;
		}
	}
	memcpy(result, rhsVectors, sizeof(double)*nRhs*dimension);

	LapackInterface::factorizeAndSolveLinearEquationRealSymmetricMatrix(dimension, nRhs, matrixWork, result);
	delete [] matrixWork;

}

// Factrize and solve a linear equation with real symmetric positive definite matrix
void Util::factorizeAndSolveLinearEquationRealSymmetricPositiveDefiniteMatrix( const int dimension, const int nRhs, const double* const matrix, const double* const rhsVectors, double* result ){

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	if( dimension < 1 ){
		ptrOutputFiles->writeErrorMessage("Dimension of linear equation is less than 1" );
	}
	if( nRhs < 1 ){
		ptrOutputFiles->writeErrorMessage("Number of the right-hand-side vectors is less than 1" );
	}

	double* matrixWork = new double[dimension * dimension];
	int icount(0);
	// Column major
	for( int col = 0; col < dimension; ++col ){
		for( int row = 0; row < col; ++row ){
			const int index = col * dimension + row;
			matrixWork[index] = 0.0;
		}
		for( int row = col; row < dimension; ++row ){
			const int index = col * dimension + row;
			matrixWork[index] = matrix[icount];
			++icount;
		}
	}
	memcpy(result, rhsVectors, sizeof(double)*nRhs*dimension);

	LapackInterface::factorizeAndSolveLinearEquationRealSymmetricPositiveDefiniteMatrix(dimension, nRhs, matrixWork, result);
	delete [] matrixWork;

}

// Factorize a real square matrix
void Util::factorizeRealSquareMatrix( const int dimension, const double* const matrix, double* factorizedMatrix, int* ipiv ){

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	if( dimension < 1 ){
		ptrOutputFiles->writeErrorMessage("Dimension of linear equation is less than 1" );
	}

	// Copy
	memcpy(factorizedMatrix, matrix, sizeof(double) * dimension * dimension);

	LapackInterface::factorizeRealSquareMatrix( dimension, factorizedMatrix, ipiv );

}

// Factorize a real symmetric matrix
void Util::factorizeRealSymmetricMatrix( const int dimension, const double* const matrix, double* factorizedMatrix, int* ipiv ){

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	if( dimension < 1 ){
		ptrOutputFiles->writeErrorMessage("Dimension of linear equation is less than 1" );
	}

	int icount(0);
	// Column major
	for( int col = 0; col < dimension; ++col ){
		for( int row = 0; row < col; ++row ){
			const int index = col * dimension + row;
			factorizedMatrix[index] = 0.0;
		}
		for( int row = col; row < dimension; ++row ){
			const int index = col * dimension + row;
			factorizedMatrix[index] = matrix[icount];
			++icount;
		}
	}
	LapackInterface::factorizeRealSymmetricMatrix( dimension, factorizedMatrix, ipiv );

}

// Factorize a real symmetric positive definite matrix
void Util::factorizeRealSymmetricPositiveDefiniteMatrix( const int dimension, const double* const matrix, double* factorizedMatrix ){

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	if( dimension < 1 ){
		ptrOutputFiles->writeErrorMessage("Dimension of linear equation is less than 1" );
	}

	int icount(0);
	// Column major
	for( int col = 0; col < dimension; ++col ){
		for( int row = 0; row < col; ++row ){
			const int index = col * dimension + row;
			factorizedMatrix[index] = 0.0;
		}
		for( int row = col; row < dimension; ++row ){
			const int index = col * dimension + row;
			factorizedMatrix[index] = matrix[icount];
			++icount;
		}
	}
	LapackInterface::factorizeRealSymmetricPositiveDefiniteMatrix( dimension, factorizedMatrix );

}

// Solve a linear equation with real square matrix
void Util::solveLinearEquationRealSquareMatrix( const int dimension, const int nRhs, const int* const ipivInt, double* factorizedMatrix, const double* const rhsVectors, double* result ){

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	if( dimension < 1 ){
		ptrOutputFiles->writeErrorMessage("Dimension of linear equation is less than 1" );
	}
	if( nRhs < 1 ){
		ptrOutputFiles->writeErrorMessage("Number of the right-hand-side vectors is less than 1" );
	}

	memcpy(result, rhsVectors, sizeof(double)*nRhs*dimension);
	LapackInterface::solveLinearEquationRealSquareMatrix(dimension, nRhs, ipivInt, factorizedMatrix, result );

}

// Solve a linear equation with real symmetric matrix
void Util::solveLinearEquationRealSymmetricMatrix( const int dimension, const int nRhs, const int* const ipivInt, double* factorizedMatrix, const double* const rhsVectors, double* result ){

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	if( dimension < 1 ){
		ptrOutputFiles->writeErrorMessage("Dimension of linear equation is less than 1" );
	}
	if( nRhs < 1 ){
		ptrOutputFiles->writeErrorMessage("Number of the right-hand-side vectors is less than 1" );
	}

	memcpy(result, rhsVectors, sizeof(double)*nRhs*dimension);
	LapackInterface::solveLinearEquationRealSymmetricMatrix(dimension, nRhs, ipivInt, factorizedMatrix, result );

}

// Solve a linear equation with real symmetric positive definite matrix
void Util::solveLinearEquationRealSymmetricPositiveDefiniteMatrix( const int dimension, const int nRhs, double* factorizedMatrix, const double* const rhsVectors, double* result ){

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	if( dimension < 1 ){
		ptrOutputFiles->writeErrorMessage("Dimension of linear equation is less than 1" );
	}
	if( nRhs < 1 ){
		ptrOutputFiles->writeErrorMessage("Number of the right-hand-side vectors is less than 1" );
	}

	memcpy(result, rhsVectors, sizeof(double)*nRhs*dimension);
	LapackInterface::solveLinearEquationRealSymmetricPositiveDefiniteMatrix(dimension, nRhs, factorizedMatrix, result );

}

std::string Util::extractFileNameWithoutExtension( const std::string& fileNameFull ){

	std::string::size_type pos;

    if( ( pos = fileNameFull.find_last_of(".") ) == std::string::npos){
        return fileNameFull;
    }
 
    return fileNameFull.substr(0, pos);

}

// Extract extension of filename
std::string Util::extractExtensionOfFileName(const std::string& fileNameFull) {

	std::string::size_type pos;

	if ((pos = fileNameFull.find_last_of(".")) == std::string::npos) {
		return fileNameFull;
	}

	return fileNameFull.substr(pos, fileNameFull.length());

}

// The second degree Lagrange interpolation
double Util::interpolation2ndOrderLagrange( const double x1, const double x2, const double x3, 
									 const double y1, const double y2, const double y3, const double x ){
	
	if( fabs(x1 - x2) < CommonParameters::EPS ){
		OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
		ptrOutputFiles->writeErrorMessage("x1 and x2 are nearly identical in the second degree Lagrange interpolation");
	}
	if( fabs(x2 - x3) < CommonParameters::EPS ){
		OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
		ptrOutputFiles->writeErrorMessage("x1 and x2 are nearly identical in the second degree Lagrange interpolation");
	}

	const double term1 = y1 * (x - x2) * (x - x3) / (x1 - x2) / (x1 - x3);
	const double term2 = y2 * (x - x1) * (x - x3) / (x2 - x1) / (x2 - x3);
	const double term3 = y3 * (x - x1) * (x - x2) / (x3 - x1) / (x3 - x2);

	return term1 + term2 + term3;

}

// Convert a decimal number to a binary-coded form
int Util::decimalToBinary(int dec ){

    int bin = 0;
    for( int i = 0; dec > 0 ; ++i){
        bin += ( dec % 2 ) * static_cast<int>( pow(10,i) );
        dec /= 2;
    }
    return bin;

}

// Convert string to integer
int Util::stringToInt( const std::string& sbuf ){

	std::istringstream oss(sbuf);
	int ret;
	oss >> ret;
	return ret;

}

// Convert string to double
double Util::stringToDouble( const std::string& sbuf ){

	std::istringstream oss(sbuf);
	double ret;
	oss >> ret;
	return ret;

}
