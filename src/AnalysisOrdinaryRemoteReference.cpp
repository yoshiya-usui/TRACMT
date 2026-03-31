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
#include "AnalysisOrdinaryRemoteReference.h"
#include "Control.h"
#include "OutputFiles.h"
#include "DoubleDenseSquareMatrix.h"
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <iomanip>
#ifdef _MERSENNE_TWISTER_ORIGINAL
#else
#include <random>
#endif

#include "Util.h"

#ifdef _MERSENNE_TWISTER_ORIGINAL
#include "mt64.h"
#endif

// Default constructer
AnalysisOrdinaryRemoteReference::AnalysisOrdinaryRemoteReference()
{
}

// Destructer
AnalysisOrdinaryRemoteReference::~AnalysisOrdinaryRemoteReference()
{
}

// Calculate partial derivatives of responses for robust bootstrap
void AnalysisOrdinaryRemoteReference::calculatePartialDerivativesOfResponses(const int numSegments, const double threshould,
	std::complex<double>** ftval, const std::complex<double> resp0, const std::complex<double> resp1, const double scale,
	const std::complex<double>* const residuals, const double* const weights, const int iOut, double** derivativesRegardingResps, 
	double* derivativesRegardingScale) const {

	const std::complex<double> resp[2] = {resp0, resp1};

	const Control* const ptrControl = Control::getInstance();
	const int out = ptrControl->getChannelIndex(CommonParameters::OUTPUT, iOut);
	const int numOfInputVariables = ptrControl->getNumInputVariables();
	assert(numOfInputVariables == 2);
	const int in0 = ptrControl->getChannelIndex(CommonParameters::INPUT, 0);
	const int in1 = ptrControl->getChannelIndex(CommonParameters::INPUT, 1);
	const int numOfReferenceVariables = ptrControl->getNumRemoteReferenceVariables();
	assert(numOfReferenceVariables == 2);
	const int rr0 = ptrControl->getChannelIndex(CommonParameters::REMOTE_REFERENCE, 0);
	const int rr1 = ptrControl->getChannelIndex(CommonParameters::REMOTE_REFERENCE, 1);

	const std::complex<double> czero = std::complex<double>(0.0, 0.0);

	std::complex<double> PMatrix[2][2] = { czero, czero, czero, czero };
	for (int iSeg = 0; iSeg < numSegments; ++iSeg) {
		PMatrix[0][0] += ftval[in0][iSeg] * std::conj(ftval[rr0][iSeg]) * weights[iSeg];
		PMatrix[1][0] += ftval[in1][iSeg] * std::conj(ftval[rr0][iSeg]) * weights[iSeg];
		PMatrix[0][1] += ftval[in0][iSeg] * std::conj(ftval[rr1][iSeg]) * weights[iSeg];
		PMatrix[1][1] += ftval[in1][iSeg] * std::conj(ftval[rr1][iSeg]) * weights[iSeg];
	}
#ifdef _DEBUG_WRITE
	std::cout << "[";
	for (int row = 0; row < numOfInputVariables; ++row) {
		for (int col = 0; col < numOfReferenceVariables; ++col) {
			std::cout << PMatrix[row][col].real() << "+" << PMatrix[row][col].imag() << "im ";
		}
		if (row + 1 < numOfInputVariables) {
			std::cout << ";";
		}
	}
	std::cout << "]" << std::endl;
#endif	
	const std::complex<double> det = PMatrix[0][0] * PMatrix[1][1] - PMatrix[1][0] * PMatrix[0][1];
	const std::complex<double> PInvMatrix[2][2] = { PMatrix[1][1] / det, -PMatrix[0][1] / det, -PMatrix[1][0] / det, PMatrix[0][0] / det };
	const std::complex<double> PInvTMatrix[2][2] = { PInvMatrix[0][0], PInvMatrix[1][0], PInvMatrix[0][1], PInvMatrix[1][1] };
	std::complex<double> PInvTIMatrix[2][2] = { czero, czero, czero, czero };
	for (int irow = 0; irow < numOfInputVariables; ++irow) {
		for (int icol = 0; icol < numOfReferenceVariables; ++icol) {
			// PInvTIMatrix = PInvTMatrix because Iqxq = 1
			PInvTIMatrix[irow][icol] = PInvTMatrix[irow][icol];
		}
	}
#ifdef _DEBUG_WRITE
	std::cout << "[";
	for( int row = 0; row < numOfReferenceVariables; ++row ){
		for( int col = 0; col < numOfReferenceVariables; ++col ){
			std::cout << PInvTIMatrix[row][col].real() << "+" << PInvTIMatrix[row][col].imag() <<"im ";
		}
		if( row+1 < numOfReferenceVariables ){
			std::cout << ";";
		}
	}
	std::cout << "]" << std::endl;
#endif	
	std::complex<double>** PInvTPInvMatrix = new std::complex<double>*[numOfInputVariables * numOfReferenceVariables];
	for (int irow = 0; irow < numOfInputVariables * numOfReferenceVariables; ++irow) {
		PInvTPInvMatrix[irow] = new std::complex<double>[numOfInputVariables * numOfReferenceVariables];
		for (int icol = 0; icol < numOfInputVariables * numOfReferenceVariables; ++icol) {
			PInvTPInvMatrix[irow][icol] = czero;// Zero clear
		}
	}
	for (int irow = 0; irow < numOfReferenceVariables; ++irow) {
		for (int icol = 0; icol < numOfInputVariables; ++icol) {
			const std::complex<double> factor = PInvTMatrix[irow][icol];
			for (int irow2 = 0; irow2 < numOfReferenceVariables; ++irow2) {
				for (int icol2 = 0; icol2 < numOfInputVariables; ++icol2) {
					const int irowOut = irow2 + irow * numOfReferenceVariables;
					const int icolOut = icol2 + icol * numOfInputVariables;
					PInvTPInvMatrix[irowOut][icolOut] = factor * PInvMatrix[irow2][icol2];
				}
			}
		}
	}
#ifdef _DEBUG_WRITE
	std::cout << "[";
	for( int row = 0; row < numOfInputVariables * numOfReferenceVariables; ++row ){
		for( int col = 0; col < numOfInputVariables * numOfReferenceVariables; ++col ){
			std::cout << PInvTPInvMatrix[row][col].real() << "+" << PInvTPInvMatrix[row][col].imag() <<"im ";
		}
		if( row+1 < numOfInputVariables * numOfReferenceVariables ){
			std::cout << ";";
		}
	}
	std::cout << "]" << std::endl;
#endif
	std::complex<double> QMatrix[2] = { czero, czero };
	for (int iSeg = 0; iSeg < numSegments; ++iSeg) {
		QMatrix[0] += ftval[out][iSeg] * std::conj(ftval[rr0][iSeg]) * weights[iSeg];
		QMatrix[1] += ftval[out][iSeg] * std::conj(ftval[rr1][iSeg]) * weights[iSeg];
	}
#ifdef _DEBUG_WRITE
	std::cout << "[";
	for( int col = 0; col < numOfReferenceVariables; ++col ){
		std::cout << QMatrix[col].real() << "+" << QMatrix[col].imag() <<"im ";
	}
	std::cout << "]" << std::endl;
	std::cout << "[";
	for( int col = 0; col < numOfReferenceVariables; ++col ){
		std::cout << resp[col].real() << "+" << resp[col].imag() <<"im ";
	}
	std::cout << "]" << std::endl;
#endif
	std::complex<double>** IQMatrix = new std::complex<double>*[numOfInputVariables];
	for (int irow = 0; irow < numOfInputVariables; ++irow) {
		IQMatrix[irow] = new std::complex<double>[numOfInputVariables * numOfReferenceVariables];
		for (int icol = 0; icol < numOfInputVariables * numOfReferenceVariables; ++icol) {
			IQMatrix[irow][icol] = czero;// Zero clear
		}
	}
	for (int irow = 0; irow < numOfInputVariables; ++irow) {
		for (int icol = 0; icol < numOfInputVariables; ++icol) {
			const double factor = irow == icol ? 1.0 : 0.0;
			for (int icol2 = 0; icol2 < numOfReferenceVariables; ++icol2) {
				const int irowOut = irow;
				const int icolOut = icol2 + icol * numOfReferenceVariables;
				IQMatrix[irowOut][icolOut] = factor * QMatrix[icol2];
			}
		}
	}
#ifdef _DEBUG_WRITE
	std::cout << "[";
	for( int row = 0; row < numOfInputVariables; ++row ){
		for( int col = 0; col < numOfInputVariables * numOfReferenceVariables; ++col ){
			std::cout << IQMatrix[row][col].real() << "+" << IQMatrix[row][col].imag() <<"im ";
		}
		if( row+1 < numOfInputVariables){
			std::cout << ";";
		}
	}
	std::cout << "]" << std::endl;
#endif
	std::complex<double>** IQPInvTPInvMatrix = new std::complex<double>*[numOfInputVariables];
	for (int irow = 0; irow < numOfInputVariables; ++irow) {
		IQPInvTPInvMatrix[irow] = new std::complex<double>[numOfInputVariables * numOfReferenceVariables];
		for (int icol = 0; icol < numOfInputVariables * numOfReferenceVariables; ++icol) {
			IQPInvTPInvMatrix[irow][icol] = czero;// Zero clear
		}
	}
	for (int irow = 0; irow < numOfInputVariables; ++irow) {
		for (int icol = 0; icol < numOfInputVariables * numOfReferenceVariables; ++icol) {
			std::complex<double> value = czero;
			for (int i = 0; i < numOfInputVariables * numOfReferenceVariables; ++i) {
				value += IQMatrix[irow][i] * PInvTPInvMatrix[i][icol];
			}
			IQPInvTPInvMatrix[irow][icol] = value;
		}
	}
#ifdef _DEBUG_WRITE
#else
	for (int irow = 0; irow < numOfInputVariables * numOfReferenceVariables; ++irow) {
		delete[] PInvTPInvMatrix[irow];
	}
	delete[] PInvTPInvMatrix;
#endif
	for (int irow = 0; irow < numOfInputVariables; ++irow) {
		delete[] IQMatrix[irow];
	}
	delete[] IQMatrix;
#ifdef _DEBUG_WRITE
	std::cout << "[";
	for( int row = 0; row < numOfInputVariables; ++row ){
		for( int col = 0; col < numOfInputVariables * numOfReferenceVariables; ++col ){
			std::cout << IQPInvTPInvMatrix[row][col].real() << "+" << IQPInvTPInvMatrix[row][col].imag() <<"im ";
		}
		if( row+1 < numOfInputVariables){
			std::cout << ";";
		}
	}
	std::cout << "]" << std::endl;
#endif
	std::complex<double>** matrixVeceh1 = new std::complex<double>*[numOfReferenceVariables];
	for (int irow = 0; irow < numOfReferenceVariables; ++irow) {
		matrixVeceh1[irow] = new std::complex<double>[2 * numOfInputVariables];
		for (int icol = 0; icol < 2 * numOfInputVariables; ++icol) {
			matrixVeceh1[irow][icol] = czero;// Zero clear
		}
	}
	std::complex<double>** matrixVechh1 = new std::complex<double>*[numOfInputVariables * numOfReferenceVariables];
	for (int irow = 0; irow < numOfInputVariables * numOfReferenceVariables; ++irow) {
		matrixVechh1[irow] = new std::complex<double>[2 * numOfInputVariables];
		for (int icol = 0; icol < 2 * numOfInputVariables; ++icol) {
			matrixVechh1[irow][icol] = czero;// Zero clear
		}
	}
	std::complex<double>* veceh = new std::complex<double>[numOfReferenceVariables];
	std::complex<double>* vechh = new std::complex<double>[numOfInputVariables * numOfReferenceVariables];
	std::complex<double>* hSigmaMatrix = new std::complex<double>[numOfInputVariables];
	double* vecSigma = new double[2 * numOfInputVariables];
	std::complex<double>* sumVeceh3 = new std::complex<double>[numOfReferenceVariables];
	for (int irow = 0; irow < numOfReferenceVariables; ++irow) {
		sumVeceh3[irow] = czero;// Zero clear
	}
	std::complex<double>* sumVechh3 = new std::complex<double>[numOfInputVariables * numOfReferenceVariables];
	for (int irow = 0; irow < numOfInputVariables * numOfReferenceVariables; ++irow) {
		sumVechh3[irow] = czero;// Zero clear
	}
#ifdef _DEBUG_WRITE
	std::cout << "[";
#endif
	for (int iSeg = 0; iSeg < numSegments; ++iSeg) {
		const double val = std::abs(residuals[iSeg]) / scale;
		const double diff = RobustWeightTukeysBiweights::calculateSecondDerivativeOfLossFunction(val, threshould) - RobustWeightTukeysBiweights::calculateWeights(val, threshould);
		double factor1(0.0);
		if (std::norm(residuals[iSeg]) < CommonParameters::EPS) {
			const double sc = scale * threshould;
			factor1 = 4.0 * ( std::norm(residuals[iSeg]) / pow(sc, 4) - 1.0 / pow(sc, 2) );
		}
		else {
			factor1 = diff / std::norm(residuals[iSeg]);
		}
		const double factor2 = diff / scale;
		for (int irr = 0; irr < numOfReferenceVariables; ++irr) {
			const int rr = ptrControl->getChannelIndex(CommonParameters::REMOTE_REFERENCE, irr);
			veceh[irr] = ftval[out][iSeg] * std::conj(ftval[rr][iSeg]) * factor1;
			sumVeceh3[irr] += ftval[out][iSeg] * std::conj(ftval[rr][iSeg]) * factor2;
		}
		for (int irr = 0; irr < numOfReferenceVariables; ++irr) {
			const int offset = numOfInputVariables * irr;
			const int rr = ptrControl->getChannelIndex(CommonParameters::REMOTE_REFERENCE, irr);
			for (int iInp = 0; iInp < numOfInputVariables; ++iInp) {
				const int index = ptrControl->getChannelIndex(CommonParameters::INPUT, iInp);
				vechh[iInp + offset] = ftval[index][iSeg] * std::conj(ftval[rr][iSeg]) * factor1;
				sumVechh3[iInp + offset] += ftval[index][iSeg] * std::conj(ftval[rr][iSeg]) * factor2;
			}
		}
#ifdef _DEBUG_WRITE
		std::cout << "[";
		for( int row = 0; row < numOfReferenceVariables; ++row ){
			std::cout << veceh[row].real() << "+" << veceh[row].imag() <<"im ";
			if( row+1 < numOfReferenceVariables ){
				std::cout << ",";
			}
		}
		std::cout << "]" << std::endl;
		std::cout << "[";
		for( int row = 0; row < numOfInputVariables * numOfReferenceVariables; ++row ){
			std::cout << vechh[row].real() << "+" << vechh[row].imag() <<"im ";
			if( row+1 < numOfInputVariables * numOfReferenceVariables ){
				std::cout << ",";
			}
		}
		std::cout << "]" << std::endl;
#endif
		calculateVectorForPartialDerivatives(iSeg, ftval, residuals, hSigmaMatrix, vecSigma);
		for (int irow = 0; irow < numOfReferenceVariables; ++irow) {
			for (int icol = 0; icol < 2 * numOfInputVariables; ++icol) {
				matrixVeceh1[irow][icol] += veceh[irow] * vecSigma[icol];
			}
		}
		for (int irow = 0; irow < numOfInputVariables * numOfReferenceVariables; ++irow) {
			for (int icol = 0; icol < 2 * numOfInputVariables; ++icol) {
				matrixVechh1[irow][icol] += vechh[irow] * vecSigma[icol];
			}
		}
#ifdef _DEBUG_WRITE
		std::cout << "[";
		for( int row = 0; row < numOfReferenceVariables; ++row ){
			for( int col = 0; col < 2 * numOfInputVariables; ++col ){
				std::cout << matrixVeceh1[row][col].real() << "+" << matrixVeceh1[row][col].imag() <<"im ";
			}
			if( row+1 < numOfReferenceVariables ){
				std::cout << ";";
			}
		}
		std::cout << "]" << std::endl;
		std::cout << "[";
		for( int row = 0; row < numOfInputVariables * numOfReferenceVariables; ++row ){
			for( int col = 0; col < 2 * numOfInputVariables; ++col ){
				std::cout << matrixVechh1[row][col].real() << "+" << matrixVechh1[row][col].imag() <<"im ";
			}
			if( row+1 < numOfInputVariables * numOfReferenceVariables ){
				std::cout << ";";
			}
		}
		std::cout << "]" << std::endl;
#endif
	}
#ifdef _DEBUG_WRITE
	std::cout << "]" << std::endl;
	std::cout << "[";
	for( int row = 0; row < numOfReferenceVariables; ++row ){
		std::cout << sumVeceh3[row].real() << "+" << sumVeceh3[row].imag() <<"im ";
		if( row+1 < numOfReferenceVariables ){
			std::cout << ",";
		}
	}
	std::cout << "]" << std::endl;
	std::cout << "[";
	for( int row = 0; row < numOfInputVariables * numOfReferenceVariables; ++row ){
		std::cout << sumVechh3[row].real() << "+" << sumVechh3[row].imag() <<"im ";
		if( row+1 < numOfInputVariables * numOfReferenceVariables ){
			std::cout << ",";
		}
	}
	std::cout << "]" << std::endl;
	std::cout << "[";
	for( int row = 0; row < numOfReferenceVariables; ++row ){
		for( int col = 0; col < 2 * numOfInputVariables; ++col ){
			const int iq = col / ( 2 * numOfInputVariables);
			const int isign = ( col / numOfInputVariables) % 2;
			const int index = 2 * ( col % numOfInputVariables) + isign + 2 * numOfInputVariables * iq;
			std::cout << -matrixVeceh1[row][index].real() << "+" << -matrixVeceh1[row][index].imag() <<"im ";
		}
		if( row+1 < numOfReferenceVariables ){
			std::cout << ";";
		}
	}
	std::cout << "]" << std::endl;
	std::cout << "[";
	for( int row = 0; row < numOfInputVariables * numOfReferenceVariables; ++row ){
		for( int col = 0; col < 2 * numOfInputVariables; ++col ){
			const int iq = col / ( 2 * numOfInputVariables);
			const int isign = ( col / numOfInputVariables) % 2;
			const int index = 2 * ( col % numOfInputVariables) + isign + 2 * numOfInputVariables * iq;
			std::cout << -matrixVechh1[row][index].real() << "+" << -matrixVechh1[row][index].imag() <<"im ";
		}
		if( row+1 < numOfInputVariables * numOfReferenceVariables ){
			std::cout << ";";
		}
	}
	std::cout << "]" << std::endl;
	std::cout << "[";
	for( int irow = 0; irow < numOfInputVariables * numOfReferenceVariables; ++irow ){
		for( int icol = 0; icol < 2 * numOfInputVariables; ++icol ){
			std::complex<double> value = czero;
			const int iq = icol / ( 2 * numOfInputVariables);
			const int isign = ( icol / numOfInputVariables) % 2;
			const int index = 2 * ( icol % numOfInputVariables) + isign + 2 * numOfInputVariables * iq;
			for( int i = 0; i < numOfInputVariables * numOfReferenceVariables; ++i ){
				value += PInvTPInvMatrix[irow][i] * matrixVechh1[i][index];
			}
			std::cout << value.real() << "+" << value.imag() <<"im ";
		}
		if( irow+1 < numOfInputVariables * numOfReferenceVariables ){
			std::cout << ";";
		}
	}
	std::cout << "]" << std::endl;
	for( int irow = 0; irow < numOfInputVariables * numOfReferenceVariables; ++irow ){
		delete [] PInvTPInvMatrix[irow];
	}
	delete [] PInvTPInvMatrix;
#endif
	delete[] veceh;
	delete[] vechh;
	delete[] hSigmaMatrix;
	delete[] vecSigma;

	// Partial derivatives regarding responses
	std::complex<double>** complexDerivatives = new std::complex<double>*[numOfInputVariables];
	for (int irow = 0; irow < numOfInputVariables; ++irow) {
		complexDerivatives[irow] = new std::complex<double>[2 * numOfInputVariables];
	}
	for (int irow = 0; irow < numOfInputVariables; ++irow) {
		for (int icol = 0; icol < 2 * numOfInputVariables; ++icol) {
			std::complex<double> value = czero;
			// The 1st term
			for (int i = 0; i < numOfReferenceVariables; ++i) {
				value -= PInvTIMatrix[irow][i] * matrixVeceh1[i][icol];
			}
			// The 2nd term
			for (int i = 0; i < numOfInputVariables * numOfReferenceVariables; ++i) {
				value += IQPInvTPInvMatrix[irow][i] * matrixVechh1[i][icol];
			}
			complexDerivatives[irow][icol] = value;
		}
	}
	for (int irow = 0; irow < numOfReferenceVariables; ++irow) {
		delete[] matrixVeceh1[irow];
	}
	delete[] matrixVeceh1;
	for (int irow = 0; irow < numOfInputVariables * numOfReferenceVariables; ++irow) {
		delete[] matrixVechh1[irow];
	}
	delete[] matrixVechh1;
#ifdef _DEBUG_WRITE
	std::cout << "[";
	for( int row = 0; row < numOfReferenceVariables; ++row ){
		for( int col = 0; col < 2 * numOfReferenceVariables; ++col ){
			std::cout << complexDerivatives[row][col].real() << "+" << complexDerivatives[row][col].imag() <<"im ";
		}
		if( row+1 < numOfReferenceVariables ){
			std::cout << ";";
		}
	}
	std::cout << "]" << std::endl;
#endif

	// Order: Re(Z11),Im(Z11),Re(Z12),Im(Z12)
	for (int irow = 0; irow < numOfReferenceVariables; ++irow) {
		const int ir = irow;
		const int index = 2 * ir;
		for (int icol = 0; icol < 2 * numOfInputVariables; ++icol) {
			if (icol < numOfInputVariables) {
				// Real part
				const int ir2 = icol;
				const int index2 = 2 * ir2;
				derivativesRegardingResps[index][index2] = complexDerivatives[irow][icol].real();
				derivativesRegardingResps[index + 1][index2] = complexDerivatives[irow][icol].imag();
			}
			else {
				// Imaginary part
				const int ir2 = icol - numOfInputVariables;
				const int index2 = 2 * ir2 + 1;
				derivativesRegardingResps[index][index2] = complexDerivatives[irow][icol].real();
				derivativesRegardingResps[index + 1][index2] = complexDerivatives[irow][icol].imag();
			}
		}
	}
#ifdef _DEBUG_WRITE
	std::cout << "drr1" << std::endl;
	std::cout << "[";
	for (int row = 0; row < 2 * numOfInputVariables; ++row) {
		for (int col = 0; col < 2 * numOfInputVariables; ++col) {
			std::cout << derivativesRegardingResps[row][col] << " ";
		}
		if (row + 1 < 2 * numOfInputVariables) {
			std::cout << ";";
		}
	}
	std::cout << "]" << std::endl;
	std::complex<double>** derivP1 = new std::complex<double>*[numOfInputVariables * numOfReferenceVariables];
	for (int irow = 0; irow < numOfInputVariables * numOfReferenceVariables; ++irow) {
		derivP1[irow] = new std::complex<double>[2 * numOfInputVariables];
	}
	std::complex<double>** derivInvP1 = new std::complex<double>*[numOfInputVariables * numOfReferenceVariables];
	for (int irow = 0; irow < numOfInputVariables * numOfReferenceVariables; ++irow) {
		derivInvP1[irow] = new std::complex<double>[2 * numOfInputVariables];
	}
	std::complex<double>** derivQ1 = new std::complex<double>*[numOfReferenceVariables];
	for (int irow = 0; irow < numOfReferenceVariables; ++irow) {
		derivQ1[irow] = new std::complex<double>[2 * numOfInputVariables];
	}
	double** derivResp1 = new double* [2 * numOfInputVariables];
	for (int irow = 0; irow < 2 * numOfInputVariables; ++irow) {
		derivResp1[irow] = new double[2 * numOfInputVariables];
	}
	for (int icol = 0; icol < numOfInputVariables; ++icol) {
		for (int isign = 0; isign < 2; ++isign) {
			std::complex<double>* residualsOrg = new std::complex<double>[numSegments];
			double* weightsOrg = new double[numSegments];
			std::complex<double>* respOrg = new std::complex<double>[numOfInputVariables];
			std::complex<double>* residualsMod = new std::complex<double>[numSegments];
			double* weightsMod = new double[numSegments];
			std::complex<double>* respMod = new std::complex<double>[numOfInputVariables];
			for (int icol2 = 0; icol2 < numOfInputVariables; ++icol2) {
				respOrg[icol2] = resp[icol2];
				respMod[icol2] = resp[icol2];
			}
			double dresp(0.0);
			if (isign == 0) {
				dresp = resp[icol].real() * 0.00001;
				respMod[icol] += std::complex<double>(dresp, 0.0);
			}
			else {
				dresp = resp[icol].imag() * 0.00001;
				respMod[icol] += std::complex<double>(0.0, dresp);
			}
			// Calculate complex residuals
			for (int iSeg = 0; iSeg < numSegments; ++iSeg) {
				const std::complex<double> outSyn = respOrg[0] * ftval[in0][iSeg] + respOrg[1] * ftval[in1][iSeg];
				residualsOrg[iSeg] = ftval[out][iSeg] - outSyn;
			}
			for (int iSeg = 0; iSeg < numSegments; ++iSeg) {
				const std::complex<double> outSyn = respMod[0] * ftval[in0][iSeg] + respMod[1] * ftval[in1][iSeg];
				residualsMod[iSeg] = ftval[out][iSeg] - outSyn;
			}
			const int index = 2 * icol + isign;
			// Calculate original weights
			for (int iSeg = 0; iSeg < numSegments; ++iSeg) {
				const double x = std::abs(residualsOrg[iSeg]) / scale;
				weightsOrg[iSeg] = RobustWeightTukeysBiweights::calculateWeights(x, threshould);
			}
			for (int iSeg = 0; iSeg < numSegments; ++iSeg) {
				const double x = std::abs(residualsMod[iSeg]) / scale;
				weightsMod[iSeg] = RobustWeightTukeysBiweights::calculateWeights(x, threshould);
			}
			std::complex<double> PMatrixOrg[2][2] = { czero, czero, czero, czero };
			for (int iSeg = 0; iSeg < numSegments; ++iSeg) {
				PMatrixOrg[0][0] += ftval[in0][iSeg] * std::conj(ftval[rr0][iSeg]) * weightsOrg[iSeg];
				PMatrixOrg[1][0] += ftval[in1][iSeg] * std::conj(ftval[rr0][iSeg]) * weightsOrg[iSeg];
				PMatrixOrg[0][1] += ftval[in0][iSeg] * std::conj(ftval[rr1][iSeg]) * weightsOrg[iSeg];
				PMatrixOrg[1][1] += ftval[in1][iSeg] * std::conj(ftval[rr1][iSeg]) * weightsOrg[iSeg];
			}
			std::complex<double> PMatrixMod[2][2] = { czero, czero, czero, czero };
			for (int iSeg = 0; iSeg < numSegments; ++iSeg) {
				PMatrixMod[0][0] += ftval[in0][iSeg] * std::conj(ftval[rr0][iSeg]) * weightsMod[iSeg];
				PMatrixMod[1][0] += ftval[in1][iSeg] * std::conj(ftval[rr0][iSeg]) * weightsMod[iSeg];
				PMatrixMod[0][1] += ftval[in0][iSeg] * std::conj(ftval[rr1][iSeg]) * weightsMod[iSeg];
				PMatrixMod[1][1] += ftval[in1][iSeg] * std::conj(ftval[rr1][iSeg]) * weightsMod[iSeg];
			}
			derivP1[0][index] = (PMatrixMod[0][0] - PMatrixOrg[0][0]) / dresp;
			derivP1[1][index] = (PMatrixMod[1][0] - PMatrixOrg[1][0]) / dresp;
			derivP1[2][index] = (PMatrixMod[0][1] - PMatrixOrg[0][1]) / dresp;
			derivP1[3][index] = (PMatrixMod[1][1] - PMatrixOrg[1][1]) / dresp;

			const std::complex<double> detOrg = PMatrixOrg[0][0] * PMatrixOrg[1][1] - PMatrixOrg[1][0] * PMatrixOrg[0][1];
			const std::complex<double> PInvMatrixOrg[2][2] = { PMatrixOrg[1][1] / detOrg, -PMatrixOrg[0][1] / detOrg, -PMatrixOrg[1][0] / detOrg, PMatrixOrg[0][0] / detOrg };
			const std::complex<double> detMod = PMatrixMod[0][0] * PMatrixMod[1][1] - PMatrixMod[1][0] * PMatrixMod[0][1];
			const std::complex<double> PInvMatrixMod[2][2] = { PMatrixMod[1][1] / detMod, -PMatrixMod[0][1] / detMod, -PMatrixMod[1][0] / detMod, PMatrixMod[0][0] / detMod };
			derivInvP1[0][index] = (PInvMatrixMod[0][0] - PInvMatrixOrg[0][0]) / dresp;
			derivInvP1[1][index] = (PInvMatrixMod[1][0] - PInvMatrixOrg[1][0]) / dresp;
			derivInvP1[2][index] = (PInvMatrixMod[0][1] - PInvMatrixOrg[0][1]) / dresp;
			derivInvP1[3][index] = (PInvMatrixMod[1][1] - PInvMatrixOrg[1][1]) / dresp;

			std::complex<double>* QMatrixOrg = new std::complex<double>[numOfReferenceVariables];
			std::complex<double>* QMatrixMod = new std::complex<double>[numOfReferenceVariables];
			for (int icol = 0; icol < numOfReferenceVariables; ++icol) {
				QMatrixOrg[icol] = czero;
				QMatrixMod[icol] = czero;
			}
			for (int irr = 0; irr < numOfReferenceVariables; ++irr) {
				const int rr = ptrControl->getChannelIndex(CommonParameters::REMOTE_REFERENCE, irr);
				for (int iSeg = 0; iSeg < numSegments; ++iSeg) {
					QMatrixOrg[irr] += ftval[out][iSeg] * std::conj(ftval[rr][iSeg]) * weightsOrg[iSeg];
					QMatrixMod[irr] += ftval[out][iSeg] * std::conj(ftval[rr][iSeg]) * weightsMod[iSeg];
				}
			}
			for (int irr = 0; irr < numOfReferenceVariables; ++irr) {
				derivQ1[irr][index] = (QMatrixMod[irr] - QMatrixOrg[irr]) / dresp;
			}
			for (int icol2 = 0; icol2 < numOfInputVariables; ++icol2) {
				std::complex<double> valueOrg = czero;
				std::complex<double> valueMod = czero;
				for (int i = 0; i < numOfReferenceVariables; ++i) {
					valueOrg += QMatrixOrg[i] * PInvMatrixOrg[i][icol2];
					valueMod += QMatrixMod[i] * PInvMatrixMod[i][icol2];
				}
				const int index2 = 2 * icol2;
				derivResp1[index2    ][index] = (valueMod.real() - valueOrg.real()) / dresp;
				derivResp1[index2 + 1][index] = (valueMod.imag() - valueOrg.imag()) / dresp;
			}
			delete[] QMatrixOrg;
			delete[] QMatrixMod;
			delete[] residualsOrg;
			delete[] weightsOrg;
			delete[] respOrg;
			delete[] residualsMod;
			delete[] weightsMod;
			delete[] respMod;
		}
	}
	std::cout << "[";
	for( int row = 0; row < numOfInputVariables * numOfReferenceVariables; ++row ){
		for( int col = 0; col < 2 * numOfInputVariables; ++col ){
			std::cout << derivP1[row][col].real() << "+" << derivP1[row][col].imag() <<"im ";
		}
		if( row+1 < numOfInputVariables * numOfReferenceVariables ){
			std::cout << ";";
		}
	}
	std::cout << "]" << std::endl;
	std::cout << "[";
	for( int row = 0; row < numOfInputVariables * numOfReferenceVariables; ++row ){
		for( int col = 0; col < 2 * numOfInputVariables; ++col ){
			std::cout << derivInvP1[row][col].real() << "+" << derivInvP1[row][col].imag() <<"im ";
		}
		if( row+1 < numOfInputVariables * numOfReferenceVariables ){
			std::cout << ";";
		}
	}
	std::cout << "]" << std::endl;
	std::cout << "[";
	for( int row = 0; row < numOfReferenceVariables; ++row ){
		for( int col = 0; col < 2 * numOfReferenceVariables; ++col ){
			std::cout << derivQ1[row][col].real() << "+" << derivQ1[row][col].imag() <<"im ";
		}
		if( row+1 < numOfReferenceVariables ){
			std::cout << ";";
		}
	}
	std::cout << "]" << std::endl;
	std::cout << "drr2" << std::endl;
	std::cout << "[";
	for (int row = 0; row < 2 * numOfInputVariables; ++row) {
		for (int col = 0; col < 2 * numOfInputVariables; ++col) {
			std::cout << derivResp1[row][col] << " ";
		}
		if (row + 1 < 2 * numOfInputVariables) {
			std::cout << ";";
		}
	}
	std::cout << "]" << std::endl;
	for (int irow = 0; irow < numOfInputVariables * numOfReferenceVariables; ++irow) {
		delete[] derivP1[irow];
	}
	delete[] derivP1;
	for (int irow = 0; irow < numOfInputVariables * numOfReferenceVariables; ++irow) {
		delete[] derivInvP1[irow];
	}
	delete[] derivInvP1;
	for (int irow = 0; irow < numOfReferenceVariables; ++irow) {
		delete[] derivQ1[irow];
	}
	delete[] derivQ1;
	for (int irow = 0; irow < 2 * numOfInputVariables; ++irow) {
		delete[] derivResp1[irow];
	}
	delete[] derivResp1;
#endif

	// Partial derivatives regarding scale
	for (int irow = 0; irow < numOfInputVariables; ++irow) {
		for (int icol = 0; icol < 2 * numOfInputVariables; ++icol) {
			complexDerivatives[irow][icol] = czero;// Zero clear
		}
	}
	for (int irowL = 0; irowL < numOfInputVariables; ++irowL) {
		std::complex<double> value = czero;
		// The 1st term
		for (int i = 0; i < numOfReferenceVariables; ++i) {
			value -= PInvTIMatrix[irowL][i] * sumVeceh3[i];
		}
		// The 2nd term
		for (int i = 0; i < numOfInputVariables * numOfReferenceVariables; ++i) {
			value += IQPInvTPInvMatrix[irowL][i] * sumVechh3[i];
		}
		complexDerivatives[irowL][1] = value;
	}
	delete[] sumVeceh3;
	delete[] sumVechh3;
	for (int irow = 0; irow < numOfInputVariables; ++irow) {
		const int index = 2 * irow;
		derivativesRegardingScale[index   ] = complexDerivatives[irow][1].real();
		derivativesRegardingScale[index + 1] = complexDerivatives[irow][1].imag();
	}
#ifdef _DEBUG_WRITE
	std::cout << "drs1" << std::endl;
	std::cout << "[";
	for (int row = 0; row < 2 * numOfInputVariables; ++row) {
		std::cout << derivativesRegardingScale[row] << " ";
		if (row + 1 < 2 * numOfInputVariables) {
			std::cout << ",";
		}
	}
	std::cout << "]" << std::endl;
	std::complex<double>* derivP3 = new std::complex<double>[numOfInputVariables * numOfReferenceVariables];
	std::complex<double>* derivInvP3 = new std::complex<double>[numOfInputVariables * numOfReferenceVariables];
	std::complex<double>* derivQ3 = new std::complex<double>[numOfReferenceVariables ];
	double* derivResp3 = new double[2 * numOfInputVariables];
	{
		std::complex<double>* residualsOrg = new std::complex<double>[numSegments];
		double* weightsOrg = new double[numSegments];
		double* weightsMod = new double[numSegments];
		const double scaleOrg = scale;
		const double dscale = scale * 0.001;
		const double scaleMod = scale + dscale;
		// Calculate complex residuals
		for (int iSeg = 0; iSeg < numSegments; ++iSeg) {
			const std::complex<double> outSyn = resp[0] * ftval[in0][iSeg] + resp[1] * ftval[in1][iSeg];
			residualsOrg[iSeg] = ftval[out][iSeg] - outSyn;
		}
		// Calculate original weights
		for (int iSeg = 0; iSeg < numSegments; ++iSeg) {
			const double x = std::abs(residualsOrg[iSeg]) / scaleOrg;
			weightsOrg[iSeg] = RobustWeightTukeysBiweights::calculateWeights(x, threshould);
		}
		for (int iSeg = 0; iSeg < numSegments; ++iSeg) {
			const double x = std::abs(residualsOrg[iSeg]) / scaleMod;
			weightsMod[iSeg] = RobustWeightTukeysBiweights::calculateWeights(x, threshould);
		}
		std::complex<double> PMatrixOrg[2][2] = { czero, czero, czero, czero };
		for (int iSeg = 0; iSeg < numSegments; ++iSeg) {
			PMatrixOrg[0][0] += ftval[in0][iSeg] * std::conj(ftval[rr0][iSeg]) * weightsOrg[iSeg];
			PMatrixOrg[1][0] += ftval[in1][iSeg] * std::conj(ftval[rr0][iSeg]) * weightsOrg[iSeg];
			PMatrixOrg[0][1] += ftval[in0][iSeg] * std::conj(ftval[rr1][iSeg]) * weightsOrg[iSeg];
			PMatrixOrg[1][1] += ftval[in1][iSeg] * std::conj(ftval[rr1][iSeg]) * weightsOrg[iSeg];
		}
		std::complex<double> PMatrixMod[2][2] = { czero, czero, czero, czero };
		for (int iSeg = 0; iSeg < numSegments; ++iSeg) {
			PMatrixMod[0][0] += ftval[in0][iSeg] * std::conj(ftval[rr0][iSeg]) * weightsMod[iSeg];
			PMatrixMod[1][0] += ftval[in1][iSeg] * std::conj(ftval[rr0][iSeg]) * weightsMod[iSeg];
			PMatrixMod[0][1] += ftval[in0][iSeg] * std::conj(ftval[rr1][iSeg]) * weightsMod[iSeg];
			PMatrixMod[1][1] += ftval[in1][iSeg] * std::conj(ftval[rr1][iSeg]) * weightsMod[iSeg];
		}
		derivP3[0] = (PMatrixMod[0][0] - PMatrixOrg[0][0]) / dscale;
		derivP3[1] = (PMatrixMod[1][0] - PMatrixOrg[1][0]) / dscale;
		derivP3[2] = (PMatrixMod[0][1] - PMatrixOrg[0][1]) / dscale;
		derivP3[3] = (PMatrixMod[1][1] - PMatrixOrg[1][1]) / dscale;
		const std::complex<double> detOrg = PMatrixOrg[0][0] * PMatrixOrg[1][1] - PMatrixOrg[1][0] * PMatrixOrg[0][1];
		const std::complex<double> PInvMatrixOrg[2][2] = { PMatrixOrg[1][1] / detOrg, -PMatrixOrg[0][1] / detOrg, -PMatrixOrg[1][0] / detOrg, PMatrixOrg[0][0] / detOrg };
		const std::complex<double> detMod = PMatrixMod[0][0] * PMatrixMod[1][1] - PMatrixMod[1][0] * PMatrixMod[0][1];
		const std::complex<double> PInvMatrixMod[2][2] = { PMatrixMod[1][1] / detMod, -PMatrixMod[0][1] / detMod, -PMatrixMod[1][0] / detMod, PMatrixMod[0][0] / detMod };
		derivInvP3[0] = (PInvMatrixMod[0][0] - PInvMatrixOrg[0][0]) / dscale;
		derivInvP3[1] = (PInvMatrixMod[1][0] - PInvMatrixOrg[1][0]) / dscale;
		derivInvP3[2] = (PInvMatrixMod[0][1] - PInvMatrixOrg[0][1]) / dscale;
		derivInvP3[3] = (PInvMatrixMod[1][1] - PInvMatrixOrg[1][1]) / dscale;
		std::complex<double>* QMatrixOrg = new std::complex<double>[numOfReferenceVariables];
		std::complex<double>* QMatrixMod = new std::complex<double>[numOfReferenceVariables];
		for (int icol = 0; icol < numOfReferenceVariables; ++icol) {
			QMatrixOrg[icol] = czero;
			QMatrixMod[icol] = czero;
		}
		for (int irr = 0; irr < numOfReferenceVariables; ++irr) {
			const int rr = ptrControl->getChannelIndex(CommonParameters::REMOTE_REFERENCE, irr);
			for (int iSeg = 0; iSeg < numSegments; ++iSeg) {
				QMatrixOrg[irr] += ftval[out][iSeg] * std::conj(ftval[rr][iSeg]) * weightsOrg[iSeg];
				QMatrixMod[irr] += ftval[out][iSeg] * std::conj(ftval[rr][iSeg]) * weightsMod[iSeg];
			}
		}
		for (int irr = 0; irr < numOfReferenceVariables; ++irr) {
			derivQ3[irr] = (QMatrixMod[irr] - QMatrixOrg[irr]) / dscale;
		}
		for (int icol2 = 0; icol2 < numOfInputVariables; ++icol2) {
			std::complex<double> valueOrg = czero;
			std::complex<double> valueMod = czero;
			for (int i = 0; i < numOfReferenceVariables; ++i) {
				valueOrg += QMatrixOrg[i] * PInvMatrixOrg[i][icol2];
				valueMod += QMatrixMod[i] * PInvMatrixMod[i][icol2];
			}
			const int index2 = 2 * icol2;
			derivResp3[index2    ] = (valueMod.real() - valueOrg.real()) / dscale;
			derivResp3[index2 + 1] = (valueMod.imag() - valueOrg.imag()) / dscale;
		}
		delete[] QMatrixOrg;
		delete[] QMatrixMod;
		delete[] residualsOrg;
		delete[] weightsOrg;
		delete[] weightsMod;
	}
	std::cout << "[";
	for( int row = 0; row < numOfInputVariables * numOfReferenceVariables; ++row ){
		std::cout << derivP3[row].real() << "+" << derivP3[row].imag() <<"im ";
		if( row+1 < numOfInputVariables * numOfReferenceVariables ){
			std::cout << ",";
		}
	}
	std::cout << "]" << std::endl;
	std::cout << "[";
	for( int row = 0; row < numOfInputVariables * numOfReferenceVariables; ++row ){
		std::cout << derivInvP3[row].real() << "+" << derivInvP3[row].imag() <<"im ";
		if( row+1 < numOfInputVariables * numOfReferenceVariables ){
			std::cout << ",";
		}
	}
	std::cout << "]" << std::endl;
	std::cout << "[";
	for( int row = 0; row < numOfReferenceVariables ; ++row ){
		std::cout << derivQ3[row].real() << "+" << derivQ3[row].imag() <<"im ";
		if( row+1 < numOfReferenceVariables ){
			std::cout << ",";
		}
	}
	std::cout << "]" << std::endl;
	std::cout << "drs3" << std::endl;
	std::cout << "[";
	for (int row = 0; row < 2 * numOfInputVariables; ++row) {
		std::cout << derivResp3[row] << " ";
		if (row + 1 < 2 * numOfInputVariables) {
			std::cout << ",";
		}
	}
	std::cout << "]" << std::endl;
	delete[] derivP3;
	delete[] derivInvP3;
	delete[] derivQ3;
	delete[] derivResp3;
#endif
	for (int irow = 0; irow < numOfInputVariables; ++irow) {
		delete[] IQPInvTPInvMatrix[irow];
	}
	delete[] IQPInvTPInvMatrix;
	for (int irow = 0; irow < numOfInputVariables; ++irow) {
		delete[] complexDerivatives[irow];
	}
	delete[] complexDerivatives;
}

// Calculate partial derivatives of scale
void AnalysisOrdinaryRemoteReference::calculatePartialDerivativesOfScale(const int numSegments, const double threshould, const double paramB,
	std::complex<double>** ftval, const std::complex<double> resp0, const std::complex<double> resp1, const double scale, 
	const std::complex<double>* const residuals, const double* const weights, const int iOut, double* derivativesRegardingResps, 
	double& derivativesRegardingScale) const {

	const std::complex<double> resp[2] = { resp0, resp1 };

	const Control* const ptrControl = Control::getInstance();
	const int out = ptrControl->getChannelIndex(CommonParameters::OUTPUT, iOut);
	const int numOfInputVariables = ptrControl->getNumInputVariables();
	assert(numOfInputVariables == 2);
	const int in0 = ptrControl->getChannelIndex(CommonParameters::INPUT, 0);
	const int in1 = ptrControl->getChannelIndex(CommonParameters::INPUT, 1);

	const std::complex<double> czero = std::complex<double>(0.0, 0.0);

	double* value1 = new double[2 * numOfInputVariables];
	for (int icol = 0; icol < 2 * numOfInputVariables; ++icol) {
		value1[icol] = 0.0;// Zero clear
	}
	double value3(0.0);
	double* vecSigma = new double[2 * numOfInputVariables];
	std::complex<double>* hSigmaMatrix = new std::complex<double>[numOfInputVariables];
	for (int iSeg = 0; iSeg < numSegments; ++iSeg) {
		calculateVectorForPartialDerivatives(iSeg, ftval, residuals, hSigmaMatrix, vecSigma);
		for (int icol = 0; icol < 2 * numOfInputVariables; ++icol) {
			double temp = vecSigma[icol];
			if (std::abs(residuals[iSeg]) < CommonParameters::EPS) {
				// Nothing to do
			}
			else {
				temp *= weights[iSeg];
			}
			value1[icol] -= temp;
		}
		if (std::abs(residuals[iSeg]) < CommonParameters::EPS) {
			// Nothing to do
		}
		else {
			const double val = std::abs(residuals[iSeg]) / scale;
			const double temp = 2.0 * pow(scale, 2) * RobustWeightTukeysBiweights::calculateLossFunction(val, threshould)
				- std::norm(residuals[iSeg]) * weights[iSeg];
			value3 += temp / scale;
		}
	}
	delete[] vecSigma;
	delete[] hSigmaMatrix;
#ifdef _DEBUG_WRITE
	std::cout << "[";
	for( int row = 0; row < 2 * numOfInputVariables; ++row ){
		const int isign = row / numOfInputVariables;
		const int index = 2 * ( row % numOfInputVariables ) + isign;
		std::cout << value1[index] <<" ";
		if( row+1 < 2 * numOfInputVariables){
			std::cout << ",";
		}
	}
	std::cout << "]" << std::endl;
#endif
	// Order: Re(Z11),Im(Z11),Re(Z12),Im(Z12)
	const double factor = 1.0 / (2.0 * scale * static_cast<double>(numSegments) * paramB);
	for (int icol = 0; icol < 2 * numOfInputVariables; ++icol) {
		if (icol < numOfInputVariables) {
			// Real part
			const int ir2 = icol;
			const int index2 = 2 * ir2;
			derivativesRegardingResps[index2] = factor * value1[icol];
		}
		else {
			// Imaginary part
			const int ir2 = icol - numOfInputVariables;
			const int index2 = 2 * ir2 + 1;
			derivativesRegardingResps[index2] = factor * value1[icol];
		}
	}
	derivativesRegardingScale = factor * value3;
	delete[] value1;
#ifdef _DEBUG_WRITE
	std::cout << "dsr1" << std::endl;
	std::cout << "[";
	for (int row = 0; row < 2 * numOfInputVariables; ++row) {
		std::cout << derivativesRegardingResps[row] << " ";
		if (row + 1 < 2 * numOfInputVariables) {
			std::cout << ",";
		}
	}
	std::cout << "]" << std::endl;
	std::cout << "dss1" << std::endl;
	std::cout << derivativesRegardingScale << std::endl;
	double** derivSumU1 = new double* [numSegments];
	for (int iSeg = 0; iSeg < numSegments; ++iSeg) {
		derivSumU1[iSeg] = new double[2 * numOfInputVariables];
	}
	double* derivSquaredScale1 = new double[2 * numOfInputVariables];
	double* derivScale1 = new double[2 * numOfInputVariables];
	for (int icol = 0; icol < numOfInputVariables; ++icol) {
		for (int isign = 0; isign < 2; ++isign) {
			std::complex<double>* residualsOrg = new std::complex<double>[numSegments];
			std::complex<double>* respOrg = new std::complex<double>[numOfInputVariables];
			std::complex<double>* residualsMod = new std::complex<double>[numSegments];
			std::complex<double>* respMod = new std::complex<double>[numOfInputVariables];
			for (int icol2 = 0; icol2 < numOfInputVariables; ++icol2) {
				respOrg[icol2] = resp[icol2];
				respMod[icol2] = resp[icol2];
			}
			double dresp(0.0);
			if (isign == 0) {
				dresp = resp[icol].real() * 0.00001;
				respMod[icol] += std::complex<double>(dresp, 0.0);
			}
			else {
				dresp = resp[icol].imag() * 0.00001;
				respMod[icol] += std::complex<double>(0.0, dresp);
			}
			// Calculate complex residuals
			for (int iSeg = 0; iSeg < numSegments; ++iSeg) {
				const std::complex<double> outSyn = respOrg[0] * ftval[in0][iSeg] + respOrg[1] * ftval[in1][iSeg];
				residualsOrg[iSeg] = ftval[out][iSeg] - outSyn;
			}
			for (int iSeg = 0; iSeg < numSegments; ++iSeg) {
				const std::complex<double> outSyn = respMod[0] * ftval[in0][iSeg] + respMod[1] * ftval[in1][iSeg];
				residualsMod[iSeg] = ftval[out][iSeg] - outSyn;
			}
			const int index = 2 * icol + isign;
			for (int iSeg = 0; iSeg < numSegments; ++iSeg) {
				const double valOrg = std::abs(residualsOrg[iSeg]) / scale;
				const double valMod = std::abs(residualsMod[iSeg]) / scale;
				const double uOrg = RobustWeightTukeysBiweights::calculateLossFunction(valOrg, threshould) * pow(scale, 2);
				const double uMod = RobustWeightTukeysBiweights::calculateLossFunction(valMod, threshould) * pow(scale, 2);
				derivSumU1[iSeg][index] = (uMod - uOrg) / dresp;
			}
			double squareScaleOrg(0.0);
			for (int iSeg = 0; iSeg < numSegments; ++iSeg) {
				if (std::abs(residualsOrg[iSeg]) < CommonParameters::EPS) {
					// rho''(0) / 2 = 1 / 2
					squareScaleOrg += 0.5 * std::norm(residualsOrg[iSeg]);
				}
				else {
					const double val = std::abs(residualsOrg[iSeg]) / scale;
					squareScaleOrg += RobustWeightTukeysBiweights::calculateLossFunction(val, threshould) * pow(scale, 2);
				}
			}
			double squareScaleMod(0.0);
			for (int iSeg = 0; iSeg < numSegments; ++iSeg) {
				if (std::abs(residualsMod[iSeg]) < CommonParameters::EPS) {
					// rho''(0) / 2 = 1 / 2
					squareScaleMod += 0.5 * std::norm(residualsMod[iSeg]);
				}
				else {
					const double val = std::abs(residualsMod[iSeg]) / scale;
					squareScaleMod += RobustWeightTukeysBiweights::calculateLossFunction(val, threshould) * pow(scale, 2);
				}
			}
			derivSquaredScale1[index] = (squareScaleMod - squareScaleOrg) / dresp;
			squareScaleOrg /= static_cast<double>(numSegments);
			squareScaleOrg /= paramB;
			const double scaleOrg = sqrt(squareScaleOrg);
			squareScaleMod /= static_cast<double>(numSegments);
			squareScaleMod /= paramB;
			const double scaleMod = sqrt(squareScaleMod);
			derivScale1[index] = (scaleMod - scaleOrg) / dresp;
			delete[] residualsOrg;
			delete[] respOrg;
			delete[] residualsMod;
			delete[] respMod;
		}
	}
	std::cout << "[";
	for( int row = 0; row < numSegments; ++row ){
		for( int col = 0; col < 2 * numOfInputVariables; ++col ){
			std::cout << derivSumU1[row][col] <<" ";
		}
		if( row+1 < numSegments ){
			std::cout << ";";
		}
	}
	std::cout << "]" << std::endl;
	std::cout << "[";
	for( int row = 0; row < 2 * numOfInputVariables; ++row ){
		std::cout << derivSquaredScale1[row] <<" ";
		if( row+1 < 2 * numOfInputVariables){
			std::cout << ",";
		}
	}
	std::cout << "]" << std::endl;
	std::cout << "dsr2" << std::endl;
	std::cout << "[";
	for (int row = 0; row < 2 * numOfInputVariables; ++row) {
		std::cout << derivScale1[row] << " ";
		if (row + 1 < 2 * numOfInputVariables) {
			std::cout << ",";
		}
	}
	std::cout << "]" << std::endl;
	for (int iSeg = 0; iSeg < numSegments; ++iSeg) {
		delete[] derivSumU1[iSeg];
	}
	delete[] derivSumU1;
	delete[] derivSquaredScale1;
	delete[] derivScale1;
	double derivScale3(0.0);
	{
		const double scaleOrg = scale;
		const double dscale = scale * 0.001;
		const double scaleMod = scale + dscale;
		// Calculate complex residuals
		std::complex<double>* residuals = new std::complex<double>[numSegments];
		for (int iSeg = 0; iSeg < numSegments; ++iSeg) {
			const std::complex<double> outSyn = resp[0] * ftval[in0][iSeg] + resp[1] * ftval[in1][iSeg];
			residuals[iSeg] = ftval[out][iSeg] - outSyn;
		}
		double squareScaleOrg(0.0);
		for (int iSeg = 0; iSeg < numSegments; ++iSeg) {
			if (std::abs(residuals[iSeg]) < CommonParameters::EPS) {
				// rho''(0) / 2 = 1 / 2
				squareScaleOrg += 0.5 * std::norm(residuals[iSeg]);
			}
			else {
				const double val = std::abs(residuals[iSeg]) / scaleOrg;
				squareScaleOrg += RobustWeightTukeysBiweights::calculateLossFunction(val, threshould) * pow(scaleOrg, 2);
			}
		}
		double squareScaleMod(0.0);
		for (int iSeg = 0; iSeg < numSegments; ++iSeg) {
			if (std::abs(residuals[iSeg]) < CommonParameters::EPS) {
				// rho''(0) / 2 = 1 / 2
				squareScaleMod += 0.5 * std::norm(residuals[iSeg]);
			}
			else {
				const double val = std::abs(residuals[iSeg]) / scaleMod;
				squareScaleMod += RobustWeightTukeysBiweights::calculateLossFunction(val, threshould) * pow(scaleMod, 2);
			}
		}
		squareScaleOrg /= static_cast<double>(numSegments);
		squareScaleOrg /= paramB;
		const double scaleNewOrg = sqrt(squareScaleOrg);
		squareScaleMod /= static_cast<double>(numSegments);
		squareScaleMod /= paramB;
		const double scaleNewMod = sqrt(squareScaleMod);
		derivScale3 = (scaleNewMod - scaleNewOrg) / dscale;
		delete[] residuals;
	}
	std::cout << "dss1" << std::endl;
	std::cout << derivScale3 << std::endl;
#endif
}

// Calculate response functions
void AnalysisOrdinaryRemoteReference::calculateResponseFunctions(const int iSegLen, const int freqDegree, const double timeLength, const double freq,
	const int numSegmentsTotal, std::complex<double>** ftval, const std::vector< std::pair<std::string, std::string> >& times,
	std::ofstream& ofsResp, std::ofstream& ofsRhoaPhs) {

	const Control* const ptrControl = Control::getInstance();
	const int numOutputVariables = ptrControl->getNumOutputVariables();
	std::complex<double>* resp0 = new std::complex<double>[numOutputVariables];
	std::complex<double>* resp1 = new std::complex<double>[numOutputVariables];
	double** hatDiagonals = new double* [numOutputVariables];
	for (int iOut = 0; iOut < numOutputVariables; ++iOut) {
		hatDiagonals[iOut] = new double[numSegmentsTotal];
	}

	ofsResp << std::setprecision(10) << std::scientific << freq;
	ofsResp << "," << std::setprecision(10) << std::scientific << 1.0 / freq;
	if (ptrControl->doesOutputApparentResistivityAndPhase()) {
		ofsRhoaPhs << std::setprecision(10) << std::scientific << freq;
		ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << 1.0 / freq;
	}

	// Estimate response functions
	calculateResponseFunctionsAux(iSegLen, freqDegree, timeLength, freq, numSegmentsTotal, ftval, times, ofsResp, ofsRhoaPhs, false, resp0, resp1, hatDiagonals);

	if (ptrControl->getErrorEstimationMethod() == Control::SUBSET_DELETION_JACKKNIFE) {
		subsetDeletionJackknife(iSegLen, freqDegree, timeLength, freq, numSegmentsTotal, ftval, times, ofsResp, ofsRhoaPhs, resp0, resp1, hatDiagonals);
	}
	else if (ptrControl->getErrorEstimationMethod() == Control::STRICT_BOOTSTRAP) {
		strictBootstrap(iSegLen, freqDegree, timeLength, freq, numSegmentsTotal, ftval, times, ofsResp, ofsRhoaPhs, resp0, resp1);
	}

	delete[] resp0;
	delete[] resp1;
	for (int iOut = 0; iOut < numOutputVariables; ++iOut) {
		delete[] hatDiagonals[iOut];
	}
	delete[] hatDiagonals;

}

void AnalysisOrdinaryRemoteReference::calculateResponseFunctionsAux(const int iSegLen, const int freqDegree, const double timeLength, const double freq,
	const int numSegmentsTotal, std::complex<double>** ftval, const std::vector< std::pair<std::string, std::string> >& times,
	std::ofstream& ofsResp, std::ofstream& ofsRhoaPhs, const bool forJackknife,	std::complex<double>* respOut0, std::complex<double>* respOut1, 
	double** hatDiagsOut) const {

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();

	ptrOutputFiles->writeLogMessage("Calculate response functions by ordinary remote reference");
	ptrOutputFiles->writeCvgMessage("================================================================================");
	ptrOutputFiles->writeCvgMessage("Now Frequency(Hz): " + Util::toString(freq) + ", Period(s): " + Util::toString(1.0 / freq));
	ptrOutputFiles->writeCvgMessage("================================================================================");
	const Control* const ptrControl = Control::getInstance();
	if (ptrControl->getErrorEstimationMethod() == Control::ROBUST_BOOTSTRAP) {
		if (getPointerToRobustWeight(0) == NULL || getPointerToRobustWeight(0)->getNameOfRobustWeight() != "Tukey's biweights") {
			ptrOutputFiles->writeErrorMessage("The first M-estimator should be Tukey's biweights for robust bootstap method");
		}
		if (getPointerToRobustWeight(1) != NULL) {
			ptrOutputFiles->writeErrorMessage("The second M-estimator should be null for robust bootstap method");
		}
	}

	const int in0 = ptrControl->getChannelIndex(CommonParameters::INPUT, 0);
	const int in1 = ptrControl->getChannelIndex(CommonParameters::INPUT, 1);
	const int numRemoteReferenceVariables = ptrControl->getNumRemoteReferenceVariables();
	assert(numRemoteReferenceVariables >= 2);
	const int rr0 = ptrControl->getChannelIndex(CommonParameters::REMOTE_REFERENCE, 0);
	const int rr1 = ptrControl->getChannelIndex(CommonParameters::REMOTE_REFERENCE, 1);
	const int numOutputVariables = ptrControl->getNumOutputVariables();

	double* unitWeights = new double[numSegmentsTotal];
	for (int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg) {
		unitWeights[iSeg] = 1.0;
	}
	double** weights = new double* [numOutputVariables];
	for (int iOutVar = 0; iOutVar < numOutputVariables; ++iOutVar) {
		weights[iOutVar] = new double[numSegmentsTotal];
		memcpy(weights[iOutVar], unitWeights, sizeof(double) * numSegmentsTotal);
	}
	double* errorMagnitudes0 = new double [numOutputVariables];
	double* errorMagnitudes1 = new double [numOutputVariables];
	double* scales = new double[numOutputVariables];
	for (int iOutVar = 0; iOutVar < numOutputVariables; ++iOutVar) {
		ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
		ptrOutputFiles->writeCvgAndLogMessage("Calculate response functions for output variable " + Util::toString(iOutVar));
		ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
		// Calculate response functions by the ordinary least square method
		double coherence(0.0);
		std::complex<double> resp0(0.0, 0.0);
		std::complex<double> resp1(0.0, 0.0);
		std::complex<double>* residuals = new std::complex<double>[numSegmentsTotal];
		const int out = ptrControl->getChannelIndex(CommonParameters::OUTPUT, iOutVar);
		ptrOutputFiles->writeCvgAndLogMessage("Calculate response functions by the ordinary least square method");
		calculateResponseFunctionByWLSRemoteReference(ftval[out], ftval[in0], ftval[in1], ftval[rr0], ftval[rr1],
			numSegmentsTotal, weights[iOutVar], residuals, resp0, resp1, coherence);
		std::vector<std::string> titles;
		std::vector<double>* outputValues = new std::vector<double>[numSegmentsTotal];
		const bool outputResidual = !forJackknife && ptrControl->getOutputLevel() >= 2;
		if (outputResidual) {
			titles.push_back("OLS");
			for (int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg) {
				outputValues[iSeg].push_back(residuals[iSeg].real());
				outputValues[iSeg].push_back(residuals[iSeg].imag());
				outputValues[iSeg].push_back(weights[iOutVar][iSeg]);
			}
		}
		double absMaxResidual(0.0);
		for (int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg) {
			const double absResidual = std::abs(residuals[iSeg]);
			if (absResidual > absMaxResidual) {
				absMaxResidual = absResidual;
			}
		}
		bool noRobust(false);
		scales[iOutVar] = RobustWeight::calculateScaleByMADN(numSegmentsTotal, residuals);
		if (absMaxResidual < CommonParameters::EPS) {
			// Such a case where the remote reference field is equivalent to the input field
			ptrOutputFiles->writeCvgMessage("Robust method is not performed because residuals is nearly zero");
			noRobust = true;
		}else{
			// Calculate response functions by regression using the first M-estimator
			calculateResponseFunctionsByIRWLSRemoteReference(0, ftval[out], ftval[in0], ftval[in1], ftval[rr0], ftval[rr1],
				numSegmentsTotal, false, scales[iOutVar], unitWeights, weights[iOutVar], residuals, resp0, resp1, coherence, titles, outputValues);
			// Calculate response functions by regression using the second M-estimator
			calculateResponseFunctionsByIRWLSRemoteReference(1, ftval[out], ftval[in0], ftval[in1], ftval[rr0], ftval[rr1],
				numSegmentsTotal, true, scales[iOutVar], unitWeights, weights[iOutVar], residuals, resp0, resp1, coherence, titles, outputValues);
		}
		if (!forJackknife && ptrControl->getOutputLevel() > 0) {
			// Output spectral density functions to cvg file
			const int out = ptrControl->getChannelIndex(CommonParameters::OUTPUT, iOutVar);
			outputSpectralDensityFunctionsToCvgFile(numSegmentsTotal, timeLength, ftval[out], ftval[in0], ftval[in1], weights[iOutVar]);
		}
		respOut0[iOutVar] = resp0;
		respOut1[iOutVar] = resp1;
		if (!forJackknife) {
			const double maxHatDiag = calculateDiagonalComponentsOfHatMatrix(numSegmentsTotal, in0, in1, ftval, weights[iOutVar], hatDiagsOut[iOutVar]);
			// Output results
			ofsResp << "," << std::setprecision(10) << std::scientific << resp0.real();
			ofsResp << "," << std::setprecision(10) << std::scientific << resp0.imag();
			ofsResp << "," << std::setprecision(10) << std::scientific << resp1.real();
			ofsResp << "," << std::setprecision(10) << std::scientific << resp1.imag();
			ofsResp << "," << std::setprecision(10) << std::scientific << coherence;
			if (ptrControl->doesOutputApparentResistivityAndPhase()) {
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcApparentResistivity(freq, resp0);
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcPhase(resp0);
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcApparentResistivity(freq, resp1);
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcPhase(resp1);
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << coherence;
			}
			if (outputResidual) {
				std::ostringstream oss;
				oss << "segm" << iSegLen << "_index" << freqDegree << "_output" << iOutVar << "_residuals.csv";
				writeResiduals(oss.str(), numSegmentsTotal, times, titles, outputValues);
			}
			if (ptrControl->getErrorEstimationMethod() == Control::PARAMETRIC) {
				parametricErrorEstimation(numSegmentsTotal, weights[iOutVar], ftval, residuals, scales[iOutVar], noRobust,errorMagnitudes0[iOutVar], errorMagnitudes1[iOutVar]);
			}
		}
		// Release memory
		delete[] residuals;
		delete[] outputValues;
	}

	const int typeOfErrorEstimationMethod = ptrControl->getErrorEstimationMethod();
	switch (typeOfErrorEstimationMethod) {
	case Control::FIXED_WEIGHTS_JACKKNIFE:
		// Fixed-weights jackknife
		fixedWeightsJackknife(freq, numSegmentsTotal, weights, ftval, ofsResp, ofsRhoaPhs, respOut0, respOut1, hatDiagsOut);
		break;
	case Control::FIXED_WEIGHTS_BOOTSTRAP:
		// Fixed-weights bootstrap
		// Go through
	case Control::ROBUST_BOOTSTRAP:
		// Robust bootstrap
		robustBootstrap(freq, numSegmentsTotal, weights, scales, ftval, ofsResp, ofsRhoaPhs, respOut0, respOut1);
		break;
	case Control::SUBSET_DELETION_JACKKNIFE:
		// Go through
	case Control::STRICT_BOOTSTRAP:
		break;
	case Control::PARAMETRIC:
		for (int iOutVar = 0; iOutVar < numOutputVariables; ++iOutVar) {
			const double dResp0 = errorMagnitudes0[iOutVar];
			const double dResp1 = errorMagnitudes1[iOutVar];
			ofsResp << "," << std::setprecision(10) << std::scientific << dResp0;
			ofsResp << "," << std::setprecision(10) << std::scientific << dResp1;
			if (ptrControl->doesOutputApparentResistivityAndPhase()) {
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcApparentResistivityError(freq, respOut0[iOutVar], dResp0);
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcPhaseError(respOut0[iOutVar], dResp0);
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcApparentResistivityError(freq, respOut1[iOutVar], dResp1);
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcPhaseError(respOut1[iOutVar], dResp1);
			}
		}
		ofsResp << std::endl;
		ofsResp.flush();
		if (ptrControl->doesOutputApparentResistivityAndPhase()) {
			ofsRhoaPhs << std::endl;
			ofsRhoaPhs.flush();
		}
		break;
	default:
		ptrOutputFiles->writeErrorMessage("Unsupported error estimation method : " + Util::toString(typeOfErrorEstimationMethod));
		break;
	}

	// Release memory
	delete[] unitWeights;
	for (int iOutVar = 0; iOutVar < numOutputVariables; ++iOutVar) {
		delete[] weights[iOutVar];
	}
	delete[] weights;
	delete[] errorMagnitudes0;
	delete[] errorMagnitudes1;
	delete[] scales;

}

// Calculate a vector for partial derivatives
void AnalysisOrdinaryRemoteReference::calculateVectorForPartialDerivatives(const int iSeg, std::complex<double>** ftval, const std::complex<double>* const residuals,
	std::complex<double>* hSigmaMatrix, double* vector) const {

	const Control* const ptrControl = Control::getInstance();
	const int numOfInputVariables = ptrControl->getNumInputVariables();

	for (int irow = 0; irow < numOfInputVariables; ++irow) {
		const int in = ptrControl->getChannelIndex(CommonParameters::INPUT, irow);
		hSigmaMatrix[irow] = std::conj(ftval[in][iSeg]) * residuals[iSeg];
	}
#ifdef _DEBUG_WRITE
	std::cout << "[";
	for( int row = 0; row < numOfInputVariables; ++row ){
		std::cout << hSigmaMatrix[row].real() << "+" << hSigmaMatrix[row].imag() <<"im ";
	}
	std::cout << "]" << std::endl;
#endif

	for (int irow = 0; irow < numOfInputVariables; ++irow) {
		const int index = irow;
		vector[index] = hSigmaMatrix[irow].real();
	}
	for (int irow = 0; irow < numOfInputVariables; ++irow) {
		const int index = irow + numOfInputVariables;
		vector[index] = hSigmaMatrix[irow].imag();
	}
#ifdef _DEBUG_WRITE
	std::cout << "[";
	for( int row = 0; row < 2 * numOfInputVariables; ++row ){
		std::cout << vector[row] <<" ";
		if( row+1 < 2 * numOfInputVariables){
			std::cout << ",";
		}
	}
	std::cout << "]" << std::endl;
#endif
}

// Perform fixed-weights jackknife
void AnalysisOrdinaryRemoteReference::fixedWeightsJackknife(const double freq, const int numSegmentsTotal, double** weightsOrg, std::complex<double>** ftval,
	std::ofstream& ofsResp, std::ofstream& ofsRhoaPhs, const std::complex<double>* const resp0, const std::complex<double>* const resp1, double** hatDiagonals) const{

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeLogMessage("Fixed-weights jackknife is performed to estimate errors");

	const Control* const ptrControl = Control::getInstance();
	const int numOutputVariables = ptrControl->getNumOutputVariables();
	const int in0 = ptrControl->getChannelIndex(CommonParameters::INPUT, 0);
	const int in1 = ptrControl->getChannelIndex(CommonParameters::INPUT, 1);
	const int numRemoteReferenceVariables = ptrControl->getNumRemoteReferenceVariables();
	assert(numRemoteReferenceVariables >= 2);
	const int rr0 = ptrControl->getChannelIndex(CommonParameters::REMOTE_REFERENCE, 0);
	const int rr1 = ptrControl->getChannelIndex(CommonParameters::REMOTE_REFERENCE, 1);

	double** weights = new double* [numOutputVariables];
	for (int iOut = 0; iOut < numOutputVariables; ++iOut) {
		weights[iOut] = new double[numSegmentsTotal];
		memcpy(weights[iOut], weightsOrg[iOut], sizeof(double) * numSegmentsTotal);
	}

	std::complex<double>** pseudoResp0 = new std::complex<double>*[numSegmentsTotal];
	std::complex<double>** pseudoResp1 = new std::complex<double>*[numSegmentsTotal];
	for (int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg) {
		pseudoResp0[iSeg] = new std::complex<double>[numOutputVariables];
		pseudoResp1[iSeg] = new std::complex<double>[numOutputVariables];
		for (int iOut = 0; iOut < numOutputVariables; ++iOut) {
			weights[iOut][iSeg] = 0.0;// Replace
			const int out = ptrControl->getChannelIndex(CommonParameters::OUTPUT, iOut);
			std::complex<double> temp0(0.0, 0.0);
			std::complex<double> temp1(0.0, 0.0);
			calculateResponseFunctionByWLSRemoteReferenceAux(numSegmentsTotal,	ftval[out], ftval[in0], ftval[in1], ftval[rr0], ftval[rr1], weights[iOut], temp0, temp1);
			const double hatMatrixDiagonal = hatDiagonals[iOut][iSeg];
			double factor = static_cast<double>(numSegmentsTotal) * (1.0 - hatMatrixDiagonal);
			if (hatMatrixDiagonal > 1.0) {
				factor = 0.0;
			}
			else if (hatMatrixDiagonal < 0.0) {
				factor = static_cast<double>(numSegmentsTotal);
			}
			pseudoResp0[iSeg][iOut] = resp0[iOut] + factor * (resp0[iOut] - temp0);
			pseudoResp1[iSeg][iOut] = resp1[iOut] + factor * (resp1[iOut] - temp1);
			weights[iOut][iSeg] = weightsOrg[iOut][iSeg];// Restore
		}
	}
	// Release memory
	for (int iOut = 0; iOut < numOutputVariables; ++iOut) {
		delete[] weights[iOut];
	}
	delete[] weights;

	// Calculate & output error bars
	for (int iOut = 0; iOut < numOutputVariables; ++iOut) {
		if (numSegmentsTotal > 2) {
			std::complex<double> avgResp0(0.0, 0.0);
			std::complex<double> avgResp1(0.0, 0.0);
			for (int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg) {
				avgResp0 += pseudoResp0[iSeg][iOut];
				avgResp1 += pseudoResp1[iSeg][iOut];
			}
			const double factor = 1.0 / static_cast<double>(numSegmentsTotal);
			avgResp0 *= factor;
			avgResp1 *= factor;
			double variance0(0.0);
			double variance1(0.0);
			for (int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg) {
				variance0 += std::norm(pseudoResp0[iSeg][iOut] - avgResp0);
				variance1 += std::norm(pseudoResp1[iSeg][iOut] - avgResp1);
			}
			const double factor2 = factor / static_cast<double>(2 * numSegmentsTotal - 4);
			variance0 *= factor2;
			variance1 *= factor2;
			const double dResp0 = sqrt(variance0);
			const double dResp1 = sqrt(variance1);
			ofsResp << "," << std::setprecision(10) << std::scientific << dResp0;
			ofsResp << "," << std::setprecision(10) << std::scientific << dResp1;
			if (ptrControl->doesOutputApparentResistivityAndPhase()) {
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcApparentResistivityError(freq, resp0[iOut], dResp0);
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcPhaseError(resp0[iOut], dResp0);
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcApparentResistivityError(freq, resp1[iOut], dResp1);
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcPhaseError(resp1[iOut], dResp1);
			}
		}
		else {
			ofsResp << "," << std::setprecision(10) << std::scientific << 1.0e10;
			ofsResp << "," << std::setprecision(10) << std::scientific << 1.0e10;
			if (ptrControl->doesOutputApparentResistivityAndPhase()) {
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << 1.0e10;
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << 180.0;
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << 1.0e10;
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << 180.0;
			}
		}
	}
	ofsResp << std::endl;
	ofsResp.flush();
	if (ptrControl->doesOutputApparentResistivityAndPhase()) {
		ofsRhoaPhs << std::endl;
		ofsRhoaPhs.flush();
	}

	// Release memory
	for (int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg) {
		delete[] pseudoResp0[iSeg];
		delete[] pseudoResp1[iSeg];
	}
	delete[] pseudoResp0;
	delete[] pseudoResp1;

}

// Perform fixed-weights bootstrap
void AnalysisOrdinaryRemoteReference::fixedWeightsBootstrap(const double freq, const int numSegmentsTotal, double** weights, std::complex<double>** ftval,
	std::ofstream& ofsResp, std::ofstream& ofsRhoaPhs, const std::complex<double>* const resp0, const std::complex<double>* const resp1) const {

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeLogMessage("Fixed-weights bootstrap is performed to estimate errors");

	const Control* const ptrControl = Control::getInstance();
	const int numOutputVariables = ptrControl->getNumOutputVariables();
	const int in0 = ptrControl->getChannelIndex(CommonParameters::INPUT, 0);
	const int in1 = ptrControl->getChannelIndex(CommonParameters::INPUT, 1);
	const int numRemoteReferenceVariables = ptrControl->getNumRemoteReferenceVariables();
	assert(numRemoteReferenceVariables >= 2);
	const int rr0 = ptrControl->getChannelIndex(CommonParameters::REMOTE_REFERENCE, 0);
	const int rr1 = ptrControl->getChannelIndex(CommonParameters::REMOTE_REFERENCE, 1);

	int* segmentIndexes = new int[numSegmentsTotal];
	const int numOfSamples = ptrControl->getNumRepetitionsOfBootstrap();
	std::complex<double>** resp0Sample = new std::complex<double>*[numOfSamples];
	std::complex<double>** resp1Sample = new std::complex<double>*[numOfSamples];
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
	for (int iSample = 0; iSample < numOfSamples; ++iSample) {
		// Make bootstrap samples
		for (int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg) {
#ifdef _RAND
			segmentIndexes[iSeg] = (rand() / RAND_MAX) * (numSegmentsTotal - 1);
#else
#ifdef _MERSENNE_TWISTER_ORIGINAL
			segmentIndexes[iSeg] = static_cast<int>(genrand64_real1() * numSegmentsTotal);
#else
			segmentIndexes[iSeg] = uniformDistibution(gen);
#endif
#endif
		}
		resp0Sample[iSample] = new std::complex<double>[numOutputVariables];
		resp1Sample[iSample] = new std::complex<double>[numOutputVariables];
		for (int iOut = 0; iOut < numOutputVariables; ++iOut) {
			const int out = ptrControl->getChannelIndex(CommonParameters::OUTPUT, iOut);
			std::complex<double> temp0(0.0, 0.0);
			std::complex<double> temp1(0.0, 0.0);
			calculateResponseFunctionByWLSRemoteReferenceForBootstrap(numSegmentsTotal, segmentIndexes,
				ftval[out], ftval[in0], ftval[in1], ftval[rr0], ftval[rr1], weights[iOut], temp0, temp1);
			resp0Sample[iSample][iOut] = temp0;
			resp1Sample[iSample][iOut] = temp1;
		}
	}
	delete[] segmentIndexes;

	if (ptrControl->getOutputLevel() >= 4) {// Output one-step estimates
		for (int iOut = 0; iOut < numOutputVariables; ++iOut) {
			for (int iSample = 0; iSample < numOfSamples; ++iSample) {
				ptrOutputFiles->writeCvgMessage("Dataset " + Util::toString(iSample));
				ptrOutputFiles->writeCvgMessage("One-step estimates of response functions:");
				std::ostringstream msg;
				msg << "(" << std::setw(12) << std::setprecision(4) << std::scientific << resp0Sample[iSample][iOut].real() << ","
					<< std::setw(12) << std::setprecision(4) << std::scientific << resp0Sample[iSample][iOut].imag() << "), ";
				msg << "(" << std::setw(12) << std::setprecision(4) << std::scientific << resp1Sample[iSample][iOut].real() << ","
					<< std::setw(12) << std::setprecision(4) << std::scientific << resp1Sample[iSample][iOut].imag() << ")";
				ptrOutputFiles->writeCvgMessage(msg.str());
				ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
			}
		}
	}

	// Calculate error of response functions
	assert(numOfSamples > 2);
	for (int iOut = 0; iOut < numOutputVariables; ++iOut) {
		// Calculate average
		std::complex<double> average0 = std::complex<double>(0.0, 0.0);
		std::complex<double> average1 = std::complex<double>(0.0, 0.0);
		for (int iSample = 0; iSample < numOfSamples; ++iSample) {
			average0 += resp0Sample[iSample][iOut];
			average1 += resp1Sample[iSample][iOut];
		}
		average0 /= static_cast<double>(numOfSamples);
		average1 /= static_cast<double>(numOfSamples);
		// Calculate variance
		double variance0(0.0);
		double variance1(0.0);
		for (int iSample = 0; iSample < numOfSamples; ++iSample) {
			variance0 += std::norm(resp0Sample[iSample][iOut] - average0);
			variance1 += std::norm(resp1Sample[iSample][iOut] - average1);
		}
		variance0 /= static_cast<double>(2 * numOfSamples - 4);
		variance1 /= static_cast<double>(2 * numOfSamples - 4);
		const double dResp0 = sqrt(variance0);
		const double dResp1 = sqrt(variance1);
		ofsResp << "," << std::setprecision(10) << std::scientific << dResp0;
		ofsResp << "," << std::setprecision(10) << std::scientific << dResp1;
		if (ptrControl->doesOutputApparentResistivityAndPhase()) {
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcApparentResistivityError(freq, resp0[iOut], dResp0);
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcPhaseError(resp0[iOut], dResp0);
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcApparentResistivityError(freq, resp1[iOut], dResp1);
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcPhaseError(resp1[iOut], dResp1);
		}
	}
	ofsResp << std::endl;
	ofsResp.flush();
	if (ptrControl->doesOutputApparentResistivityAndPhase()) {
		ofsRhoaPhs << std::endl;
		ofsRhoaPhs.flush();
	}
	// Release memory
	for (int iSample = 0; iSample < numOfSamples; ++iSample) {
		delete[] resp0Sample[iSample];
		delete[] resp1Sample[iSample];
	}
	delete[] resp0Sample;
	delete[] resp1Sample;

}

// Estimate errors by parametric method
void AnalysisOrdinaryRemoteReference::parametricErrorEstimation(const int numSegmentsTotal, const double* const weights, std::complex<double>** ftval,
	const std::complex<double>* const residuals, const double scale, const bool noRobust, double& error0, double& error1) const {

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeLogMessage("Parametric error estimation is performed");

	if (numSegmentsTotal <= 2) {
		error0 = 1.0e10;
		error1 = 1.0e10;
		return;
	}

	int index(1);
	const RobustWeight* const robustWeight2nd = getPointerToRobustWeight(1);
	if (robustWeight2nd == NULL || robustWeight2nd->getNumIterationMax() < 1) {
		// The second M-estimater was not active
		index = 0;
		const RobustWeight* const robustWeight1st = getPointerToRobustWeight(0);
		if (robustWeight1st == NULL || robustWeight1st->getNumIterationMax() < 1) {
			// The 1st M-estimater was not active
			index = -1;
		}
	}

	double numerator(0.0);
	if (noRobust || index < 0) {
		double variance(0.0);
		for (int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg) {
			variance += std::norm(residuals[iSeg]);
		}
		variance /= static_cast<double>(2 * numSegmentsTotal - 4);
		numerator = variance;
	}
	else {
		for (int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg) {
			const double weightedAbsResidual = weights[iSeg] * std::abs(residuals[iSeg]);
			numerator += pow(weightedAbsResidual, 2);
		}
		numerator /= static_cast<double>(numSegmentsTotal - 2);
	}

	double denominator = 0.0;
	if (noRobust || index < 0) {
		denominator = 1.0;
	}
	else {
		const RobustWeight* const robustWeight = getPointerToRobustWeight(index);
		double sumOf2ndOrderDerivativeOfLossFunction = robustWeight->calculateSumOf2ndOrderDerivativeOfLossFunction(numSegmentsTotal, residuals, scale);
		sumOf2ndOrderDerivativeOfLossFunction /= static_cast<double>(numSegmentsTotal);
		denominator = pow(sumOf2ndOrderDerivativeOfLossFunction, 2);
		if (denominator < 1.0e-10) {
			denominator = 1.0e-10;
		}
	}

	const Control* const ptrControl = Control::getInstance();
	const int numOutputVariables = ptrControl->getNumOutputVariables();
	const int in0 = ptrControl->getChannelIndex(CommonParameters::INPUT, 0);
	const int in1 = ptrControl->getChannelIndex(CommonParameters::INPUT, 1);
	const int numRemoteReferenceVariables = ptrControl->getNumRemoteReferenceVariables();
	assert(numRemoteReferenceVariables >= 2);
	const int rr0 = ptrControl->getChannelIndex(CommonParameters::REMOTE_REFERENCE, 0);
	const int rr1 = ptrControl->getChannelIndex(CommonParameters::REMOTE_REFERENCE, 1);

	double BrxBrx(0.0);
	double BryBry(0.0);
	std::complex<double> BrxBry(0.0, 0.0);
	std::complex<double> BrxBx(0.0, 0.0);
	std::complex<double> BryBy(0.0, 0.0);
	std::complex<double> BrxBy(0.0, 0.0);
	std::complex<double> BryBx(0.0, 0.0);
	for (int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg) {
		BrxBrx += std::norm(ftval[rr0][iSeg]);
		BryBry += std::norm(ftval[rr1][iSeg]);
		BrxBry += std::conj(ftval[rr0][iSeg]) * ftval[rr1][iSeg];
		BrxBx  += std::conj(ftval[rr0][iSeg]) * ftval[in0][iSeg];
		BrxBy  += std::conj(ftval[rr0][iSeg]) * ftval[in1][iSeg];
		BryBx  += std::conj(ftval[rr1][iSeg]) * ftval[in0][iSeg];
		BryBy  += std::conj(ftval[rr1][iSeg]) * ftval[in1][iSeg];
	}
	const std::complex<double>   BrBr[2][2] = { {BrxBrx, BrxBry}, {std::conj(BrxBry), BryBry} };
	const std::complex<double>    BrB[2][2] = { {BrxBx,  BrxBy},  {BryBx,  BryBy} };
	const std::complex<double>    BBr[2][2] = { {std::conj(BrxBx), std::conj(BryBx)}, {std::conj(BrxBy), std::conj(BryBy)} };
	const std::complex<double> invBrB[2][2] = { {BryBy, -BrxBy},  {-BryBx, BrxBx} };
	const std::complex<double> invBBr[2][2] = { {std::conj(BryBy), -std::conj(BryBx)}, {-std::conj(BrxBy), std::conj(BrxBx)} };
	const std::complex<double> determinantOfBrB = BrxBx * BryBy - BrxBy * BryBx;
	const std::complex<double> determinantOfBBr = std::conj(BrxBx) * std::conj(BryBy) - std::conj(BryBx) * std::conj(BrxBy);
	const std::complex<double> invBrBBrBr[2][2] = {
		{ BryBy * BrxBrx - BrxBy * std::conj(BrxBry),  BryBy * BrxBry - BrxBy * BryBry},
		{-BryBx * BrxBrx + BrxBx * std::conj(BrxBry), -BryBx * BrxBry + BrxBx * BryBry}
	};
	const std::complex<double> invBrBBrBrinvBBrDiagonalXX = ( invBrBBrBr[0][0] * std::conj(BryBy) - invBrBBrBr[0][1] * std::conj(BrxBy)) / determinantOfBrB / determinantOfBBr;
	const std::complex<double> invBrBBrBrinvBBrDiagonalYY = (-invBrBBrBr[1][0] * std::conj(BryBx) + invBrBBrBr[1][1] * std::conj(BrxBx)) / determinantOfBrB / determinantOfBBr;
	error0 = sqrt( numerator / denominator * std::abs(invBrBBrBrinvBBrDiagonalXX) );
	error1 = sqrt( numerator / denominator * std::abs(invBrBBrBrinvBBrDiagonalYY) );

}

// Estimate error by robust bootstrap
void AnalysisOrdinaryRemoteReference::robustBootstrap(const double freq, const int numSegmentsTotal, double** weightsOrg, double* scalesOrg, std::complex<double>** ftval,
	std::ofstream& ofsResp, std::ofstream& ofsRhoaPhs, const std::complex<double>* const resp0Org, const std::complex<double>* const resp1Org) const{

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	const Control* const ptrControl = Control::getInstance();

	bool fixedWeights(true);
	const int typeOfErrorEstimationMethod = ptrControl->getErrorEstimationMethod();
	switch (typeOfErrorEstimationMethod) {
	case Control::FIXED_WEIGHTS_BOOTSTRAP:
		fixedWeights = true;
		ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
		ptrOutputFiles->writeCvgAndLogMessage("Estimate errors by fixed-weights bootstrap");
		ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
		break;
	case Control::ROBUST_BOOTSTRAP:
		fixedWeights = false;
		ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
		ptrOutputFiles->writeCvgAndLogMessage("Estimate errors by robust bootstrap");
		ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
		break;
	default:
		ptrOutputFiles->writeErrorMessage("Unsupported error estimation method : " + Util::toString(typeOfErrorEstimationMethod));
		break;
	}

	const int numOfOutputVariables = ptrControl->getNumOutputVariables();
	const int numOfInputVariables = ptrControl->getNumInputVariables();
	const int numOfReferenceVariables = ptrControl->getNumRemoteReferenceVariables();

	assert(numOfInputVariables == 2);
	const int in0 = ptrControl->getChannelIndex(CommonParameters::INPUT, 0);
	const int in1 = ptrControl->getChannelIndex(CommonParameters::INPUT, 1);
	assert(numOfReferenceVariables == 2);
	const int rr0 = ptrControl->getChannelIndex(CommonParameters::REMOTE_REFERENCE, 0);
	const int rr1 = ptrControl->getChannelIndex(CommonParameters::REMOTE_REFERENCE, 1);

	double paramB(0.0);
	double paramC(0.0);
	RobustWeightTukeysBiweights::calculateParams(2, numSegmentsTotal, paramB, paramC);

	for (int iOut = 0; iOut < numOfOutputVariables; ++iOut) {
		const int out = ptrControl->getChannelIndex(CommonParameters::OUTPUT, iOut);
		// Calculate complex residuals
		std::complex<double>* complexResidualsOrg = new std::complex<double>[numSegmentsTotal];
		for (int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg) {
			const std::complex<double> outSyn = resp0Org[iOut] * ftval[in0][iSeg] + resp1Org[iOut] * ftval[in1][iSeg];
			complexResidualsOrg[iSeg] = ftval[out][iSeg] - outSyn;
		}
		// Calculate original weights
		double* termsForScaleOrg = new double[numSegmentsTotal];
		for (int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg) {
			const double val = std::abs(complexResidualsOrg[iSeg]) / scalesOrg[iOut];
			if (fabs(val) < CommonParameters::EPS) {
				// rho''(0) / 2 = 1 / 2
				termsForScaleOrg[iSeg] = 0.5 * std::norm(complexResidualsOrg[iSeg]);
			}
			else {
				termsForScaleOrg[iSeg] = RobustWeightTukeysBiweights::calculateLossFunction(val, paramC) * pow(scalesOrg[iOut], 2);
			}
		}
		// Calculate partial derivatives for robust bootstrap
		const int dofOfMatrixForCorrection = 2 * numOfInputVariables + 1;
		DoubleDenseSquareMatrix matrixForCorrection;
		if (!fixedWeights) {
			matrixForCorrection.setDegreeOfEquation(dofOfMatrixForCorrection);
			matrixForCorrection.zeroClearMatrix();
			for (int irow = 0; irow < dofOfMatrixForCorrection; ++irow) {
				matrixForCorrection.setValue(irow, irow, 1.0);// Unit matrix
			}
#ifdef _DEBUG_WRITE
			matrixForCorrection.debugWriteMatrix();
#endif
			double** derivativesRespsResps = new double* [2 * numOfInputVariables];
			for (int irow = 0; irow < 2 * numOfInputVariables; ++irow) {
				derivativesRespsResps[irow] = new double[2 * numOfInputVariables];
			}
			double* derivativesRespsScale = new double[2 * numOfInputVariables];
			calculatePartialDerivativesOfResponses(numSegmentsTotal, paramC, ftval, resp0Org[iOut], resp1Org[iOut], scalesOrg[iOut],
				complexResidualsOrg, weightsOrg[iOut], iOut, derivativesRespsResps, derivativesRespsScale);
			for (int irow = 0; irow < 2 * numOfInputVariables; ++irow) {
				int icol(0);
				for (; icol < 2 * numOfInputVariables; ++icol) {
					matrixForCorrection.addValue(irow, icol, -derivativesRespsResps[irow][icol]);
				}
				matrixForCorrection.addValue(irow, icol, -derivativesRespsScale[irow]);
			}
#ifdef _DEBUG_WRITE
			matrixForCorrection.debugWriteMatrix();
#endif
			for (int irow = 0; irow < 2 * numOfReferenceVariables; ++irow) {
				delete[] derivativesRespsResps[irow];
			}
			delete[] derivativesRespsResps;
			delete[] derivativesRespsScale;

			double* derivativesScaleResps = new double[2 * numOfInputVariables];
			double derivativesScaleScale(0.0);
			calculatePartialDerivativesOfScale(numSegmentsTotal, paramC, paramB, ftval, resp0Org[iOut], resp1Org[iOut], scalesOrg[iOut],
				complexResidualsOrg, weightsOrg[iOut], iOut, derivativesScaleResps, derivativesScaleScale);
			{
				const int irow = dofOfMatrixForCorrection - 1;
				int icol(0);
				for (; icol < 2 * numOfInputVariables; ++icol) {
					matrixForCorrection.addValue(irow, icol, -derivativesScaleResps[icol]);
				}
				matrixForCorrection.addValue(irow, icol, -derivativesScaleScale);
			}
#ifdef _DEBUG_WRITE
			matrixForCorrection.debugWriteMatrix();
#endif
			delete[] derivativesScaleResps;
			// Factorize matrix
			matrixForCorrection.factorizeMatrix();
		}

		// Bootstrap
		int* segmentIndexes = new int[numSegmentsTotal];
		const int numOfSamples = ptrControl->getNumRepetitionsOfBootstrap();
		std::complex<double>* resp0Sample = new std::complex<double>[numOfSamples];
		std::complex<double>* resp1Sample = new std::complex<double>[numOfSamples];
		double* vectorForCorrection = new double[dofOfMatrixForCorrection];
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
		for (int iSample = 0; iSample < numOfSamples; ++iSample) {
			// Make bootstrap samples
			for (int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg) {
#ifdef _RAND
				segmentIndexes[iSeg] = (rand() / RAND_MAX) * (numSegmentsTotal - 1);
#else
#ifdef _MERSENNE_TWISTER_ORIGINAL
				segmentIndexes[iSeg] = static_cast<int>(genrand64_real1() * numSegmentsTotal);
#else
				segmentIndexes[iSeg] = uniformDistibution(gen);
#endif
#endif
			}
			// Calculate one-step estimates of response functions
			std::complex<double> resp0(0.0, 0.0);
			std::complex<double> resp1(0.0, 0.0);
			calculateResponseFunctionByWLSRemoteReferenceForBootstrap(numSegmentsTotal, segmentIndexes,
				ftval[out], ftval[in0], ftval[in1], ftval[rr0], ftval[rr1], weightsOrg[iOut], resp0, resp1);
			// Calculate one-step estimates of the scale
			double squareScale(0.0);
#ifdef _USE_OMP
			int icount(0);
#pragma omp parallel for default(shared) private(icount) reduction(+ : squareScale) 
			for (icount = 0; icount < numSegmentsTotal; ++icount) {
				squareScale += termsForScaleOrg[segmentIndexes[icount]];
			}
#else
			for (int icount = 0; icount < numSegmentsTotal; ++icount) {
				const int iSeg = segmentIndexes[icount];
				squareScale += termsForScaleOrg[iSeg];
			}
#endif
			squareScale /= static_cast<double>(numSegmentsTotal);
			squareScale /= paramB;
			const double scale = sqrt(squareScale);
			if (ptrControl->getOutputLevel() >= 4) {// Output one-step estimates
				ptrOutputFiles->writeCvgMessage("Dataset " + Util::toString(iSample));
				ptrOutputFiles->writeCvgMessage("One-step estimates of response functions:");
				std::ostringstream msg;
				msg << "(" << std::setw(12) << std::setprecision(4) << std::scientific << resp0.real() << ","
					<< std::setw(12) << std::setprecision(4) << std::scientific << resp0.imag() << "), ";
				msg << "(" << std::setw(12) << std::setprecision(4) << std::scientific << resp1.real() << ","
					<< std::setw(12) << std::setprecision(4) << std::scientific << resp1.imag() << ")";
				ptrOutputFiles->writeCvgMessage(msg.str());
				if (!fixedWeights) {
					ptrOutputFiles->writeCvgMessage("One-step estimate of scale: " + Util::toString(scale));
				}
			}
			if (!fixedWeights) {
				// Calculate differences between one-step estimates and original estimates
				// Order: Re(Z11),Im(Z11),Re(Z12),Im(Z12)
				vectorForCorrection[0] = resp0.real() - resp0Org[iOut].real();
				vectorForCorrection[1] = resp0.imag() - resp0Org[iOut].imag();
				vectorForCorrection[2] = resp1.real() - resp1Org[iOut].real();
				vectorForCorrection[3] = resp1.imag() - resp1Org[iOut].imag();
				vectorForCorrection[4] = scale - scalesOrg[iOut];
#if _DEBUG_WRITE
				std::cout << "[";
				for (int row = 0; row < dofOfMatrixForCorrection; ++row) {
					std::cout << vectorForCorrection[row] << " ";
					if (row + 1 < dofOfMatrixForCorrection) {
						std::cout << ",";
					}
				}
				std::cout << "]" << std::endl;
#endif
				// Solve linear equation
				matrixForCorrection.solveLinearEquation(vectorForCorrection, vectorForCorrection);
#if _DEBUG_WRITE
				std::cout << "[";
				for (int row = 0; row < dofOfMatrixForCorrection; ++row) {
					std::cout << vectorForCorrection[row] << " ";
					if (row + 1 < dofOfMatrixForCorrection) {
						std::cout << ",";
					}
				}
				std::cout << "]" << std::endl;
#endif
				// Calculate corrected version of the output one-step estimates
				// Order: Re(Z11),Im(Z11),Re(Z12),Im(Z12)
				resp0 = resp0Org[iOut] + std::complex<double>(vectorForCorrection[0], vectorForCorrection[1]);
				resp1 = resp1Org[iOut] + std::complex<double>(vectorForCorrection[2], vectorForCorrection[3]);
				if (ptrControl->getOutputLevel() >= 4) {
					//Output corrected version of the output one-step estimates
					ptrOutputFiles->writeCvgMessage("Corrected one-step estimates of response functions:");
					std::ostringstream msg;
					msg << "(" << std::setw(12) << std::setprecision(4) << std::scientific << resp0.real() << ","
						<< std::setw(12) << std::setprecision(4) << std::scientific << resp0.imag() << "), ";
					msg << "(" << std::setw(12) << std::setprecision(4) << std::scientific << resp1.real() << ","
						<< std::setw(12) << std::setprecision(4) << std::scientific << resp1.imag() << ")";
					ptrOutputFiles->writeCvgMessage(msg.str());
					ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
				}
			}
			resp0Sample[iSample] = resp0;
			resp1Sample[iSample] = resp1;
		}
		// Calculate error of response functions
		assert(numOfSamples > 2);
		// Calculate average
		std::complex<double> average0 = std::complex<double>(0.0, 0.0);
		std::complex<double> average1 = std::complex<double>(0.0, 0.0);
		for (int iSample = 0; iSample < numOfSamples; ++iSample) {
			average0 += resp0Sample[iSample];
			average1 += resp1Sample[iSample];
		}
		average0 /= static_cast<double>(numOfSamples);
		average1 /= static_cast<double>(numOfSamples);
		// Calculate variance
		double variance0(0.0);
		double variance1(0.0);
		for (int iSample = 0; iSample < numOfSamples; ++iSample) {
			variance0 += std::norm(resp0Sample[iSample] - average0);
			variance1 += std::norm(resp1Sample[iSample] - average1);
		}
		variance0 /= static_cast<double>(2 * numOfSamples - 4);
		variance1 /= static_cast<double>(2 * numOfSamples - 4);
		const double dResp0 = sqrt(variance0);
		const double dResp1 = sqrt(variance1);
		ofsResp << "," << std::setprecision(10) << std::scientific << dResp0;
		ofsResp << "," << std::setprecision(10) << std::scientific << dResp1;
		if (ptrControl->doesOutputApparentResistivityAndPhase()) {
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcApparentResistivityError(freq, resp0Org[iOut], dResp0);
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcPhaseError(resp0Org[iOut], dResp0);
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcApparentResistivityError(freq, resp1Org[iOut], dResp1);
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcPhaseError(resp1Org[iOut], dResp1);
		}
		// Release memory
		delete[] complexResidualsOrg;
		delete[] termsForScaleOrg;
		delete[] segmentIndexes;
		delete[] resp0Sample;
		delete[] resp1Sample;
		delete[] vectorForCorrection;
	}

	ofsResp << std::endl;
	ofsResp.flush();
	if (ptrControl->doesOutputApparentResistivityAndPhase()) {
		ofsRhoaPhs << std::endl;
		ofsRhoaPhs.flush();
	}

}

// Estimate error by strict bootstrap
void AnalysisOrdinaryRemoteReference::strictBootstrap(const int iSegLen, const int freqDegree, const double timeLength, const double freq,
	const int numSegmentsTotal, std::complex<double>** ftval, const std::vector< std::pair<std::string, std::string> >& times,
	std::ofstream& ofsResp, std::ofstream& ofsRhoaPhs, const std::complex<double>* const resp0, const std::complex<double>* const resp1) const {

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeLogMessage("Strict bootstrap is performed to estimate errors");

	const Control* const ptrControl = Control::getInstance();
	const int numOutputVariables = ptrControl->getNumOutputVariables();

	// Copy Fourier transformed values
	const int numChannels = ptrControl->getNumberOfChannels();
	std::complex<double>** ftvalForBootstrap = new std::complex<double>*[numChannels];
	for (int iChan = 0; iChan < numChannels; ++iChan) {
		ftvalForBootstrap[iChan] = new std::complex<double>[numSegmentsTotal];
	}
	const int numOfSamples = ptrControl->getNumRepetitionsOfBootstrap();
	std::complex<double>** resp0Sample = new std::complex<double>*[numOfSamples];
	std::complex<double>** resp1Sample = new std::complex<double>*[numOfSamples];
	int* segmentIndexes = new int[numSegmentsTotal];
	ptrOutputFiles->stopToWriteCvgMessage();
	ptrOutputFiles->stopToWriteLogMessage();
	ptrOutputFiles->stopToWriteWarningMessage();
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
	for (int iSample = 0; iSample < numOfSamples; ++iSample) {
		// Make bootstrap samples
		for (int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg) {
#ifdef _RAND
			segmentIndexes[iSeg] = (rand() / RAND_MAX) * (numSegmentsTotal - 1);
#else
#ifdef _MERSENNE_TWISTER_ORIGINAL
			segmentIndexes[iSeg] = static_cast<int>(genrand64_real1() * numSegmentsTotal);
#else
			segmentIndexes[iSeg] = uniformDistibution(gen);
#endif
#endif
		}
		// Copy data
		for (int iChan = 0; iChan < numChannels; ++iChan) {
			for (int icount = 0; icount < numSegmentsTotal; ++icount) {
				const int iSeg = segmentIndexes[icount];
				ftvalForBootstrap[iChan][icount] = ftval[iChan][iSeg];
			}
		}
		resp0Sample[iSample] = new std::complex<double>[numOutputVariables];
		resp1Sample[iSample] = new std::complex<double>[numOutputVariables];
		calculateResponseFunctionsAux(iSegLen, freqDegree, timeLength, freq, numSegmentsTotal, ftvalForBootstrap, times,
			ofsResp, ofsRhoaPhs, true, resp0Sample[iSample], resp1Sample[iSample], NULL);
	}
	ptrOutputFiles->restartToWriteCvgMessage();
	ptrOutputFiles->restartToWriteLogMessage();
	ptrOutputFiles->restartToWriteWarningMessage();

	if (ptrControl->getOutputLevel() >= 4) {// Output estimates of final response functions
		for (int iSample = 0; iSample < numOfSamples; ++iSample) {
			ptrOutputFiles->writeCvgMessage("Dataset " + Util::toString(iSample));
			ptrOutputFiles->writeCvgMessage("Estimates of final response functions:");
			std::ostringstream msg;
			for (int iOut = 0; iOut < numOutputVariables; ++iOut) {
				msg << "(" << std::setw(12) << std::setprecision(4) << std::scientific << resp0Sample[iSample][iOut].real() << ","
					<< std::setw(12) << std::setprecision(4) << std::scientific << resp0Sample[iSample][iOut].imag() << "), ";
				msg << "(" << std::setw(12) << std::setprecision(4) << std::scientific << resp1Sample[iSample][iOut].real() << ","
					<< std::setw(12) << std::setprecision(4) << std::scientific << resp1Sample[iSample][iOut].imag() << ")";
				if (iOut + 1 < numOutputVariables) {
					msg << std::endl;
				}
			}
			ptrOutputFiles->writeCvgMessage(msg.str());
			ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
		}
	}

	// Calculate error of response functions
	assert(numOfSamples > 2);
	for (int iOut = 0; iOut < numOutputVariables; ++iOut) {
		// Calculate average
		std::complex<double> average0 = std::complex<double>(0.0, 0.0);
		std::complex<double> average1 = std::complex<double>(0.0, 0.0);
		for (int iSample = 0; iSample < numOfSamples; ++iSample) {
			average0 += resp0Sample[iSample][iOut];
			average1 += resp1Sample[iSample][iOut];
		}
		average0 /= static_cast<double>(numOfSamples);
		average1 /= static_cast<double>(numOfSamples);
		// Calculate variance
		double variance0(0.0);
		double variance1(0.0);
		for (int iSample = 0; iSample < numOfSamples; ++iSample) {
			variance0 += std::norm(resp0Sample[iSample][iOut] - average0);
			variance1 += std::norm(resp1Sample[iSample][iOut] - average1);
		}
		variance0 /= static_cast<double>(2 * numOfSamples - 4);
		variance1 /= static_cast<double>(2 * numOfSamples - 4);
		const double dResp0 = sqrt(variance0);
		const double dResp1 = sqrt(variance1);
		ofsResp << "," << std::setprecision(10) << std::scientific << dResp0;
		ofsResp << "," << std::setprecision(10) << std::scientific << dResp1;
		if (ptrControl->doesOutputApparentResistivityAndPhase()) {
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcApparentResistivityError(freq, resp0[iOut], dResp0);
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcPhaseError(resp0[iOut], dResp0);
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcApparentResistivityError(freq, resp1[iOut], dResp1);
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcPhaseError(resp1[iOut], dResp1);
		}
	}
	ofsResp << std::endl;
	ofsResp.flush();
	if (ptrControl->doesOutputApparentResistivityAndPhase()) {
		ofsRhoaPhs << std::endl;
		ofsRhoaPhs.flush();
	}

	// Release memory
	for (int iChan = 0; iChan < numChannels; ++iChan) {
		delete[] ftvalForBootstrap[iChan];
	}
	delete[] ftvalForBootstrap;
	for (int iSample = 0; iSample < numOfSamples; ++iSample) {
		delete[] resp0Sample[iSample];
		delete[] resp1Sample[iSample];
	}
	delete[] resp0Sample;
	delete[] resp1Sample;
	delete[] segmentIndexes;

}

// Perform subset deletion jackknife
void AnalysisOrdinaryRemoteReference::subsetDeletionJackknife(const int iSegLen, const int freqDegree, const double timeLength, const double freq,
	const int numSegmentsTotal, std::complex<double>** ftval, const std::vector< std::pair<std::string, std::string> >& times,
	std::ofstream& ofsResp, std::ofstream& ofsRhoaPhs, const std::complex<double>* const resp0, const std::complex<double>* const resp1, double** hatDiagonals) const {

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeLogMessage("Subset deletion jackknife is performed to estimate errors");

	const Control* const ptrControl = Control::getInstance();
	int numOmittedData = static_cast<int>(0.01 * ptrControl->getPercentageOfOmmitedDataSubsetDeletionJackknife() * static_cast<double>(numSegmentsTotal));
	if (numOmittedData < 1) {
		numOmittedData = 1;
	}
	const int numOfSubsets = numSegmentsTotal / numOmittedData;

	ptrOutputFiles->writeLogMessage("Number of ommited data : " + Util::toString(numOmittedData));
	ptrOutputFiles->writeLogMessage("Number of subsets : " + Util::toString(numOfSubsets));
	ptrOutputFiles->stopToWriteCvgMessage();
	ptrOutputFiles->stopToWriteLogMessage();
	ptrOutputFiles->stopToWriteWarningMessage();

	const int numOutputVariables = ptrControl->getNumOutputVariables();

	// Copy Fourier transformed values
	const int numChannels = ptrControl->getNumberOfChannels();
	std::complex<double>** ftvalForJackknife = new std::complex<double>*[numChannels];
	for (int iChan = 0; iChan < numChannels; ++iChan) {
		ftvalForJackknife[iChan] = new std::complex<double>[numSegmentsTotal - numOmittedData];
	}
	std::complex<double>** pseudoResp0 = new std::complex<double>*[numOfSubsets];
	std::complex<double>** pseudoResp1 = new std::complex<double>*[numOfSubsets];
	for (int iSubset = 0; iSubset < numOfSubsets; ++iSubset) {
		double* averageHatDiags = new double[numOutputVariables];
		for (int iOut = 0; iOut < numOutputVariables; ++iOut) {
			// Zero clear
			averageHatDiags[iOut] = 0.0;
		}
		const int iSegOmitStart = iSubset * numOmittedData;
		const int iSegOmitEnd = iSegOmitStart + numOmittedData;
		assert(iSegOmitEnd <= numSegmentsTotal);
		for (int iSeg = iSegOmitStart; iSeg < iSegOmitEnd; ++iSeg) {
			for (int iOut = 0; iOut < numOutputVariables; ++iOut) {
				averageHatDiags[iOut] += hatDiagonals[iOut][iSeg];
			}
		}
		for (int iOut = 0; iOut < numOutputVariables; ++iOut) {
			averageHatDiags[iOut] /= static_cast<double>(numOmittedData);
		}
		for (int iChan = 0; iChan < numChannels; ++iChan) {
			int icount(0);
			for (int iSeg = 0; iSeg < iSegOmitStart; ++iSeg, ++icount) {
				ftvalForJackknife[iChan][icount] = ftval[iChan][iSeg];
			}
			for (int iSeg = iSegOmitEnd; iSeg < numSegmentsTotal; ++iSeg, ++icount) {
				ftvalForJackknife[iChan][icount] = ftval[iChan][iSeg];
			}
			assert(icount == numSegmentsTotal - numOmittedData);
		}
		std::complex<double>* temp0 = new std::complex<double>[numOutputVariables];
		std::complex<double>* temp1 = new std::complex<double>[numOutputVariables];
		calculateResponseFunctionsAux(iSegLen, freqDegree, timeLength, freq, numSegmentsTotal - numOmittedData,
			ftvalForJackknife, times, ofsResp, ofsRhoaPhs, true, temp0, temp1, NULL);
		pseudoResp0[iSubset] = new std::complex<double>[numOutputVariables];
		pseudoResp1[iSubset] = new std::complex<double>[numOutputVariables];
		for (int iOut = 0; iOut < numOutputVariables; ++iOut) {
			double factor = static_cast<double>(numOfSubsets) * (1.0 - averageHatDiags[iOut]);
			if (averageHatDiags[iOut] > 1.0) {
				factor = 0.0;
			}
			else if (averageHatDiags[iOut] < 0.0) {
				factor = static_cast<double>(numOfSubsets);
			}
			pseudoResp0[iSubset][iOut] = resp0[iOut] + factor * (resp0[iOut] - temp0[iOut]);
			pseudoResp1[iSubset][iOut] = resp1[iOut] + factor * (resp1[iOut] - temp1[iOut]);
		}
		delete[] averageHatDiags;
		delete[] temp0;
		delete[] temp1;
	}
	ptrOutputFiles->restartToWriteCvgMessage();
	ptrOutputFiles->restartToWriteLogMessage();
	ptrOutputFiles->restartToWriteWarningMessage();
	for (int iChan = 0; iChan < numChannels; ++iChan) {
		delete[] ftvalForJackknife[iChan];
	}
	delete[] ftvalForJackknife;

	// Calculate & output error bars
	for (int iOut = 0; iOut < numOutputVariables; ++iOut) {
		if (numOfSubsets > 2) {
			std::complex<double> avgResp0(0.0, 0.0);
			std::complex<double> avgResp1(0.0, 0.0);
			for (int iSubset = 0; iSubset < numOfSubsets; ++iSubset) {
				avgResp0 += pseudoResp0[iSubset][iOut];
				avgResp1 += pseudoResp1[iSubset][iOut];
			}
			const double factor = 1.0 / static_cast<double>(numOfSubsets);
			avgResp0 *= factor;
			avgResp1 *= factor;
			double variance0(0.0);
			double variance1(0.0);
			for (int iSubset = 0; iSubset < numOfSubsets; ++iSubset) {
				variance0 += std::norm(pseudoResp0[iSubset][iOut] - avgResp0);
				variance1 += std::norm(pseudoResp1[iSubset][iOut] - avgResp1);
			}
			const double factor2 = factor / static_cast<double>(2 * numOfSubsets - 4);
			variance0 *= factor2;
			variance1 *= factor2;
			const double dResp0 = sqrt(variance0);
			const double dResp1 = sqrt(variance1);
			ofsResp << "," << std::setprecision(10) << std::scientific << dResp0;
			ofsResp << "," << std::setprecision(10) << std::scientific << dResp1;
			if (ptrControl->doesOutputApparentResistivityAndPhase()) {
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcApparentResistivityError(freq, resp0[iOut], dResp0);
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcPhaseError(resp0[iOut], dResp0);
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcApparentResistivityError(freq, resp1[iOut], dResp1);
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcPhaseError(resp1[iOut], dResp1);
			}
		}
		else {
			ofsResp << "," << std::setprecision(10) << std::scientific << 1.0e10;
			ofsResp << "," << std::setprecision(10) << std::scientific << 1.0e10;
			if (ptrControl->doesOutputApparentResistivityAndPhase()) {
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << 1.0e10;
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << 180.0;
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << 1.0e10;
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << 180.0;
			}
		}
	}
	ofsResp << std::endl;
	ofsResp.flush();
	if (ptrControl->doesOutputApparentResistivityAndPhase()) {
		ofsRhoaPhs << std::endl;
		ofsRhoaPhs.flush();
	}

	for (int iSubset = 0; iSubset < numOfSubsets; ++iSubset) {
		delete[] pseudoResp0[iSubset];
		delete[] pseudoResp1[iSubset];
	}
	delete[] pseudoResp0;
	delete[] pseudoResp1;

}

// Write residuals
void AnalysisOrdinaryRemoteReference::writeResiduals(const std::string& fileName, const int numSegmentsTotal,
	const std::vector< std::pair<std::string, std::string> >& times, const std::vector<std::string>& titles, const std::vector<double>* outputValues) const {

	std::ofstream ofs;
	ofs.open(fileName.c_str(), std::ios::out);
	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	if (ofs.fail()) {
		ptrOutputFiles->writeErrorMessage("File open error : " + fileName);
	}
	ofs << "index";
	ofs << ",start_time,end_time";
	for (std::vector<std::string>::const_iterator itr = titles.begin(); itr != titles.end(); ++itr) {
		ofs << ",residual_real_" << *itr;
		ofs << ",residual_imag_" << *itr;
		ofs << ",weight_" << *itr;
	}
	ofs << std::endl;

	int index(0);
	for (int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg) {
		ofs << index;
		const std::string timeStart = times[iSeg].first;
		const std::string timeEnd = times[iSeg].second;
		ofs << "," << timeStart << "," << timeEnd;
		for (std::vector<double>::const_iterator itr = outputValues[index].begin(); itr != outputValues[index].end(); ++itr) {
			ofs << "," << std::setprecision(10) << std::scientific << *itr;
		}
		ofs << std::endl;
		++index;
	}
	ofs.close();

}