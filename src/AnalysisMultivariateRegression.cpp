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
#include "AnalysisMultivariateRegression.h"
#include "Control.h"
#include "OutputFiles.h"
#include "DoubleDenseSquareMatrix.h"
#include <algorithm>

#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <iomanip>
#include <set>
#include <vector>
#ifdef _MERSENNE_TWISTER_ORIGINAL
#else
#include <random>
#endif

#include "Util.h"

#ifdef _MERSENNE_TWISTER_ORIGINAL
#include "mt64.h"
#endif

#ifdef _USE_OMP
#include <omp.h>
#endif

// Default constructer
AnalysisMultivariateRegression::AnalysisMultivariateRegression():
	m_numOfOutputAndInputVariables(0),
	m_responseFunctions(NULL)
{
}

// Destructer
AnalysisMultivariateRegression::~AnalysisMultivariateRegression()
{
	for( int i = 0; i < m_numOfOutputAndInputVariables; ++i ){
		delete [] m_responseFunctions[i];
	}	
	delete m_responseFunctions;
}

// Calculate complex residuals
void AnalysisMultivariateRegression::calculateComplexResiduals( const int numSegmentsTotal, std::complex<double>** ftval, 
	std::complex<double>** resp, std::complex<double>** complexResiduals ) const{

	const Control* const ptrControl = Control::getInstance();
	const int numOfOutputVariables = ptrControl->getNumOutputVariables();
	const int numOfInputVariables = ptrControl->getNumInputVariables();
	const int rr0 = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, 0 );
	const int rr1 = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, 1 );

	int iVar(0);
	for( int iOut = 0; iOut < numOfOutputVariables; ++iOut ){
		const int index = ptrControl->getChannelIndex( CommonParameters::OUTPUT, iOut );
		for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
			complexResiduals[iVar][iSeg] = ftval[index][iSeg]
				- resp[iVar][0] * ftval[rr0][iSeg] - resp[iVar][1] * ftval[rr1][iSeg];
		}
		++iVar;
	}
	for( int iInp = 0; iInp < numOfInputVariables; ++iInp ){
		const int index = ptrControl->getChannelIndex( CommonParameters::INPUT, iInp );
		for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
			complexResiduals[iVar][iSeg] = ftval[index][iSeg]
				- resp[iVar][0] * ftval[rr0][iSeg] - resp[iVar][1] * ftval[rr1][iSeg];
		}
		++iVar;
	}
	assert( iVar == numOfOutputVariables + numOfInputVariables );

}

// Calculate partial derivatives of responses for robust bootstrap
void AnalysisMultivariateRegression::calculatePartialDerivativesOfResponses( const int numSegments, const double paramC,
		std::complex<double>** ftval, std::complex<double>** resp, const double* const variancesWithoutScale, const double scale, 
		std::complex<double>** complexResiduals, const double* const MD, const double* const weights,
		double** derivativesRegardingResps, double** derivativesRegardingVariancesWithoutScale, double* derivativesRegardingScale ) const{

	const Control* const ptrControl = Control::getInstance();
	const int numOfOutputVariables = ptrControl->getNumOutputVariables();
	const int numOfInputVariables = ptrControl->getNumInputVariables();
	const int numOfOutputAndInputVariables = numOfOutputVariables + numOfInputVariables;
	const int numOfReferenceVariables = ptrControl->getNumRemoteReferenceVariables();
	assert(numOfReferenceVariables == 2); 
	const int rr0 = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, 0 );
	const int rr1 = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, 1 );

	const std::complex<double> czero = std::complex<double>( 0.0, 0.0 );

	std::complex<double> PMatrix[2][2] = { czero, czero, czero, czero };
	for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
		PMatrix[0][0] += ftval[rr0][iSeg] * std::conj(ftval[rr0][iSeg]) * weights[iSeg];
		PMatrix[1][0] += ftval[rr1][iSeg] * std::conj(ftval[rr0][iSeg]) * weights[iSeg];
		PMatrix[0][1] += ftval[rr0][iSeg] * std::conj(ftval[rr1][iSeg]) * weights[iSeg];
		PMatrix[1][1] += ftval[rr1][iSeg] * std::conj(ftval[rr1][iSeg]) * weights[iSeg];
	}
	const std::complex<double> det = PMatrix[0][0] * PMatrix[1][1] - PMatrix[1][0] * PMatrix[0][1];
	const std::complex<double> PInvMatrix[2][2] = { PMatrix[1][1]/det, -PMatrix[0][1]/det, -PMatrix[1][0]/det, PMatrix[0][0]/det };
	const std::complex<double> PInvTMatrix[2][2] = { PInvMatrix[0][0], PInvMatrix[1][0], PInvMatrix[0][1], PInvMatrix[1][1] };
	std::complex<double>** PInvTIMatrix = new std::complex<double>*[numOfReferenceVariables * numOfOutputAndInputVariables];
	for( int irow = 0; irow < numOfReferenceVariables * numOfOutputAndInputVariables; ++irow ){
		PInvTIMatrix[irow] = new std::complex<double>[numOfReferenceVariables * numOfOutputAndInputVariables];
		for( int icol = 0; icol < numOfReferenceVariables * numOfOutputAndInputVariables; ++icol ){
			PInvTIMatrix[irow][icol] = czero;
		}
	}
	for( int irow = 0; irow < numOfReferenceVariables; ++irow ){
		for( int icol = 0; icol < numOfReferenceVariables; ++icol ){
			const std::complex<double> value = PInvTMatrix[irow][icol];
			for( int iq = 0; iq < numOfOutputAndInputVariables; ++iq ){
				const int irowOut = iq + irow * numOfOutputAndInputVariables;
				const int icolOut = iq + icol * numOfOutputAndInputVariables;
				PInvTIMatrix[irowOut][icolOut] = value;
			}
		}
	}
//#ifdef _DEBUG_WRITE
//	std::cout << "[";
//	for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
//		std::cout << MD[iSeg] <<" ";
//		if( iSeg+1 < numSegments ){
//			std::cout << ",";
//		}
//	}
//	std::cout << "]" << std::endl;
//	std::cout << "[";
//	for( int row = 0; row < numOfOutputAndInputVariables; ++row ){
//		std::cout << variancesWithoutScale[row] <<" ";
//		if( row+1 < numOfOutputAndInputVariables ){
//			std::cout << ",";
//		}
//	}
//	std::cout << "]" << std::endl;
//	std::cout << "[";
//	for( int row = 0; row < numOfReferenceVariables * numOfOutputAndInputVariables; ++row ){
//		for( int col = 0; col < numOfReferenceVariables * numOfOutputAndInputVariables; ++col ){
//			std::cout << PInvTIMatrix[row][col].real() << "+" << PInvTIMatrix[row][col].imag() <<"im ";
//		}
//		if( row+1 < numOfReferenceVariables * numOfOutputAndInputVariables ){
//			std::cout << ";";
//		}
//	}
//	std::cout << "]" << std::endl;
//#endif
//
	std::complex<double>** PInvTPInvMatrix = new std::complex<double>*[numOfReferenceVariables * numOfReferenceVariables];
	for( int irow = 0; irow < numOfReferenceVariables * numOfReferenceVariables; ++irow ){
		PInvTPInvMatrix[irow] = new std::complex<double>[numOfReferenceVariables * numOfReferenceVariables];
		for( int icol = 0; icol < numOfReferenceVariables * numOfReferenceVariables; ++icol ){
			PInvTPInvMatrix[irow][icol] = czero;
		}
	}
	for( int irow = 0; irow <numOfReferenceVariables; ++irow ){
		for( int icol = 0; icol < numOfReferenceVariables; ++icol ){
			const std::complex<double> factor = PInvTMatrix[irow][icol];
			for( int irow2 = 0; irow2 < numOfReferenceVariables; ++irow2 ){
				for( int icol2 = 0; icol2 < numOfReferenceVariables; ++icol2 ){
					const int irowOut = irow2 + irow * numOfReferenceVariables;
					const int icolOut = icol2 + icol * numOfReferenceVariables;
					PInvTPInvMatrix[irowOut][icolOut] = factor * PInvMatrix[irow2][icol2];
				}
			}
		}
	}
//#ifdef _DEBUG_WRITE
//	std::cout << "[";
//	for( int row = 0; row < numOfReferenceVariables * numOfReferenceVariables; ++row ){
//		for( int col = 0; col < numOfReferenceVariables * numOfReferenceVariables; ++col ){
//			std::cout << PInvTPInvMatrix[row][col].real() << "+" << PInvTPInvMatrix[row][col].imag() <<"im ";
//		}
//		if( row+1 < numOfReferenceVariables * numOfReferenceVariables ){
//			std::cout << ";";
//		}
//	}
//	std::cout << "]" << std::endl;
//#endif

	std::complex<double>** QMatrix = new std::complex<double>*[numOfOutputAndInputVariables];
	for( int irow = 0; irow < numOfOutputAndInputVariables; ++irow ){
		QMatrix[irow] = new std::complex<double>[numOfReferenceVariables];
		for( int icol = 0; icol < numOfReferenceVariables; ++icol ){
			QMatrix[irow][icol] = czero;
		}
	}
	for( int irr = 0; irr < numOfReferenceVariables; ++irr ){
		const int rr = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, irr );
		int iVar(0);
		for( int iOut = 0; iOut < numOfOutputVariables; ++iOut ){
			const int index = ptrControl->getChannelIndex( CommonParameters::OUTPUT, iOut );
			for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
				QMatrix[iVar][irr] += ftval[index][iSeg] * std::conj(ftval[rr][iSeg]) * weights[iSeg];
			}
			assert(index == iVar);
			++iVar;
		}
		for( int iInp = 0; iInp < numOfInputVariables; ++iInp ){
			const int index = ptrControl->getChannelIndex( CommonParameters::INPUT, iInp );
			for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
				QMatrix[iVar][irr] += ftval[index][iSeg] * std::conj(ftval[rr][iSeg]) * weights[iSeg];
			}
			assert(index == iVar);
			++iVar;
		}
	}

//#ifdef _DEBUG_WRITE
//	std::cout << "[";
//	for( int row = 0; row < numOfOutputAndInputVariables; ++row ){
//		for( int col = 0; col < numOfReferenceVariables; ++col ){
//			std::cout << QMatrix[row][col].real() << "+" << QMatrix[row][col].imag() <<"im ";
//		}
//		if( row+1 < numOfOutputAndInputVariables ){
//			std::cout << ";";
//		}
//	}
//	std::cout << "]" << std::endl;
//	std::cout << "[";
//	for( int row = 0; row < numOfOutputAndInputVariables; ++row ){
//		for( int col = 0; col < numOfReferenceVariables; ++col ){
//			std::cout << resp[row][col].real() << "+" << resp[row][col].imag() <<"im ";
//		}
//		if( row+1 < numOfOutputAndInputVariables ){
//			std::cout << ";";
//		}
//	}
//	std::cout << "]" << std::endl;
//#endif
//
	std::complex<double>** IQMatrix = new std::complex<double>*[numOfReferenceVariables * numOfOutputAndInputVariables];
	for( int irow = 0; irow < numOfReferenceVariables * numOfOutputAndInputVariables; ++irow ){
		IQMatrix[irow] = new std::complex<double>[numOfReferenceVariables * numOfReferenceVariables];
		for( int icol = 0; icol < numOfReferenceVariables * numOfReferenceVariables; ++icol ){
			IQMatrix[irow][icol] = czero;
		}
	}
	for( int irow = 0; irow < numOfReferenceVariables; ++irow ){
		for( int icol = 0; icol < numOfReferenceVariables; ++icol ){
			const double factor = irow == icol ? 1.0 : 0.0;
			for( int irow2 = 0; irow2 < numOfOutputAndInputVariables; ++irow2 ){\
				for( int icol2 = 0; icol2 < numOfReferenceVariables; ++icol2 ){
					const int irowOut = irow2 + irow * numOfOutputAndInputVariables;
					const int icolOut = icol2 + icol * numOfReferenceVariables;
					IQMatrix[irowOut][icolOut] = factor * QMatrix[irow2][icol2];
				}
			}
		}
	}
	for( int irow = 0; irow < numOfOutputAndInputVariables; ++irow ){
		delete [] QMatrix[irow];
	}
	delete [] QMatrix;

//#ifdef _DEBUG_WRITE
//	std::cout << "[";
//	for( int row = 0; row < numOfReferenceVariables * numOfOutputAndInputVariables; ++row ){
//		for( int col = 0; col < numOfReferenceVariables * numOfReferenceVariables; ++col ){
//			std::cout << IQMatrix[row][col].real() << "+" << IQMatrix[row][col].imag() <<"im ";
//		}
//		if( row+1 < numOfReferenceVariables * numOfOutputAndInputVariables ){
//			std::cout << ";";
//		}
//	}
//	std::cout << "]" << std::endl;
//#endif

	std::complex<double>** IQPInvTPInvMatrix = new std::complex<double>*[numOfReferenceVariables * numOfOutputAndInputVariables];
	for( int irow = 0; irow < numOfReferenceVariables * numOfOutputAndInputVariables; ++irow ){
		IQPInvTPInvMatrix[irow] = new std::complex<double>[numOfReferenceVariables * numOfReferenceVariables];
		for( int icol = 0; icol < numOfReferenceVariables * numOfReferenceVariables; ++icol ){
			IQPInvTPInvMatrix[irow][icol] = czero;
		}
	}
	for( int irow = 0; irow < numOfReferenceVariables * numOfOutputAndInputVariables; ++irow ){
		for( int icol = 0; icol < numOfReferenceVariables * numOfReferenceVariables; ++icol ){
			std::complex<double> value = czero;
			for( int i = 0; i < numOfReferenceVariables * numOfReferenceVariables; ++i ){
				value += IQMatrix[irow][i] * PInvTPInvMatrix[i][icol];
			}
			IQPInvTPInvMatrix[irow][icol] = value;
		}
	}
#ifdef _DEBUG_WRITE
#else
	for( int irow = 0; irow < numOfReferenceVariables * numOfReferenceVariables; ++irow ){
		delete [] PInvTPInvMatrix[irow];
	}
	delete [] PInvTPInvMatrix;
#endif
	for( int irow = 0; irow < numOfReferenceVariables * numOfOutputAndInputVariables; ++irow ){
		delete [] IQMatrix[irow];
	}
	delete [] IQMatrix;

//#ifdef _DEBUG_WRITE
//	std::cout << "[";
//	for( int row = 0; row < numOfReferenceVariables * numOfOutputAndInputVariables; ++row ){
//		for( int col = 0; col < numOfReferenceVariables * numOfReferenceVariables; ++col ){
//			std::cout << IQPInvTPInvMatrix[row][col].real() << "+" << IQPInvTPInvMatrix[row][col].imag() <<"im ";
//		}
//		if( row+1 < numOfReferenceVariables * numOfOutputAndInputVariables ){
//			std::cout << ";";
//		}
//	}
//	std::cout << "]" << std::endl;
//#endif

	std::complex<double>** matrixVeceh1 = new std::complex<double>*[numOfReferenceVariables * numOfOutputAndInputVariables];
	for( int irow = 0; irow < numOfReferenceVariables * numOfOutputAndInputVariables; ++irow ){
		matrixVeceh1[irow] = new std::complex<double>[2 * numOfReferenceVariables * numOfOutputAndInputVariables];
		for( int icol = 0; icol < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++icol ){
			matrixVeceh1[irow][icol] = czero;// Zero clear
		}
	}
	std::complex<double>** matrixVechh1 = new std::complex<double>*[numOfReferenceVariables * numOfReferenceVariables];
	for( int irow = 0; irow < numOfReferenceVariables * numOfReferenceVariables; ++irow ){
		matrixVechh1[irow] = new std::complex<double>[2 * numOfReferenceVariables * numOfOutputAndInputVariables];
		for( int icol = 0; icol < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++icol ){
			matrixVechh1[irow][icol] = czero;// Zero clear
		}
	}
	std::complex<double>** matrixVeceh2 = new std::complex<double>*[numOfReferenceVariables * numOfOutputAndInputVariables];
	for( int irow = 0; irow < numOfReferenceVariables * numOfOutputAndInputVariables; ++irow ){
		matrixVeceh2[irow] = new std::complex<double>[numOfOutputAndInputVariables];
		for( int icol = 0; icol < numOfOutputAndInputVariables; ++icol ){
			matrixVeceh2[irow][icol] = czero;// Zero clear
		}
	}
	std::complex<double>** matrixVechh2 = new std::complex<double>*[numOfReferenceVariables * numOfReferenceVariables];
	for( int irow = 0; irow < numOfReferenceVariables * numOfReferenceVariables; ++irow ){
		matrixVechh2[irow] = new std::complex<double>[numOfOutputAndInputVariables];
		for( int icol = 0; icol < numOfOutputAndInputVariables; ++icol ){
			matrixVechh2[irow][icol] = czero;// Zero clear
		}
	}
	std::complex<double>* veceh = new std::complex<double>[numOfReferenceVariables * numOfOutputAndInputVariables];
	std::complex<double>* vechh = new std::complex<double>[numOfReferenceVariables * numOfReferenceVariables];
	std::complex<double>** hSigmaMatrix = new std::complex<double>*[numOfReferenceVariables];
	for( int irow = 0; irow < numOfReferenceVariables; ++irow ){
		hSigmaMatrix[irow] = new std::complex<double>[numOfOutputAndInputVariables];
	}
	double* vecSigma = new double[2 * numOfReferenceVariables * numOfOutputAndInputVariables];
	std::complex<double>* sumVeceh3 = new std::complex<double>[numOfReferenceVariables * numOfOutputAndInputVariables];
	for( int irow = 0; irow < numOfReferenceVariables * numOfOutputAndInputVariables; ++irow ){
		sumVeceh3[irow] = czero;// Zero clear
	}
	std::complex<double>* sumVechh3 = new std::complex<double>[numOfReferenceVariables * numOfReferenceVariables];
	for( int irow = 0; irow < numOfReferenceVariables * numOfReferenceVariables; ++irow ){
		sumVechh3[irow] = czero;// Zero clear
	}
#ifdef _DEBUG_WRITE
//	std::cout << "[";
#endif
	for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
		const double val = MD[iSeg] / scale;
		const double diff = RobustWeightTukeysBiweights::calculateSecondDerivativeOfLossFunction(val, paramC) - RobustWeightTukeysBiweights::calculateWeights(val, paramC);
		double factor1(0.0);
		if( fabs(MD[iSeg]) < CommonParameters::EPS ){
			const double sc = scale * paramC;
			factor1 = 4.0 * ( pow(MD[iSeg],2) / pow(sc,4) - 1.0 / pow(sc,2) );
		}else{
			factor1 = diff / pow(MD[iSeg], 2);
		}
		const double factor2 = diff / scale;
		for( int irr = 0; irr < numOfReferenceVariables; ++irr ){
			const int offset = numOfOutputAndInputVariables * irr;
			const int rr = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, irr );
			int iVar(0);
			for( int iOut = 0; iOut < numOfOutputVariables; ++iOut ){
				const int index = ptrControl->getChannelIndex( CommonParameters::OUTPUT, iOut );
				veceh[iVar + offset] = ftval[index][iSeg] * std::conj(ftval[rr][iSeg]) * factor1;
				sumVeceh3[iVar + offset] += ftval[index][iSeg] * std::conj(ftval[rr][iSeg]) * factor2;
				assert(iVar == index);
				++iVar;
			}
			for( int iInp = 0; iInp < numOfInputVariables; ++iInp ){
				const int index = ptrControl->getChannelIndex( CommonParameters::INPUT, iInp );
				veceh[iVar + offset] = ftval[index][iSeg] * std::conj(ftval[rr][iSeg]) * factor1;
				sumVeceh3[iVar + offset] += ftval[index][iSeg] * std::conj(ftval[rr][iSeg]) * factor2;
				assert(iVar == index);
				++iVar;
			}
		}
		for( int irr = 0; irr < numOfReferenceVariables; ++irr ){
			const int offset = numOfReferenceVariables * irr;
			const int rr = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, irr );
			for( int irr2 = 0; irr2 < numOfReferenceVariables; ++irr2 ){
				const int index = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, irr2 );
				vechh[irr2 + offset] = ftval[index][iSeg] * std::conj(ftval[rr][iSeg]) * factor1;
				sumVechh3[irr2 + offset] += ftval[index][iSeg] * std::conj(ftval[rr][iSeg]) * factor2;
			}
		}
#ifdef _DEBUG_WRITE
		//std::cout << "[";
		//for( int row = 0; row < numOfReferenceVariables * numOfOutputAndInputVariables; ++row ){
		//	std::cout << veceh[row].real() << "+" << veceh[row].imag() <<"im ";
		//	if( row+1 < numOfReferenceVariables * numOfOutputAndInputVariables ){
		//		std::cout << ",";
		//	}
		//}
		//std::cout << "]" << std::endl;
		//std::cout << "[";
		//for( int row = 0; row < numOfReferenceVariables * numOfReferenceVariables; ++row ){
		//	std::cout << vechh[row].real() << "+" << vechh[row].imag() <<"im ";
		//	if( row+1 < numOfReferenceVariables * numOfReferenceVariables ){
		//		std::cout << ",";
		//	}
		//}
		//std::cout << "]" << std::endl;
#endif
		calculateVectorForPartialDerivatives( iSeg, ftval, variancesWithoutScale, complexResiduals, hSigmaMatrix, vecSigma );
//#ifdef _DEBUG_WRITE
//		for( int col = 0; col < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++col ){
//			const int iq = col / ( 2 * numOfReferenceVariables );
//			const int isign = ( col / numOfReferenceVariables ) % 2;
//			const int index = 2 * ( col % numOfReferenceVariables ) + isign + 2 * numOfReferenceVariables * iq;
//			const double deriv = - 1.0 / MD[iSeg] *  vecSigma[index];
//			std::cout << deriv <<" ";
//		}
//		if( iSeg+1 < numSegments ){
//			std::cout << ";";
//		}
//#endif
		for( int irow = 0; irow < numOfReferenceVariables * numOfOutputAndInputVariables; ++irow ){
			for( int icol = 0; icol < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++icol ){
				matrixVeceh1[irow][icol] += veceh[irow] * vecSigma[icol];
			}
		}
		for( int irow = 0; irow < numOfReferenceVariables * numOfReferenceVariables; ++irow ){
			for( int icol = 0; icol < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++icol ){
				matrixVechh1[irow][icol] += vechh[irow] * vecSigma[icol];
			}
		}
//#ifdef _DEBUG_WRITE
//		std::cout << "[";
//		for( int row = 0; row < numOfReferenceVariables * numOfOutputAndInputVariables; ++row ){
//			for( int col = 0; col < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++col ){
//				std::cout << matrixVeceh1[row][col].real() << "+" << matrixVeceh1[row][col].imag() <<"im ";
//			}
//			if( row+1 < numOfReferenceVariables * numOfOutputAndInputVariables ){
//				std::cout << ";";
//			}
//		}
//		std::cout << "]" << std::endl;
//		std::cout << "[";
//		for( int row = 0; row < numOfReferenceVariables * numOfReferenceVariables; ++row ){
//			for( int col = 0; col < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++col ){
//				std::cout << matrixVechh1[row][col].real() << "+" << matrixVechh1[row][col].imag() <<"im ";
//			}
//			if( row+1 < numOfReferenceVariables * numOfReferenceVariables ){
//				std::cout << ";";
//			}
//		}
//		std::cout << "]" << std::endl;
//#endif
		for( int irow = 0; irow < numOfReferenceVariables * numOfOutputAndInputVariables; ++irow ){
			for( int icol = 0; icol < numOfOutputAndInputVariables; ++icol ){
				matrixVeceh2[irow][icol] += veceh[irow] * 0.5 * std::norm(complexResiduals[icol][iSeg]) / pow(variancesWithoutScale[icol], 2);
			}
		}
		for( int irow = 0; irow < numOfReferenceVariables * numOfReferenceVariables; ++irow ){
			for( int icol = 0; icol < numOfOutputAndInputVariables; ++icol ){
				matrixVechh2[irow][icol] += vechh[irow] * 0.5 * std::norm(complexResiduals[icol][iSeg]) / pow(variancesWithoutScale[icol], 2);
			}
		}
//#ifdef _DEBUG_WRITE
//		std::cout << "[";
//		for( int row = 0; row < numOfOutputAndInputVariables; ++row ){
//			std::cout << complexResiduals[row][iSeg].real() << "+" << complexResiduals[row][iSeg].imag() <<"im ";
//			if( row+1 < numOfOutputAndInputVariables ){
//				std::cout << ",";
//			}
//		}
//		std::cout << "]" << std::endl;
//		std::cout << "[";
//		for( int row = 0; row < numOfReferenceVariables * numOfOutputAndInputVariables; ++row ){
//			for( int col = 0; col < numOfOutputAndInputVariables; ++col ){
//				std::cout << matrixVeceh2[row][col].real() << "+" << matrixVeceh2[row][col].imag() <<"im ";
//			}
//			if( row+1 < numOfReferenceVariables * numOfOutputAndInputVariables ){
//				std::cout << ";";
//			}
//		}
//		std::cout << "]" << std::endl;
//		std::cout << "[";
//		for( int row = 0; row < numOfReferenceVariables * numOfReferenceVariables; ++row ){
//			for( int col = 0; col < numOfOutputAndInputVariables; ++col ){
//				std::cout << matrixVechh2[row][col].real() << "+" << matrixVechh2[row][col].imag() <<"im ";
//			}
//			if( row+1 < numOfReferenceVariables * numOfReferenceVariables ){
//				std::cout << ";";
//			}
//		}
//		std::cout << "]" << std::endl;
//#endif
	}
//#ifdef _DEBUG_WRITE
//	std::cout << "]" << std::endl;
//	std::cout << "[";
//	for( int row = 0; row < numOfReferenceVariables * numOfOutputAndInputVariables; ++row ){
//		std::cout << sumVeceh3[row].real() << "+" << sumVeceh3[row].imag() <<"im ";
//		if( row+1 < numOfReferenceVariables * numOfOutputAndInputVariables ){
//			std::cout << ",";
//		}
//	}
//	std::cout << "]" << std::endl;
//	std::cout << "[";
//	for( int row = 0; row < numOfReferenceVariables * numOfReferenceVariables; ++row ){
//		std::cout << sumVechh3[row].real() << "+" << sumVechh3[row].imag() <<"im ";
//		if( row+1 < numOfReferenceVariables * numOfReferenceVariables ){
//			std::cout << ",";
//		}
//	}
//	std::cout << "]" << std::endl;
//	std::cout << "[";
//	for( int row = 0; row < numOfReferenceVariables * numOfOutputAndInputVariables; ++row ){
//		for( int col = 0; col < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++col ){
//			const int iq = col / ( 2 * numOfReferenceVariables );
//			const int isign = ( col / numOfReferenceVariables ) % 2;
//			const int index = 2 * ( col % numOfReferenceVariables ) + isign + 2 * numOfReferenceVariables * iq;
//			std::cout << -matrixVeceh1[row][index].real() << "+" << -matrixVeceh1[row][index].imag() <<"im ";
//		}
//		if( row+1 < numOfReferenceVariables * numOfOutputAndInputVariables ){
//			std::cout << ";";
//		}
//	}
//	std::cout << "]" << std::endl;
//	std::cout << "[";
//	for( int row = 0; row < numOfReferenceVariables * numOfReferenceVariables; ++row ){
//		for( int col = 0; col < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++col ){
//			const int iq = col / ( 2 * numOfReferenceVariables );
//			const int isign = ( col / numOfReferenceVariables ) % 2;
//			const int index = 2 * ( col % numOfReferenceVariables ) + isign + 2 * numOfReferenceVariables * iq;
//			std::cout << -matrixVechh1[row][index].real() << "+" << -matrixVechh1[row][index].imag() <<"im ";
//		}
//		if( row+1 < numOfReferenceVariables * numOfReferenceVariables ){
//			std::cout << ";";
//		}
//	}
//	std::cout << "]" << std::endl;
//	std::cout << "[";
//	for( int irow = 0; irow < numOfReferenceVariables * numOfReferenceVariables; ++irow ){
//		for( int icol = 0; icol < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++icol ){
//			std::complex<double> value = czero;
//			const int iq = icol / ( 2 * numOfReferenceVariables );
//			const int isign = ( icol / numOfReferenceVariables ) % 2;
//			const int index = 2 * ( icol % numOfReferenceVariables ) + isign + 2 * numOfReferenceVariables * iq;
//			for( int i = 0; i < numOfReferenceVariables * numOfReferenceVariables; ++i ){
//				value += PInvTPInvMatrix[irow][i] * matrixVechh1[i][index];
//			}
//			std::cout << value.real() << "+" << value.imag() <<"im ";
//		}
//		if( irow+1 < numOfReferenceVariables * numOfReferenceVariables ){
//			std::cout << ";";
//		}
//	}
//	std::cout << "]" << std::endl;
//	for( int irow = 0; irow < numOfReferenceVariables * numOfReferenceVariables; ++irow ){
//		delete [] PInvTPInvMatrix[irow];
//	}
//	delete [] PInvTPInvMatrix;
//#endif
	delete [] veceh;
	delete [] vechh;
	for( int irow = 0; irow < numOfReferenceVariables; ++irow ){
		delete [] hSigmaMatrix[irow];
	}
	delete [] hSigmaMatrix;
	delete [] vecSigma;

	// Partial derivatives regarding responses
	std::complex<double>** complexDerivatives = new std::complex<double>*[numOfReferenceVariables * numOfOutputAndInputVariables];
	for( int irow = 0; irow < numOfReferenceVariables * numOfOutputAndInputVariables; ++irow ){
		complexDerivatives[irow] = new std::complex<double>[2 * numOfReferenceVariables * numOfOutputAndInputVariables];
	}
	for( int irow = 0; irow < numOfReferenceVariables * numOfOutputAndInputVariables; ++irow ){
		for( int icol = 0; icol < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++icol ){
			std::complex<double> value = czero;
			// The 1st term
			for( int i = 0; i < numOfReferenceVariables * numOfOutputAndInputVariables; ++i ){
				value -= PInvTIMatrix[irow][i] * matrixVeceh1[i][icol];
			}
			// The 2nd term
			for( int i = 0; i < numOfReferenceVariables * numOfReferenceVariables; ++i ){
				value += IQPInvTPInvMatrix[irow][i] * matrixVechh1[i][icol];
			}
			complexDerivatives[irow][icol] = value;
		}
	}
	for( int irow = 0; irow < numOfReferenceVariables * numOfOutputAndInputVariables; ++irow ){
		delete [] matrixVeceh1[irow];
	}
	delete [] matrixVeceh1;
	for( int irow = 0; irow < numOfReferenceVariables * numOfReferenceVariables; ++irow ){
		delete [] matrixVechh1[irow];
	}
	delete [] matrixVechh1;
//#ifdef _DEBUG_WRITE
//	std::cout << "[";
//	for( int row = 0; row < numOfReferenceVariables * numOfOutputAndInputVariables; ++row ){
//		for( int col = 0; col < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++col ){
//			std::cout << complexDerivatives[row][col].real() << "+" << complexDerivatives[row][col].imag() <<"im ";
//		}
//		if( row+1 < numOfReferenceVariables * numOfOutputAndInputVariables ){
//			std::cout << ";";
//		}
//	}
//	std::cout << "]" << std::endl;
//#endif

	// Order: Re(Z11),Im(Z11),Re(Z12),Im(Z12),Re(Z21),Im(Z21),Re(Z22),Im(Z22),...,Re(Zq1),Im(Zq1),Re(Zq2),Im(Zq2)
	for( int irow = 0; irow < numOfReferenceVariables * numOfOutputAndInputVariables; ++irow ){
		const int ir = irow / numOfOutputAndInputVariables;
		const int iq = irow % numOfOutputAndInputVariables;
		const int index = 2 * ir + iq * 2 * numOfReferenceVariables;
		for( int icol = 0; icol < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++icol ){
			const int iq2 = icol / ( 2 * numOfReferenceVariables );
			const int offset = iq2 * 2 * numOfReferenceVariables;
			const int amari = icol % ( 2 * numOfReferenceVariables );
			if( amari < numOfReferenceVariables ){
				// Real part
				const int ir2 = amari;
				const int index2 = 2 * ir2 + offset;
				derivativesRegardingResps[index  ][index2] = complexDerivatives[irow][icol].real();
				derivativesRegardingResps[index+1][index2] = complexDerivatives[irow][icol].imag();
			}else{
				// Imaginary part
				const int ir2 = amari - numOfReferenceVariables;
				const int index2 = 2 * ir2 + 1 + offset;
				derivativesRegardingResps[index  ][index2] = complexDerivatives[irow][icol].real();
				derivativesRegardingResps[index+1][index2] = complexDerivatives[irow][icol].imag();
			}
		}
	}
#ifdef _DEBUG_WRITE
	std::cout << "drr1" << std::endl;
	std::cout << "[";
	for( int row = 0; row < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++row ){
		for( int col = 0; col < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++col ){
			std::cout << derivativesRegardingResps[row][col] <<" ";
		}
		if( row+1 < 2 * numOfReferenceVariables * numOfOutputAndInputVariables ){
			std::cout << ";";
		}
	}
	std::cout << "]" << std::endl;
#endif
#ifdef _DEBUG_20230128
	std::cout << "derivatives" << std::endl;
	for( int row = 0; row < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++row ){
		for( int col = 0; col < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++col ){
			std::cout << derivativesRegardingResps[row][col] << ",";
		}
		std::cout << std::endl;
	}
#endif
#ifdef _DEBUG_WRITE
	double** derivMD1 = new double*[numSegments];
	for( int irow = 0; irow < numSegments; ++irow ){
		derivMD1[irow] = new double[2 * numOfReferenceVariables * numOfOutputAndInputVariables];
	}
	std::complex<double>** derivP1 = new std::complex<double>*[numOfReferenceVariables * numOfReferenceVariables];
	for( int irow = 0; irow < numOfReferenceVariables * numOfReferenceVariables; ++irow ){
		derivP1[irow] = new std::complex<double>[2 * numOfReferenceVariables * numOfOutputAndInputVariables];
	}
	std::complex<double>** derivInvP1 = new std::complex<double>*[numOfReferenceVariables * numOfReferenceVariables];
	for( int irow = 0; irow < numOfReferenceVariables * numOfReferenceVariables; ++irow ){
		derivInvP1[irow] = new std::complex<double>[2 * numOfReferenceVariables * numOfOutputAndInputVariables];
	}
	std::complex<double>** derivQ1 = new std::complex<double>*[numOfReferenceVariables * numOfOutputAndInputVariables];
	for( int irow = 0; irow < numOfReferenceVariables * numOfOutputAndInputVariables; ++irow ){
		derivQ1[irow] = new std::complex<double>[2 * numOfReferenceVariables * numOfOutputAndInputVariables];
	}
	double** derivResp1 = new double*[2 * numOfReferenceVariables * numOfOutputAndInputVariables];
	for( int irow = 0; irow < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++irow ){
		derivResp1[irow] = new double[2 * numOfReferenceVariables * numOfOutputAndInputVariables];
	}
	for( int irow = 0; irow < numOfOutputAndInputVariables; ++irow ){
		for( int icol = 0; icol < numOfReferenceVariables; ++icol ){
			for( int isign = 0; isign < 2; ++isign ){
				std::complex<double>** complexResidualsOrg = new std::complex<double>*[numOfOutputAndInputVariables];
				for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
					complexResidualsOrg[iVar] = new std::complex<double>[numSegments]; 
				}
				double* MDOrg = new double[numSegments];
				double* weightsOrg = new double[numSegments];
				std::complex<double>** respOrg = new std::complex<double>*[numOfOutputAndInputVariables];
				for( int i = 0; i < numOfOutputAndInputVariables; ++i ){
					respOrg[i] = new std::complex<double>[numOfReferenceVariables];
				}
				std::complex<double>** complexResidualsMod = new std::complex<double>*[numOfOutputAndInputVariables];
				for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
					complexResidualsMod[iVar] = new std::complex<double>[numSegments]; 
				}
				double* MDMod = new double[numSegments];
				double* weightsMod = new double[numSegments];
				std::complex<double>** respMod = new std::complex<double>*[numOfOutputAndInputVariables];
				for( int i = 0; i < numOfOutputAndInputVariables; ++i ){
					respMod[i] = new std::complex<double>[numOfReferenceVariables];
				}
				for( int irow2 = 0; irow2 < numOfOutputAndInputVariables; ++irow2 ){
					for( int icol2 = 0; icol2 < numOfReferenceVariables; ++icol2 ){
						respOrg[irow2][icol2] = resp[irow2][icol2];
						respMod[irow2][icol2] = resp[irow2][icol2];
					}
				}
				double dresp(0.0);
				if( isign == 0 ){
					dresp = resp[irow][icol].real() * 0.001;
					respMod[irow][icol] += std::complex<double>( dresp, 0.0 );
				}else{
					dresp = resp[irow][icol].imag() * 0.001;
					respMod[irow][icol] += std::complex<double>( 0.0, dresp );
				}
				// Calculate complex residuals
				calculateComplexResiduals( numSegments, ftval, respOrg, complexResidualsOrg );
				calculateComplexResiduals( numSegments, ftval, respMod, complexResidualsMod );
				// Calculate Mahalanobis distance
				calculateMD( numSegments, numOfOutputAndInputVariables, complexResidualsOrg,
					variancesWithoutScale, MDOrg );
				calculateMD( numSegments, numOfOutputAndInputVariables, complexResidualsMod,
					variancesWithoutScale, MDMod );
				const int index = 2 * icol + isign + irow * 2 * numOfReferenceVariables;
				for( int irow2 = 0; irow2 < numSegments; ++irow2 ){
					const double deriv = ( MDMod[irow2] - MDOrg[irow2] ) / dresp;
					derivMD1[irow2][index] = deriv;
				}

				// Calculate original weights
				for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
					const double val = MDOrg[iSeg] / scale;
					weightsOrg[iSeg] = RobustWeightTukeysBiweights::calculateWeights(val, paramC);
				}
				for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
					const double val = MDMod[iSeg] / scale;
					weightsMod[iSeg] = RobustWeightTukeysBiweights::calculateWeights(val, paramC);
				}
				std::complex<double> PMatrixOrg[2][2] = { czero, czero, czero, czero };
				for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
					PMatrixOrg[0][0] += ftval[rr0][iSeg] * std::conj(ftval[rr0][iSeg]) * weightsOrg[iSeg];
					PMatrixOrg[1][0] += ftval[rr1][iSeg] * std::conj(ftval[rr0][iSeg]) * weightsOrg[iSeg];
					PMatrixOrg[0][1] += ftval[rr0][iSeg] * std::conj(ftval[rr1][iSeg]) * weightsOrg[iSeg];
					PMatrixOrg[1][1] += ftval[rr1][iSeg] * std::conj(ftval[rr1][iSeg]) * weightsOrg[iSeg];
				}
				std::complex<double> PMatrixMod[2][2] = { czero, czero, czero, czero };
				for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
					PMatrixMod[0][0] += ftval[rr0][iSeg] * std::conj(ftval[rr0][iSeg]) * weightsMod[iSeg];
					PMatrixMod[1][0] += ftval[rr1][iSeg] * std::conj(ftval[rr0][iSeg]) * weightsMod[iSeg];
					PMatrixMod[0][1] += ftval[rr0][iSeg] * std::conj(ftval[rr1][iSeg]) * weightsMod[iSeg];
					PMatrixMod[1][1] += ftval[rr1][iSeg] * std::conj(ftval[rr1][iSeg]) * weightsMod[iSeg];
				}
				derivP1[0][index] = ( PMatrixMod[0][0] - PMatrixOrg[0][0] ) / dresp;
				derivP1[1][index] = ( PMatrixMod[1][0] - PMatrixOrg[1][0] ) / dresp;
				derivP1[2][index] = ( PMatrixMod[0][1] - PMatrixOrg[0][1] ) / dresp;
				derivP1[3][index] = ( PMatrixMod[1][1] - PMatrixOrg[1][1] ) / dresp;

				const std::complex<double> detOrg = PMatrixOrg[0][0] * PMatrixOrg[1][1] - PMatrixOrg[1][0] * PMatrixOrg[0][1];
				const std::complex<double> PInvMatrixOrg[2][2] = { PMatrixOrg[1][1]/detOrg, -PMatrixOrg[0][1]/detOrg, -PMatrixOrg[1][0]/detOrg, PMatrixOrg[0][0]/detOrg };
				const std::complex<double> detMod = PMatrixMod[0][0] * PMatrixMod[1][1] - PMatrixMod[1][0] * PMatrixMod[0][1];
				const std::complex<double> PInvMatrixMod[2][2] = { PMatrixMod[1][1]/detMod, -PMatrixMod[0][1]/detMod, -PMatrixMod[1][0]/detMod, PMatrixMod[0][0]/detMod };
				derivInvP1[0][index] = ( PInvMatrixMod[0][0] - PInvMatrixOrg[0][0] ) / dresp;
				derivInvP1[1][index] = ( PInvMatrixMod[1][0] - PInvMatrixOrg[1][0] ) / dresp;
				derivInvP1[2][index] = ( PInvMatrixMod[0][1] - PInvMatrixOrg[0][1] ) / dresp;
				derivInvP1[3][index] = ( PInvMatrixMod[1][1] - PInvMatrixOrg[1][1] ) / dresp;

				std::complex<double>** QMatrixOrg = new std::complex<double>*[numOfOutputAndInputVariables];
				std::complex<double>** QMatrixMod = new std::complex<double>*[numOfOutputAndInputVariables];
				for( int irow = 0; irow < numOfOutputAndInputVariables; ++irow ){
					QMatrixOrg[irow] = new std::complex<double>[numOfReferenceVariables];
					QMatrixMod[irow] = new std::complex<double>[numOfReferenceVariables];
					for( int icol = 0; icol < numOfReferenceVariables; ++icol ){
						QMatrixOrg[irow][icol] = czero;
						QMatrixMod[irow][icol] = czero;
					}
				}
				for( int irr = 0; irr < numOfReferenceVariables; ++irr ){
					const int rr = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, irr );
					int iVar(0);
					for( int iOut = 0; iOut < numOfOutputVariables; ++iOut ){
						const int index = ptrControl->getChannelIndex( CommonParameters::OUTPUT, iOut );
						for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
							QMatrixOrg[iVar][irr] += ftval[index][iSeg] * std::conj(ftval[rr][iSeg]) * weightsOrg[iSeg];
							QMatrixMod[iVar][irr] += ftval[index][iSeg] * std::conj(ftval[rr][iSeg]) * weightsMod[iSeg];
						}
						assert(index == iVar);
						++iVar;
					}
					for( int iInp = 0; iInp < numOfInputVariables; ++iInp ){
						const int index = ptrControl->getChannelIndex( CommonParameters::INPUT, iInp );
						for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
							QMatrixOrg[iVar][irr] += ftval[index][iSeg] * std::conj(ftval[rr][iSeg]) * weightsOrg[iSeg];
							QMatrixMod[iVar][irr] += ftval[index][iSeg] * std::conj(ftval[rr][iSeg]) * weightsMod[iSeg];
						}
						assert(index == iVar);
						++iVar;
					}
				}
				for( int irr = 0; irr < numOfReferenceVariables; ++irr ){
					for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
						const int irow = iVar + numOfOutputAndInputVariables * irr;
						derivQ1[irow][index] = ( QMatrixMod[iVar][irr] - QMatrixOrg[iVar][irr] ) / dresp;
					}
				}
				for( int irow2 = 0; irow2 < numOfOutputAndInputVariables; ++irow2 ){
					for( int icol2 = 0; icol2 < numOfReferenceVariables; ++icol2 ){
						std::complex<double> valueOrg = czero;
						std::complex<double> valueMod = czero;
						for( int i = 0; i < numOfReferenceVariables; ++i ){
							valueOrg += QMatrixOrg[irow2][i] * PInvMatrixOrg[i][icol2];
							valueMod += QMatrixMod[irow2][i] * PInvMatrixMod[i][icol2];
						}
						const int index2 = 2 * icol2 + irow2 * 2 * numOfReferenceVariables;
						derivResp1[index2  ][index] = ( valueMod.real() - valueOrg.real() ) / dresp;
						derivResp1[index2+1][index] = ( valueMod.imag() - valueOrg.imag() )/ dresp;
					}
				}
				for( int irow2 = 0; irow2 < numOfOutputAndInputVariables; ++irow2 ){
					delete [] QMatrixOrg[irow2];
				}
				delete [] QMatrixOrg;
				for( int irow2 = 0; irow2 < numOfOutputAndInputVariables; ++irow2 ){
					delete [] QMatrixMod[irow2];
				}
				delete [] QMatrixMod;
				for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
					delete [] complexResidualsOrg[iVar]; 
				}
				delete [] complexResidualsOrg;
				delete [] MDOrg;
				delete [] weightsOrg;
				for( int i = 0; i < numOfOutputAndInputVariables; ++i ){
					delete [] respOrg[i];
				}
				delete [] respOrg;
				for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
					delete [] complexResidualsMod[iVar]; 
				}
				delete [] complexResidualsMod;
				delete [] MDMod;
				delete [] weightsMod;
				for( int i = 0; i < numOfOutputAndInputVariables; ++i ){
					delete [] respMod[i];
				}
				delete [] respMod;
			}
		}
	}
	//std::cout << "[";
	//for( int row = 0; row < numSegments; ++row ){
	//	for( int col = 0; col < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++col ){
	//		std::cout << derivMD[row][col] <<" ";
	//	}
	//	if( row+1 < numSegments ){
	//		std::cout << ";";
	//	}
	//}
	//std::cout << "]" << std::endl;
	//std::cout << "[";
	//for( int row = 0; row < numOfReferenceVariables * numOfReferenceVariables; ++row ){
	//	for( int col = 0; col < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++col ){
	//		std::cout << derivP[row][col].real() << "+" << derivP[row][col].imag() <<"im ";
	//	}
	//	if( row+1 < numOfReferenceVariables * numOfReferenceVariables ){
	//		std::cout << ";";
	//	}
	//}
	//std::cout << "]" << std::endl;
	//std::cout << "[";
	//for( int row = 0; row < numOfReferenceVariables * numOfReferenceVariables; ++row ){
	//	for( int col = 0; col < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++col ){
	//		std::cout << derivInvP[row][col].real() << "+" << derivInvP[row][col].imag() <<"im ";
	//	}
	//	if( row+1 < numOfReferenceVariables * numOfReferenceVariables ){
	//		std::cout << ";";
	//	}
	//}
	//std::cout << "]" << std::endl;
	//std::cout << "[";
	//for( int row = 0; row < numOfReferenceVariables * numOfOutputAndInputVariables; ++row ){
	//	for( int col = 0; col < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++col ){
	//		std::cout << derivQ[row][col].real() << "+" << derivQ[row][col].imag() <<"im ";
	//	}
	//	if( row+1 < numOfReferenceVariables * numOfOutputAndInputVariables ){
	//		std::cout << ";";
	//	}
	//}
	//std::cout << "]" << std::endl;
	std::cout << "drr2" << std::endl;
	std::cout << "[";
	for( int row = 0; row < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++row ){
		for( int col = 0; col < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++col ){
			std::cout << derivResp1[row][col] <<" ";
		}
		if( row+1 < 2 * numOfReferenceVariables * numOfOutputAndInputVariables ){
			std::cout << ";";
		}
	}
	std::cout << "]" << std::endl;
	for( int irow = 0; irow < numSegments; ++irow ){
		delete [] derivMD1[irow];
	}
	delete [] derivMD1;
	for( int irow = 0; irow < numOfReferenceVariables*numOfReferenceVariables; ++irow ){
		delete [] derivP1[irow];
	}
	delete [] derivP1;
	for( int irow = 0; irow < numOfReferenceVariables*numOfReferenceVariables; ++irow ){
		delete [] derivInvP1[irow];
	}
	delete [] derivInvP1;
	for( int irow = 0; irow < numOfReferenceVariables * numOfOutputAndInputVariables; ++irow ){
		delete [] derivQ1[irow];
	}
	delete [] derivQ1;
	for( int irow = 0; irow < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++irow ){
		delete [] derivResp1[irow];
	}
	delete [] derivResp1;
#endif
#ifdef _DEBUG_20230128
	////double** derivMD1 = new double*[numSegments];
	////for( int irow = 0; irow < numSegments; ++irow ){
	////	derivMD1[irow] = new double[2 * numOfReferenceVariables * numOfOutputAndInputVariables];
	////}
	////std::complex<double>** derivP1 = new std::complex<double>*[numOfReferenceVariables * numOfReferenceVariables];
	////for( int irow = 0; irow < numOfReferenceVariables * numOfReferenceVariables; ++irow ){
	////	derivP1[irow] = new std::complex<double>[2 * numOfReferenceVariables * numOfOutputAndInputVariables];
	////}
	////std::complex<double>** derivInvP1 = new std::complex<double>*[numOfReferenceVariables * numOfReferenceVariables];
	////for( int irow = 0; irow < numOfReferenceVariables * numOfReferenceVariables; ++irow ){
	////	derivInvP1[irow] = new std::complex<double>[2 * numOfReferenceVariables * numOfOutputAndInputVariables];
	////}
	////std::complex<double>** derivQ1 = new std::complex<double>*[numOfReferenceVariables * numOfOutputAndInputVariables];
	////for( int irow = 0; irow < numOfReferenceVariables * numOfOutputAndInputVariables; ++irow ){
	////	derivQ1[irow] = new std::complex<double>[2 * numOfReferenceVariables * numOfOutputAndInputVariables];
	////}
	////double** derivResp1 = new double*[2 * numOfReferenceVariables * numOfOutputAndInputVariables];
	////for( int irow = 0; irow < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++irow ){
	////	derivResp1[irow] = new double[2 * numOfReferenceVariables * numOfOutputAndInputVariables];
	////}
	//const double dresps[11] = { -0.010, -0.008, -0.006, -0.004, -0.002, 0.0, 0.002, 0.004, 0.006, 0.008, 0.010 };
	//for( int idresp = 0; idresp < 11; ++idresp ){
	//	std::complex<double>** complexResidualsOrg = new std::complex<double>*[numOfOutputAndInputVariables];
	//	for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
	//		complexResidualsOrg[iVar] = new std::complex<double>[numSegments]; 
	//	}
	//	double* MDOrg = new double[numSegments];
	//	double* weightsOrg = new double[numSegments];
	//	std::complex<double>** respOrg = new std::complex<double>*[numOfOutputAndInputVariables];
	//	for( int i = 0; i < numOfOutputAndInputVariables; ++i ){
	//		respOrg[i] = new std::complex<double>[numOfReferenceVariables];
	//	}
	//	std::complex<double>** complexResidualsMod = new std::complex<double>*[numOfOutputAndInputVariables];
	//	for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
	//		complexResidualsMod[iVar] = new std::complex<double>[numSegments]; 
	//	}
	//	double* MDMod = new double[numSegments];
	//	double* weightsMod = new double[numSegments];
	//	std::complex<double>** respMod = new std::complex<double>*[numOfOutputAndInputVariables];
	//	for( int i = 0; i < numOfOutputAndInputVariables; ++i ){
	//		respMod[i] = new std::complex<double>[numOfReferenceVariables];
	//	}
	//	for( int irow2 = 0; irow2 < numOfOutputAndInputVariables; ++irow2 ){
	//		for( int icol2 = 0; icol2 < numOfReferenceVariables; ++icol2 ){
	//			respOrg[irow2][icol2] = resp[irow2][icol2];
	//			respMod[irow2][icol2] = resp[irow2][icol2];
	//		}
	//	}
	//	const double dresp = dresps[idresp];
	//	respMod[3][0] += std::complex<double>( 0.0, dresp );
	//	// Calculate complex residuals
	//	calculateComplexResiduals( numSegments, ftval, respOrg, complexResidualsOrg );
	//	calculateComplexResiduals( numSegments, ftval, respMod, complexResidualsMod );
	//	// Calculate Mahalanobis distance
	//	calculateMD( numSegments, numOfOutputAndInputVariables, complexResidualsOrg,
	//		variancesWithoutScale, MDOrg );
	//	calculateMD( numSegments, numOfOutputAndInputVariables, complexResidualsMod,
	//		variancesWithoutScale, MDMod );
	//	//const int index = 2 * icol + isign + irow * 2 * numOfReferenceVariables;
	//	//for( int irow2 = 0; irow2 < numSegments; ++irow2 ){
	//	//	const double deriv = ( MDMod[irow2] - MDOrg[irow2] ) / dresp;
	//	//	derivMD1[irow2][index] = deriv;
	//	//}

	//	// Calculate original weights
	//	for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
	//		const double val = MDOrg[iSeg] / scale;
	//		weightsOrg[iSeg] = UtilRobust::calculateTukeysBiweights(val, paramC);
	//	}
	//	for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
	//		const double val = MDMod[iSeg] / scale;
	//		weightsMod[iSeg] = UtilRobust::calculateTukeysBiweights(val, paramC);
	//	}
	//	std::complex<double> PMatrixOrg[2][2] = { czero, czero, czero, czero };
	//	for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
	//		PMatrixOrg[0][0] += ftval[rr0][iSeg] * std::conj(ftval[rr0][iSeg]) * weightsOrg[iSeg];
	//		PMatrixOrg[1][0] += ftval[rr1][iSeg] * std::conj(ftval[rr0][iSeg]) * weightsOrg[iSeg];
	//		PMatrixOrg[0][1] += ftval[rr0][iSeg] * std::conj(ftval[rr1][iSeg]) * weightsOrg[iSeg];
	//		PMatrixOrg[1][1] += ftval[rr1][iSeg] * std::conj(ftval[rr1][iSeg]) * weightsOrg[iSeg];
	//	}
	//	std::complex<double> PMatrixMod[2][2] = { czero, czero, czero, czero };
	//	for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
	//		PMatrixMod[0][0] += ftval[rr0][iSeg] * std::conj(ftval[rr0][iSeg]) * weightsMod[iSeg];
	//		PMatrixMod[1][0] += ftval[rr1][iSeg] * std::conj(ftval[rr0][iSeg]) * weightsMod[iSeg];
	//		PMatrixMod[0][1] += ftval[rr0][iSeg] * std::conj(ftval[rr1][iSeg]) * weightsMod[iSeg];
	//		PMatrixMod[1][1] += ftval[rr1][iSeg] * std::conj(ftval[rr1][iSeg]) * weightsMod[iSeg];
	//	}
	//	//derivP1[0][index] = ( PMatrixMod[0][0] - PMatrixOrg[0][0] ) / dresp;
	//	//derivP1[1][index] = ( PMatrixMod[1][0] - PMatrixOrg[1][0] ) / dresp;
	//	//derivP1[2][index] = ( PMatrixMod[0][1] - PMatrixOrg[0][1] ) / dresp;
	//	//derivP1[3][index] = ( PMatrixMod[1][1] - PMatrixOrg[1][1] ) / dresp;

	//	const std::complex<double> detOrg = PMatrixOrg[0][0] * PMatrixOrg[1][1] - PMatrixOrg[1][0] * PMatrixOrg[0][1];
	//	const std::complex<double> PInvMatrixOrg[2][2] = { PMatrixOrg[1][1]/detOrg, -PMatrixOrg[0][1]/detOrg, -PMatrixOrg[1][0]/detOrg, PMatrixOrg[0][0]/detOrg };
	//	const std::complex<double> detMod = PMatrixMod[0][0] * PMatrixMod[1][1] - PMatrixMod[1][0] * PMatrixMod[0][1];
	//	const std::complex<double> PInvMatrixMod[2][2] = { PMatrixMod[1][1]/detMod, -PMatrixMod[0][1]/detMod, -PMatrixMod[1][0]/detMod, PMatrixMod[0][0]/detMod };
	//	//derivInvP1[0][index] = ( PInvMatrixMod[0][0] - PInvMatrixOrg[0][0] ) / dresp;
	//	//derivInvP1[1][index] = ( PInvMatrixMod[1][0] - PInvMatrixOrg[1][0] ) / dresp;
	//	//derivInvP1[2][index] = ( PInvMatrixMod[0][1] - PInvMatrixOrg[0][1] ) / dresp;
	//	//derivInvP1[3][index] = ( PInvMatrixMod[1][1] - PInvMatrixOrg[1][1] ) / dresp;

	//	std::complex<double>** QMatrixOrg = new std::complex<double>*[numOfOutputAndInputVariables];
	//	std::complex<double>** QMatrixMod = new std::complex<double>*[numOfOutputAndInputVariables];
	//	for( int irow = 0; irow < numOfOutputAndInputVariables; ++irow ){
	//		QMatrixOrg[irow] = new std::complex<double>[numOfReferenceVariables];
	//		QMatrixMod[irow] = new std::complex<double>[numOfReferenceVariables];
	//		for( int icol = 0; icol < numOfReferenceVariables; ++icol ){
	//			QMatrixOrg[irow][icol] = czero;
	//			QMatrixMod[irow][icol] = czero;
	//		}
	//	}
	//	for( int irr = 0; irr < numOfReferenceVariables; ++irr ){
	//		const int rr = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, irr );
	//		int iVar(0);
	//		for( int iOut = 0; iOut < numOfOutputVariables; ++iOut ){
	//			const int index = ptrControl->getChannelIndex( CommonParameters::OUTPUT, iOut );
	//			for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
	//				QMatrixOrg[iVar][irr] += ftval[index][iSeg] * std::conj(ftval[rr][iSeg]) * weightsOrg[iSeg];
	//				QMatrixMod[iVar][irr] += ftval[index][iSeg] * std::conj(ftval[rr][iSeg]) * weightsMod[iSeg];
	//			}
	//			assert(index == iVar);
	//			++iVar;
	//		}
	//		for( int iInp = 0; iInp < numOfInputVariables; ++iInp ){
	//			const int index = ptrControl->getChannelIndex( CommonParameters::INPUT, iInp );
	//			for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
	//				QMatrixOrg[iVar][irr] += ftval[index][iSeg] * std::conj(ftval[rr][iSeg]) * weightsOrg[iSeg];
	//				QMatrixMod[iVar][irr] += ftval[index][iSeg] * std::conj(ftval[rr][iSeg]) * weightsMod[iSeg];
	//			}
	//			assert(index == iVar);
	//			++iVar;
	//		}
	//	}
	//	//for( int irr = 0; irr < numOfReferenceVariables; ++irr ){
	//	//	for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
	//	//		const int irow = iVar + numOfOutputAndInputVariables * irr;
	//	//		derivQ1[irow][index] = ( QMatrixMod[iVar][irr] - QMatrixOrg[iVar][irr] ) / dresp;
	//	//	}
	//	//}
	//	std::cout << dresp;
	//	for( int irow2 = 0; irow2 < numOfOutputAndInputVariables; ++irow2 ){
	//		for( int icol2 = 0; icol2 < numOfReferenceVariables; ++icol2 ){
	//			std::complex<double> valueOrg = czero;
	//			std::complex<double> valueMod = czero;
	//			for( int i = 0; i < numOfReferenceVariables; ++i ){
	//				valueOrg += QMatrixOrg[irow2][i] * PInvMatrixOrg[i][icol2];
	//				valueMod += QMatrixMod[irow2][i] * PInvMatrixMod[i][icol2];
	//			}
	//			std::cout << ", " << valueMod.real() << ", " << valueMod.imag();
	//		}
	//	}
	//	std::cout << std::endl;
	//	for( int irow2 = 0; irow2 < numOfOutputAndInputVariables; ++irow2 ){
	//		delete [] QMatrixOrg[irow2];
	//	}
	//	delete [] QMatrixOrg;
	//	for( int irow2 = 0; irow2 < numOfOutputAndInputVariables; ++irow2 ){
	//		delete [] QMatrixMod[irow2];
	//	}
	//	delete [] QMatrixMod;
	//	for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
	//		delete [] complexResidualsOrg[iVar]; 
	//	}
	//	delete [] complexResidualsOrg;
	//	delete [] MDOrg;
	//	delete [] weightsOrg;
	//	for( int i = 0; i < numOfOutputAndInputVariables; ++i ){
	//		delete [] respOrg[i];
	//	}
	//	delete [] respOrg;
	//	for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
	//		delete [] complexResidualsMod[iVar]; 
	//	}
	//	delete [] complexResidualsMod;
	//	delete [] MDMod;
	//	delete [] weightsMod;
	//	for( int i = 0; i < numOfOutputAndInputVariables; ++i ){
	//		delete [] respMod[i];
	//	}
	//	delete [] respMod;
	//}
	////for( int irow = 0; irow < numSegments; ++irow ){
	////	delete [] derivMD1[irow];
	////}
	////delete [] derivMD1;
	////for( int irow = 0; irow < numOfReferenceVariables*numOfReferenceVariables; ++irow ){
	////	delete [] derivP1[irow];
	////}
	////delete [] derivP1;
	////for( int irow = 0; irow < numOfReferenceVariables*numOfReferenceVariables; ++irow ){
	////	delete [] derivInvP1[irow];
	////}
	////delete [] derivInvP1;
	////for( int irow = 0; irow < numOfReferenceVariables * numOfOutputAndInputVariables; ++irow ){
	////	delete [] derivQ1[irow];
	////}
	////delete [] derivQ1;
	////for( int irow = 0; irow < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++irow ){
	////	delete [] derivResp1[irow];
	////}
	////delete [] derivResp1;
	//double** derivMD1 = new double*[numSegments];
	//for( int irow = 0; irow < numSegments; ++irow ){
	//	derivMD1[irow] = new double[2 * numOfReferenceVariables * numOfOutputAndInputVariables];
	//}
	//std::complex<double>** derivP1 = new std::complex<double>*[numOfReferenceVariables * numOfReferenceVariables];
	//for( int irow = 0; irow < numOfReferenceVariables * numOfReferenceVariables; ++irow ){
	//	derivP1[irow] = new std::complex<double>[2 * numOfReferenceVariables * numOfOutputAndInputVariables];
	//}
	//std::complex<double>** derivInvP1 = new std::complex<double>*[numOfReferenceVariables * numOfReferenceVariables];
	//for( int irow = 0; irow < numOfReferenceVariables * numOfReferenceVariables; ++irow ){
	//	derivInvP1[irow] = new std::complex<double>[2 * numOfReferenceVariables * numOfOutputAndInputVariables];
	//}
	//std::complex<double>** derivQ1 = new std::complex<double>*[numOfReferenceVariables * numOfOutputAndInputVariables];
	//for( int irow = 0; irow < numOfReferenceVariables * numOfOutputAndInputVariables; ++irow ){
	//	derivQ1[irow] = new std::complex<double>[2 * numOfReferenceVariables * numOfOutputAndInputVariables];
	//}
	//double** derivResp1 = new double*[2 * numOfReferenceVariables * numOfOutputAndInputVariables];
	//for( int irow = 0; irow < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++irow ){
	//	derivResp1[irow] = new double[2 * numOfReferenceVariables * numOfOutputAndInputVariables];
	//}
	double** predictors = new double*[2 * numOfReferenceVariables * numOfOutputAndInputVariables];
	for( int irow = 0; irow < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++irow ){
		predictors[irow] = new double[2 * numOfReferenceVariables * numOfOutputAndInputVariables];
	}
	double** actuals = new double*[2 * numOfReferenceVariables * numOfOutputAndInputVariables];
	for( int irow = 0; irow < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++irow ){
		actuals[irow] = new double[2 * numOfReferenceVariables * numOfOutputAndInputVariables];
	}
	const double epsilons[4] = { 0.001, 0.01, 0.1, 1.0 };
	for( int idresp = 0; idresp < 4; ++idresp ){
		for( int irow = 0; irow < numOfOutputAndInputVariables; ++irow ){
			for( int icol = 0; icol < numOfReferenceVariables; ++icol ){
				for( int isign = 0; isign < 2; ++isign ){
					std::complex<double>** complexResidualsOrg = new std::complex<double>*[numOfOutputAndInputVariables];
					for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
						complexResidualsOrg[iVar] = new std::complex<double>[numSegments]; 
					}
					double* MDOrg = new double[numSegments];
					double* weightsOrg = new double[numSegments];
					std::complex<double>** respOrg = new std::complex<double>*[numOfOutputAndInputVariables];
					for( int i = 0; i < numOfOutputAndInputVariables; ++i ){
						respOrg[i] = new std::complex<double>[numOfReferenceVariables];
					}
					std::complex<double>** complexResidualsMod = new std::complex<double>*[numOfOutputAndInputVariables];
					for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
						complexResidualsMod[iVar] = new std::complex<double>[numSegments]; 
					}
					double* MDMod = new double[numSegments];
					double* weightsMod = new double[numSegments];
					std::complex<double>** respMod = new std::complex<double>*[numOfOutputAndInputVariables];
					for( int i = 0; i < numOfOutputAndInputVariables; ++i ){
						respMod[i] = new std::complex<double>[numOfReferenceVariables];
					}
					for( int irow2 = 0; irow2 < numOfOutputAndInputVariables; ++irow2 ){
						for( int icol2 = 0; icol2 < numOfReferenceVariables; ++icol2 ){
							respOrg[irow2][icol2] = resp[irow2][icol2];
							respMod[irow2][icol2] = resp[irow2][icol2];
						}
					}
					//double dresp(0.0);
					const double dresp = std::abs(resp[irow][icol]) * epsilons[idresp];
					if( isign == 0 ){
						//dresp = resp[irow][icol].real() * epsilons[idresp];
						respMod[irow][icol] += std::complex<double>( dresp, 0.0 );
					}else{
						//dresp = resp[irow][icol].imag() * epsilons[idresp];
						respMod[irow][icol] += std::complex<double>( 0.0, dresp );
					}
					// Calculate complex residuals
					calculateComplexResiduals( numSegments, ftval, respOrg, complexResidualsOrg );
					calculateComplexResiduals( numSegments, ftval, respMod, complexResidualsMod );
					// Calculate Mahalanobis distance
					calculateMD( numSegments, numOfOutputAndInputVariables, complexResidualsOrg,
						variancesWithoutScale, MDOrg );
					calculateMD( numSegments, numOfOutputAndInputVariables, complexResidualsMod,
						variancesWithoutScale, MDMod );
					const int index = 2 * icol + isign + irow * 2 * numOfReferenceVariables;
					//for( int irow2 = 0; irow2 < numSegments; ++irow2 ){
					//	const double deriv = ( MDMod[irow2] - MDOrg[irow2] ) / dresp;
					//	derivMD1[irow2][index] = deriv;
					//}

					// Calculate original weights
					for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
						const double val = MDOrg[iSeg] / scale;
						weightsOrg[iSeg] = UtilRobust::calculateTukeysBiweights(val, paramC);
					}
					for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
						const double val = MDMod[iSeg] / scale;
						weightsMod[iSeg] = UtilRobust::calculateTukeysBiweights(val, paramC);
					}
					std::complex<double> PMatrixOrg[2][2] = { czero, czero, czero, czero };
					for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
						PMatrixOrg[0][0] += ftval[rr0][iSeg] * std::conj(ftval[rr0][iSeg]) * weightsOrg[iSeg];
						PMatrixOrg[1][0] += ftval[rr1][iSeg] * std::conj(ftval[rr0][iSeg]) * weightsOrg[iSeg];
						PMatrixOrg[0][1] += ftval[rr0][iSeg] * std::conj(ftval[rr1][iSeg]) * weightsOrg[iSeg];
						PMatrixOrg[1][1] += ftval[rr1][iSeg] * std::conj(ftval[rr1][iSeg]) * weightsOrg[iSeg];
					}
					std::complex<double> PMatrixMod[2][2] = { czero, czero, czero, czero };
					for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
						PMatrixMod[0][0] += ftval[rr0][iSeg] * std::conj(ftval[rr0][iSeg]) * weightsMod[iSeg];
						PMatrixMod[1][0] += ftval[rr1][iSeg] * std::conj(ftval[rr0][iSeg]) * weightsMod[iSeg];
						PMatrixMod[0][1] += ftval[rr0][iSeg] * std::conj(ftval[rr1][iSeg]) * weightsMod[iSeg];
						PMatrixMod[1][1] += ftval[rr1][iSeg] * std::conj(ftval[rr1][iSeg]) * weightsMod[iSeg];
					}
					//derivP1[0][index] = ( PMatrixMod[0][0] - PMatrixOrg[0][0] ) / dresp;
					//derivP1[1][index] = ( PMatrixMod[1][0] - PMatrixOrg[1][0] ) / dresp;
					//derivP1[2][index] = ( PMatrixMod[0][1] - PMatrixOrg[0][1] ) / dresp;
					//derivP1[3][index] = ( PMatrixMod[1][1] - PMatrixOrg[1][1] ) / dresp;

					const std::complex<double> detOrg = PMatrixOrg[0][0] * PMatrixOrg[1][1] - PMatrixOrg[1][0] * PMatrixOrg[0][1];
					const std::complex<double> PInvMatrixOrg[2][2] = { PMatrixOrg[1][1]/detOrg, -PMatrixOrg[0][1]/detOrg, -PMatrixOrg[1][0]/detOrg, PMatrixOrg[0][0]/detOrg };
					const std::complex<double> detMod = PMatrixMod[0][0] * PMatrixMod[1][1] - PMatrixMod[1][0] * PMatrixMod[0][1];
					const std::complex<double> PInvMatrixMod[2][2] = { PMatrixMod[1][1]/detMod, -PMatrixMod[0][1]/detMod, -PMatrixMod[1][0]/detMod, PMatrixMod[0][0]/detMod };
					//derivInvP1[0][index] = ( PInvMatrixMod[0][0] - PInvMatrixOrg[0][0] ) / dresp;
					//derivInvP1[1][index] = ( PInvMatrixMod[1][0] - PInvMatrixOrg[1][0] ) / dresp;
					//derivInvP1[2][index] = ( PInvMatrixMod[0][1] - PInvMatrixOrg[0][1] ) / dresp;
					//derivInvP1[3][index] = ( PInvMatrixMod[1][1] - PInvMatrixOrg[1][1] ) / dresp;

					std::complex<double>** QMatrixOrg = new std::complex<double>*[numOfOutputAndInputVariables];
					std::complex<double>** QMatrixMod = new std::complex<double>*[numOfOutputAndInputVariables];
					for( int irow = 0; irow < numOfOutputAndInputVariables; ++irow ){
						QMatrixOrg[irow] = new std::complex<double>[numOfReferenceVariables];
						QMatrixMod[irow] = new std::complex<double>[numOfReferenceVariables];
						for( int icol = 0; icol < numOfReferenceVariables; ++icol ){
							QMatrixOrg[irow][icol] = czero;
							QMatrixMod[irow][icol] = czero;
						}
					}
					for( int irr = 0; irr < numOfReferenceVariables; ++irr ){
						const int rr = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, irr );
						int iVar(0);
						for( int iOut = 0; iOut < numOfOutputVariables; ++iOut ){
							const int index = ptrControl->getChannelIndex( CommonParameters::OUTPUT, iOut );
							for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
								QMatrixOrg[iVar][irr] += ftval[index][iSeg] * std::conj(ftval[rr][iSeg]) * weightsOrg[iSeg];
								QMatrixMod[iVar][irr] += ftval[index][iSeg] * std::conj(ftval[rr][iSeg]) * weightsMod[iSeg];
							}
							assert(index == iVar);
							++iVar;
						}
						for( int iInp = 0; iInp < numOfInputVariables; ++iInp ){
							const int index = ptrControl->getChannelIndex( CommonParameters::INPUT, iInp );
							for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
								QMatrixOrg[iVar][irr] += ftval[index][iSeg] * std::conj(ftval[rr][iSeg]) * weightsOrg[iSeg];
								QMatrixMod[iVar][irr] += ftval[index][iSeg] * std::conj(ftval[rr][iSeg]) * weightsMod[iSeg];
							}
							assert(index == iVar);
							++iVar;
						}
					}
					//for( int irr = 0; irr < numOfReferenceVariables; ++irr ){
					//	for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
					//		const int irow = iVar + numOfOutputAndInputVariables * irr;
					//		derivQ1[irow][index] = ( QMatrixMod[iVar][irr] - QMatrixOrg[iVar][irr] ) / dresp;
					//	}
					//}
					for( int irow2 = 0; irow2 < numOfOutputAndInputVariables; ++irow2 ){
						for( int icol2 = 0; icol2 < numOfReferenceVariables; ++icol2 ){
							std::complex<double> valueOrg = czero;
							std::complex<double> valueMod = czero;
							for( int i = 0; i < numOfReferenceVariables; ++i ){
								valueOrg += QMatrixOrg[irow2][i] * PInvMatrixOrg[i][icol2];
								valueMod += QMatrixMod[irow2][i] * PInvMatrixMod[i][icol2];
							}
							const int index2 = 2 * icol2 + irow2 * 2 * numOfReferenceVariables;
							//derivResp1[index2  ][index] = ( valueMod.real() - valueOrg.real() ) / dresp;
							//derivResp1[index2+1][index] = ( valueMod.imag() - valueOrg.imag() ) / dresp;
							predictors[index2  ][index] = valueOrg.real() + derivativesRegardingResps[index2  ][index] * dresp;
							predictors[index2+1][index] = valueOrg.imag() + derivativesRegardingResps[index2+1][index] * dresp;
							actuals[index2  ][index] = valueMod.real();
							actuals[index2+1][index] = valueMod.imag();
						}
					}
					for( int irow2 = 0; irow2 < numOfOutputAndInputVariables; ++irow2 ){
						delete [] QMatrixOrg[irow2];
					}
					delete [] QMatrixOrg;
					for( int irow2 = 0; irow2 < numOfOutputAndInputVariables; ++irow2 ){
						delete [] QMatrixMod[irow2];
					}
					delete [] QMatrixMod;
					for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
						delete [] complexResidualsOrg[iVar]; 
					}
					delete [] complexResidualsOrg;
					delete [] MDOrg;
					delete [] weightsOrg;
					for( int i = 0; i < numOfOutputAndInputVariables; ++i ){
						delete [] respOrg[i];
					}
					delete [] respOrg;
					for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
						delete [] complexResidualsMod[iVar]; 
					}
					delete [] complexResidualsMod;
					delete [] MDMod;
					delete [] weightsMod;
					for( int i = 0; i < numOfOutputAndInputVariables; ++i ){
						delete [] respMod[i];
					}
					delete [] respMod;
				}
			}
		}
		std::cout << "epsilon " << epsilons[idresp] << std::endl;
		std::cout << "predictors" << std::endl;
		for( int row = 0; row < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++row ){
			for( int col = 0; col < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++col ){
				std::cout << predictors[row][col] <<",";
			}
			std::cout << std::endl;
		}
		std::cout << "actuals" << std::endl;
		for( int row = 0; row < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++row ){
			for( int col = 0; col < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++col ){
				std::cout << actuals[row][col] <<",";
			}
			std::cout << std::endl;
		}
	}
	//for( int irow = 0; irow < numSegments; ++irow ){
	//	delete [] derivMD1[irow];
	//}
	//delete [] derivMD1;
	//for( int irow = 0; irow < numOfReferenceVariables*numOfReferenceVariables; ++irow ){
	//	delete [] derivP1[irow];
	//}
	//delete [] derivP1;
	//for( int irow = 0; irow < numOfReferenceVariables*numOfReferenceVariables; ++irow ){
	//	delete [] derivInvP1[irow];
	//}
	//delete [] derivInvP1;
	//for( int irow = 0; irow < numOfReferenceVariables * numOfOutputAndInputVariables; ++irow ){
	//	delete [] derivQ1[irow];
	//}
	//delete [] derivQ1;
	//for( int irow = 0; irow < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++irow ){
	//	delete [] derivResp1[irow];
	//}
	//delete [] derivResp1;
	for( int irow = 0; irow < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++irow ){
		delete [] predictors[irow];
	}
	delete [] predictors;
	for( int irow = 0; irow < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++irow ){
		delete [] actuals[irow];
	}
	delete [] actuals;
#endif

	// Partial derivatives regarding variances without scale
	for( int irow = 0; irow < numOfReferenceVariables * numOfOutputAndInputVariables; ++irow ){
		for( int icol = 0; icol < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++icol ){
			complexDerivatives[irow][icol] = czero;// Zero clear
		}
	}
	for( int irowL = 0; irowL < numOfReferenceVariables * numOfOutputAndInputVariables; ++irowL ){
		for( int icolR = 0; icolR < numOfOutputAndInputVariables; ++icolR ){
			std::complex<double> value = czero;
			// The 1st term
			for( int i = 0; i < numOfReferenceVariables * numOfOutputAndInputVariables; ++i ){
				value -= PInvTIMatrix[irowL][i] * matrixVeceh2[i][icolR];
			}
			// The 2nd term
			for( int i = 0; i < numOfReferenceVariables * numOfReferenceVariables; ++i ){
				value += IQPInvTPInvMatrix[irowL][i] * matrixVechh2[i][icolR];
			}
			complexDerivatives[irowL][icolR] = value;
		}
	}
	for( int irow = 0; irow < numOfReferenceVariables * numOfOutputAndInputVariables; ++irow ){
		delete [] matrixVeceh2[irow];
	}
	delete [] matrixVeceh2;
	for( int irow = 0; irow < numOfReferenceVariables * numOfReferenceVariables; ++irow ){
		delete [] matrixVechh2[irow];
	}
	delete [] matrixVechh2;
	for( int irow = 0; irow < numOfReferenceVariables * numOfOutputAndInputVariables; ++irow ){
		const int ir = irow / numOfOutputAndInputVariables;
		const int iq = irow % numOfOutputAndInputVariables;
		const int index = 2 * ir + iq * 2 * numOfReferenceVariables;
		for( int icol = 0; icol < numOfOutputAndInputVariables; ++icol ){
			derivativesRegardingVariancesWithoutScale[index  ][icol] = complexDerivatives[irow][icol].real();
			derivativesRegardingVariancesWithoutScale[index+1][icol] = complexDerivatives[irow][icol].imag();
		}
	}
#ifdef _DEBUG_WRITE
	std::cout << "drv1" << std::endl;
	std::cout << "[";
	for( int row = 0; row < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++row ){
		for( int col = 0; col < numOfOutputAndInputVariables; ++col ){
			std::cout << derivativesRegardingVariancesWithoutScale[row][col] <<" ";
		}
		if( row+1 < 2 * numOfReferenceVariables * numOfOutputAndInputVariables ){
			std::cout << ";";
		}
	}
	std::cout << "]" << std::endl;
#endif
#ifdef _DEBUG_WRITE
	double** derivMD2 = new double*[numSegments];
	for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
		derivMD2[iSeg] = new double[numOfOutputAndInputVariables];
	}
	std::complex<double>** derivP2 = new std::complex<double>*[numOfReferenceVariables * numOfReferenceVariables];
	for( int irow = 0; irow < numOfReferenceVariables * numOfReferenceVariables; ++irow ){
		derivP2[irow] = new std::complex<double>[numOfOutputAndInputVariables];
	}
	std::complex<double>** derivInvP2 = new std::complex<double>*[numOfReferenceVariables * numOfReferenceVariables];
	for( int irow = 0; irow < numOfReferenceVariables * numOfReferenceVariables; ++irow ){
		derivInvP2[irow] = new std::complex<double>[numOfOutputAndInputVariables];
	}
	std::complex<double>** derivQ2 = new std::complex<double>*[numOfReferenceVariables * numOfOutputAndInputVariables];
	for( int irow = 0; irow < numOfReferenceVariables * numOfOutputAndInputVariables; ++irow ){
		derivQ2[irow] = new std::complex<double>[numOfOutputAndInputVariables];
	}
	double** derivResp2 = new double*[2 * numOfReferenceVariables * numOfOutputAndInputVariables];
	for( int irow = 0; irow < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++irow ){
		derivResp2[irow] = new double[numOfOutputAndInputVariables];
	}
	for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
		std::complex<double>** complexResidualsOrg = new std::complex<double>*[numOfOutputAndInputVariables];
		for( int iVar2 = 0; iVar2 < numOfOutputAndInputVariables; ++iVar2 ){
			complexResidualsOrg[iVar2] = new std::complex<double>[numSegments]; 
		}
		double* MDOrg = new double[numSegments];
		double* weightsOrg = new double[numSegments];
		double* variancesWithoutScaleOrg = new double[numOfOutputAndInputVariables];
		double* MDMod = new double[numSegments];
		double* weightsMod = new double[numSegments];
		double* variancesWithoutScaleMod = new double[numOfOutputAndInputVariables];
		for( int iVar2 = 0; iVar2 < numOfOutputAndInputVariables; ++iVar2 ){
			variancesWithoutScaleOrg[iVar2] = variancesWithoutScale[iVar2];
			variancesWithoutScaleMod[iVar2] = variancesWithoutScale[iVar2];
		}
		const double dvar = variancesWithoutScale[iVar] * 0.001;
		variancesWithoutScaleMod[iVar] += dvar;
		// Calculate complex residuals
		calculateComplexResiduals( numSegments, ftval, resp, complexResidualsOrg );
		// Calculate Mahalanobis distance
		calculateMD( numSegments, numOfOutputAndInputVariables, complexResidualsOrg,
			variancesWithoutScaleOrg, MDOrg );
		calculateMD( numSegments, numOfOutputAndInputVariables, complexResidualsOrg,
			variancesWithoutScaleMod, MDMod );
		for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
			const double deriv = ( MDMod[iSeg] - MDOrg[iSeg] ) / dvar;
			derivMD2[iSeg][iVar] = deriv;
		}

		// Calculate original weights
		for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
			const double val = MDOrg[iSeg] / scale;
			weightsOrg[iSeg] = RobustWeightTukeysBiweights::calculateWeights(val, paramC);
		}
		for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
			const double val = MDMod[iSeg] / scale;
			weightsMod[iSeg] = RobustWeightTukeysBiweights::calculateWeights(val, paramC);
		}
		std::complex<double> PMatrixOrg[2][2] = { czero, czero, czero, czero };
		for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
			PMatrixOrg[0][0] += ftval[rr0][iSeg] * std::conj(ftval[rr0][iSeg]) * weightsOrg[iSeg];
			PMatrixOrg[1][0] += ftval[rr1][iSeg] * std::conj(ftval[rr0][iSeg]) * weightsOrg[iSeg];
			PMatrixOrg[0][1] += ftval[rr0][iSeg] * std::conj(ftval[rr1][iSeg]) * weightsOrg[iSeg];
			PMatrixOrg[1][1] += ftval[rr1][iSeg] * std::conj(ftval[rr1][iSeg]) * weightsOrg[iSeg];
		}
		std::complex<double> PMatrixMod[2][2] = { czero, czero, czero, czero };
		for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
			PMatrixMod[0][0] += ftval[rr0][iSeg] * std::conj(ftval[rr0][iSeg]) * weightsMod[iSeg];
			PMatrixMod[1][0] += ftval[rr1][iSeg] * std::conj(ftval[rr0][iSeg]) * weightsMod[iSeg];
			PMatrixMod[0][1] += ftval[rr0][iSeg] * std::conj(ftval[rr1][iSeg]) * weightsMod[iSeg];
			PMatrixMod[1][1] += ftval[rr1][iSeg] * std::conj(ftval[rr1][iSeg]) * weightsMod[iSeg];
		}
		derivP2[0][iVar] = ( PMatrixMod[0][0] - PMatrixOrg[0][0] ) / dvar;
		derivP2[1][iVar] = ( PMatrixMod[1][0] - PMatrixOrg[1][0] ) / dvar;
		derivP2[2][iVar] = ( PMatrixMod[0][1] - PMatrixOrg[0][1] ) / dvar;
		derivP2[3][iVar] = ( PMatrixMod[1][1] - PMatrixOrg[1][1] ) / dvar;

		const std::complex<double> detOrg = PMatrixOrg[0][0] * PMatrixOrg[1][1] - PMatrixOrg[1][0] * PMatrixOrg[0][1];
		const std::complex<double> PInvMatrixOrg[2][2] = { PMatrixOrg[1][1]/detOrg, -PMatrixOrg[0][1]/detOrg, -PMatrixOrg[1][0]/detOrg, PMatrixOrg[0][0]/detOrg };
		const std::complex<double> detMod = PMatrixMod[0][0] * PMatrixMod[1][1] - PMatrixMod[1][0] * PMatrixMod[0][1];
		const std::complex<double> PInvMatrixMod[2][2] = { PMatrixMod[1][1]/detMod, -PMatrixMod[0][1]/detMod, -PMatrixMod[1][0]/detMod, PMatrixMod[0][0]/detMod };
		derivInvP2[0][iVar] = ( PInvMatrixMod[0][0] - PInvMatrixOrg[0][0] ) / dvar;
		derivInvP2[1][iVar] = ( PInvMatrixMod[1][0] - PInvMatrixOrg[1][0] ) / dvar;
		derivInvP2[2][iVar] = ( PInvMatrixMod[0][1] - PInvMatrixOrg[0][1] ) / dvar;
		derivInvP2[3][iVar] = ( PInvMatrixMod[1][1] - PInvMatrixOrg[1][1] ) / dvar;

		std::complex<double>** QMatrixOrg = new std::complex<double>*[numOfOutputAndInputVariables];
		std::complex<double>** QMatrixMod = new std::complex<double>*[numOfOutputAndInputVariables];
		for( int irow = 0; irow < numOfOutputAndInputVariables; ++irow ){
			QMatrixOrg[irow] = new std::complex<double>[numOfReferenceVariables];
			QMatrixMod[irow] = new std::complex<double>[numOfReferenceVariables];
			for( int icol = 0; icol < numOfReferenceVariables; ++icol ){
				QMatrixOrg[irow][icol] = czero;
				QMatrixMod[irow][icol] = czero;
			}
		}
		for( int irr = 0; irr < numOfReferenceVariables; ++irr ){
			const int rr = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, irr );
			int iVar(0);
			for( int iOut = 0; iOut < numOfOutputVariables; ++iOut ){
				const int index = ptrControl->getChannelIndex( CommonParameters::OUTPUT, iOut );
				for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
					QMatrixOrg[iVar][irr] += ftval[index][iSeg] * std::conj(ftval[rr][iSeg]) * weightsOrg[iSeg];
					QMatrixMod[iVar][irr] += ftval[index][iSeg] * std::conj(ftval[rr][iSeg]) * weightsMod[iSeg];
				}
				assert(index == iVar);
				++iVar;
			}
			for( int iInp = 0; iInp < numOfInputVariables; ++iInp ){
				const int index = ptrControl->getChannelIndex( CommonParameters::INPUT, iInp );
				for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
					QMatrixOrg[iVar][irr] += ftval[index][iSeg] * std::conj(ftval[rr][iSeg]) * weightsOrg[iSeg];
					QMatrixMod[iVar][irr] += ftval[index][iSeg] * std::conj(ftval[rr][iSeg]) * weightsMod[iSeg];
				}
				assert(index == iVar);
				++iVar;
			}
		}
		for( int irr = 0; irr < numOfReferenceVariables; ++irr ){
			for( int iVar2 = 0; iVar2 < numOfOutputAndInputVariables; ++iVar2 ){
				const int irow = iVar + numOfOutputAndInputVariables * irr;
				derivQ2[irow][iVar2] = ( QMatrixMod[iVar2][irr] - QMatrixOrg[iVar2][irr] ) / dvar;
			}
		}
		for( int irow2 = 0; irow2 < numOfOutputAndInputVariables; ++irow2 ){
			for( int icol2 = 0; icol2 < numOfReferenceVariables; ++icol2 ){
				std::complex<double> valueOrg = czero;
				std::complex<double> valueMod = czero;
				for( int i = 0; i < numOfReferenceVariables; ++i ){
					valueOrg += QMatrixOrg[irow2][i] * PInvMatrixOrg[i][icol2];
					valueMod += QMatrixMod[irow2][i] * PInvMatrixMod[i][icol2];
				}
				const int index2 = 2 * icol2 + irow2 * 2 * numOfReferenceVariables;
				derivResp2[index2  ][iVar] = ( valueMod.real() - valueOrg.real() ) / dvar;
				derivResp2[index2+1][iVar] = ( valueMod.imag() - valueOrg.imag() )/ dvar;
			}
		}
		for( int irow2 = 0; irow2 < numOfOutputAndInputVariables; ++irow2 ){
			delete [] QMatrixOrg[irow2];
		}
		delete [] QMatrixOrg;
		for( int irow2 = 0; irow2 < numOfOutputAndInputVariables; ++irow2 ){
			delete [] QMatrixMod[irow2];
		}
		delete [] QMatrixMod;

		for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
			delete [] complexResidualsOrg[iVar]; 
		}
		delete [] complexResidualsOrg;
		delete [] MDOrg;
		delete [] weightsOrg;
		delete [] variancesWithoutScaleOrg;
		delete [] MDMod;
		delete [] weightsMod;
		delete [] variancesWithoutScaleMod;
	}
	//std::cout << "[";
	//for( int row = 0; row < numSegments; ++row ){
	//	for( int col = 0; col < numOfOutputAndInputVariables; ++col ){
	//		std::cout << derivMD[row][col] <<" ";
	//	}
	//	if( row+1 < numSegments ){
	//		std::cout << ";";
	//	}
	//}
	//std::cout << "]" << std::endl;
	//std::cout << "[";
	//for( int row = 0; row < numOfReferenceVariables * numOfReferenceVariables; ++row ){
	//	for( int col = 0; col < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++col ){
	//		std::cout << derivP[row][col].real() << "+" << derivP[row][col].imag() <<"im ";
	//	}
	//	if( row+1 < numOfReferenceVariables * numOfReferenceVariables ){
	//		std::cout << ";";
	//	}
	//}
	//std::cout << "]" << std::endl;
	//std::cout << "[";
	//for( int row = 0; row < numOfReferenceVariables * numOfReferenceVariables; ++row ){
	//	for( int col = 0; col < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++col ){
	//		std::cout << derivInvP[row][col].real() << "+" << derivInvP[row][col].imag() <<"im ";
	//	}
	//	if( row+1 < numOfReferenceVariables * numOfReferenceVariables ){
	//		std::cout << ";";
	//	}
	//}
	//std::cout << "]" << std::endl;
	//std::cout << "[";
	//for( int row = 0; row < numOfReferenceVariables * numOfOutputAndInputVariables; ++row ){
	//	for( int col = 0; col < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++col ){
	//		std::cout << derivQ[row][col].real() << "+" << derivQ[row][col].imag() <<"im ";
	//	}
	//	if( row+1 < numOfReferenceVariables * numOfOutputAndInputVariables ){
	//		std::cout << ";";
	//	}
	//}
	//std::cout << "]" << std::endl;
	std::cout << "drv2" << std::endl;
	std::cout << "[";
	for( int row = 0; row < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++row ){
		for( int col = 0; col < numOfOutputAndInputVariables; ++col ){
			std::cout << derivResp2[row][col] <<" ";
		}
		if( row+1 < 2 * numOfReferenceVariables * numOfOutputAndInputVariables ){
			std::cout << ";";
		}
	}
	std::cout << "]" << std::endl;
	for( int irow = 0; irow < numSegments; ++irow ){
		delete [] derivMD2[irow];
	}
	delete [] derivMD2;
	for( int irow = 0; irow < numOfReferenceVariables*numOfReferenceVariables; ++irow ){
		delete [] derivP2[irow];
	}
	delete [] derivP2;
	for( int irow = 0; irow < numOfReferenceVariables*numOfReferenceVariables; ++irow ){
		delete [] derivInvP2[irow];
	}
	delete [] derivInvP2;
	for( int irow = 0; irow < numOfReferenceVariables * numOfOutputAndInputVariables; ++irow ){
		delete [] derivQ2[irow];
	}
	delete [] derivQ2;
	for( int irow = 0; irow < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++irow ){
		delete [] derivResp2[irow];
	}
	delete [] derivResp2;
#endif

	// Partial derivatives regarding scale
	for( int irow = 0; irow < numOfReferenceVariables * numOfOutputAndInputVariables; ++irow ){
		for( int icol = 0; icol < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++icol ){
			complexDerivatives[irow][icol] = czero;// Zero clear
		}
	}
	for( int irowL = 0; irowL < numOfReferenceVariables * numOfOutputAndInputVariables; ++irowL ){
		std::complex<double> value = czero;
		// The 1st term
		for( int i = 0; i < numOfReferenceVariables * numOfOutputAndInputVariables; ++i ){
			value -= PInvTIMatrix[irowL][i] * sumVeceh3[i];
		}
		// The 2nd term
		for( int i = 0; i < numOfReferenceVariables * numOfReferenceVariables; ++i ){
			value += IQPInvTPInvMatrix[irowL][i] * sumVechh3[i];
		}
		complexDerivatives[irowL][1] = value;
	}
	delete [ ]sumVeceh3;
	delete [] sumVechh3;
	for( int irow = 0; irow < numOfReferenceVariables * numOfOutputAndInputVariables; ++irow ){
		const int ir = irow / numOfOutputAndInputVariables;
		const int iq = irow % numOfOutputAndInputVariables;
		const int index = 2 * ir + iq * 2 * numOfReferenceVariables;
		derivativesRegardingScale[index  ] = complexDerivatives[irow][1].real();
		derivativesRegardingScale[index+1] = complexDerivatives[irow][1].imag();
	}
#ifdef _DEBUG_WRITE
	std::cout << "drs1" << std::endl;
	std::cout << "[";
	for( int row = 0; row < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++row ){
		std::cout << derivativesRegardingScale[row] <<" ";
		if( row+1 < 2 * numOfReferenceVariables * numOfOutputAndInputVariables ){
			std::cout << ",";
		}
	}
	std::cout << "]" << std::endl;
#endif
#ifdef _DEBUG_WRITE
	std::complex<double>* derivP3 = new std::complex<double>[numOfReferenceVariables * numOfReferenceVariables];
	std::complex<double>* derivInvP3 = new std::complex<double>[numOfReferenceVariables * numOfReferenceVariables];
	std::complex<double>* derivQ3 = new std::complex<double>[numOfReferenceVariables * numOfOutputAndInputVariables];
	double* derivResp3 = new double[2 * numOfReferenceVariables * numOfOutputAndInputVariables];
	{
		std::complex<double>** complexResidualsOrg = new std::complex<double>*[numOfOutputAndInputVariables];
		for( int iVar2 = 0; iVar2 < numOfOutputAndInputVariables; ++iVar2 ){
			complexResidualsOrg[iVar2] = new std::complex<double>[numSegments]; 
		}
		double* MDOrg = new double[numSegments];
		double* weightsOrg = new double[numSegments];
		double* weightsMod = new double[numSegments];
		const double scaleOrg = scale;
		const double dscale = scale * 0.001;
		const double scaleMod = scale + dscale;
		// Calculate complex residuals
		calculateComplexResiduals( numSegments, ftval, resp, complexResidualsOrg );
		// Calculate Mahalanobis distance
		calculateMD( numSegments, numOfOutputAndInputVariables, complexResidualsOrg,
			variancesWithoutScale, MDOrg );
		// Calculate original weights
		for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
			const double val = MDOrg[iSeg] / scaleOrg;
			weightsOrg[iSeg] = RobustWeightTukeysBiweights::calculateWeights(val, paramC);
		}
		for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
			const double val = MDOrg[iSeg] / scaleMod;
			weightsMod[iSeg] = RobustWeightTukeysBiweights::calculateWeights(val, paramC);
		}
		std::complex<double> PMatrixOrg[2][2] = { czero, czero, czero, czero };
		for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
			PMatrixOrg[0][0] += ftval[rr0][iSeg] * std::conj(ftval[rr0][iSeg]) * weightsOrg[iSeg];
			PMatrixOrg[1][0] += ftval[rr1][iSeg] * std::conj(ftval[rr0][iSeg]) * weightsOrg[iSeg];
			PMatrixOrg[0][1] += ftval[rr0][iSeg] * std::conj(ftval[rr1][iSeg]) * weightsOrg[iSeg];
			PMatrixOrg[1][1] += ftval[rr1][iSeg] * std::conj(ftval[rr1][iSeg]) * weightsOrg[iSeg];
		}
		std::complex<double> PMatrixMod[2][2] = { czero, czero, czero, czero };
		for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
			PMatrixMod[0][0] += ftval[rr0][iSeg] * std::conj(ftval[rr0][iSeg]) * weightsMod[iSeg];
			PMatrixMod[1][0] += ftval[rr1][iSeg] * std::conj(ftval[rr0][iSeg]) * weightsMod[iSeg];
			PMatrixMod[0][1] += ftval[rr0][iSeg] * std::conj(ftval[rr1][iSeg]) * weightsMod[iSeg];
			PMatrixMod[1][1] += ftval[rr1][iSeg] * std::conj(ftval[rr1][iSeg]) * weightsMod[iSeg];
		}
		derivP3[0] = ( PMatrixMod[0][0] - PMatrixOrg[0][0] ) / dscale;
		derivP3[1] = ( PMatrixMod[1][0] - PMatrixOrg[1][0] ) / dscale;
		derivP3[2] = ( PMatrixMod[0][1] - PMatrixOrg[0][1] ) / dscale;
		derivP3[3] = ( PMatrixMod[1][1] - PMatrixOrg[1][1] ) / dscale;
		const std::complex<double> detOrg = PMatrixOrg[0][0] * PMatrixOrg[1][1] - PMatrixOrg[1][0] * PMatrixOrg[0][1];
		const std::complex<double> PInvMatrixOrg[2][2] = { PMatrixOrg[1][1]/detOrg, -PMatrixOrg[0][1]/detOrg, -PMatrixOrg[1][0]/detOrg, PMatrixOrg[0][0]/detOrg };
		const std::complex<double> detMod = PMatrixMod[0][0] * PMatrixMod[1][1] - PMatrixMod[1][0] * PMatrixMod[0][1];
		const std::complex<double> PInvMatrixMod[2][2] = { PMatrixMod[1][1]/detMod, -PMatrixMod[0][1]/detMod, -PMatrixMod[1][0]/detMod, PMatrixMod[0][0]/detMod };
		derivInvP3[0] = ( PInvMatrixMod[0][0] - PInvMatrixOrg[0][0] ) / dscale;
		derivInvP3[1] = ( PInvMatrixMod[1][0] - PInvMatrixOrg[1][0] ) / dscale;
		derivInvP3[2] = ( PInvMatrixMod[0][1] - PInvMatrixOrg[0][1] ) / dscale;
		derivInvP3[3] = ( PInvMatrixMod[1][1] - PInvMatrixOrg[1][1] ) / dscale;
		std::complex<double>** QMatrixOrg = new std::complex<double>*[numOfOutputAndInputVariables];
		std::complex<double>** QMatrixMod = new std::complex<double>*[numOfOutputAndInputVariables];
		for( int irow = 0; irow < numOfOutputAndInputVariables; ++irow ){
			QMatrixOrg[irow] = new std::complex<double>[numOfReferenceVariables];
			QMatrixMod[irow] = new std::complex<double>[numOfReferenceVariables];
			for( int icol = 0; icol < numOfReferenceVariables; ++icol ){
				QMatrixOrg[irow][icol] = czero;
				QMatrixMod[irow][icol] = czero;
			}
		}
		for( int irr = 0; irr < numOfReferenceVariables; ++irr ){
			const int rr = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, irr );
			int iVar(0);
			for( int iOut = 0; iOut < numOfOutputVariables; ++iOut ){
				const int index = ptrControl->getChannelIndex( CommonParameters::OUTPUT, iOut );
				for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
					QMatrixOrg[iVar][irr] += ftval[index][iSeg] * std::conj(ftval[rr][iSeg]) * weightsOrg[iSeg];
					QMatrixMod[iVar][irr] += ftval[index][iSeg] * std::conj(ftval[rr][iSeg]) * weightsMod[iSeg];
				}
				assert(index == iVar);
				++iVar;
			}
			for( int iInp = 0; iInp < numOfInputVariables; ++iInp ){
				const int index = ptrControl->getChannelIndex( CommonParameters::INPUT, iInp );
				for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
					QMatrixOrg[iVar][irr] += ftval[index][iSeg] * std::conj(ftval[rr][iSeg]) * weightsOrg[iSeg];
					QMatrixMod[iVar][irr] += ftval[index][iSeg] * std::conj(ftval[rr][iSeg]) * weightsMod[iSeg];
				}
				assert(index == iVar);
				++iVar;
			}
		}
		for( int irr = 0; irr < numOfReferenceVariables; ++irr ){
			for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
				const int irow = iVar + numOfOutputAndInputVariables * irr;
				derivQ3[irow] = ( QMatrixMod[iVar][irr] - QMatrixOrg[iVar][irr] ) / dscale;
			}
		}
		for( int irow2 = 0; irow2 < numOfOutputAndInputVariables; ++irow2 ){
			for( int icol2 = 0; icol2 < numOfReferenceVariables; ++icol2 ){
				std::complex<double> valueOrg = czero;
				std::complex<double> valueMod = czero;
				for( int i = 0; i < numOfReferenceVariables; ++i ){
					valueOrg += QMatrixOrg[irow2][i] * PInvMatrixOrg[i][icol2];
					valueMod += QMatrixMod[irow2][i] * PInvMatrixMod[i][icol2];
				}
				const int index2 = 2 * icol2 + irow2 * 2 * numOfReferenceVariables;
				derivResp3[index2  ] = ( valueMod.real() - valueOrg.real() ) / dscale;
				derivResp3[index2+1] = ( valueMod.imag() - valueOrg.imag() )/ dscale;
			}
		}
		for( int irow2 = 0; irow2 < numOfOutputAndInputVariables; ++irow2 ){
			delete [] QMatrixOrg[irow2];
		}
		delete [] QMatrixOrg;
		for( int irow2 = 0; irow2 < numOfOutputAndInputVariables; ++irow2 ){
			delete [] QMatrixMod[irow2];
		}
		delete [] QMatrixMod;

		for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
			delete [] complexResidualsOrg[iVar]; 
		}
		delete [] complexResidualsOrg;
		delete [] MDOrg;
		delete [] weightsOrg;
		delete [] weightsMod;
	}
	//std::cout << "[";
	//for( int row = 0; row < numOfReferenceVariables * numOfReferenceVariables; ++row ){
	//	std::cout << derivP[row].real() << "+" << derivP[row].imag() <<"im ";
	//	if( row+1 < numOfReferenceVariables * numOfReferenceVariables ){
	//		std::cout << ",";
	//	}
	//}
	//std::cout << "]" << std::endl;
	//std::cout << "[";
	//for( int row = 0; row < numOfReferenceVariables * numOfReferenceVariables; ++row ){
	//	std::cout << derivInvP[row].real() << "+" << derivInvP[row].imag() <<"im ";
	//	if( row+1 < numOfReferenceVariables * numOfReferenceVariables ){
	//		std::cout << ",";
	//	}
	//}
	//std::cout << "]" << std::endl;
	//std::cout << "[";
	//for( int row = 0; row < numOfReferenceVariables * numOfOutputAndInputVariables; ++row ){
	//	std::cout << derivQ[row].real() << "+" << derivQ[row].imag() <<"im ";
	//	if( row+1 < numOfReferenceVariables * numOfOutputAndInputVariables ){
	//		std::cout << ",";
	//	}
	//}
	//std::cout << "]" << std::endl;
	std::cout << "drs2" << std::endl;
	std::cout << "[";
	for( int row = 0; row < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++row ){
		std::cout << derivResp3[row] <<" ";
		if( row+1 < 2 * numOfReferenceVariables * numOfOutputAndInputVariables ){
			std::cout << ",";
		}
	}
	std::cout << "]" << std::endl;
	delete [] derivP3;
	delete [] derivInvP3;
	delete [] derivQ3;
	delete [] derivResp3;
#endif

	for( int irow = 0; irow < 2 * numOfOutputAndInputVariables; ++irow ){
		delete [] PInvTIMatrix[irow];
	}
	delete [] PInvTIMatrix;
	for( int irow = 0; irow < numOfReferenceVariables * numOfOutputAndInputVariables; ++irow ){
		delete [] IQPInvTPInvMatrix[irow];
	}
	delete [] IQPInvTPInvMatrix;
	for( int irow = 0; irow < numOfReferenceVariables * numOfOutputAndInputVariables; ++irow ){
		delete [] complexDerivatives[irow];
	}
	delete [] complexDerivatives;

}

// Calculate partial derivatives of scale
void AnalysisMultivariateRegression::calculatePartialDerivativesOfScale( const int numSegments, const double paramB, const double paramC,
	std::complex<double>** ftval, std::complex<double>** resp, const double* const variancesWithoutScale, const double scale, 
	std::complex<double>** complexResiduals, const double* const MD, const double* const weights,
	double* derivativesRegardingResps, double* derivativesRegardingVariancesWithoutScale, double& derivativesRegardingScale ) const{

	const Control* const ptrControl = Control::getInstance();
	const int numOfOutputVariables = ptrControl->getNumOutputVariables();
	const int numOfInputVariables = ptrControl->getNumInputVariables();
	const int numOfOutputAndInputVariables = numOfOutputVariables + numOfInputVariables;
	const int numOfReferenceVariables = ptrControl->getNumRemoteReferenceVariables();
	const int degreeOfFreedom = 2 * numOfOutputAndInputVariables;
	assert(numOfReferenceVariables == 2); 
	const int rr0 = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, 0 );
	const int rr1 = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, 1 );

	const std::complex<double> czero = std::complex<double>( 0.0, 0.0 );

	double* value1 = new double[2 * numOfReferenceVariables * numOfOutputAndInputVariables];
	for( int icol = 0; icol < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++icol ){
		value1[icol] = 0.0;// Zero clear
	}
	double* value2 = new double[numOfOutputAndInputVariables];
	for( int icol = 0; icol < numOfOutputAndInputVariables; ++icol ){
		value2[icol] = 0.0;// Zero clear
	}
	double value3(0.0);
	double* vecSigma = new double[2 * numOfReferenceVariables * numOfOutputAndInputVariables];
	std::complex<double>** hSigmaMatrix = new std::complex<double>*[numOfReferenceVariables];
	for( int irow = 0; irow < numOfReferenceVariables; ++irow ){
		hSigmaMatrix[irow] = new std::complex<double>[numOfOutputAndInputVariables];
	}
//#ifdef _DEBUG_WRITE
//	std::cout << "[";
//#endif
	for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
		calculateVectorForPartialDerivatives( iSeg, ftval, variancesWithoutScale, complexResiduals, hSigmaMatrix, vecSigma);
		for( int icol = 0; icol < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++icol ){
			double temp = vecSigma[icol];
			if( fabs(MD[iSeg]) < CommonParameters::EPS ){
				// Nothing to do
			}else{
				temp *= weights[iSeg];
			}
			value1[icol] -= temp;
//#ifdef _DEBUG_WRITE
//			const int iq = icol / ( 2 * numOfReferenceVariables );
//			const int isign = ( icol / numOfReferenceVariables ) % 2;
//			const int index = 2 * ( icol % numOfReferenceVariables ) + isign + 2 * numOfReferenceVariables * iq;
//			std::cout << -temp <<" ";
//#endif
		}
//#ifdef _DEBUG_WRITE
//		if( iSeg + 1 < numSegments){
//			std::cout << ";";
//		}
//#endif
		for( int icol = 0; icol < numOfOutputAndInputVariables; ++icol ){
			double temp = 0.5 * std::norm(complexResiduals[icol][iSeg]) / pow(variancesWithoutScale[icol], 2);
			if( fabs(MD[iSeg]) < CommonParameters::EPS ){
				// Nothing to do
			}else{
				temp *= weights[iSeg];
			}
			value2[icol] -= temp;
		}
		const double val = MD[iSeg] / scale;
		if( fabs(val) < CommonParameters::EPS ){
			// Nothing to do
		}else{
			const double temp = 2.0 * pow(scale, 2) * RobustWeightTukeysBiweights::calculateLossFunction(val, paramC)
				- pow(MD[iSeg], 2) * weights[iSeg];
			value3 += temp / scale;
		}
	}
//#ifdef _DEBUG_WRITE
//	std::cout << "]" << std::endl;
//#endif
	delete [] vecSigma;
	for( int irow = 0; irow < numOfReferenceVariables; ++irow ){
		delete [] hSigmaMatrix[irow];
	}
	delete [] hSigmaMatrix;
	
//#ifdef _DEBUG_WRITE
//	std::cout << "[";
//	for( int row = 0; row < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++row ){
//		const int iq = row / ( 2 * numOfReferenceVariables );
//		const int isign = ( row / numOfReferenceVariables ) % 2;
//		const int index = 2 * ( row % numOfReferenceVariables ) + isign + 2 * numOfReferenceVariables * iq;
//		std::cout << value1[index] <<" ";
//		if( row+1 < 2 * numOfReferenceVariables * numOfOutputAndInputVariables ){
//			std::cout << ",";
//		}
//	}
//	std::cout << "]" << std::endl;
//#endif

	// Order: Re(Z11),Im(Z11),Re(Z12),Im(Z12),Re(Z21),Im(Z21),Re(Z22),Im(Z22),...,Re(Zq1),Im(Zq1),Re(Zq2),Im(Zq2)
	const double factor = 1.0 / ( 2.0 * scale * static_cast<double>(numSegments) * paramB );
	for( int icol = 0; icol < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++icol ){
		const int iq2 = icol / ( 2 * numOfReferenceVariables );
		const int offset = iq2 * 2 * numOfReferenceVariables;
		const int amari = icol % ( 2 * numOfReferenceVariables );
		if( amari < numOfReferenceVariables ){
			// Real part
			const int ir2 = amari;
			const int index2 = 2 * ir2 + offset;
			derivativesRegardingResps[index2] = factor * value1[icol];
		}else{
			// Imaginary part
			const int ir2 = amari - numOfReferenceVariables;
			const int index2 = 2 * ir2 + 1 + offset;
			derivativesRegardingResps[index2] = factor * value1[icol];
		}
	}
	for( int icol = 0; icol < numOfOutputAndInputVariables; ++icol ){
		derivativesRegardingVariancesWithoutScale[icol] = factor * value2[icol];
	}
	derivativesRegardingScale = factor * value3;

	delete [] value1;
	delete [] value2;

#ifdef _DEBUG_WRITE
	std::cout << "dsr1" << std::endl;
	std::cout << "[";
	for( int row = 0; row < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++row ){
		std::cout << derivativesRegardingResps[row] <<" ";
		if( row+1 < 2 * numOfReferenceVariables * numOfOutputAndInputVariables ){
			std::cout << ",";
		}
	}
	std::cout << "]" << std::endl;
	std::cout << "dsv1" << std::endl;
	std::cout << "[";
	for( int row = 0; row < numOfOutputAndInputVariables; ++row ){
		std::cout << derivativesRegardingVariancesWithoutScale[row] <<" ";
		if( row+1 < numOfOutputAndInputVariables ){
			std::cout << ",";
		}
	}
	std::cout << "]" << std::endl;
	std::cout << "dss1" << std::endl;
	std::cout << derivativesRegardingScale << std::endl;
#endif
#ifdef _DEBUG_WRITE
	const double valq = static_cast<double>(numOfOutputAndInputVariables);
	const double val2q = static_cast<double>(degreeOfFreedom);
	double** derivSumU1 = new double*[numSegments];
	for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
		derivSumU1[iSeg] = new double[2 * numOfReferenceVariables * numOfOutputAndInputVariables];
	}
	double* derivSquaredScale1 = new double[2 * numOfReferenceVariables * numOfOutputAndInputVariables];
	double* derivScale1 = new double[2 * numOfReferenceVariables * numOfOutputAndInputVariables];
	for( int irow = 0; irow < numOfOutputAndInputVariables; ++irow ){
		for( int icol = 0; icol < numOfReferenceVariables; ++icol ){
			for( int isign = 0; isign < 2; ++isign ){
				std::complex<double>** complexResidualsOrg = new std::complex<double>*[numOfOutputAndInputVariables];
				for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
					complexResidualsOrg[iVar] = new std::complex<double>[numSegments]; 
				}
				double* MDOrg = new double[numSegments];
				std::complex<double>** respOrg = new std::complex<double>*[numOfOutputAndInputVariables];
				for( int i = 0; i < numOfOutputAndInputVariables; ++i ){
					respOrg[i] = new std::complex<double>[numOfReferenceVariables];
				}
				std::complex<double>** complexResidualsMod = new std::complex<double>*[numOfOutputAndInputVariables];
				for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
					complexResidualsMod[iVar] = new std::complex<double>[numSegments]; 
				}
				double* MDMod = new double[numSegments];
				std::complex<double>** respMod = new std::complex<double>*[numOfOutputAndInputVariables];
				for( int i = 0; i < numOfOutputAndInputVariables; ++i ){
					respMod[i] = new std::complex<double>[numOfReferenceVariables];
				}

				for( int irow2 = 0; irow2 < numOfOutputAndInputVariables; ++irow2 ){
					for( int icol2 = 0; icol2 < numOfReferenceVariables; ++icol2 ){
						respOrg[irow2][icol2] = resp[irow2][icol2];
						respMod[irow2][icol2] = resp[irow2][icol2];
					}
				}
				double dresp(0.0);
				if( isign == 0 ){
					dresp = resp[irow][icol].real() * 0.001;
					respMod[irow][icol] += std::complex<double>( dresp, 0.0 );
				}else{
					dresp = resp[irow][icol].imag() * 0.001;
					respMod[irow][icol] += std::complex<double>( 0.0, dresp );
				}

				// Calculate complex residuals
				calculateComplexResiduals( numSegments, ftval, respOrg, complexResidualsOrg );
				calculateComplexResiduals( numSegments, ftval, respMod, complexResidualsMod );
				// Calculate Mahalanobis distance
				calculateMD( numSegments, numOfOutputAndInputVariables, complexResidualsOrg,
					variancesWithoutScale, MDOrg );
				calculateMD( numSegments, numOfOutputAndInputVariables, complexResidualsMod,
					variancesWithoutScale, MDMod );

				const int index = 2 * icol + isign + irow * 2 * numOfReferenceVariables;

				for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
					const double valOrg = MDOrg[iSeg] / scale;
					const double valMod = MDMod[iSeg] / scale;
					const double uOrg = RobustWeightTukeysBiweights::calculateLossFunction(valOrg, paramC) * pow(scale, 2);
					const double uMod = RobustWeightTukeysBiweights::calculateLossFunction(valMod, paramC) * pow(scale, 2);
					derivSumU1[iSeg][index] = ( uMod - uOrg ) / dresp;
				}

				double squareScaleOrg(0.0);
				for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
					if( fabs(MDOrg[iSeg]) < CommonParameters::EPS ){
						// rho''(0) / 2 = 1 / 2
						squareScaleOrg += 0.5 * pow(MDOrg[iSeg], 2);
					}else{
						const double val = MDOrg[iSeg] / scale;
						squareScaleOrg += RobustWeightTukeysBiweights::calculateLossFunction(val, paramC) * pow(scale, 2);
					}
				}
				double squareScaleMod(0.0);
				for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
					const double val = MDMod[iSeg] / scale;
					if( fabs(val) < CommonParameters::EPS ){
						// rho''(0) / 2 = 1 / 2
						squareScaleMod += 0.5 * pow(MDMod[iSeg], 2);
					}else{
						squareScaleMod += RobustWeightTukeysBiweights::calculateLossFunction(val, paramC) * pow(scale, 2);
					}
				}

				derivSquaredScale1[index] = ( squareScaleMod - squareScaleOrg ) / dresp;

				squareScaleOrg /= static_cast<double>(numSegments);
				squareScaleOrg /= paramB;
				const double scaleOrg = sqrt(squareScaleOrg);
				squareScaleMod /= static_cast<double>(numSegments);
				squareScaleMod /= paramB;
				const double scaleMod = sqrt(squareScaleMod);

				derivScale1[index] = ( scaleMod - scaleOrg ) / dresp;
				
				for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
					delete [] complexResidualsOrg[iVar]; 
				}
				delete [] complexResidualsOrg;
				delete [] MDOrg;
				for( int i = 0; i < numOfOutputAndInputVariables; ++i ){
					delete [] respOrg[i];
				}
				delete [] respOrg;
				for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
					delete [] complexResidualsMod[iVar]; 
				}
				delete [] complexResidualsMod;
				delete [] MDMod;
				for( int i = 0; i < numOfOutputAndInputVariables; ++i ){
					delete [] respMod[i];
				}
				delete [] respMod;
			}
		}
	}
	//std::cout << "[";
	//for( int row = 0; row < numSegments; ++row ){
	//	for( int col = 0; col < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++col ){
	//		std::cout << derivSumU1[row][col] <<" ";
	//	}
	//	if( row+1 < numSegments ){
	//		std::cout << ";";
	//	}
	//}
	//std::cout << "]" << std::endl;
	//std::cout << "[";
	//for( int row = 0; row < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++row ){
	//	std::cout << derivSquaredScale1[row] <<" ";
	//	if( row+1 < 2 * numOfReferenceVariables * numOfOutputAndInputVariables ){
	//		std::cout << ",";
	//	}
	//}
	//std::cout << "]" << std::endl;
	std::cout << "dsr2" << std::endl;
	std::cout << "[";
	for( int row = 0; row < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++row ){
		std::cout << derivScale1[row] <<" ";
		if( row+1 < 2 * numOfReferenceVariables * numOfOutputAndInputVariables ){
			std::cout << ",";
		}
	}
	std::cout << "]" << std::endl;
	for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
		delete [] derivSumU1[iSeg];
	}
	delete [] derivSumU1;
	delete [] derivSquaredScale1;
	delete [] derivScale1;

	double* derivScale2 = new double[numOfOutputAndInputVariables];
	for( int icol = 0; icol < numOfOutputAndInputVariables; ++icol ){
		std::complex<double>** complexResidualsOrg = new std::complex<double>*[numOfOutputAndInputVariables];
		for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
			complexResidualsOrg[iVar] = new std::complex<double>[numSegments]; 
		}
		double* MDOrg = new double[numSegments];
		double* variancesWithoutScaleOrg = new double[numOfOutputAndInputVariables];
		double* MDMod = new double[numSegments];
		double* variancesWithoutScaleMod = new double[numOfOutputAndInputVariables];
		for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
			variancesWithoutScaleOrg[iVar] = variancesWithoutScale[iVar];
			variancesWithoutScaleMod[iVar] = variancesWithoutScale[iVar];
		}
		const double dvar = variancesWithoutScale[icol] * 0.001;
		variancesWithoutScaleMod[icol] += dvar;
		// Calculate complex residuals
		calculateComplexResiduals( numSegments, ftval, resp, complexResidualsOrg );
		// Calculate Mahalanobis distance
		calculateMD( numSegments, numOfOutputAndInputVariables, complexResidualsOrg,
			variancesWithoutScaleOrg, MDOrg );
		calculateMD( numSegments, numOfOutputAndInputVariables, complexResidualsOrg,
			variancesWithoutScaleMod, MDMod );
				
		double squareScaleOrg(0.0);
		for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
			const double val = MDOrg[iSeg] / scale;
			if( fabs(val) < CommonParameters::EPS ){
				// rho''(0) / 2 = 1 / 2
				squareScaleOrg += 0.5 * pow(MDOrg[iSeg], 2);
			}else{
				squareScaleOrg += RobustWeightTukeysBiweights::calculateLossFunction(val, paramC) * pow(scale, 2);
			}
		}
		double squareScaleMod(0.0);
		for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
			const double val = MDMod[iSeg] / scale;
			if( fabs(val) < CommonParameters::EPS ){
				// rho''(0) / 2 = 1 / 2
				squareScaleMod += 0.5 * pow(MDMod[iSeg], 2);
			}else{
				squareScaleMod += RobustWeightTukeysBiweights::calculateLossFunction(val, paramC) * pow(scale, 2);
			}
		}

		squareScaleOrg /= static_cast<double>(numSegments);
		squareScaleOrg /= paramB;
		const double scaleOrg = sqrt(squareScaleOrg);
		squareScaleMod /= static_cast<double>(numSegments);
		squareScaleMod /= paramB;
		const double scaleMod = sqrt(squareScaleMod);

		derivScale2[icol] = ( scaleMod - scaleOrg ) / dvar;

		for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
			delete [] complexResidualsOrg[iVar]; 
		}
		delete [] complexResidualsOrg;
		delete [] MDOrg;
		delete [] variancesWithoutScaleOrg;
		delete [] MDMod;
		delete [] variancesWithoutScaleMod;
	}
	std::cout << "dsv2" << std::endl;
	std::cout << "[";
	for( int row = 0; row < numOfOutputAndInputVariables; ++row ){
		std::cout << derivScale2[row] << " ";
		if( row+1 < numOfOutputAndInputVariables ){
			std::cout << ",";
		}
	}
	std::cout << "]" << std::endl;

	double derivScale3(0.0);
	{
		std::complex<double>** complexResiduals = new std::complex<double>*[numOfOutputAndInputVariables];
		for( int iVar2 = 0; iVar2 < numOfOutputAndInputVariables; ++iVar2 ){
			complexResiduals[iVar2] = new std::complex<double>[numSegments]; 
		}
		double* MD = new double[numSegments];
		const double scaleOrg = scale;
		const double dscale = scale * 0.001;
		const double scaleMod = scale + dscale;
		// Calculate complex residuals
		calculateComplexResiduals( numSegments, ftval, resp, complexResiduals );
		// Calculate Mahalanobis distance
		calculateMD( numSegments, numOfOutputAndInputVariables, complexResiduals,
			variancesWithoutScale, MD );

		double squareScaleOrg(0.0);
		for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
			const double val = MD[iSeg] / scaleOrg;
			if( fabs(val) < CommonParameters::EPS ){
				// rho''(0) / 2 = 1 / 2
				squareScaleOrg += 0.5 * pow(MD[iSeg], 2);
			}else{
				squareScaleOrg += RobustWeightTukeysBiweights::calculateLossFunction(val, paramC) * pow(scaleOrg, 2);
			}
		}
		double squareScaleMod(0.0);
		for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
			const double val = MD[iSeg] / scaleMod;
			if( fabs(val) < CommonParameters::EPS ){
				// rho''(0) / 2 = 1 / 2
				squareScaleMod += 0.5 * pow(MD[iSeg], 2);
			}else{
				squareScaleMod += RobustWeightTukeysBiweights::calculateLossFunction(val, paramC) * pow(scaleMod, 2);
			}
		}

		squareScaleOrg /= static_cast<double>(numSegments);
		squareScaleOrg /= paramB;
		const double scaleNewOrg = sqrt(squareScaleOrg);
		squareScaleMod /= static_cast<double>(numSegments);
		squareScaleMod /= paramB;
		const double scaleNewMod = sqrt(squareScaleMod);

		derivScale3 = ( scaleNewMod - scaleNewOrg ) / dscale;

		for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
			delete [] complexResiduals[iVar]; 
		}
		delete [] complexResiduals;
		delete [] MD;
	}
	std::cout << "dss2" << std::endl;
	std::cout << derivScale3 << std::endl;
#endif

}

// Calculate partial derivatives of variances without scale for robust bootstrap
void AnalysisMultivariateRegression::calculatePartialDerivativesOfVariancesWithoutScale( const int numSegments, const double paramC,
	std::complex<double>** ftval, std::complex<double>** resp, const double* const variancesWithoutScale, const double scale, const double determinant,
	std::complex<double>** complexResiduals, const double* const MD, const double* const weights,
	double** derivativesRegardingResps, double** derivativesRegardingVariancesWithoutScale, double* derivativesRegardingScale ) const{

	const Control* const ptrControl = Control::getInstance();
	const int numOfOutputVariables = ptrControl->getNumOutputVariables();
	const int numOfInputVariables = ptrControl->getNumInputVariables();
	const int numOfOutputAndInputVariables = numOfOutputVariables + numOfInputVariables;
	const int numOfReferenceVariables = ptrControl->getNumRemoteReferenceVariables();
	const int degreeOfFreedom = 2 * numOfOutputAndInputVariables;
	assert(numOfReferenceVariables == 2); 
	const int rr0 = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, 0 );
	const int rr1 = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, 1 );
	const double valq = static_cast<double>(numOfOutputAndInputVariables);

	const std::complex<double> czero = std::complex<double>( 0.0, 0.0 );

	// Variances with scale
	double* variances = new double[numOfOutputAndInputVariables];
	for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
		variances[iVar] = variancesWithoutScale[iVar] * pow(determinant, 1.0/static_cast<double>(degreeOfFreedom));
	}

	double denominator(0.0);
	for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
		denominator += weights[iSeg] * pow(MD[iSeg] / scale, 2);
	}
	double* term1 = new double[2 * numOfReferenceVariables * numOfOutputAndInputVariables];
	double* term2 = new double[2 * numOfReferenceVariables * numOfOutputAndInputVariables];
	double* term3 = new double[numOfOutputAndInputVariables];
	double* term4 = new double[numOfOutputAndInputVariables];
	double** dSigma1 = new double*[numOfOutputAndInputVariables];
	for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
		dSigma1[iVar] = new double[2 * numOfReferenceVariables * numOfOutputAndInputVariables];
	}
	double** dSigma2 = new double*[numOfOutputAndInputVariables];
	for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
		dSigma2[iVar] = new double[numOfOutputAndInputVariables];
	}
	double* dSigma3 = new double[numOfOutputAndInputVariables];
	double* vecSigma = new double[2 * numOfReferenceVariables * numOfOutputAndInputVariables];
	std::complex<double>** hSigmaMatrix = new std::complex<double>*[numOfReferenceVariables];
	for( int irow = 0; irow < numOfReferenceVariables; ++irow ){
		hSigmaMatrix[irow] = new std::complex<double>[numOfOutputAndInputVariables];
	}
	for( int irow = 0; irow < numOfOutputAndInputVariables; ++irow ){
		const double sc = scale * paramC;
		std::complex<double> termFirst = czero;
		for( int icol = 0; icol < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++icol ){
			term1[icol] = 0.0;// Zero clear
			term2[icol] = 0.0;// Zero clear
		}
		for( int icol = 0; icol < numOfOutputAndInputVariables; ++icol ){
			term3[icol] = 0.0;// Zero clear
			term4[icol] = 0.0;// Zero clear
		}
		double term5(0.0);
		double term6(0.0);
		double numerator (0.0);
		for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
			const double val = MD[iSeg] / scale;
			const double diff = RobustWeightTukeysBiweights::calculateSecondDerivativeOfLossFunction(val, paramC) - RobustWeightTukeysBiweights::calculateWeights(val, paramC);
			double factor1(0.0);
			if( fabs(MD[iSeg]) < CommonParameters::EPS ){
				const double sc = scale * paramC;
				factor1 = 4.0 * ( pow(MD[iSeg],2) / pow(sc,4) - 1.0 / pow(sc,2) );
			}else{
				factor1 = diff / pow(MD[iSeg], 2);
			}
			numerator += std::norm(complexResiduals[irow][iSeg]) * weights[iSeg];
			calculateVectorForPartialDerivatives( iSeg, ftval, variancesWithoutScale, complexResiduals, hSigmaMatrix, vecSigma );
			for( int icol = 0; icol < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++icol ){
				term1[icol] += std::norm(complexResiduals[irow][iSeg]) * factor1 * vecSigma[icol];
				term2[icol] += ( diff + 2.0 * weights[iSeg] ) / pow(scale, 2) * vecSigma[icol];
			}
			for( int icol = 0; icol < numOfOutputAndInputVariables; ++icol ){
				const double temp = std::norm(complexResiduals[icol][iSeg]) / pow(variancesWithoutScale[icol], 2);
				term3[icol] += std::norm(complexResiduals[irow][iSeg]) * 0.5 * factor1 * temp;
				term4[icol] += ( diff * 0.5 + weights[iSeg] ) / pow(scale, 2) * temp;
			}
			term5 += std::norm(complexResiduals[irow][iSeg]) * diff / scale; 
			term6 += ( diff + 2.0 * weights[iSeg] ) / scale * pow(MD[iSeg] / scale, 2);
		}
		for( int icol = 0; icol < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++icol ){
			dSigma1[irow][icol] = - valq / denominator * term1[icol]	+ valq * numerator / pow(denominator,2) * term2[icol];
		}
		for( int icol = 0; icol < numOfOutputAndInputVariables; ++icol ){
			dSigma2[irow][icol] = - valq / denominator * term3[icol]	+ valq * numerator / pow(denominator,2) * term4[icol];
		}
		dSigma3[irow] = - valq / denominator * term5 + valq * numerator / pow(denominator,2) * term6;
	}
	for( int irow = 0; irow < numOfReferenceVariables; ++irow ){
		delete [] hSigmaMatrix[irow];
	}
	delete [] hSigmaMatrix;
	delete [] vecSigma;
	delete [] term1;
	delete [] term2;
	delete [] term3;
	delete [] term4;

//#ifdef _DEBUG_WRITE
//	std::cout << "[";
//	for( int row = 0; row < numOfOutputAndInputVariables; ++row ){
//		for( int col = 0; col < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++col ){
//			const int iq = col / ( 2 * numOfReferenceVariables );
//			const int isign = ( col / numOfReferenceVariables ) % 2;
//			const int index = 2 * ( col % numOfReferenceVariables ) + isign + 2 * numOfReferenceVariables * iq;
//			std::cout << dSigma1[row][index] <<" ";
//		}
//		if( row+1 < numOfOutputAndInputVariables ){
//			std::cout << ";";
//		}
//	}
//	std::cout << "]" << std::endl;
//	std::cout << "[";
//	for( int row = 0; row < numOfOutputAndInputVariables; ++row ){
//		for( int col = 0; col < numOfOutputAndInputVariables; ++col ){
//			std::cout << dSigma2[row][col] <<" ";
//		}
//		if( row+1 < numOfOutputAndInputVariables ){
//			std::cout << ";";
//		}
//	}
//	std::cout << "]" << std::endl;
//	std::cout << "[";
//	for( int row = 0; row < numOfOutputAndInputVariables; ++row ){
//		std::cout << dSigma3[row];
//		if( row+1 < numOfOutputAndInputVariables ){
//			std::cout << ",";
//		}
//	}
//	std::cout << "]" << std::endl;
//#endif

	const double val2q = static_cast<double>(degreeOfFreedom);
	for( int irow = 0; irow < numOfOutputAndInputVariables; ++irow ){
		// Partial derivatives regarding responses
		for( int icol = 0; icol < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++icol ){
			const double value1 = pow(determinant, -1.0/val2q) * dSigma1[irow][icol];
			double temp(0.0);
			for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
				temp += 2.0 / variances[iVar] * dSigma1[iVar][icol];
			}
			const double value2 = variancesWithoutScale[irow] / val2q * temp;
			const int iq2 = icol / ( 2 * numOfReferenceVariables );
			const int offset = iq2 * 2 * numOfReferenceVariables;
			const int amari = icol % ( 2 * numOfReferenceVariables );
			if( amari < numOfReferenceVariables ){
				// Real part
				const int ir2 = amari;
				const int index2 = 2 * ir2 + offset;
				derivativesRegardingResps[irow][index2] = value1 - value2;
			}else{
				// Imaginary part
				const int ir2 = amari - numOfReferenceVariables;
				const int index2 = 2 * ir2 + 1 + offset;
				derivativesRegardingResps[irow][index2] = value1 - value2;
			}
		}
		// Partial derivatives regarding variances without scale
		for( int icol = 0; icol < numOfOutputAndInputVariables; ++icol ){
			const double value1 = pow(determinant, -1.0/val2q) * dSigma2[irow][icol];
			double temp(0.0);
			for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
				temp += 2.0 / variances[iVar] * dSigma2[iVar][icol];
			}
			const double value2 = variancesWithoutScale[irow] / val2q * temp;
			derivativesRegardingVariancesWithoutScale[irow][icol] = value1 - value2;
		}
		// Partial derivatives regarding scale
		const double value1 = pow(determinant, -1.0/val2q) * dSigma3[irow];
		double temp(0.0);
		for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
			temp += 2.0 / variances[iVar] * dSigma3[iVar];
		}
		const double value2 = variancesWithoutScale[irow] / val2q * temp;
		derivativesRegardingScale[irow] = value1 - value2;
	}
	for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
		delete [] dSigma1[iVar];
	}
	delete [] dSigma1;
	for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
		delete [] dSigma2[iVar];
	}
	delete [] dSigma2;
	delete [] dSigma3;
	delete [] variances;

#ifdef _DEBUG_WRITE
	std::cout << "dvr1" << std::endl;
	std::cout << "[";
	for( int row = 0; row < numOfOutputAndInputVariables; ++row ){
		for( int col = 0; col < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++col ){
			std::cout << derivativesRegardingResps[row][col] <<" ";
		}
		if( row+1 < numOfOutputAndInputVariables ){
			std::cout << ";";
		}
	}
	std::cout << "]" << std::endl;
	std::cout << "dvv1" << std::endl;
	std::cout << "[";
	for( int row = 0; row < numOfOutputAndInputVariables; ++row ){
		for( int col = 0; col < numOfOutputAndInputVariables; ++col ){
			std::cout << derivativesRegardingVariancesWithoutScale[row][col] <<" ";
		}
		if( row+1 < numOfOutputAndInputVariables ){
			std::cout << ";";
		}
	}
	std::cout << "]" << std::endl;
	std::cout << "dvs1" << std::endl;
	std::cout << "[";
	for( int row = 0; row < numOfOutputAndInputVariables; ++row ){
		std::cout << derivativesRegardingScale[row];
		if( row+1 < numOfOutputAndInputVariables ){
			std::cout << ",";
		}
	}
	std::cout << "]" << std::endl;
#endif
#ifdef _DEBUG_WRITE
	double** derivVar1 = new double*[numOfOutputAndInputVariables];
	for( int irow = 0; irow < numOfOutputAndInputVariables; ++irow ){
		derivVar1[irow] = new double[2 * numOfReferenceVariables * numOfOutputAndInputVariables];
	}
	for( int irow = 0; irow < numOfOutputAndInputVariables; ++irow ){
		for( int icol = 0; icol < numOfReferenceVariables; ++icol ){
			for( int isign = 0; isign < 2; ++isign ){
				std::complex<double>** complexResidualsOrg = new std::complex<double>*[numOfOutputAndInputVariables];
				for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
					complexResidualsOrg[iVar] = new std::complex<double>[numSegments]; 
				}
				double* MDOrg = new double[numSegments];
				double* weightsOrg = new double[numSegments];
				std::complex<double>** respOrg = new std::complex<double>*[numOfOutputAndInputVariables];
				for( int i = 0; i < numOfOutputAndInputVariables; ++i ){
					respOrg[i] = new std::complex<double>[numOfReferenceVariables];
				}
				double* variancesWithoutScaleOrg = new double[numOfOutputAndInputVariables];
				std::complex<double>** complexResidualsMod = new std::complex<double>*[numOfOutputAndInputVariables];
				for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
					complexResidualsMod[iVar] = new std::complex<double>[numSegments]; 
				}
				double* MDMod = new double[numSegments];
				double* weightsMod = new double[numSegments];
				std::complex<double>** respMod = new std::complex<double>*[numOfOutputAndInputVariables];
				for( int i = 0; i < numOfOutputAndInputVariables; ++i ){
					respMod[i] = new std::complex<double>[numOfReferenceVariables];
				}
				double* variancesWithoutScaleMod = new double[numOfOutputAndInputVariables];

				for( int irow2 = 0; irow2 < numOfOutputAndInputVariables; ++irow2 ){
					for( int icol2 = 0; icol2 < numOfReferenceVariables; ++icol2 ){
						respOrg[irow2][icol2] = resp[irow2][icol2];
						respMod[irow2][icol2] = resp[irow2][icol2];
					}
				}
				double dresp(0.0);
				if( isign == 0 ){
					dresp = resp[irow][icol].real() * 0.001;
					respMod[irow][icol] += std::complex<double>( dresp, 0.0 );
				}else{
					dresp = resp[irow][icol].imag() * 0.001;
					respMod[irow][icol] += std::complex<double>( 0.0, dresp );
				}

				// Calculate complex residuals
				calculateComplexResiduals( numSegments, ftval, respOrg, complexResidualsOrg );
				calculateComplexResiduals( numSegments, ftval, respMod, complexResidualsMod );
				// Calculate Mahalanobis distance
				calculateMD( numSegments, numOfOutputAndInputVariables, complexResidualsOrg,
					variancesWithoutScale, MDOrg );
				calculateMD( numSegments, numOfOutputAndInputVariables, complexResidualsMod,
					variancesWithoutScale, MDMod );

				// Calculate original weights
				for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
					const double val = MDOrg[iSeg] / scale;
					weightsOrg[iSeg] = RobustWeightTukeysBiweights::calculateWeights(val, paramC);
				}
				for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
					const double val = MDMod[iSeg] / scale;
					weightsMod[iSeg] = RobustWeightTukeysBiweights::calculateWeights(val, paramC);
				}
				{
					double determinantOrg(1.0);
					for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
						// Make real residual vector from complex residual vector
						double numerator (0.0);
						double denominator(0.0);
						for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
							const double val = MDOrg[iSeg] / scale;
							numerator += std::norm(complexResidualsOrg[iVar][iSeg]) * weightsOrg[iSeg];
							denominator += pow(val, 2) * weightsOrg[iSeg];
						}
						variancesWithoutScaleOrg[iVar] = valq * numerator / denominator;
						determinantOrg *= variancesWithoutScaleOrg[iVar];// Real part
						determinantOrg *= variancesWithoutScaleOrg[iVar];// Imaginary part
					}
					const double factor = pow(determinantOrg, -1.0/val2q);
					for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
						// Make determinant of covariance matrix one
						variancesWithoutScaleOrg[iVar] *= factor;
					}
				}
				{
					double determinantMod(1.0);
					for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
						// Make real residual vector from complex residual vector
						double numerator (0.0);
						double denominator(0.0);
						for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
							const double val = MDMod[iSeg] / scale;
							numerator += std::norm(complexResidualsMod[iVar][iSeg]) * weightsMod[iSeg];
							denominator += pow(val, 2) * weightsMod[iSeg];
						}
						variancesWithoutScaleMod[iVar] = valq * numerator / denominator;
						determinantMod *= variancesWithoutScaleMod[iVar];// Real part
						determinantMod *= variancesWithoutScaleMod[iVar];// Imaginary part
					}
					const double factor = pow(determinantMod, -1.0/val2q);
					for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
						// Make determinant of covariance matrix one
						variancesWithoutScaleMod[iVar] *= factor;
					}
				}

				const int index = 2 * icol + isign + irow * 2 * numOfReferenceVariables;

				for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
					derivVar1[iVar][index] = ( variancesWithoutScaleMod[iVar] - variancesWithoutScaleOrg[iVar] ) / dresp;
				}
				
				for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
					delete [] complexResidualsOrg[iVar]; 
				}
				delete [] complexResidualsOrg;
				delete [] MDOrg;
				delete [] weightsOrg;
				for( int i = 0; i < numOfOutputAndInputVariables; ++i ){
					delete [] respOrg[i];
				}
				delete [] respOrg;
				delete [] variancesWithoutScaleOrg;
				for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
					delete [] complexResidualsMod[iVar]; 
				}
				delete [] complexResidualsMod;
				delete [] MDMod;
				delete [] weightsMod;
				for( int i = 0; i < numOfOutputAndInputVariables; ++i ){
					delete [] respMod[i];
				}
				delete [] respMod;
				delete [] variancesWithoutScaleMod;
			}
		}
	}
	std::cout << "dvr2" << std::endl;
	std::cout << "[";
	for( int row = 0; row < numOfOutputAndInputVariables; ++row ){
		for( int col = 0; col < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++col ){
			std::cout << derivVar1[row][col] <<" ";
		}
		if( row+1 < numOfOutputAndInputVariables ){
			std::cout << ";";
		}
	}
	std::cout << "]" << std::endl;
	for( int irow = 0; irow < numOfOutputAndInputVariables; ++irow ){
		delete [] derivVar1[irow];
	}
	delete [] derivVar1;

	double** derivVar2 = new double*[numOfOutputAndInputVariables];
	for( int irow = 0; irow < numOfOutputAndInputVariables; ++irow ){
		derivVar2[irow] = new double[numOfOutputAndInputVariables];
	}
	for( int icol = 0; icol < numOfOutputAndInputVariables; ++icol ){
		std::complex<double>** complexResidualsOrg = new std::complex<double>*[numOfOutputAndInputVariables];
		for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
			complexResidualsOrg[iVar] = new std::complex<double>[numSegments]; 
		}
		double* MDOrg = new double[numSegments];
		double* weightsOrg = new double[numSegments];
		double* variancesWithoutScaleOrg = new double[numOfOutputAndInputVariables];
		double* MDMod = new double[numSegments];
		double* weightsMod = new double[numSegments];
		double* variancesWithoutScaleMod = new double[numOfOutputAndInputVariables];
		for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
			variancesWithoutScaleOrg[iVar] = variancesWithoutScale[iVar];
			variancesWithoutScaleMod[iVar] = variancesWithoutScale[iVar];
		}
		const double dvar = variancesWithoutScale[icol] * 0.001;
		variancesWithoutScaleMod[icol] += dvar;
		// Calculate complex residuals
		calculateComplexResiduals( numSegments, ftval, resp, complexResidualsOrg );
		// Calculate Mahalanobis distance
		calculateMD( numSegments, numOfOutputAndInputVariables, complexResidualsOrg,
			variancesWithoutScaleOrg, MDOrg );
		calculateMD( numSegments, numOfOutputAndInputVariables, complexResidualsOrg,
			variancesWithoutScaleMod, MDMod );

		// Calculate original weights
		for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
			const double val = MDOrg[iSeg] / scale;
			weightsOrg[iSeg] = RobustWeightTukeysBiweights::calculateWeights(val, paramC);
		}
		for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
			const double val = MDMod[iSeg] / scale;
			weightsMod[iSeg] = RobustWeightTukeysBiweights::calculateWeights(val, paramC);
		}

		{
			double determinantOrg(1.0);
			for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
				// Make real residual vector from complex residual vector
				double numerator (0.0);
				double denominator(0.0);
				for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
					const double val = MDOrg[iSeg] / scale;
					numerator += std::norm(complexResidualsOrg[iVar][iSeg]) * weightsOrg[iSeg];
					denominator += pow(val, 2) * weightsOrg[iSeg];
				}
				variancesWithoutScaleOrg[iVar] = valq * numerator / denominator;
				determinantOrg *= variancesWithoutScaleOrg[iVar];// Real part
				determinantOrg *= variancesWithoutScaleOrg[iVar];// Imaginary part
			}
			const double factor = pow(determinantOrg, -1.0/val2q);
			for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
				// Make determinant of covariance matrix one
				variancesWithoutScaleOrg[iVar] *= factor;
			}
		}
		{
			double determinantMod(1.0);
			for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
				// Make real residual vector from complex residual vector
				double numerator (0.0);
				double denominator(0.0);
				for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
					const double val = MDMod[iSeg] / scale;
					numerator += std::norm(complexResidualsOrg[iVar][iSeg]) * weightsMod[iSeg];
					denominator += pow(val, 2) * weightsMod[iSeg];
				}
				variancesWithoutScaleMod[iVar] = valq * numerator / denominator;
				determinantMod *= variancesWithoutScaleMod[iVar];// Real part
				determinantMod *= variancesWithoutScaleMod[iVar];// Imaginary part
			}
			const double factor = pow(determinantMod, -1.0/val2q);
			for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
				// Make determinant of covariance matrix one
				variancesWithoutScaleMod[iVar] *= factor;
			}
		}

		for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
			derivVar2[iVar][icol] = ( variancesWithoutScaleMod[iVar] - variancesWithoutScaleOrg[iVar] ) / dvar;
		}

		for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
			delete [] complexResidualsOrg[iVar]; 
		}
		delete [] complexResidualsOrg;
		delete [] MDOrg;
		delete [] weightsOrg;
		delete [] variancesWithoutScaleOrg;
		delete [] MDMod;
		delete [] weightsMod;
		delete [] variancesWithoutScaleMod;
	}
	std::cout << "dvv2" << std::endl;
	std::cout << "[";
	for( int row = 0; row < numOfOutputAndInputVariables; ++row ){
		for( int col = 0; col < numOfOutputAndInputVariables; ++col ){
			std::cout << derivVar2[row][col] <<" ";
		}
		if( row+1 < numOfOutputAndInputVariables ){
			std::cout << ";";
		}
	}
	std::cout << "]" << std::endl;
	for( int row = 0; row < numOfOutputAndInputVariables; ++row ){
		delete [] derivVar2[row];
	}
	delete [] derivVar2;

	double* derivVar3 = new double[numOfOutputAndInputVariables];
	{
		std::complex<double>** complexResidualsOrg = new std::complex<double>*[numOfOutputAndInputVariables];
		for( int iVar2 = 0; iVar2 < numOfOutputAndInputVariables; ++iVar2 ){
			complexResidualsOrg[iVar2] = new std::complex<double>[numSegments]; 
		}
		double* MDOrg = new double[numSegments];
		double* weightsOrg = new double[numSegments];
		double* variancesWithoutScaleOrg = new double[numOfOutputAndInputVariables];
		double* weightsMod = new double[numSegments];
		double* variancesWithoutScaleMod = new double[numOfOutputAndInputVariables];
		const double scaleOrg = scale;
		const double dscale = scale * 0.001;
		const double scaleMod = scale + dscale;
		// Calculate complex residuals
		calculateComplexResiduals( numSegments, ftval, resp, complexResidualsOrg );
		// Calculate Mahalanobis distance
		calculateMD( numSegments, numOfOutputAndInputVariables, complexResidualsOrg,
			variancesWithoutScale, MDOrg );
		// Calculate original weights
		for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
			const double val = MDOrg[iSeg] / scaleOrg;
			weightsOrg[iSeg] = RobustWeightTukeysBiweights::calculateWeights(val, paramC);
		}
		for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
			const double val = MDOrg[iSeg] / scaleMod;
			weightsMod[iSeg] = RobustWeightTukeysBiweights::calculateWeights(val, paramC);
		}

		{
			double determinantOrg(1.0);
			for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
				// Make real residual vector from complex residual vector
				double numerator (0.0);
				double denominator(0.0);
				for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
					numerator += std::norm(complexResidualsOrg[iVar][iSeg]) * weightsOrg[iSeg];
					const double val = MDOrg[iSeg] / scaleOrg;
					denominator += pow(val, 2) * weightsOrg[iSeg];
				}
				variancesWithoutScaleOrg[iVar] = valq * numerator / denominator;
				determinantOrg *= variancesWithoutScaleOrg[iVar];// Real part
				determinantOrg *= variancesWithoutScaleOrg[iVar];// Imaginary part
			}
			const double factor = pow(determinantOrg, -1.0/val2q);
			for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
				// Make determinant of covariance matrix one
				variancesWithoutScaleOrg[iVar] *= factor;
			}
		}
		{
			double determinantMod(1.0);
			for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
				// Make real residual vector from complex residual vector
				double numerator (0.0);
				double denominator(0.0);
				for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
					numerator += std::norm(complexResidualsOrg[iVar][iSeg]) * weightsMod[iSeg];
					const double val = MDOrg[iSeg] / scaleMod;
					denominator += pow(val, 2) * weightsMod[iSeg];
				}
				variancesWithoutScaleMod[iVar] = valq * numerator / denominator;
				determinantMod *= variancesWithoutScaleMod[iVar];// Real part
				determinantMod *= variancesWithoutScaleMod[iVar];// Imaginary part
			}
			const double factor = pow(determinantMod, -1.0/val2q);
			for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
				// Make determinant of covariance matrix one
				variancesWithoutScaleMod[iVar] *= factor;
			}
		}

		for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
			derivVar3[iVar] = ( variancesWithoutScaleMod[iVar] - variancesWithoutScaleOrg[iVar] ) / dscale;
		}
		for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
			delete [] complexResidualsOrg[iVar]; 
		}
		delete [] complexResidualsOrg;
		delete [] MDOrg;
		delete [] weightsOrg;
		delete [] variancesWithoutScaleOrg;
		delete [] weightsMod;
		delete [] variancesWithoutScaleMod;
	}
	std::cout << "dvs2" << std::endl;
	std::cout << "[";
	for( int row = 0; row < numOfOutputAndInputVariables; ++row ){
		std::cout << derivVar3[row] <<" ";
		if( row+1 < numOfOutputAndInputVariables ){
			std::cout << ",";
		}
	}
	std::cout << "]" << std::endl;
	delete [] derivVar3;
#endif

}

// Calculate response functions
void AnalysisMultivariateRegression::calculateResponseFunctions( const int iSegLen,const int freqDegree, const double timeLength, const double freq, 
	const int numSegmentsTotal, std::complex<double>** ftval, const std::vector< std::pair<std::string, std::string> >& times,
	std::ofstream& ofsResp, std::ofstream& ofsRhoaPhs ){

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeLogMessage("Calculate response functions by multivariate regression");
	ptrOutputFiles->writeCvgMessage("================================================================================");
	ptrOutputFiles->writeCvgMessage("Now Frequency(Hz): " + Util::toString(freq) + ", Period(s): " + Util::toString(1.0/freq));
	ptrOutputFiles->writeCvgMessage("================================================================================");

	const Control* const ptrControl = Control::getInstance();
	const int numOfOutputVariables = ptrControl->getNumOutputVariables();
	const int numOfInputVariables = ptrControl->getNumInputVariables();
	const int numOfOutputAndInputVariables = numOfOutputVariables + numOfInputVariables;
	const int degreeOfFreedom = 2 * numOfOutputAndInputVariables;
	const int numOfReferenceVariables = ptrControl->getNumRemoteReferenceVariables();
	assert(numOfReferenceVariables == 2); 

	// Determine candidates
	std::vector< std::pair<int,int> > candidatesOfPairs;
	std::vector< std::complex<double>** > resp;
	determineCandidates( freq, numSegmentsTotal, ftval, times, candidatesOfPairs, resp );
	const int numOfCancidates = static_cast<int>( resp.size() );
	ptrOutputFiles->writeCvgAndLogMessage("Number of candidates: " + Util::toString(numOfCancidates));
	if( ptrControl->getOutputLevel() >= 3 ){
		ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
		ptrOutputFiles->writeCvgMessage("Initial response functions of all candidates:");
		ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
		// Output candidates
		for( int iCan = 0; iCan < numOfCancidates; ++iCan ){
			std::ostringstream msg;
			msg << "Candidate: " << std::setw(10) << iCan << ", ";
			if( m_frequencies.empty() ){
				const int iSeg1 = candidatesOfPairs[iCan].first;
				const int iSeg2 = candidatesOfPairs[iCan].second;
				msg << "Segment pair: ("  << std::setw(10) << iSeg1 << "," << std::setw(10) << iSeg2 << "), ";
			}
			msg << "Responses: ";
			int iVar = 0;
			for( int iOut = 0; iOut < numOfOutputVariables; ++iOut ){
				msg << "("  << std::setw(12) << std::setprecision(4) << std::scientific << resp[iCan][iVar][0].real() << "," 
							<< std::setw(12) << std::setprecision(4) << std::scientific << resp[iCan][iVar][0].imag() << "), ";
				msg << "("  << std::setw(12) << std::setprecision(4) << std::scientific << resp[iCan][iVar][1].real() << "," 
							<< std::setw(12) << std::setprecision(4) << std::scientific << resp[iCan][iVar][1].imag() << "); ";
				++iVar;
			}
			for( int iInp = 0; iInp < numOfInputVariables; ++iInp ){
				msg << "("  << std::setw(12) << std::setprecision(4) << std::scientific << resp[iCan][iVar][0].real() << "," 
							<< std::setw(12) << std::setprecision(4) << std::scientific << resp[iCan][iVar][0].imag() << "), ";
				msg << "("  << std::setw(12) << std::setprecision(4) << std::scientific << resp[iCan][iVar][1].real() << "," 
							<< std::setw(12) << std::setprecision(4) << std::scientific << resp[iCan][iVar][1].imag() << ")";
				++iVar;
				if( iInp + 1 < numOfInputVariables ){
					msg << "; ";
				}
			}
			assert(iVar == numOfOutputAndInputVariables);
			ptrOutputFiles->writeCvgMessage(msg.str());
		}
	}

	double** variancesWithoutScale = new double*[numOfCancidates];
	for( int iCan = 0; iCan < numOfCancidates; ++iCan ){
		variancesWithoutScale[iCan] = new double[numOfOutputAndInputVariables];
	}
	double* scales = new double[numOfCancidates];
	double* determinants = new double[numOfCancidates];

	// Perform I-steps for all candidates
	ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
	ptrOutputFiles->writeCvgAndLogMessage("Improve all candidates");
	double paramB(0.0);
	double paramC(0.0);
	RobustWeightTukeysBiweights::calculateParams(degreeOfFreedom, numSegmentsTotal, paramB, paramC);
	ptrOutputFiles->writeCvgMessage("Parameter b: " + Util::toString(paramB) + ", Parameter c: " + Util::toString(paramC));

	const Control::ParamsForRobustMultivariateRegression params = ptrControl->getParamsForRobustMultivariateRegression();
	double* coherences = new double[numOfOutputVariables];
	for( int iCan = 0; iCan < numOfCancidates; ++iCan ){
		const int numOfMaxIterations = params.numOfMaxIterationsOfFirstIstep;
		const double convergenceCriteria = params.convergenceCriteriaOfFirstIstep;
		if( ptrControl->getOutputLevel() >= 3 ){
			ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
			ptrOutputFiles->writeCvgMessage("Candidate: " + Util::toString(iCan));
			ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
		}
		improveCandidate( numSegmentsTotal, numOfMaxIterations, convergenceCriteria, true, paramB, paramC,
			ftval, resp[iCan], variancesWithoutScale[iCan], scales[iCan], determinants[iCan], coherences );
	}

	if( ptrControl->getOutputLevel() >= 3 ){
		ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
		ptrOutputFiles->writeCvgMessage("Improved response functions of all candidates:");
		ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
		// Output resultant response functions
		for( int iCan = 0; iCan < numOfCancidates; ++iCan ){
			std::ostringstream msg;
			msg << "Candidate: " << std::setw(10) << iCan << ", "; 
			if( m_frequencies.empty() ){
				const int iSeg1 = candidatesOfPairs[iCan].first;
				const int iSeg2 = candidatesOfPairs[iCan].second;
				msg << "Segment pair: ("  << std::setw(10) << iSeg1 << "," << std::setw(10) << iSeg2 << "), ";
			}
			msg << "Responses: ";
			int iVar(0);
			for( int iOut = 0; iOut < numOfOutputVariables; ++iOut ){
				msg << "("  << std::setw(12) << std::setprecision(4) << std::scientific << resp[iCan][iVar][0].real() << "," 
							<< std::setw(12) << std::setprecision(4) << std::scientific << resp[iCan][iVar][0].imag() << "), ";
				msg << "("  << std::setw(12) << std::setprecision(4) << std::scientific << resp[iCan][iVar][1].real() << "," 
							<< std::setw(12) << std::setprecision(4) << std::scientific << resp[iCan][iVar][1].imag() << "); ";
				++iVar;
			}
			for( int iInp = 0; iInp < numOfInputVariables; ++iInp ){
				msg << "("  << std::setw(12) << std::setprecision(4) << std::scientific << resp[iCan][iVar][0].real() << "," 
							<< std::setw(12) << std::setprecision(4) << std::scientific << resp[iCan][iVar][0].imag() << "), ";
				msg << "("  << std::setw(12) << std::setprecision(4) << std::scientific << resp[iCan][iVar][1].real() << "," 
							<< std::setw(12) << std::setprecision(4) << std::scientific << resp[iCan][iVar][1].imag() << ")";
				++iVar;
				if( iInp + 1 < numOfInputVariables ){
					msg << "; ";
				}
			}
			assert(iVar == numOfOutputAndInputVariables);
			ptrOutputFiles->writeCvgMessage(msg.str());
		}

		// Output scales and variances without scale
		ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
		ptrOutputFiles->writeCvgMessage("Scales and variances without scale of all candidates:");
		ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
		for( int iCan = 0; iCan < numOfCancidates; ++iCan ){
			std::ostringstream msg;
			msg << "Candidate: " << std::setw(10) << iCan << ", "; 
			if( m_frequencies.empty() ){
				const int iSeg1 = candidatesOfPairs[iCan].first;
				const int iSeg2 = candidatesOfPairs[iCan].second;
				msg << "Segment pair: ("  << std::setw(10) << iSeg1 << "," << std::setw(10) << iSeg2 << "), ";
			}
			msg << "Scale:" << std::setw(12) << std::setprecision(4) << std::scientific << scales[iCan] << ", Variances without scale:";
			for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
				msg << std::setw(12) << std::setprecision(4) << std::scientific << variancesWithoutScale[iCan][iVar];
				if( iVar + 1 < numOfOutputAndInputVariables ){
					msg << ", ";
				}
			}
			ptrOutputFiles->writeCvgMessage(msg.str());
		}
	}

	// Detetermine the best improved candidates, to which further I-step is applied
	ptrOutputFiles->writeLogMessage("Detetermine the best improved candidates");
	const int rr0 = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, 0 );
	const int rr1 = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, 1 );
	std::complex<double>** complexResiduals = new std::complex<double>*[numOfOutputAndInputVariables];
	for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
		complexResiduals[iVar] = new std::complex<double>[numSegmentsTotal]; 
	}
	double* MD = new double[numSegmentsTotal];

	int numOfBestCandidates = params.numOfMaxCandidatesOfSecondIstep;
	std::vector<int> bestCandidates;
	bestCandidates.reserve(numOfBestCandidates);
	double maxScaleOfBestCandidates(-1.0);
	int candidateWithMaxScale(-1);
	for( int iCan = 0; iCan < numOfCancidates; ++iCan ){
		assert( bestCandidates.size() <= numOfBestCandidates );
		if( bestCandidates.size() >= numOfBestCandidates ){
			// Calculate complex residuals
			calculateComplexResiduals( numSegmentsTotal, ftval, resp[iCan], complexResiduals );
			// Calculate Mahalanobis distance
			calculateMD( numSegmentsTotal, numOfOutputAndInputVariables, complexResiduals, variancesWithoutScale[iCan], MD );
			double averageOfLossFunction(0.0);
			for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
				const double val = MD[iSeg] / maxScaleOfBestCandidates;
				averageOfLossFunction += RobustWeightTukeysBiweights::calculateLossFunction(val, paramC);
			}
			averageOfLossFunction /= static_cast<double>(numSegmentsTotal);
			if( averageOfLossFunction < paramB ){
				// Adopt this candidates as one of the best candidates
				std::vector<int>::iterator itrCan = std::find(bestCandidates.begin(), bestCandidates.end(), candidateWithMaxScale);
				assert( itrCan != bestCandidates.end() );
				*itrCan = iCan;
				maxScaleOfBestCandidates = -1.0;
				candidateWithMaxScale = -1;
				for( std::vector<int>::const_iterator itr = bestCandidates.begin(); itr != bestCandidates.end(); ++itr ){
					if( scales[*itr] > maxScaleOfBestCandidates ){
						maxScaleOfBestCandidates = scales[*itr];
						candidateWithMaxScale= *itr;
					}
				}
			}
		}else{
			bestCandidates.push_back(iCan);
			if( scales[iCan] > maxScaleOfBestCandidates ){
				maxScaleOfBestCandidates = scales[iCan];
				candidateWithMaxScale = iCan;
			}
		}
	}
	for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
		delete [] complexResiduals[iVar];
	}
	delete [] complexResiduals;
	delete [] MD;

	std::sort( bestCandidates.begin(), bestCandidates.end() );
	if( ptrControl->getOutputLevel() >= 3 ){
		ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
		ptrOutputFiles->writeCvgMessage("Response functions of the best " + Util::toString(bestCandidates.size()) + " candidates:");
		ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
		for( std::vector<int>::const_iterator itrCanBest = bestCandidates.begin(); itrCanBest != bestCandidates.end(); ++itrCanBest ){
			const int iCan = *itrCanBest;	
			std::ostringstream msg;
			msg << "Candidate: " << std::setw(10) << iCan << ", "; 
			if( m_frequencies.empty() ){
				const int iSeg1 = candidatesOfPairs[iCan].first;
				const int iSeg2 = candidatesOfPairs[iCan].second;
				msg << "Segment pair: ("  << std::setw(10) << iSeg1 << "," << std::setw(10) << iSeg2 << "), ";
			}
			msg << "Responses: ";
			int iVar = 0;
			for( int iOut = 0; iOut < numOfOutputVariables; ++iOut ){
				msg << "("  << std::setw(12) << std::setprecision(4) << std::scientific << resp[iCan][iVar][0].real() << "," 
							<< std::setw(12) << std::setprecision(4) << std::scientific << resp[iCan][iVar][0].imag() << "), ";
				msg << "("  << std::setw(12) << std::setprecision(4) << std::scientific << resp[iCan][iVar][1].real() << "," 
							<< std::setw(12) << std::setprecision(4) << std::scientific << resp[iCan][iVar][1].imag() << "); ";
				++iVar;
			}
			for( int iInp = 0; iInp < numOfInputVariables; ++iInp ){
				msg << "("  << std::setw(12) << std::setprecision(4) << std::scientific << resp[iCan][iVar][0].real() << "," 
							<< std::setw(12) << std::setprecision(4) << std::scientific << resp[iCan][iVar][0].imag() << "), ";
				msg << "("  << std::setw(12) << std::setprecision(4) << std::scientific << resp[iCan][iVar][1].real() << "," 
							<< std::setw(12) << std::setprecision(4) << std::scientific << resp[iCan][iVar][1].imag() << ")";
				++iVar;
				if( iInp + 1 < numOfInputVariables ){
					msg << "; ";
				}
			}
			assert(iVar == numOfOutputAndInputVariables);
			ptrOutputFiles->writeCvgMessage(msg.str());
		}
		ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
		ptrOutputFiles->writeCvgMessage("Scales and variances without scale of the best " + Util::toString(bestCandidates.size()) + " candidates:");
		ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
		for( std::vector<int>::const_iterator itrCanBest = bestCandidates.begin(); itrCanBest != bestCandidates.end(); ++itrCanBest ){
			const int iCan = *itrCanBest;	
			std::ostringstream msg;
			msg << "Candidate: " << std::setw(10) << iCan << ", "; 
			if( m_frequencies.empty() ){
				const int iSeg1 = candidatesOfPairs[iCan].first;
				const int iSeg2 = candidatesOfPairs[iCan].second;
				msg << "Segment pair: ("  << std::setw(10) << iSeg1 << "," << std::setw(10) << iSeg2 << "), ";
			}
			msg << "Scale:" << std::setw(12) << std::setprecision(4) << std::scientific << scales[iCan] << ", Variances without scale:";
			for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
				msg << std::setw(12) << std::setprecision(4) << std::scientific << variancesWithoutScale[iCan][iVar];
				if( iVar + 1 < numOfOutputAndInputVariables ){
					msg << ", ";
				}
			}
			ptrOutputFiles->writeCvgMessage(msg.str());
		}
	}

	// Apply further I-step to the best improved candidates
	ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
	ptrOutputFiles->writeCvgAndLogMessage("Perform further improvements to the best " + Util::toString(bestCandidates.size()) + " candidates");
	double smallestScale(1.0e20);
	int candidateWithSmallestScale(-1);
	for( std::vector<int>::const_iterator itrCanBest = bestCandidates.begin(); itrCanBest != bestCandidates.end(); ++itrCanBest ){
		// Improve candidate by I-step
		const int numOfMaxIterations = params.numOfMaxIterationsOfSecondIstep;
		const double convergenceCriteria = params.convergenceCriteriaOfSecondIstep;
		if( ptrControl->getOutputLevel() >= 3 ){
			ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
			ptrOutputFiles->writeCvgMessage("Candidate: " + Util::toString(*itrCanBest));
			ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
		}
		improveCandidate( numSegmentsTotal, numOfMaxIterations, convergenceCriteria, false, paramB, paramC,
			ftval, resp[*itrCanBest], variancesWithoutScale[*itrCanBest], scales[*itrCanBest], determinants[*itrCanBest], coherences );
		if( scales[*itrCanBest] < smallestScale) {
			smallestScale = scales[*itrCanBest];
			candidateWithSmallestScale= *itrCanBest;
		}
	}

	// Calculate complex residuals
	const bool outputResiduals = ptrControl->getOutputLevel() >= 2;
	if( outputResiduals ){
		const int iCan = candidateWithSmallestScale;
		// Calculate complex residuals
		std::complex<double>** complexResiduals = new std::complex<double>*[numOfOutputAndInputVariables];
		for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
			complexResiduals[iVar] = new std::complex<double>[numSegmentsTotal]; 
		}
		calculateComplexResiduals( numSegmentsTotal, ftval, resp[iCan], complexResiduals );
		// Calculate Mahalanobis distance
		double* MD = new double[numSegmentsTotal];
		calculateMD( numSegmentsTotal, numOfOutputAndInputVariables, complexResiduals,
			variancesWithoutScale[iCan], MD );
		double* weights = new double[numSegmentsTotal];
		for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
			const double val = MD[iSeg] / scales[iCan];
			weights[iSeg] = RobustWeightTukeysBiweights::calculateWeights(val, paramC);
		}
		std::ostringstream oss;
		oss << "segm" << iSegLen << "_index" << freqDegree << "_residuals.csv";
		writeResiduals( oss.str().c_str(), numSegmentsTotal, numOfOutputAndInputVariables, times, complexResiduals,
			MD, weights );
		for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
			delete [] complexResiduals[iVar];
		}
		delete [] complexResiduals;
		delete [] MD;
		delete [] weights;
	}

	// Select the response with the smallest scale
	ptrOutputFiles->writeLogMessage("Select the response with the smallest scale");
	{
		const int iCan = candidateWithSmallestScale;
		ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
		ptrOutputFiles->writeCvgMessage("Best estimator of frequency " + Util::toString(freq) + "(Hz)");
		if( ptrControl->getOutputLevel() >= 3 ){
			ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
			ptrOutputFiles->writeCvgMessage("Candidate: " + Util::toString(iCan));
			if( m_frequencies.empty() ){
				const int iSeg1 = candidatesOfPairs[iCan].first;
				const int iSeg2 = candidatesOfPairs[iCan].second;
				ptrOutputFiles->writeCvgMessage("Segment pair: (" + Util::toString(iSeg1) + "," + Util::toString(iSeg2) + ")");
			}
			ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
		}
	}
	{
		const int iCan = candidateWithSmallestScale;
		ptrOutputFiles->writeCvgMessage("Response functions:");
		std::ostringstream msg;
		int iVar = 0;
		for( int iOut = 0; iOut < numOfOutputVariables; ++iOut ){
			msg << "("  << std::setw(12) << std::setprecision(4) << std::scientific << resp[iCan][iVar][0].real() << "," 
					   << std::setw(12) << std::setprecision(4) << std::scientific << resp[iCan][iVar][0].imag() << "), ";
			msg << "("  << std::setw(12) << std::setprecision(4) << std::scientific << resp[iCan][iVar][1].real() << "," 
					   << std::setw(12) << std::setprecision(4) << std::scientific << resp[iCan][iVar][1].imag() << ")" << std::endl;
			++iVar;
		}
		for( int iInp = 0; iInp < numOfInputVariables; ++iInp ){
			msg << "("  << std::setw(12) << std::setprecision(4) << std::scientific << resp[iCan][iVar][0].real() << "," 
					   << std::setw(12) << std::setprecision(4) << std::scientific << resp[iCan][iVar][0].imag() << "), ";
			msg << "("  << std::setw(12) << std::setprecision(4) << std::scientific << resp[iCan][iVar][1].real() << "," 
					   << std::setw(12) << std::setprecision(4) << std::scientific << resp[iCan][iVar][1].imag() << ")";
			++iVar;
			if( iInp + 1 < numOfInputVariables ){
				msg << std::endl;
			}
		}
		assert(iVar == numOfOutputAndInputVariables);
		ptrOutputFiles->writeCvgMessage(msg.str());
	}
	{
		const int iCan = candidateWithSmallestScale;
		ptrOutputFiles->writeCvgMessage("Scale: " + Util::toString(scales[iCan]));
		std::ostringstream msg;
		msg << "Variances without scale:";
		for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
			msg << std::setw(12) << std::setprecision(4) << std::scientific << variancesWithoutScale[iCan][iVar];
			if( iVar + 1 < numOfOutputAndInputVariables ){
				msg << ", ";
			}
		}
		ptrOutputFiles->writeCvgMessage(msg.str());
	}

	if( ptrControl->getOutputLevel() > 0 ){
		// Output spectral density functions to cvg file
		std::complex<double>** complexResiduals = new std::complex<double>*[numOfOutputAndInputVariables];
		for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
			complexResiduals[iVar] = new std::complex<double>[numSegmentsTotal]; 
		}
		double* MD = new double[numSegmentsTotal];
		double* weights = new double[numSegmentsTotal];
		const int iCan = candidateWithSmallestScale;
		// Calculate complex residuals
		calculateComplexResiduals( numSegmentsTotal, ftval, resp[iCan], complexResiduals );
		// Calculate Mahalanobis distance
		calculateMD( numSegmentsTotal, numOfOutputAndInputVariables, complexResiduals,
			variancesWithoutScale[iCan], MD );
		// Calculate weights
		for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
			const double val = MD[iSeg] / scales[iCan];
			weights[iSeg] = RobustWeightTukeysBiweights::calculateWeights(val, paramC);
		}
		for( int iInVar = 0; iInVar < numOfInputVariables; ++iInVar ){
			ptrOutputFiles->writeCvgMessage("Spectral density functions for input variable " + Util::toString(iInVar));
			const int inp = ptrControl->getChannelIndex( CommonParameters::INPUT, iInVar );
			outputSpectralDensityFunctionsToCvgFile( numSegmentsTotal, timeLength, ftval[inp], ftval[rr0], ftval[rr1], weights );
		}
		for( int iOutVar = 0; iOutVar < numOfOutputVariables; ++iOutVar ){
			ptrOutputFiles->writeCvgMessage("Spectral density functions for output variable " + Util::toString(iOutVar));
			const int out = ptrControl->getChannelIndex( CommonParameters::OUTPUT, iOutVar );
			outputSpectralDensityFunctionsToCvgFile( numSegmentsTotal, timeLength, ftval[out], ftval[rr0], ftval[rr1], weights );
		}
		// Delete dynamic arrays
		for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
			delete [] complexResiduals[iVar];
		}
		delete [] complexResiduals;
		delete [] MD;
		delete [] weights;
	}

	// Output response functions
	std::complex<double>* resp0 = new std::complex<double>[numOfOutputVariables];
	std::complex<double>* resp1 = new std::complex<double>[numOfOutputVariables];
	assert(numOfInputVariables == 2);
	const int in0 = ptrControl->getChannelIndex( CommonParameters::INPUT, 0 );
	const int in1 = ptrControl->getChannelIndex( CommonParameters::INPUT, 1 );
	const std::complex<double> Txx = resp[candidateWithSmallestScale][in0][0];
	const std::complex<double> Txy = resp[candidateWithSmallestScale][in0][1];
	const std::complex<double> Tyx = resp[candidateWithSmallestScale][in1][0];
	const std::complex<double> Tyy = resp[candidateWithSmallestScale][in1][1];
	const std::complex<double> det = Txx * Tyy - Txy * Tyx;
	if( std::abs(det) < CommonParameters::EPS ){
		ptrOutputFiles->writeErrorMessage("Determinant is too small: " + Util::toString(std::abs(det)));
	}
	ofsResp << std::setprecision(10) << std::scientific << freq;
	ofsResp << "," << std::setprecision(10) << std::scientific << 1.0/freq;
	if( ptrControl->doesOutputApparentResistivityAndPhase() ){
		ofsRhoaPhs << std::setprecision(10) << std::scientific << freq;
		ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << 1.0/freq;
	}
	for( int iOut = 0; iOut < numOfOutputVariables; ++iOut ){
		const int index = ptrControl->getChannelIndex( CommonParameters::OUTPUT, iOut );
		const std::complex<double> U_x = resp[candidateWithSmallestScale][index][0];
		const std::complex<double> U_y = resp[candidateWithSmallestScale][index][1];
		resp0[iOut] = ( U_x * Tyy - U_y * Tyx ) / det;
		resp1[iOut] = ( U_y * Txx - U_x * Txy ) / det;
		ofsResp << "," << std::setprecision(10) << std::scientific << resp0[iOut].real();
		ofsResp << "," << std::setprecision(10) << std::scientific << resp0[iOut].imag();
		ofsResp << "," << std::setprecision(10) << std::scientific << resp1[iOut].real();
		ofsResp << "," << std::setprecision(10) << std::scientific << resp1[iOut].imag();
		ofsResp << "," << std::setprecision(10) << std::scientific << coherences[iOut];
		if( ptrControl->doesOutputApparentResistivityAndPhase() ){
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcApparentResistivity(freq, resp0[iOut]);
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcPhase(resp0[iOut]);
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcApparentResistivity(freq, resp1[iOut]);
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcPhase(resp1[iOut]);
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << coherences[iOut];
		}
	}
	delete [] coherences;

	// Estimate errors
	double** respErr = new double*[numOfOutputVariables];
	for( int iOut = 0; iOut < numOfOutputVariables; ++iOut ){
		respErr[iOut] = new double[numOfReferenceVariables];
	}
	const int typeOfErrorEstimationMethod = ptrControl->getErrorEstimationMethod();
	switch(typeOfErrorEstimationMethod){
	case Control::PARAMETRIC:
		estimateErrorParametric(numSegmentsTotal, paramB, paramC, ftval, resp[candidateWithSmallestScale], variancesWithoutScale[candidateWithSmallestScale],
			scales[candidateWithSmallestScale], respErr);
		break;
	case Control::FIXED_WEIGHTS_JACKKNIFE:
		estimateErrorByFixedWeightsJackknife( numSegmentsTotal, paramB, paramC, ftval, resp0, resp1, resp[candidateWithSmallestScale],
			variancesWithoutScale[candidateWithSmallestScale], scales[candidateWithSmallestScale], respErr );
		break;
	case Control::FIXED_WEIGHTS_BOOTSTRAP:
		estimateErrorByFixedWeightsBootstrap(numSegmentsTotal, paramB, paramC, ftval, resp[candidateWithSmallestScale], variancesWithoutScale[candidateWithSmallestScale],
			scales[candidateWithSmallestScale], determinants[candidateWithSmallestScale], respErr);
		break;
	case Control::SUBSET_DELETION_JACKKNIFE:
		estimateErrorBySubsetDeletionJackknife( numSegmentsTotal, ftval, resp[candidateWithSmallestScale], variancesWithoutScale[candidateWithSmallestScale],
			scales[candidateWithSmallestScale], determinants[candidateWithSmallestScale], respErr );
		break;
	case Control::STRICT_BOOTSTRAP:
		estimateErrorByStrictBootstrap( numSegmentsTotal, paramB, paramC, ftval, resp[candidateWithSmallestScale], variancesWithoutScale[candidateWithSmallestScale],
			scales[candidateWithSmallestScale], determinants[candidateWithSmallestScale], respErr );
		break;
	default:
		ptrOutputFiles->writeErrorMessage("Unsupported error estimation method : " + Util::toString(typeOfErrorEstimationMethod));
		break;
	}

	for( int iCan = 0; iCan < numOfCancidates; ++iCan ){
		delete [] variancesWithoutScale[iCan];
	}
	delete [] variancesWithoutScale;
	delete [] scales;
	delete [] determinants;

	// Output errors
	for( int iOut = 0; iOut < numOfOutputVariables; ++iOut ){
		const double dResp0 = respErr[iOut][0];
		const double dResp1 = respErr[iOut][1];
		ofsResp << "," << std::setprecision(10) << std::scientific << dResp0;		
		ofsResp << "," << std::setprecision(10) << std::scientific << dResp1;
		if( ptrControl->doesOutputApparentResistivityAndPhase() ){
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcApparentResistivityError(freq, resp0[iOut], dResp0);
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcPhaseError(resp0[iOut], dResp0);
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcApparentResistivityError(freq, resp1[iOut], dResp1);
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcPhaseError(resp1[iOut], dResp1);
		}
		delete [] respErr[iOut];
	}
	ofsResp << std::endl;
	ofsResp.flush();
	if( ptrControl->doesOutputApparentResistivityAndPhase() ){
		ofsRhoaPhs << std::endl;
		ofsRhoaPhs.flush();
	}
	delete [] resp0;
	delete [] resp1;
	delete [] respErr;

	// Add response functions to member variables
	m_numOfOutputAndInputVariables = numOfOutputAndInputVariables;
	if( m_responseFunctions == NULL ){
		m_responseFunctions = new std::vector< std::complex<double> >*[m_numOfOutputAndInputVariables];
		for( int iVar = 0; iVar < m_numOfOutputAndInputVariables; ++iVar ){
			m_responseFunctions[iVar] = new std::vector< std::complex<double> >[numOfReferenceVariables];
		}
	}
	m_frequencies.push_back(freq);
	for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
		for( int iRR = 0; iRR < numOfReferenceVariables; ++iRR ){
			m_responseFunctions[iVar][iRR].push_back( resp[candidateWithSmallestScale][iVar][iRR] );
		}
	}

	// Delete dynamic arrays
	for( int iCan = 0; iCan < numOfCancidates; ++iCan ){
		for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
			delete [] resp[iCan][iVar];
		}	
		delete [] resp[iCan];
	}

}

// Calculate response functions from two samples
bool AnalysisMultivariateRegression::calculateResponseFunctionsFromTwoSamples(
	const std::complex<double>& out1, const std::complex<double>& x1, const std::complex<double>& y1,
	const std::complex<double>& out2, const std::complex<double>& x2, const std::complex<double>& y2,
	std::complex<double>& resp1, std::complex<double>& resp2) const {

	const std::complex<double> det = x1 * y2 - x2 * y1;
	if (abs(det) < CommonParameters::EPS) {
		return true;
	}
	resp1 = (out1 * y2 - out2 * y1) / det;
	resp2 = (out2 * x1 - out1 * x2) / det;

	return false;
}

// Calculate a vector for partial derivatives
void AnalysisMultivariateRegression::calculateVectorForPartialDerivatives( const int iSeg, std::complex<double>** ftval, 
	const double* const variancesWithoutScale, std::complex<double>** complexResiduals, 
	std::complex<double>** hSigmaMatrix, double* vector ) const{

	const Control* const ptrControl = Control::getInstance();
	const int numOfOutputVariables = ptrControl->getNumOutputVariables();
	const int numOfInputVariables = ptrControl->getNumInputVariables();
	const int numOfOutputAndInputVariables = numOfOutputVariables + numOfInputVariables;
	const int numOfReferenceVariables = ptrControl->getNumRemoteReferenceVariables();

	for( int irow = 0; irow < numOfReferenceVariables; ++irow ){
		const int rr = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, irow );
		for( int icol = 0; icol < numOfOutputAndInputVariables; ++icol ){
			hSigmaMatrix[irow][icol] = std::conj(ftval[rr][iSeg]) * complexResiduals[icol][iSeg] / variancesWithoutScale[icol];
		}
	}
//#ifdef _DEBUG_WRITE
//	std::cout << "[";
//	for( int row = 0; row < numOfReferenceVariables; ++row ){
//		for( int col = 0; col < numOfOutputAndInputVariables; ++col ){
//			std::cout << hSigmaMatrix[row][col].real() << "+" << hSigmaMatrix[row][col].imag() <<"im ";
//		}
//		if( row+1 < numOfReferenceVariables ){
//			std::cout << ";";
//		}
//	}
//	std::cout << "]" << std::endl;
//#endif

	for( int icol = 0; icol < numOfOutputAndInputVariables; ++icol ){
		const int offset = 2 * numOfReferenceVariables * icol;
		for( int irow = 0; irow < numOfReferenceVariables; ++irow ){
			const int index = irow + offset;
			vector[index] = hSigmaMatrix[irow][icol].real();
		}
		for( int irow = 0; irow < numOfReferenceVariables; ++irow ){
			const int index = irow + numOfReferenceVariables + offset;
			vector[index] = hSigmaMatrix[irow][icol].imag();
		}
	}

//#ifdef _DEBUG_WRITE
//	std::cout << "[";
//	for( int row = 0; row < 2 * numOfReferenceVariables * numOfOutputAndInputVariables; ++row ){
//		std::cout << vector[row] <<" ";
//		if( row+1 < 2 * numOfReferenceVariables * numOfOutputAndInputVariables ){
//			std::cout << ",";
//		}
//	}
//	std::cout << "]" << std::endl;
//#endif

}

// Determine candidates
void AnalysisMultivariateRegression::determineCandidates( const double freq, const int numSegmentsTotal, std::complex<double>** ftval, 
	const std::vector< std::pair<std::string, std::string> >& times, std::vector< std::pair<int,int> >& candidatesOfPairs, 
	std::vector< std::complex<double>** >& resp ) const{

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeLogMessage("Determine candicates");

	const Control* const ptrControl = Control::getInstance();
	const Control::ParamsForRobustMultivariateRegression paramsForRobustMultivariateRegression = ptrControl->getParamsForRobustMultivariateRegression();

	if( m_frequencies.empty() || paramsForRobustMultivariateRegression.selectInitialCandidatesByRandomSamplingAtEachFrequency ){
		determineCandidatesByRandomSampling( freq, numSegmentsTotal, ftval, times, candidatesOfPairs, resp );
	}
	else if( static_cast<int>(m_frequencies.size()) == 1 ){
		candidatesOfPairs.push_back( std::make_pair(-1,-1) );
		const int numOfOutputVariables = ptrControl->getNumOutputVariables();
		const int numOfInputVariables = ptrControl->getNumInputVariables();
		const int numOfOutputAndInputVariables = numOfOutputVariables + numOfInputVariables;
		const int numOfReferenceVariables = ptrControl->getNumRemoteReferenceVariables();
		// This array must not be deleted because they are used out of the function
		std::complex<double>** respWork = new std::complex<double>*[numOfOutputAndInputVariables];
		for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
			respWork[iVar] = new std::complex<double>[numOfReferenceVariables];
			for( int iRR = 0; iRR < numOfReferenceVariables; ++iRR ){
				respWork[iVar][iRR] = m_responseFunctions[iVar][iRR].back();
			}
		}
		resp.push_back(respWork);
	}
	else{
		const Control::ParamsForDecidingCandicatesForSubsequentFrequencies params = ptrControl->getParamsForDecidingCandicatesForSubsequentFrequencies();
		if(params.useResponseFunctionsOfPreviousFrequency){
			candidatesOfPairs.push_back( std::make_pair(-1,-1) );
			const int numOfOutputVariables = ptrControl->getNumOutputVariables();
			const int numOfInputVariables = ptrControl->getNumInputVariables();
			const int numOfOutputAndInputVariables = numOfOutputVariables + numOfInputVariables;
			const int numOfReferenceVariables = ptrControl->getNumRemoteReferenceVariables();
			// This array must not be deleted because they are used out of the function
			std::complex<double>** respWork = new std::complex<double>*[numOfOutputAndInputVariables];
			for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
				respWork[iVar] = new std::complex<double>[numOfReferenceVariables];
				for( int iRR = 0; iRR < numOfReferenceVariables; ++iRR ){
					respWork[iVar][iRR] = m_responseFunctions[iVar][iRR].back();
				}
			}
			resp.push_back(respWork);
		}else{
			determineCandidatesForLaterFrequencies( freq, numSegmentsTotal, resp );
		}
	}

}

// Determine candidates for later frequencies
void AnalysisMultivariateRegression::determineCandidatesForLaterFrequencies( const double freq, const int numSegmentsTotal, std::vector< std::complex<double>** >& resp ) const{

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	const Control* const ptrControl = Control::getInstance();
	const int numOfOutputVariables = ptrControl->getNumOutputVariables();
	const int numOfInputVariables = ptrControl->getNumInputVariables();
	const int numOfOutputAndInputVariables = numOfOutputVariables + numOfInputVariables;
	const int numOfReferenceVariables = ptrControl->getNumRemoteReferenceVariables();

	std::vector<double>::const_reverse_iterator ritrFreq = m_frequencies.rbegin();
	const double freqPre = *ritrFreq;
	++ritrFreq;
	const double freqPrePre = *ritrFreq;

	const Control::ParamsForDecidingCandicatesForSubsequentFrequencies params = ptrControl->getParamsForDecidingCandicatesForSubsequentFrequencies();
	std::vector<double> denominators;
	denominators.reserve(numOfOutputAndInputVariables);
	for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
		std::complex<double> respPre[2];
		std::complex<double> respPrePre[2];
		assert( numOfReferenceVariables == 2 );
		const int typeOfThreshould = params.thresholdOfDifferences[iVar].type;
		for( int iRR = 0; iRR < numOfReferenceVariables; ++iRR ){
			// Divide the impedance tensor by frequency to aviod the dependancy of impedance tensor on the frequency
			std::vector< std::complex<double> >::const_reverse_iterator citrResp = m_responseFunctions[iVar][iRR].rbegin();
			respPre[iRR] = *citrResp;
			++citrResp;
			respPrePre[iRR] = *citrResp;
			if( typeOfThreshould == Control::DIFFERENCE_OF_RESPONSES_DIVIDED_BY_SQUARE_ROOT_FREQUENCY ){
				respPre[iRR] /= sqrt(freqPre);
				respPrePre[iRR] /= sqrt(freqPrePre);
			}
		}
		const double absRespPre = sqrt( 0.5 * ( std::norm(respPre[0]) + std::norm(respPre[1]) ) );
		const double absRespPrePre = sqrt( 0.5 * ( std::norm(respPrePre[0]) + std::norm(respPrePre[1]) ) );
		denominators.push_back(std::min( absRespPre, absRespPrePre));
	}

	std::vector< std::complex<double> >** candidates = new std::vector< std::complex<double> >*[numOfOutputAndInputVariables];
	int numOfCandidates(1);
	for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
		candidates[iVar] = new std::vector< std::complex<double> >[numOfReferenceVariables];
		for( int iRR = 0; iRR < numOfReferenceVariables; ++iRR ){
			std::vector< std::complex<double> >::const_reverse_iterator citrResp = m_responseFunctions[iVar][iRR].rbegin();
			const std::complex<double> respPre = *citrResp;
			candidates[iVar][iRR].push_back(respPre);
			++citrResp;
			const std::complex<double> respPrePre = *citrResp;
			const int typeOfThreshould = params.thresholdOfDifferences[iVar].type;
			const double threshold = params.thresholdOfDifferences[iVar].threshold;
			bool largeDifference(false);
			if( typeOfThreshould == Control::DIFFERENCE_OF_RAW_RESPONSES ){
				if( fabs(respPre.real() - respPrePre.real()) > threshold * denominators[iVar] ||
					fabs(respPre.imag() - respPrePre.imag()) > threshold * denominators[iVar] ){
					largeDifference = true;
				}
			}else if( typeOfThreshould == Control::DIFFERENCE_OF_RESPONSES_DIVIDED_BY_SQUARE_ROOT_FREQUENCY ){
				const std::complex<double> respModPre = respPre / sqrt(freqPre);
				const std::complex<double> respModPrePre = respPrePre / sqrt(freqPrePre);
				if( fabs(respModPre.real() - respModPrePre.real()) > threshold * denominators[iVar] || 
					fabs(respModPre.imag() - respModPrePre.imag()) > threshold * denominators[iVar] ){
					largeDifference = true;
				}
			}else{
				ptrOutputFiles->writeErrorMessage("Unsupported type of response function: " +Util::toString(typeOfThreshould));
			}
			if( largeDifference ){
				candidates[iVar][iRR].push_back(respPrePre);
				numOfCandidates *= 2;
				ptrOutputFiles->writeCvgAndLogMessage("The differene of response function (" +
					Util::toString(iVar) + ", " + Util::toString(iRR) + ") is large");
				std::ostringstream msg;
				msg << "     Real part: " 
					<< std::setw(12) << std::setprecision(4) << std::scientific << respPrePre.real()
					<< "->" 
					<< std::setw(12) << std::setprecision(4) << std::scientific << respPre.real() << std::endl;
				msg << "Imaginary part: " 
					<< std::setw(12) << std::setprecision(4) << std::scientific << respPrePre.imag()
					<< "->" 
					<< std::setw(12) << std::setprecision(4) << std::scientific << respPre.imag();
				ptrOutputFiles->writeCvgMessage(msg.str());
			}
		}
	}

	int*** mapCandidateToIndexes = new int**[numOfCandidates];
	for( int iCan = 0; iCan < numOfCandidates; ++iCan ){
		mapCandidateToIndexes[iCan] = new int*[numOfOutputAndInputVariables];
		for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
			mapCandidateToIndexes[iCan][iVar] = new int[numOfReferenceVariables];
		}
	}
	int chunkSize = numOfCandidates;
	for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
		for( int iRR = 0; iRR < numOfReferenceVariables; ++iRR ){
			const int numCans = static_cast<int>(candidates[iVar][iRR].size());
			chunkSize /= numCans;
			assert(chunkSize >= 1);
			int counter(0);
			for( int iChunk = 0; iChunk < numOfCandidates; iChunk += chunkSize, ++counter ){
				for( int iCan = iChunk; iCan < iChunk + chunkSize; ++iCan ){
					const int index = counter % numCans;
					mapCandidateToIndexes[iCan][iVar][iRR] = index; 
				}
			}
		}
	}

#ifdef _DEBUG_WRITE
	for( int iCan = 0; iCan < numOfCandidates; ++iCan ){
		std::cout << iCan << " ";
		for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
			for( int iRR = 0; iRR < numOfReferenceVariables; ++iRR ){
				const int index = mapCandidateToIndexes[iCan][iVar][iRR];
				std::cout << index << " ";
			}
		}
		std::cout << std::endl;
	}
#endif

	resp.reserve(numOfCandidates);
	for( int iCan = 0; iCan < numOfCandidates; ++iCan ){
		// This array must not be deleted because they are used out of the function
		std::complex<double>** respWork = new std::complex<double>*[numOfOutputAndInputVariables];
		for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
			respWork[iVar] = new std::complex<double>[numOfReferenceVariables];
			for( int iRR = 0; iRR < numOfReferenceVariables; ++iRR ){
				const int index = mapCandidateToIndexes[iCan][iVar][iRR];
				respWork[iVar][iRR] = candidates[iVar][iRR][index];
			}
		}
		resp.push_back(respWork);
	}

	// Delete  dymanic arrays
	for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
		delete [] candidates[iVar];
	}
	delete [] candidates;
	for( int iCan = 0; iCan < numOfCandidates; ++iCan ){
		for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
			delete [] mapCandidateToIndexes[iCan][iVar];
		}
		delete [] mapCandidateToIndexes[iCan];
	}
	delete [] mapCandidateToIndexes;

}

// Determine candidates by random sampling
void AnalysisMultivariateRegression::determineCandidatesByRandomSampling( const double freq, 
	const int numSegmentsTotal, std::complex<double>** ftval, const std::vector< std::pair<std::string, std::string> >& times, 
	std::vector< std::pair<int,int> >& candidatesOfPairs, std::vector< std::complex<double>** >& resp ) const{

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();

	const Control* const ptrControl = Control::getInstance();
	const int numOfOutputVariables = ptrControl->getNumOutputVariables();
	const int numOfInputVariables = ptrControl->getNumInputVariables();
	const int numOfOutputAndInputVariables = numOfOutputVariables + numOfInputVariables;
	const int degreeOfFreedom = 2 * numOfOutputAndInputVariables;
	const int numOfReferenceVariables = ptrControl->getNumRemoteReferenceVariables();
	assert(numOfReferenceVariables == 2); 
	const int rr0 = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, 0 );
	const int rr1 = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, 1 );

	// Determine candidates
	const Control::ParamsForRobustMultivariateRegression params = ptrControl->getParamsForRobustMultivariateRegression();
	const int timeStart = Util::convertHHMMSSToSeconds(params.startOfTimeRange);
	const int timeEnd = Util::convertHHMMSSToSeconds(params.endOfTimeRange);
	std::vector<int> segmentsInTimeRange;
	if( timeStart <= Util::convertHHMMSSToSeconds("00:00:00") && timeEnd >= Util::convertHHMMSSToSeconds("24:00:00") ){
		// Include all segments
		for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
			segmentsInTimeRange.push_back(iSeg);
		}
	}else{
		// Include only the segments which locate in the specified time range
		for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
			const int time1 = Util::convertHHMMSSToSeconds(times[iSeg].first);
			int time2 = Util::convertHHMMSSToSeconds(times[iSeg].second);
			if( time2 < time1 ){
				// Add 24 hr
				time2 += 24 * 3600;
			}
			if( time1 >= timeStart && time2 <= timeEnd ){
				segmentsInTimeRange.push_back(iSeg);
			}
		}
	}
	std::sort(segmentsInTimeRange.begin(), segmentsInTimeRange.end());
	segmentsInTimeRange.erase( std::unique( segmentsInTimeRange.begin(), segmentsInTimeRange.end() ), segmentsInTimeRange.end() );
	const int numOfSegmentsInTimeRange = static_cast<int>(segmentsInTimeRange.size());
	const long long int numOfSegmentsInTimeRange_64 = static_cast<long long int>(numOfSegmentsInTimeRange);
	const long long int numOfCandidatesOfPairsAll = (numOfSegmentsInTimeRange_64 * numOfSegmentsInTimeRange_64 - numOfSegmentsInTimeRange_64) / 2;

	const long long int numOfCancidatesMax = static_cast<long long int>(params.numOfMaxInitialCandidates);
	std::vector< std::pair<int,int> > candidatesOfPairsInit;
	if( numOfCancidatesMax >= numOfCandidatesOfPairsAll ){
		for( int i1 = 0; i1 < numOfSegmentsInTimeRange; ++i1 ){
			const int iSeg1 = segmentsInTimeRange[i1];
			for( int i2 = i1 + 1; i2 < numOfSegmentsInTimeRange; ++i2 ){
				const int iSeg2 = segmentsInTimeRange[i2];
				candidatesOfPairsInit.push_back(std::make_pair(iSeg1, iSeg2));
			}
		}
	}else{
#ifdef _RAND
		srand(1234);
#else
#ifdef _MERSENNE_TWISTER_ORIGINAL
		init_genrand64(1234);
#else
		std::mt19937_64 gen(1234);
		std::uniform_int_distribution<int> uniformDistibution(0, numOfSegmentsInTimeRange - 1);
#endif
#endif
		int numIterationMax = 1000;
		bool converge(false);
		for( int iter = 0; iter < numIterationMax; ++iter ){
#ifdef _RAND
			const int i1 = (rand() / RAND_MAX) * (numOfSegmentsInTimeRange - 1);
			const int i2 = (rand() / RAND_MAX) * (numOfSegmentsInTimeRange - 1);
#else
#ifdef _MERSENNE_TWISTER_ORIGINAL
			const int i1 = static_cast<int>(genrand64_real1() * numOfSegmentsInTimeRange);
			const int i2 = static_cast<int>(genrand64_real1() * numOfSegmentsInTimeRange);
#else
			const int i1 = uniformDistibution(gen);
			const int i2 = uniformDistibution(gen);
#endif
#endif
			const int iSeg1 = segmentsInTimeRange[i1];
			const int iSeg2 = segmentsInTimeRange[i2];
			if( iSeg2 > iSeg1 ){
				candidatesOfPairsInit.push_back(std::make_pair(iSeg1, iSeg2));
			}else if( iSeg1 > iSeg2 ){
				candidatesOfPairsInit.push_back(std::make_pair(iSeg2, iSeg1));
			}
			std::sort(candidatesOfPairsInit.begin(), candidatesOfPairsInit.end());
			candidatesOfPairsInit.erase( std::unique( candidatesOfPairsInit.begin(), candidatesOfPairsInit.end() ), candidatesOfPairsInit.end() );
			if( static_cast<long long int>(candidatesOfPairsInit.size()) >= numOfCancidatesMax){
				converge = true;
				break;
			}
		}
		if(!converge){
			ptrOutputFiles->writeWarningMessage("Rearch maximum iteration: " + Util::toString(numIterationMax));
		}
	}

	const int numOfInitialCancidates = static_cast<int>( candidatesOfPairsInit.size() );
	ptrOutputFiles->writeCvgAndLogMessage("Number of initial candidates: " + Util::toString(numOfInitialCancidates));

	resp.reserve(numOfInitialCancidates);
	// This array must not be deleted because they are used out of the function
	std::complex<double>*** respWork = new std::complex<double>**[numOfInitialCancidates];
	for( int iCan = 0; iCan < numOfInitialCancidates; ++iCan ){
		respWork[iCan] = new std::complex<double>*[numOfOutputAndInputVariables];
		for( int i = 0; i < numOfOutputAndInputVariables; ++i ){
			respWork[iCan][i] = new std::complex<double>[numOfReferenceVariables];
		}
	}
	for( int iCan = 0; iCan < numOfInitialCancidates; ++iCan ){
		const int iSeg1 = candidatesOfPairsInit[iCan].first;
		const int iSeg2 = candidatesOfPairsInit[iCan].second;
		const int rr0 = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, 0 );
		const int rr1 = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, 1 );
		bool determinantIsTooSmall(false);
		int iVar(0);
		for( int iOut = 0; iOut < numOfOutputVariables; ++iOut ){
			const int index = ptrControl->getChannelIndex( CommonParameters::OUTPUT, iOut );
			const bool tooSmall = calculateResponseFunctionsFromTwoSamples( 
				ftval[index][iSeg1], ftval[rr0][iSeg1], ftval[rr1][iSeg1],
				ftval[index][iSeg2], ftval[rr0][iSeg2], ftval[rr1][iSeg2],
				respWork[iCan][iVar][0], respWork[iCan][iVar][1] );
			++iVar;
			if(tooSmall){
				determinantIsTooSmall = true;
			}
		}
		for( int iInp = 0; iInp < numOfInputVariables; ++iInp ){
			const int index = ptrControl->getChannelIndex( CommonParameters::INPUT, iInp );
			const bool tooSmall = calculateResponseFunctionsFromTwoSamples( 
				ftval[index][iSeg1], ftval[rr0][iSeg1], ftval[rr1][iSeg1],
				ftval[index][iSeg2], ftval[rr0][iSeg2], ftval[rr1][iSeg2],
				respWork[iCan][iVar][0], respWork[iCan][iVar][1] );
			++iVar;
			if(tooSmall){
				determinantIsTooSmall = true;
			}
		}
		assert(iVar == numOfOutputAndInputVariables);
		if(determinantIsTooSmall){
			ptrOutputFiles->writeWarningMessage("Pair of segments (" + Util::toString(iSeg1) + ", "
				+ Util::toString(iSeg2) + ") is omitted because determinant is too small");
			// Delete because this candidate is not used out of the function
			for( int i = 0; i < numOfOutputAndInputVariables; ++i ){
				delete [] respWork[iCan][i];
			}
			delete [] respWork[iCan];
		}else{
			// Add candidate
			candidatesOfPairs.push_back( candidatesOfPairsInit[iCan] );
			resp.push_back(respWork[iCan]);
		}
	}

}

// Improve candidate by I-step
void AnalysisMultivariateRegression::improveCandidate( const int numSegmentsTotal, const int numOfMaxIterations, 
	const double convergenceCriteria, const bool initialCalculation, const double paramB, const double paramC,
	std::complex<double>** ftval, std::complex<double>** resp, double* variancesWithoutScale, double& scale, double& determinant,
	double* coherencesMin ) const{

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	const Control* const ptrControl = Control::getInstance();

	if( ptrControl->getOutputLevel() < 3 ){
		ptrOutputFiles->stopToWriteCvgMessage();
	}

	const int numOfOutputVariables = ptrControl->getNumOutputVariables();
	const int numOfInputVariables = ptrControl->getNumInputVariables();
	const int numOfOutputAndInputVariables = numOfOutputVariables + numOfInputVariables;
	const int degreeOfFreedom = 2 * numOfOutputAndInputVariables;

	const int numOfReferenceVariables = ptrControl->getNumRemoteReferenceVariables();
	assert(numOfReferenceVariables == 2); 

	const int rr0 = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, 0 );
	const int rr1 = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, 1 );
	std::complex<double>** complexResiduals = new std::complex<double>*[numOfOutputAndInputVariables];
	for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
		complexResiduals[iVar] = new std::complex<double>[numSegmentsTotal]; 
	}
	// Calculate complex residuals
	calculateComplexResiduals( numSegmentsTotal, ftval, resp, complexResiduals );
	if(initialCalculation){
		ptrOutputFiles->writeCvgMessage("Initial response functions:");
		std::ostringstream msg;
		int iVar = 0;
		for( int iOut = 0; iOut < numOfOutputVariables; ++iOut ){
			msg << "("  << std::setw(12) << std::setprecision(4) << std::scientific << resp[iVar][0].real() << "," 
					   << std::setw(12) << std::setprecision(4) << std::scientific << resp[iVar][0].imag() << "), ";
			msg << "("  << std::setw(12) << std::setprecision(4) << std::scientific << resp[iVar][1].real() << "," 
					   << std::setw(12) << std::setprecision(4) << std::scientific << resp[iVar][1].imag() << ")" << std::endl;
			++iVar;
		}
		for( int iInp = 0; iInp < numOfInputVariables; ++iInp ){
			msg << "("  << std::setw(12) << std::setprecision(4) << std::scientific << resp[iVar][0].real() << "," 
					   << std::setw(12) << std::setprecision(4) << std::scientific << resp[iVar][0].imag() << "), ";
			msg << "("  << std::setw(12) << std::setprecision(4) << std::scientific << resp[iVar][1].real() << "," 
					   << std::setw(12) << std::setprecision(4) << std::scientific << resp[iVar][1].imag() << ")";
			++iVar;
			if( iInp + 1 < numOfInputVariables ){
				msg << std::endl;
			}
		}
		assert(iVar == numOfOutputAndInputVariables);
		ptrOutputFiles->writeCvgMessage(msg.str());
		ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
		// Calculate variances withoutu scale
		double* work = new double[2 * numSegmentsTotal];
		for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
			// Make real residual vector from complex residual vector
			for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
				work[2*iSeg  ] = complexResiduals[iVar][iSeg].real();
				work[2*iSeg+1] = complexResiduals[iVar][iSeg].imag();		
			}
			const double sigma = Util::calculateMADN(2 * numSegmentsTotal, work);
#ifdef _DEBUG_WRITE
			std::cout << "iVar , sigma: " << iVar << " " << sigma << std::endl;
			for( int iSample = 0; iSample < 2 * numSegmentsTotal; ++iSample ){
				std::cout << work[iSample] << std::endl;
			}
#endif
			variancesWithoutScale[iVar] = pow(sigma, 2);
		}
		delete [] work;
		determinant = 1.0;
		for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
			determinant *= variancesWithoutScale[iVar];// Real part
			determinant *= variancesWithoutScale[iVar];// Imaginary part
		}
		const double factor = pow(determinant, -1.0/static_cast<double>(degreeOfFreedom));
		// Make determinant of covariance matrix one
		for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
			variancesWithoutScale[iVar] *= factor;
		}
#ifdef _DEBUG_WRITE
		for( int i = 0; i < numOfOutputAndInputVariables; ++i ){
			std::cout << "iVar, variancesWithoutScale: " << i << " " << variancesWithoutScale[i] << std::endl;
		}
#endif
	}

	// Calculate Mahalanobis distance
	double* MD = new double[numSegmentsTotal];
	calculateMD( numSegmentsTotal, numOfOutputAndInputVariables, complexResiduals, variancesWithoutScale, MD );

	if(initialCalculation){
		scale = Util::calculateMedian(numSegmentsTotal, MD);
		scale = RobustWeightTukeysBiweights::calculateRobustScale(numSegmentsTotal, MD, scale, paramB, paramC);
	}

	// Calculate response functions
	double* weights = new double[numSegmentsTotal];
	double* Gres = new double[numOfOutputAndInputVariables];
	double* GresPre = new double[numOfOutputAndInputVariables];
	double scalePre = scale;
	for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
		// Initialize
		Gres[iVar] = 0.0;
		GresPre[iVar] = 0.0;
	}
	bool converge(false);
	for( int iter = 0; iter < numOfMaxIterations; ++iter ){
		double sumOfWeights(0.0);
		for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
			weights[iSeg] = RobustWeightTukeysBiweights::calculateWeights(MD[iSeg] / scale, paramC);
			sumOfWeights += weights[iSeg];
		}
		ptrOutputFiles->writeCvgMessage("Iteration number = " + Util::toString(iter) + ", Scale = " + Util::toString(scale) 
			+ ", Sum of weights = " + Util::toString(sumOfWeights));
		std::ostringstream msgVariance;
		msgVariance << "Variances without scale: ";
		for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
			msgVariance << std::setw(12) << std::setprecision(4) << std::scientific << variancesWithoutScale[iVar];
			if( iVar + 1 < numOfOutputAndInputVariables ){
				msgVariance << ",";
			}
		}
		ptrOutputFiles->writeCvgMessage(msgVariance.str());
		if( sumOfWeights < CommonParameters::EPS ){
			ptrOutputFiles->writeErrorMessage("Sum of weights is too small: " + Util::toString(sumOfWeights));
		}
		for( int iOut = 0; iOut < numOfOutputVariables; ++iOut ){
			coherencesMin[iOut] = 1.0e20;// Initialize
		}
		int iVar(0);
		bool diverge(false);
		for( int iOut = 0; iOut < numOfOutputVariables; ++iOut ){
			const int index = ptrControl->getChannelIndex( CommonParameters::OUTPUT, iOut );
			double coherence(0.0);
			ptrOutputFiles->writeCvgMessage("Calculate response functions for output variable " +
				Util::toString(iOut) + " by the weighted least square method");
			const double GresWork = calculateResponseFunctionByWLS( ftval[index], ftval[rr0], ftval[rr1], numSegmentsTotal,
				weights, complexResiduals[iVar], resp[iVar][0], resp[iVar][1], coherence );
			if( coherence < coherencesMin[iOut] ){
				coherencesMin[iOut] = coherence;
			}
			if( GresWork < CommonParameters::EPS ){
				ptrOutputFiles->writeCvgMessage("	Weighted residual power is too small (" +  Util::toString(GresWork) + ")");
				ptrOutputFiles->writeWarningMessage("Weighted residual power is too small (" +  Util::toString(GresWork) + ")");
				diverge = true;
			}
			Gres[iVar] = GresWork;
			++iVar;
		}
		for( int iInp = 0; iInp < numOfInputVariables; ++iInp ){
			const int index = ptrControl->getChannelIndex( CommonParameters::INPUT, iInp );
			double coherence(0.0);
			ptrOutputFiles->writeCvgMessage("Calculate response functions for input variable " +
				Util::toString(iInp) + " by the weighted least square method");
			const double GresWork = calculateResponseFunctionByWLS( ftval[index], ftval[rr0], ftval[rr1], numSegmentsTotal,
				weights, complexResiduals[iVar], resp[iVar][0], resp[iVar][1], coherence );
			for( int iOut = 0; iOut < numOfOutputVariables; ++iOut ){
				if( coherence < coherencesMin[iOut] ){
					coherencesMin[iOut] = coherence;
				}
			}
			if( GresWork < CommonParameters::EPS ){
				ptrOutputFiles->writeCvgMessage("	Weighted residual power is too small (" +  Util::toString(GresWork) + ")");
				ptrOutputFiles->writeWarningMessage("Weighted residual power is too small (" +  Util::toString(GresWork) + ")");
				diverge = true;
			}
			Gres[iVar] = GresWork;
			++iVar;
		}
		std::ostringstream msgGres;
		msgGres << "Weighted residual powers: ";
		for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
			msgGres << std::setw(12) << std::setprecision(4) << std::scientific << Gres[iVar];
			if( iVar + 1 < numOfOutputAndInputVariables ){
				msgGres << ",";
			}
		}
		ptrOutputFiles->writeCvgMessage(msgGres.str());
		assert(iVar == numOfOutputAndInputVariables);
		if(diverge){
			// Go out from the iteration
			break;
		}
		// Calculate variances without scale
		determinant = 1.0;
		for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
			// Make real residual vector from complex residual vector
			double numerator (0.0);
			double denominator(0.0);
			for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
				numerator += std::norm(complexResiduals[iVar][iSeg]) * weights[iSeg];
				const double val = MD[iSeg] / scale;
				denominator += RobustWeightTukeysBiweights::calculateTermInDenominatorOfRobustCovariance(val, paramC);
			}
			if( denominator < CommonParameters::EPS ){
				ptrOutputFiles->writeErrorMessage("Denominator of robust covariance is too small (" +  Util::toString(denominator) + ")");
			}
			variancesWithoutScale[iVar] = static_cast<double>(numOfOutputAndInputVariables) * numerator / denominator;
			determinant *= variancesWithoutScale[iVar];// Real part
			determinant *= variancesWithoutScale[iVar];// Imaginary part
		}
		const double factor = pow(determinant, -1.0/static_cast<double>(degreeOfFreedom));
		// Make determinant of covariance matrix one
		for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
			variancesWithoutScale[iVar] *= factor;
		}
		// Calculate Mahalanobis distance
		calculateMD( numSegmentsTotal, numOfOutputAndInputVariables, complexResiduals,
			variancesWithoutScale, MD );
		scale = RobustWeightTukeysBiweights::calculateRobustScale(numSegmentsTotal, MD,
			scale, paramB, paramC );
		// Convergence judgment
		converge = true;
		if( fabs(scale - scalePre) / fabs(scalePre) > convergenceCriteria ){
			converge = false;
		}
		for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
			if( fabs(Gres[iVar] - GresPre[iVar]) / fabs(GresPre[iVar]) > convergenceCriteria ){
				converge = false;
				break;
			}
		}
		if(converge){
			// Converged
			// Go out from the iteration
			ptrOutputFiles->writeCvgMessage("Iteration converged");
			break;
		}
		scalePre = scale;
		for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
			GresPre[iVar] = Gres[iVar];
		}		
	}

	if(!converge){
		ptrOutputFiles->writeCvgMessage("Iteration does not converge");
	}

	for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
		delete [] complexResiduals[iVar];
	}
	delete [] complexResiduals;
	delete [] MD;
	delete [] weights;
	delete [] Gres;
	delete [] GresPre;

	if( ptrControl->getOutputLevel() < 3 ){
		ptrOutputFiles->restartToWriteCvgMessage();
	}

}

// Calculate Mahalanobis distances
void AnalysisMultivariateRegression::calculateMD( const int numSegmentsTotal, const int numOfOutputAndInputVariables,
		std::complex<double>** complexResiduals, const double* const variancesWithoutScale, double* MD ) const{

	for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
		// Zero clear
		MD[iSeg] = 0.0;
	}
	for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
		for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
			MD[iSeg] += pow(complexResiduals[iVar][iSeg].real(), 2) / variancesWithoutScale[iVar];// Real part
			MD[iSeg] += pow(complexResiduals[iVar][iSeg].imag(), 2) / variancesWithoutScale[iVar];// Imaginary part
		}
	}
	// Calculate Mahalanobis distance
	for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
		const double work = sqrt(MD[iSeg]);
		MD[iSeg] = work; 
//#ifdef _DEBUG_WRITE
//		std::cout << "iSeg , MD: " << iSeg << " " << MD[iSeg] << std::endl;
//#endif
	}

}

// Estimate error by fixed-weights bootstrap
void AnalysisMultivariateRegression::estimateErrorByFixedWeightsBootstrap( const int numSegmentsTotal, const double paramB, const double paramC,
	std::complex<double>** ftval, std::complex<double>** respOrg, const double* const variancesWithoutScaleOrg, const double scaleOrg,
	const double determinantOrg, double** respErr ) const{

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	const Control* const ptrControl = Control::getInstance();

	bool fixedWeights(true);
	ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
	ptrOutputFiles->writeCvgAndLogMessage("Estimate errors by fixed-weights bootstrap");
	ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");

	const int numOfOutputVariables = ptrControl->getNumOutputVariables();
	const int numOfInputVariables = ptrControl->getNumInputVariables();
	const int numOfOutputAndInputVariables = numOfOutputVariables + numOfInputVariables;
	const int numOfReferenceVariables = ptrControl->getNumRemoteReferenceVariables();
	const int degreeOfFreedom = 2 * numOfOutputAndInputVariables;

	assert(numOfReferenceVariables == 2); 
	const int rr0 = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, 0 );
	const int rr1 = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, 1 );

	std::complex<double>** complexResidualsOrg = new std::complex<double>*[numOfOutputAndInputVariables];
	for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
		complexResidualsOrg[iVar] = new std::complex<double>[numSegmentsTotal]; 
	}
	double* MDOrg = new double[numSegmentsTotal];
	double* weightsOrg = new double[numSegmentsTotal];
	double* termsInDenominatorOfRobustCovariance = new double[numSegmentsTotal];
	double* termsForScaleOrg = new double[numSegmentsTotal];

	// Calculate complex residuals
	calculateComplexResiduals( numSegmentsTotal, ftval, respOrg, complexResidualsOrg );
	// Calculate Mahalanobis distance
	calculateMD( numSegmentsTotal, numOfOutputAndInputVariables, complexResidualsOrg,
		variancesWithoutScaleOrg, MDOrg );
	// Calculate original weights
	double sumOfWeights(0.0);
	for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
		const double val = MDOrg[iSeg] / scaleOrg;
		weightsOrg[iSeg] = RobustWeightTukeysBiweights::calculateWeights(val, paramC);
		sumOfWeights += weightsOrg[iSeg];
		termsInDenominatorOfRobustCovariance[iSeg] = RobustWeightTukeysBiweights::calculateTermInDenominatorOfRobustCovariance(val, paramC);
		if( fabs(val) < CommonParameters::EPS ){
			// rho''(0) / 2 = 1 / 2
			termsForScaleOrg[iSeg] = 0.5 * pow(MDOrg[iSeg], 2);
		}else{
			termsForScaleOrg[iSeg] = RobustWeightTukeysBiweights::calculateLossFunction(val, paramC) * pow(scaleOrg, 2);
		}
	}
#ifdef _DEBUG_WRITE
	for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
		std::cout << iSeg << " " << weightsOrg[iSeg] << std::endl;
	}		
	for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
		std::cout << iSeg << " " << termsForScaleOrg[iSeg] << std::endl;
	}
	{
		double* variancesWithoutScale = new double[numOfOutputAndInputVariables];
		double determinant(1.0);
		for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
			// Make real residual vector from complex residual vector
			double numerator (0.0);
			double denominator(0.0);
			for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
				numerator += std::norm(complexResidualsOrg[iVar][iSeg]) * weightsOrg[iSeg];
				denominator += termsInDenominatorOfRobustCovariance[iSeg];
			}
			if( denominator < CommonParameters::EPS ){
				ptrOutputFiles->writeErrorMessage("Denominator of robust covariance is too small (" +  Util::toString(denominator) + ")");
			}
			variancesWithoutScale[iVar] = static_cast<double>(numOfOutputAndInputVariables) * numerator / denominator;
			determinant *= variancesWithoutScale[iVar];// Real part
			determinant *= variancesWithoutScale[iVar];// Imaginary part
		}
		const double factor = pow(determinant, -1.0/static_cast<double>(degreeOfFreedom));
		for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
			// Make determinant of covariance matrix one
			variancesWithoutScale[iVar] *= factor;
			std::cout << "variance[" << iVar << "] = " << variancesWithoutScale[iVar] << std::endl;
		}
		// Calculate one-step estimates of the scale of Mahalanobis distance
		double squareScale(0.0);
		for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
			squareScale += termsForScaleOrg[iSeg];
		}
		squareScale /= static_cast<double>(numSegmentsTotal);
		squareScale /= paramB;
		const double scale = sqrt(squareScale);
		std::cout << "scale = " << scale << std::endl;
	}
#endif
	if( sumOfWeights < CommonParameters::EPS ){
		ptrOutputFiles->writeErrorMessage("Sum of weights is too small: " + Util::toString(sumOfWeights));
	}

	// Find the minimum variance
	double minimumVariance(1.0e20);
	int indexOfMinimumVariance(-1);
	for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
		if( variancesWithoutScaleOrg[iVar] < minimumVariance ){
			minimumVariance = variancesWithoutScaleOrg[iVar];
			indexOfMinimumVariance = iVar;
		}
	}

//#ifdef _DEBUG_WRITE
//	DoubleDenseSquareMatrix testMatrix;
//	testMatrix.setDegreeOfEquation(4);
//	double tempMat[16] = {1.80, 2.88, 2.05,-0.89,5.25,-2.95,-0.95,-3.80,1.58,-2.69,-2.90,-1.04, -1.11,-0.66,-0.59, 0.80};
//	double tempVec[4] = {1,2,3,4};
//	double tempRes[4] = {0.0,0.0,0.0,0.0};
//	int itemp(0);
//	for( int irow = 0; irow < 4; ++irow ){
//		for( int icol = 0; icol < 4; ++icol ){
//			testMatrix.setValue(irow, icol, tempMat[itemp]);
//			++itemp;
//		}
//	}
//	testMatrix.debugWriteMatrix();
//	testMatrix.factorizeMatrix();
//	testMatrix.solveLinearEquation(tempVec, tempRes);
//	for( int i = 0; i < 4; ++i ){
//		std::cout << tempRes[i] << std::endl;
//	}
//#endif

	// Bootstrap
	int* segmentIndexes = new int [numSegmentsTotal];
	std::complex<double>** resp = new std::complex<double>*[numOfOutputAndInputVariables];
	for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
		resp[iVar] = new std::complex<double>[numOfReferenceVariables];
	}
	const int numOfDataSet = ptrControl->getNumRepetitionsOfBootstrap();
	std::complex<double>*** respFinal = new std::complex<double>**[numOfDataSet];
	double* variancesWithoutScale = new double[numOfOutputAndInputVariables];
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
	for( int iDataSet = 0; iDataSet < numOfDataSet; ++iDataSet ){
		// Make bootstrap samples
		for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
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
#ifdef _DEBUG_WRITE
		for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
			std::cout << iSeg << " " << segmentIndexes[iSeg] << std::endl;
		}		
#endif
		// Calculate one-step estimates of response functions
		int iVar(0);
		for( int iOut = 0; iOut < numOfOutputVariables; ++iOut ){
			const int index = ptrControl->getChannelIndex( CommonParameters::OUTPUT, iOut );
			calculateResponseFunctionByWLSForBootstrap( numSegmentsTotal, segmentIndexes, ftval[index], ftval[rr0], ftval[rr1],
				weightsOrg, resp[iVar][0], resp[iVar][1] );
			++iVar;
		}
		for( int iInp = 0; iInp < numOfInputVariables; ++iInp ){
			const int index = ptrControl->getChannelIndex( CommonParameters::INPUT, iInp );
			calculateResponseFunctionByWLSForBootstrap( numSegmentsTotal, segmentIndexes, ftval[index], ftval[rr0], ftval[rr1],
				weightsOrg, resp[iVar][0], resp[iVar][1] );
			++iVar;
		}
		assert( iVar == numOfOutputAndInputVariables );
		// Calculate one-step estimates of variances without scale
		double determinant(1.0);
		for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
			// Make real residual vector from complex residual vector
			double numerator (0.0);
			double denominator(0.0);
#ifdef _USE_OMP
			int icount(0);
			#pragma omp parallel for default(shared) private(icount) reduction(+ : numerator) 
			for( icount = 0; icount < numSegmentsTotal; ++icount ){
				numerator += std::norm(complexResidualsOrg[iVar][segmentIndexes[icount]]) * weightsOrg[segmentIndexes[icount]];
			}
			#pragma omp parallel for default(shared) private(icount) reduction(+ : denominator) 
			for( icount = 0; icount < numSegmentsTotal; ++icount ){
				denominator += termsInDenominatorOfRobustCovariance[segmentIndexes[icount]];
			}
#else
			for( int icount = 0; icount < numSegmentsTotal; ++icount ){
				const int iSeg = segmentIndexes[icount];
				numerator += std::norm(complexResidualsOrg[iVar][iSeg]) * weightsOrg[iSeg];
				denominator += termsInDenominatorOfRobustCovariance[iSeg];
			}
#endif
#ifdef _DEBUG_WRITE
			std::cout << "numerator, denominator = " << numerator << "," << denominator << std::endl;
#endif
			if( denominator < CommonParameters::EPS ){
				ptrOutputFiles->writeErrorMessage("Denominator of robust covariance is too small (" +  Util::toString(denominator) + ")");
			}
			variancesWithoutScale[iVar] = static_cast<double>(numOfOutputAndInputVariables) * numerator / denominator;
			determinant *= variancesWithoutScale[iVar];// Real part
			determinant *= variancesWithoutScale[iVar];// Imaginary part
		}
		const double factor = pow(determinant, -1.0/static_cast<double>(degreeOfFreedom));
		for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
			// Make determinant of covariance matrix one
			variancesWithoutScale[iVar] *= factor;
		}
		// Calculate one-step estimates of the scale of Mahalanobis distance
		double squareScale(0.0);
#ifdef _USE_OMP
		int icount(0);
		#pragma omp parallel for default(shared) private(icount) reduction(+ : squareScale) 
		for( icount = 0; icount < numSegmentsTotal; ++icount ){
			squareScale += termsForScaleOrg[segmentIndexes[icount]];
		}
#else
		for( int icount = 0; icount < numSegmentsTotal; ++icount ){
			const int iSeg = segmentIndexes[icount];
			squareScale += termsForScaleOrg[iSeg];
		}
#endif
#ifdef _DEBUG_WRITE
		std::cout << "squareScale = " << squareScale << std::endl;
#endif
		squareScale /= static_cast<double>(numSegmentsTotal);
		squareScale /= paramB;
		const double scale = sqrt(squareScale);
		if( ptrControl->getOutputLevel() >= 4 ){// Output one-step estimates
			ptrOutputFiles->writeCvgMessage("Dataset " + Util::toString(iDataSet));
			ptrOutputFiles->writeCvgMessage("One-step estimates of response functions:");
			std::ostringstream msg;
			int iVar = 0;
			for( int iOut = 0; iOut < numOfOutputVariables; ++iOut ){
				msg << "("  << std::setw(12) << std::setprecision(4) << std::scientific << resp[iVar][0].real() << "," 
						   << std::setw(12) << std::setprecision(4) << std::scientific << resp[iVar][0].imag() << "), ";
				msg << "("  << std::setw(12) << std::setprecision(4) << std::scientific << resp[iVar][1].real() << "," 
						   << std::setw(12) << std::setprecision(4) << std::scientific << resp[iVar][1].imag() << ")" << std::endl;
				++iVar;
			}
			for( int iInp = 0; iInp < numOfInputVariables; ++iInp ){
				msg << "("  << std::setw(12) << std::setprecision(4) << std::scientific << resp[iVar][0].real() << "," 
						   << std::setw(12) << std::setprecision(4) << std::scientific << resp[iVar][0].imag() << "), ";
				msg << "("  << std::setw(12) << std::setprecision(4) << std::scientific << resp[iVar][1].real() << "," 
						   << std::setw(12) << std::setprecision(4) << std::scientific << resp[iVar][1].imag() << ")";
				++iVar;
				if( iInp + 1 < numOfInputVariables ){
					msg << std::endl;
				}
			}
			assert(iVar == numOfOutputAndInputVariables);
			ptrOutputFiles->writeCvgMessage(msg.str());
		}
		respFinal[iDataSet] = new std::complex<double>*[numOfOutputVariables];
		for (int iOut = 0; iOut < numOfOutputVariables; ++iOut) {
			respFinal[iDataSet][iOut] = new std::complex<double>[numOfReferenceVariables];
		}
		if( ptrControl->getOutputLevel() >= 4 ){
			ptrOutputFiles->writeCvgMessage("One-step estimate of scale: " + Util::toString(scale));
			std::ostringstream msg;
			msg << "One-step estimates of variances without scale:";
			for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
				msg << std::setw(12) << std::setprecision(4) << std::scientific << variancesWithoutScale[iVar];
				if( iVar + 1 < numOfOutputAndInputVariables ){
					msg << ", ";
				}
			}
			ptrOutputFiles->writeCvgMessage(msg.str());
			// Calculate estimates of final response functions
			const int in0 = ptrControl->getChannelIndex(CommonParameters::INPUT, 0);
			const int in1 = ptrControl->getChannelIndex(CommonParameters::INPUT, 1);
			const std::complex<double> Txx = resp[in0][0];
			const std::complex<double> Txy = resp[in0][1];
			const std::complex<double> Tyx = resp[in1][0];
			const std::complex<double> Tyy = resp[in1][1];
			const std::complex<double> det = Txx * Tyy - Txy * Tyx;
			if (std::abs(det) < CommonParameters::EPS) {
				ptrOutputFiles->writeErrorMessage("Determinant is too small: " + Util::toString(std::abs(det)));
			}
			for (int iOut = 0; iOut < numOfOutputVariables; ++iOut) {
				const int index = ptrControl->getChannelIndex(CommonParameters::OUTPUT, iOut);
				const std::complex<double> U_x = resp[index][0];
				const std::complex<double> U_y = resp[index][1];
				const std::complex<double> resp0 = (U_x * Tyy - U_y * Tyx) / det;
				const std::complex<double> resp1 = (U_y * Txx - U_x * Txy) / det;
				respFinal[iDataSet][iOut][0] = resp0;
				respFinal[iDataSet][iOut][1] = resp1;
			}
			ptrOutputFiles->writeCvgMessage("One-step estimates of final response functions:");
			std::ostringstream msg2;
			for (int iOut = 0; iOut < numOfOutputVariables; ++iOut) {
				msg2 << "(" << std::setw(12) << std::setprecision(4) << std::scientific << respFinal[iDataSet][iOut][0].real() << ","
					<< std::setw(12) << std::setprecision(4) << std::scientific << respFinal[iDataSet][iOut][0].imag() << "), ";
				msg2 << "(" << std::setw(12) << std::setprecision(4) << std::scientific << respFinal[iDataSet][iOut][1].real() << ","
					<< std::setw(12) << std::setprecision(4) << std::scientific << respFinal[iDataSet][iOut][1].imag() << ")";
				if (iOut + 1 < numOfOutputVariables) {
					msg2 << std::endl;
				}
			}
			ptrOutputFiles->writeCvgMessage(msg2.str());
			ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
		}

		// Calculate estimates of final response functions
		const int in0 = ptrControl->getChannelIndex(CommonParameters::INPUT, 0);
		const int in1 = ptrControl->getChannelIndex(CommonParameters::INPUT, 1);
		const std::complex<double> Txx = resp[in0][0];
		const std::complex<double> Txy = resp[in0][1];
		const std::complex<double> Tyx = resp[in1][0];
		const std::complex<double> Tyy = resp[in1][1];
		const std::complex<double> det = Txx * Tyy - Txy * Tyx;
		if( std::abs(det) < CommonParameters::EPS ){
			ptrOutputFiles->writeErrorMessage("Determinant is too small: " + Util::toString(std::abs(det)));
		}
		for( int iOut = 0; iOut < numOfOutputVariables; ++iOut ){
			const int index = ptrControl->getChannelIndex( CommonParameters::OUTPUT, iOut );
			const std::complex<double> U_x = resp[index][0];
			const std::complex<double> U_y = resp[index][1];
			const std::complex<double> resp0 = ( U_x * Tyy - U_y * Tyx ) / det;
			const std::complex<double> resp1 = ( U_y * Txx - U_x * Txy ) / det;
			respFinal[iDataSet][iOut][0] = resp0;
			respFinal[iDataSet][iOut][1] = resp1;
		}
		if( ptrControl->getOutputLevel() >= 4 ){// Output estimates of final response functions
			ptrOutputFiles->writeCvgMessage("Estimates of final response functions:");
			std::ostringstream msg;
			for( int iOut = 0; iOut < numOfOutputVariables; ++iOut ){
				msg << "("  << std::setw(12) << std::setprecision(4) << std::scientific << respFinal[iDataSet][iOut][0].real() << "," 
						   << std::setw(12) << std::setprecision(4) << std::scientific << respFinal[iDataSet][iOut][0].imag() << "), ";
				msg << "("  << std::setw(12) << std::setprecision(4) << std::scientific << respFinal[iDataSet][iOut][1].real() << "," 
						   << std::setw(12) << std::setprecision(4) << std::scientific << respFinal[iDataSet][iOut][1].imag() << ")";
				if( iOut + 1 < numOfOutputVariables ){
					msg << std::endl;
				}
			}
			ptrOutputFiles->writeCvgMessage(msg.str());
			ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
		}
	}
	delete [] MDOrg;
	for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
		delete [] complexResidualsOrg[iVar];
	}
	delete [] complexResidualsOrg;
	delete [] weightsOrg;
	delete [] termsInDenominatorOfRobustCovariance;
	delete [] termsForScaleOrg;
	delete [] segmentIndexes;
	for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
		delete [] resp[iVar];
	}
	delete [] resp;
	delete [] variancesWithoutScale;

	// Calculate error of response functions
	assert(numOfDataSet > 2);
	for( int iOut = 0; iOut < numOfOutputVariables; ++iOut ){
		for( int iRR = 0; iRR < numOfReferenceVariables; ++iRR ){
			// Calculate average
			std::complex<double> average = std::complex<double>(0.0, 0.0);
			for( int iDataSet = 0; iDataSet < numOfDataSet; ++iDataSet ){
				average += respFinal[iDataSet][iOut][iRR];
			}
			average /= static_cast<double>(numOfDataSet);
			// Calculate variance
			double variance(0.0);
			for( int iDataSet = 0; iDataSet < numOfDataSet; ++iDataSet ){
				variance += std::norm(respFinal[iDataSet][iOut][iRR] - average);
			}
			variance /= static_cast<double>(2 * numOfDataSet - 4);
			respErr[iOut][iRR] = sqrt(variance);
		}
	}

	for( int iDataSet = 0; iDataSet < numOfDataSet; ++iDataSet ){
		for( int iOut = 0; iOut < numOfOutputVariables; ++iOut ){
			delete [] respFinal[iDataSet][iOut];
		}
		delete [] respFinal[iDataSet];
	}
	delete [] respFinal;

}

// Estimate error by strict bootstrap
void AnalysisMultivariateRegression::estimateErrorByStrictBootstrap( const int numSegmentsTotal, const double paramB, const double paramC, 
	std::complex<double>** ftvalOrg, std::complex<double>** respOrg, const double* const variancesWithoutScaleOrg, const double scaleOrg,
	const double determinantOrg, double** respErr ) const{

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeLogMessage("Strict bootstrap is performed to estimate errors");

	const Control* const ptrControl = Control::getInstance();
	const int numOfOutputVariables = ptrControl->getNumOutputVariables();
	const int numOfInputVariables = ptrControl->getNumInputVariables();
	const int numOfOutputAndInputVariables = numOfOutputVariables + numOfInputVariables;
	const int numOfReferenceVariables = ptrControl->getNumRemoteReferenceVariables();
	assert(numOfReferenceVariables == 2); 
	const int numChannels = ptrControl->getNumberOfChannels();

	const Control::ParamsForRobustMultivariateRegression params = ptrControl->getParamsForRobustMultivariateRegression();
	const int numOfMaxIterations = params.numOfMaxIterationsOfSecondIstep;
	const double convergenceCriteria = params.convergenceCriteriaOfSecondIstep;

	// Copy Fourier transformed values
	std::complex<double>** ftvalForBootstrap = new std::complex<double>*[numChannels];
	for( int iChan = 0; iChan < numChannels; ++iChan ){
		ftvalForBootstrap[iChan] = new std::complex<double>[numSegmentsTotal];
	}
	std::complex<double>** resp = new std::complex<double>*[numOfOutputAndInputVariables];
	for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
		resp[iVar] = new std::complex<double>[numOfReferenceVariables];
	}
	double* variancesWithoutScale = new double[numOfOutputAndInputVariables];
	double scales = scaleOrg;
	double determinants = determinantOrg;
	double* coherences = new double[numOfOutputVariables];
	int* segmentIndexes = new int[numSegmentsTotal];

	const int numOfSamples = ptrControl->getNumRepetitionsOfBootstrap();
	std::complex<double>*** respFinal = new std::complex<double>**[numOfSamples];
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
	for( int iSample = 0; iSample < numOfSamples; ++iSample ){
		// Make bootstrap samples
		for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
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
#ifdef _DEBUG_WRITE
		for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
			std::cout << iSeg << " " << segmentIndexes[iSeg] << std::endl;
		}
#endif
		// Copy data
		for( int iChan = 0; iChan < numChannels; ++iChan ){
			for( int icount = 0; icount < numSegmentsTotal; ++icount ){
				const int iSeg = segmentIndexes[icount];
				ftvalForBootstrap[iChan][icount] = ftvalOrg[iChan][iSeg];
			}
		}
		for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
			for( int iRR = 0; iRR < numOfReferenceVariables; ++iRR ){
				resp[iVar][iRR] = respOrg[iVar][iRR];
			}
		}
		for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
			variancesWithoutScale[iVar] = variancesWithoutScaleOrg[iVar];
		}
		double scales = scaleOrg;
		double determinants = determinantOrg;
		for( int iVar = 0; iVar < numOfOutputVariables; ++iVar ){
			coherences[iVar] = 0.0;
		}
		improveCandidate( numSegmentsTotal, numOfMaxIterations, convergenceCriteria, true, paramB, paramC,
			ftvalForBootstrap, resp, variancesWithoutScale, scales, determinants, coherences );
		// Calculate estimates of final response functions
		const int in0 = ptrControl->getChannelIndex( CommonParameters::INPUT, 0 );
		const int in1 = ptrControl->getChannelIndex( CommonParameters::INPUT, 1 );
		const std::complex<double> Txx = resp[in0][0];
		const std::complex<double> Txy = resp[in0][1];
		const std::complex<double> Tyx = resp[in1][0];
		const std::complex<double> Tyy = resp[in1][1];
		const std::complex<double> det = Txx * Tyy - Txy * Tyx;
		if( std::abs(det) < CommonParameters::EPS ){
			ptrOutputFiles->writeErrorMessage("Determinant is too small: " + Util::toString(std::abs(det)));
		}
		respFinal[iSample] = new std::complex<double>*[numOfOutputVariables];
		for( int iOut = 0; iOut < numOfOutputVariables; ++iOut ){
			respFinal[iSample][iOut] = new std::complex<double>[numOfReferenceVariables];
			const int index = ptrControl->getChannelIndex( CommonParameters::OUTPUT, iOut );
			const std::complex<double> U_x = resp[index][0];
			const std::complex<double> U_y = resp[index][1];
			const std::complex<double> resp0 = ( U_x * Tyy - U_y * Tyx ) / det;
			const std::complex<double> resp1 = ( U_y * Txx - U_x * Txy ) / det;
			respFinal[iSample][iOut][0] = resp0;
			respFinal[iSample][iOut][1] = resp1;
		}
	}
	ptrOutputFiles->restartToWriteCvgMessage();
	ptrOutputFiles->restartToWriteLogMessage();
	ptrOutputFiles->restartToWriteWarningMessage();

	if( ptrControl->getOutputLevel() >= 4 ){// Output estimates of final response functions
		for( int iSample = 0; iSample < numOfSamples; ++iSample ){
			ptrOutputFiles->writeCvgMessage("Dataset " + Util::toString(iSample));
			ptrOutputFiles->writeCvgMessage("Estimates of final response functions:");
			std::ostringstream msg;
			for( int iOut = 0; iOut < numOfOutputVariables; ++iOut ){
				msg << "("  << std::setw(12) << std::setprecision(4) << std::scientific << respFinal[iSample][iOut][0].real() << "," 
						   << std::setw(12) << std::setprecision(4) << std::scientific << respFinal[iSample][iOut][0].imag() << "), ";
				msg << "("  << std::setw(12) << std::setprecision(4) << std::scientific << respFinal[iSample][iOut][1].real() << "," 
						   << std::setw(12) << std::setprecision(4) << std::scientific << respFinal[iSample][iOut][1].imag() << ")";
				if( iOut + 1 < numOfOutputVariables ){
					msg << std::endl;
				}
			}
			ptrOutputFiles->writeCvgMessage(msg.str());
			ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
		}
	}

	// Delete arrays
	for( int iChan = 0; iChan < numChannels; ++iChan ){
		delete [] ftvalForBootstrap[iChan];
	}
	delete [] ftvalForBootstrap;
	for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
		delete [] resp[iVar];
	}
	delete [] resp;
	delete [] variancesWithoutScale;
	delete [] coherences;
	delete [] segmentIndexes;

	// Calculate error of response functions
	assert(numOfSamples > 2);
	for( int iOut = 0; iOut < numOfOutputVariables; ++iOut ){
		for( int iRR = 0; iRR < numOfReferenceVariables; ++iRR ){
			// Calculate average
			std::complex<double> average = std::complex<double>(0.0, 0.0);
			for( int iSample = 0; iSample < numOfSamples; ++iSample ){
				average += respFinal[iSample][iOut][iRR];
			}
			average /= static_cast<double>(numOfSamples);
			// Calculate variance
			double variance(0.0);
			for( int iSample = 0; iSample < numOfSamples; ++iSample ){
				variance += std::norm(respFinal[iSample][iOut][iRR] - average);
			}
			variance /= static_cast<double>(2 * numOfSamples - 4);
			respErr[iOut][iRR] = sqrt(variance);
		}
	}

	// Delete arrays
	for( int iSample = 0; iSample < numOfSamples; ++iSample ){
		for( int iOut = 0; iOut < numOfOutputVariables; ++iOut ){
			delete [] respFinal[iSample][iOut];
		}
		delete [] respFinal[iSample];
	}
	delete [] respFinal;

}

// Estimate error by fixed-weights jackknife
void AnalysisMultivariateRegression::estimateErrorByFixedWeightsJackknife( const int numSegmentsTotal, const double paramB, const double paramC, 
	std::complex<double>** ftval, const std::complex<double>* const respOrg0, const std::complex<double>* const respOrg1,
	std::complex<double>** respOrg, const double* const variancesWithoutScaleOrg, const double scaleOrg, double** respErr ) const{

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
	ptrOutputFiles->writeCvgAndLogMessage("Estimate errors by fixed-weights jackknife");
	ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");

	const Control* const ptrControl = Control::getInstance();
	const int numOfOutputVariables = ptrControl->getNumOutputVariables();
	const int numOfInputVariables = ptrControl->getNumInputVariables();
	const int numOfOutputAndInputVariables = numOfOutputVariables + numOfInputVariables;
	const int numOfReferenceVariables = ptrControl->getNumRemoteReferenceVariables();

	assert(numOfReferenceVariables == 2); 
	const int rr0 = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, 0 );
	const int rr1 = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, 1 );

	std::complex<double>** complexResidualsOrg = new std::complex<double>*[numOfOutputAndInputVariables];
	for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
		complexResidualsOrg[iVar] = new std::complex<double>[numSegmentsTotal]; 
	}
	double* MDOrg = new double[numSegmentsTotal];
	double* weightsOrg = new double[numSegmentsTotal];
	// Calculate complex residuals
	calculateComplexResiduals( numSegmentsTotal, ftval, respOrg, complexResidualsOrg );
	// Calculate Mahalanobis distance
	calculateMD( numSegmentsTotal, numOfOutputAndInputVariables, complexResidualsOrg,
		variancesWithoutScaleOrg, MDOrg );
	// Calculate original weights
	double sumOfWeights(0.0);
	for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
		const double val = MDOrg[iSeg] / scaleOrg;
		weightsOrg[iSeg] = RobustWeightTukeysBiweights::calculateWeights(val, paramC);
		sumOfWeights += weightsOrg[iSeg];
	}
#ifdef _DEBUG_WRITE
	for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
		std::cout << iSeg << " " << weightsOrg[iSeg] << std::endl;
	}
#endif
	if( sumOfWeights < CommonParameters::EPS ){
		ptrOutputFiles->writeErrorMessage("Sum of weights is too small: " + Util::toString(sumOfWeights));
	}
	// Release memory
	for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
		delete [] complexResidualsOrg[iVar];
	}
	delete [] complexResidualsOrg;
	delete [] MDOrg;

	double * hatDiagonals = new double [numSegmentsTotal];
	const double maxHatDiag = calculateDiagonalComponentsOfHatMatrix( numSegmentsTotal, rr0, rr1, ftval, weightsOrg, hatDiagonals );
#ifdef _DEBUG_WRITE
	for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
		std::cout << iSeg << " " << hatDiagonals[iSeg] << std::endl;
	}
#endif

	double* weights = new double[numSegmentsTotal];
	memcpy(weights, weightsOrg, sizeof(double)*numSegmentsTotal);
	std::complex<double>* temp0 = new std::complex<double>[numOfOutputAndInputVariables];
	std::complex<double>* temp1 = new std::complex<double>[numOfOutputAndInputVariables];
	std::complex<double>* resp0 = new std::complex<double>[numOfOutputVariables];
	std::complex<double>* resp1 = new std::complex<double>[numOfOutputVariables];
	std::complex<double>** pseudoResp0 = new std::complex<double>*[numSegmentsTotal];
	std::complex<double>** pseudoResp1 = new std::complex<double>*[numSegmentsTotal];
	for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
		weights[iSeg] = 0.0;// Replace
#ifdef _DEBUG_WRITE
		for( int i = 0; i < numSegmentsTotal; ++i ){
			std::cout << i << " " << weights[i] << std::endl;
		}
#endif
		// Calculate one-step estimates of response functions
		int iVar(0);
		for( int iOut = 0; iOut < numOfOutputVariables; ++iOut ){
			const int index = ptrControl->getChannelIndex( CommonParameters::OUTPUT, iOut );
			calculateResponseFunctionByWLSAux( numSegmentsTotal, ftval[index], ftval[rr0], ftval[rr1], weights, temp0[iVar], temp1[iVar] );
			++iVar;
		}
		for( int iInp = 0; iInp < numOfInputVariables; ++iInp ){
			const int index = ptrControl->getChannelIndex( CommonParameters::INPUT, iInp );
			calculateResponseFunctionByWLSAux( numSegmentsTotal, ftval[index], ftval[rr0], ftval[rr1], weights, temp0[iVar], temp1[iVar] );
			++iVar;
		}
		weights[iSeg] = weightsOrg[iSeg];// Restore
#ifdef _DEBUG_WRITE
		for( int i = 0; i < numSegmentsTotal; ++i ){
			std::cout << i << " " << weights[i] << std::endl;
		}
#endif
		assert(numOfInputVariables == 2);
		const int in0 = ptrControl->getChannelIndex( CommonParameters::INPUT, 0 );
		const int in1 = ptrControl->getChannelIndex( CommonParameters::INPUT, 1 );
		const std::complex<double> Txx = temp0[in0];
		const std::complex<double> Txy = temp1[in0];
		const std::complex<double> Tyx = temp0[in1];
		const std::complex<double> Tyy = temp1[in1];
		const std::complex<double> det = Txx * Tyy - Txy * Tyx;
		if( std::abs(det) < CommonParameters::EPS ){
			ptrOutputFiles->writeErrorMessage("Determinant is too small: " + Util::toString(std::abs(det)));
		}
		pseudoResp0[iSeg] = new std::complex<double>[numOfOutputVariables];
		pseudoResp1[iSeg] = new std::complex<double>[numOfOutputVariables];
		for( int iOut = 0; iOut < numOfOutputVariables; ++iOut ){
			const int index = ptrControl->getChannelIndex( CommonParameters::OUTPUT, iOut );
			const std::complex<double> U_x = temp0[index];
			const std::complex<double> U_y = temp1[index];
			resp0[iOut] = ( U_x * Tyy - U_y * Tyx ) / det;
			resp1[iOut] = ( U_y * Txx - U_x * Txy ) / det;
			const double hatMatrixDiagonal = hatDiagonals[iSeg];
			double factor = static_cast<double>(numSegmentsTotal) * (1.0 - hatMatrixDiagonal);
			if( hatMatrixDiagonal > 1.0 ){
				factor= 0.0;
			}else if( hatMatrixDiagonal < 0.0 ){
				factor = static_cast<double>(numSegmentsTotal);
			}
			pseudoResp0[iSeg][iOut] = respOrg0[iOut] + factor * (respOrg0[iOut] - resp0[iOut]);
			pseudoResp1[iSeg][iOut] = respOrg1[iOut] + factor * (respOrg1[iOut] - resp1[iOut]);
		}
	}
	delete [] hatDiagonals;
	delete [] weights;
	delete [] weightsOrg;
	delete [] temp0;
	delete [] temp1;
	delete [] resp0;
	delete [] resp1;

	// Calculate & output error bars
	for( int iOut = 0; iOut < numOfOutputVariables; ++iOut ){
		if( numSegmentsTotal > 2 ){
			std::complex<double> avgResp0(0.0, 0.0);
			std::complex<double> avgResp1(0.0, 0.0);
			for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
				avgResp0 += pseudoResp0[iSeg][iOut];
				avgResp1 += pseudoResp1[iSeg][iOut];
			}
			const double factor = 1.0 / static_cast<double>(numSegmentsTotal);
			avgResp0 *= factor;
			avgResp1 *= factor;
			double variance0(0.0);
			double variance1(0.0);
			for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
				variance0 += std::norm( pseudoResp0[iSeg][iOut] - avgResp0 );
				variance1 += std::norm( pseudoResp1[iSeg][iOut] - avgResp1 );
			}
			const double factor2 = factor / static_cast<double>(2 * numSegmentsTotal - 4);	
			variance0 *= factor2;
			variance1 *= factor2;
			respErr[iOut][0] = sqrt(variance0);
			respErr[iOut][1] = sqrt(variance1);
		}else{
			respErr[iOut][0] = 1.0e10;
			respErr[iOut][1] = 1.0e10;
		}
	}
	// Release memory
	for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
		delete [] pseudoResp0[iSeg];
		delete [] pseudoResp1[iSeg];
	}
	delete [] pseudoResp0;
	delete [] pseudoResp1;

}

// Estimate error by subset deletion jackknife
void AnalysisMultivariateRegression::estimateErrorBySubsetDeletionJackknife( const int numSegmentsTotal, std::complex<double>** ftvalOrg,
	std::complex<double>** respOrg, const double* const variancesWithoutScaleOrg, const double scaleOrg, const double determinantOrg, 
	double** respErr ) const{

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeLogMessage("Subset deletion jackknife is performed to estimate errors");
	const Control* const ptrControl = Control::getInstance();
	int numOmittedData = static_cast<int>( 0.01 * ptrControl->getPercentageOfOmmitedDataSubsetDeletionJackknife() * static_cast<double>(numSegmentsTotal) );
	if( numOmittedData < 1 ){
		numOmittedData = 1;
	}
	const int numOfSubsets = numSegmentsTotal / numOmittedData;
	ptrOutputFiles->writeLogMessage("Number of ommited data : " + Util::toString(numOmittedData));
	ptrOutputFiles->writeLogMessage("Number of subsets : " + Util::toString(numOfSubsets));

	const int numOfOutputVariables = ptrControl->getNumOutputVariables();
	const int numOfInputVariables = ptrControl->getNumInputVariables();
	const int numOfOutputAndInputVariables = numOfOutputVariables + numOfInputVariables;
	const int degreeOfFreedom = 2 * numOfOutputAndInputVariables;
	const int numOfReferenceVariables = ptrControl->getNumRemoteReferenceVariables();
	assert(numOfReferenceVariables == 2); 
	const int numChannels = ptrControl->getNumberOfChannels();

	double paramB(0.0);
	double paramC(0.0);
	RobustWeightTukeysBiweights::calculateParams(degreeOfFreedom, numSegmentsTotal - numOmittedData, paramB, paramC);

	std::complex<double>** complexResidualsOrg = new std::complex<double>*[numOfOutputAndInputVariables];
	for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
		complexResidualsOrg[iVar] = new std::complex<double>[numSegmentsTotal]; 
	}
	double* MDOrg = new double[numSegmentsTotal];
	double* weightsOrg = new double[numSegmentsTotal];
	// Calculate complex residuals
	calculateComplexResiduals( numSegmentsTotal, ftvalOrg, respOrg, complexResidualsOrg );
	// Calculate Mahalanobis distance
	calculateMD( numSegmentsTotal, numOfOutputAndInputVariables, complexResidualsOrg,
		variancesWithoutScaleOrg, MDOrg );
	// Calculate original weights
	double sumOfWeights(0.0);
	for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
		const double val = MDOrg[iSeg] / scaleOrg;
		weightsOrg[iSeg] = RobustWeightTukeysBiweights::calculateWeights(val, paramC);
		sumOfWeights += weightsOrg[iSeg];
	}
#ifdef _DEBUG_WRITE
	for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
		std::cout << iSeg << " " << weightsOrg[iSeg] << std::endl;
	}
#endif
	if( sumOfWeights < CommonParameters::EPS ){
		ptrOutputFiles->writeErrorMessage("Sum of weights is too small: " + Util::toString(sumOfWeights));
	}
	// Release memory
	for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
		delete [] complexResidualsOrg[iVar];
	}
	delete [] complexResidualsOrg;
	delete [] MDOrg;
	double * hatDiagonals = new double [numSegmentsTotal];
	assert(numOfReferenceVariables == 2); 
	const int rr0 = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, 0 );
	const int rr1 = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, 1 );
	const double maxHatDiag = calculateDiagonalComponentsOfHatMatrix( numSegmentsTotal, rr0, rr1, ftvalOrg, weightsOrg, hatDiagonals );
	delete [] weightsOrg;
#ifdef _DEBUG_WRITE
	for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
		std::cout << iSeg << " " << hatDiagonals[iSeg] << std::endl;
	}
#endif

	const Control::ParamsForRobustMultivariateRegression params = ptrControl->getParamsForRobustMultivariateRegression();
	const int numOfMaxIterations = params.numOfMaxIterationsOfSecondIstep;
	const double convergenceCriteria = params.convergenceCriteriaOfSecondIstep;

#ifdef _DEBUG_WRITE
	for( int iSeg = 0;iSeg < numSegmentsTotal; ++iSeg){
		std::cout << std::setw(15) << hatDiagonals[iSeg] << std::endl;
	}
	for( int iChan = 0; iChan < ptrControl->getNumberOfChannels(); ++iChan ){
		for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
			std::cout << "iChan iSeg val: " << iChan << " " << iSeg << " " << ftvalOrg[iChan][iSeg] << std::endl;
		}
	}
#endif

	// Copy Fourier transformed values
	std::complex<double>** ftvalForJackknife = new std::complex<double>*[numChannels];
	for( int iChan = 0; iChan < numChannels; ++iChan ){
		ftvalForJackknife[iChan] = new std::complex<double>[numSegmentsTotal - numOmittedData];
	}
	std::complex<double>** resp = new std::complex<double>*[numOfOutputAndInputVariables];
	for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
		resp[iVar] = new std::complex<double>[numOfReferenceVariables];
	}
	double* variancesWithoutScale = new double[numOfOutputAndInputVariables];
	double scales = scaleOrg;
	double determinants = determinantOrg;
	double* coherences = new double[numOfOutputVariables];
	std::complex<double>** pseudoResp0 = new std::complex<double>*[numOfSubsets];
	std::complex<double>** pseudoResp1 = new std::complex<double>*[numOfSubsets];
	ptrOutputFiles->stopToWriteCvgMessage();
	ptrOutputFiles->stopToWriteLogMessage();
	ptrOutputFiles->stopToWriteWarningMessage();
	for( int iSubset = 0; iSubset < numOfSubsets; ++iSubset ){
		double averageHatDiags(0.0);
		const int iSegOmitStart = iSubset * numOmittedData;
		const int iSegOmitEnd = iSegOmitStart + numOmittedData;
		assert(iSegOmitEnd <= numSegmentsTotal);
		for( int iSeg = iSegOmitStart; iSeg < iSegOmitEnd; ++iSeg ){
			averageHatDiags += hatDiagonals[iSeg];
		}
		averageHatDiags /= static_cast<double>(numOmittedData);
#ifdef _DEBUG_WRITE
		std::cout << std::setw(20) << averageHatDiags << std::endl;
#endif
		// Copy data
		for( int iChan = 0; iChan < numChannels; ++iChan ){
			int icount(0);
			for( int iSeg = 0; iSeg < iSegOmitStart; ++iSeg, ++icount ){
				ftvalForJackknife[iChan][icount] = ftvalOrg[iChan][iSeg];
			}
			for( int iSeg = iSegOmitEnd; iSeg < numSegmentsTotal; ++iSeg, ++icount ){
				ftvalForJackknife[iChan][icount] = ftvalOrg[iChan][iSeg];
			}
			assert(icount == numSegmentsTotal - numOmittedData);
		}
#ifdef _DEBUG_WRITE
		for( int iChan = 0; iChan < numChannels; ++iChan ){
			for( int i = 0; i < numSegmentsTotal - numOmittedData; ++i ){
				std::cout << "iChan i val: " << iChan << " " << i << " " << ftvalForJackknife[iChan][i] << std::endl;
			}
		}
#endif
		for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
			for( int iRR = 0; iRR < numOfReferenceVariables; ++iRR ){
				resp[iVar][iRR] = respOrg[iVar][iRR];
			}
		}
		for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
			variancesWithoutScale[iVar] = variancesWithoutScaleOrg[iVar];
		}
		double scales = scaleOrg;
		double determinants = determinantOrg;
		for( int iVar = 0; iVar < numOfOutputVariables; ++iVar ){
			coherences[iVar] = 0.0;
		}
		improveCandidate( numSegmentsTotal - numOmittedData, numOfMaxIterations, convergenceCriteria, true, paramB, paramC,
			ftvalForJackknife, resp, variancesWithoutScale, scales, determinants, coherences );
		// Calculate estimates of final response functions
		const int in0 = ptrControl->getChannelIndex( CommonParameters::INPUT, 0 );
		const int in1 = ptrControl->getChannelIndex( CommonParameters::INPUT, 1 );
		const std::complex<double> Txx = resp[in0][0];
		const std::complex<double> Txy = resp[in0][1];
		const std::complex<double> Tyx = resp[in1][0];
		const std::complex<double> Tyy = resp[in1][1];
		const std::complex<double> det = Txx * Tyy - Txy * Tyx;
		if( std::abs(det) < CommonParameters::EPS ){
			ptrOutputFiles->writeErrorMessage("Determinant is too small: " + Util::toString(std::abs(det)));
		}
		pseudoResp0[iSubset] = new std::complex<double>[numOfOutputVariables];
		pseudoResp1[iSubset] = new std::complex<double>[numOfOutputVariables];
		for( int iOut = 0; iOut < numOfOutputVariables; ++iOut ){
			double factor = static_cast<double>(numOfSubsets) * (1.0 - averageHatDiags);
			if( averageHatDiags > 1.0 ){
				factor= 0.0;
			}else if( averageHatDiags < 0.0 ){
				factor = static_cast<double>(numOfSubsets);
			}
			const int index = ptrControl->getChannelIndex( CommonParameters::OUTPUT, iOut );
			const std::complex<double> U_x = resp[index][0];
			const std::complex<double> U_y = resp[index][1];
			const std::complex<double> resp0 = ( U_x * Tyy - U_y * Tyx ) / det;
			const std::complex<double> resp1 = ( U_y * Txx - U_x * Txy ) / det;
			pseudoResp0[iSubset][iOut] = respOrg[index][0] + factor * (respOrg[index][0] - resp0);
			pseudoResp1[iSubset][iOut] = respOrg[index][1] + factor * (respOrg[index][1] - resp1);
		}
	}
	ptrOutputFiles->restartToWriteCvgMessage();
	ptrOutputFiles->restartToWriteLogMessage();
	ptrOutputFiles->restartToWriteWarningMessage();

	// Delete arrays
	delete [] hatDiagonals;
	for( int iChan = 0; iChan < numChannels; ++iChan ){
		delete [] ftvalForJackknife[iChan];
	}
	delete [] ftvalForJackknife;
	for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
		delete [] resp[iVar];
	}
	delete [] resp;
	delete [] variancesWithoutScale;
	delete [] coherences;

	// Calculate & output error bars
	for( int iOut = 0; iOut < numOfOutputVariables; ++iOut ){
		if( numOfSubsets > 2 ){
			std::complex<double> avgResp0(0.0, 0.0);
			std::complex<double> avgResp1(0.0, 0.0);
			for( int iSubset = 0; iSubset < numOfSubsets; ++iSubset ){
				avgResp0 += pseudoResp0[iSubset][iOut];
				avgResp1 += pseudoResp1[iSubset][iOut];
			}
			const double factor = 1.0 / static_cast<double>(numOfSubsets);
			avgResp0 *= factor;
			avgResp1 *= factor;
			double variance0(0.0);
			double variance1(0.0);
			for( int iSubset = 0; iSubset < numOfSubsets; ++iSubset ){
				variance0 += std::norm( pseudoResp0[iSubset][iOut] - avgResp0 );
				variance1 += std::norm( pseudoResp1[iSubset][iOut] - avgResp1 );
			}
			const double factor2 = factor / static_cast<double>(2 * numOfSubsets - 4);	
			variance0 *= factor2;
			variance1 *= factor2;
			respErr[iOut][0] = sqrt(variance0);
			respErr[iOut][1] = sqrt(variance1);
		}else{
			respErr[iOut][0] = 1.0e10;
			respErr[iOut][1] = 1.0e10;
		}
	}
	// Release memory
	for( int iSeg = 0; iSeg < numOfSubsets; ++iSeg ){
		delete [] pseudoResp0[iSeg];
		delete [] pseudoResp1[iSeg];
	}
	delete [] pseudoResp0;
	delete [] pseudoResp1;


}

// Estimate error by a parametric approach
void AnalysisMultivariateRegression::estimateErrorParametric(const int numSegmentsTotal, const double paramB, const double paramC,
	std::complex<double>** ftval, std::complex<double>** respOrg, const double* const variancesWithoutScaleOrg, const double scaleOrg, 
	double** respErr) const {

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
	ptrOutputFiles->writeCvgAndLogMessage("Estimate errors by a parametric approach");
	ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");

	const Control* const ptrControl = Control::getInstance();
	const int numOfOutputVariables = ptrControl->getNumOutputVariables();
	const int numOfInputVariables = ptrControl->getNumInputVariables();
	const int numOfOutputAndInputVariables = numOfOutputVariables + numOfInputVariables;
	const int numOfReferenceVariables = ptrControl->getNumRemoteReferenceVariables();

	if (numSegmentsTotal <= 2) {
		for (int iOut = 0; iOut < numOfOutputVariables; ++iOut) {
			respErr[iOut][0] = 1.0e10;
			respErr[iOut][1] = 1.0e10;
		}
		return;
	}

	assert(numOfReferenceVariables == 2);
	const int rr0 = ptrControl->getChannelIndex(CommonParameters::REMOTE_REFERENCE, 0);
	const int rr1 = ptrControl->getChannelIndex(CommonParameters::REMOTE_REFERENCE, 1);

	std::complex<double>** complexResidualsOrg = new std::complex<double>*[numOfOutputAndInputVariables];
	for (int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar) {
		complexResidualsOrg[iVar] = new std::complex<double>[numSegmentsTotal];
	}
	// Calculate complex residuals
	calculateComplexResiduals(numSegmentsTotal, ftval, respOrg, complexResidualsOrg);
	double* MDOrg = new double[numSegmentsTotal];
	// Calculate Mahalanobis distance
	calculateMD(numSegmentsTotal, numOfOutputAndInputVariables, complexResidualsOrg,
		variancesWithoutScaleOrg, MDOrg);
	double sumOfSquaredDerivativeOfLossFunction(0.0);
	double sumOfSecondDerivativeOfLossFunction(0.0);
	double sumfOfBrxNorm(0.0);
	double sumfOfBryNorm(0.0);
	for (int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg) {
		const double val = MDOrg[iSeg] / scaleOrg;
		const double influenceFunction = RobustWeightTukeysBiweights::calculateWeights(val, paramC) * val;
		sumOfSquaredDerivativeOfLossFunction += pow(influenceFunction, 2);
		sumOfSecondDerivativeOfLossFunction += RobustWeightTukeysBiweights::calculateSecondDerivativeOfLossFunction(val, paramC);
		sumfOfBrxNorm += std::norm(ftval[rr0][iSeg]);
		sumfOfBryNorm += std::norm(ftval[rr1][iSeg]);
	}
	double* variance = new double[numOfOutputAndInputVariables];
	for (int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar) {
		variance[iVar] = 0.0;
		for (int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg) {
			variance[iVar] += std::norm(complexResidualsOrg[iVar][iSeg]);
		}
		variance[iVar] /= static_cast<double>(2 * numSegmentsTotal - 4);
	}
	const double numerator = sumOfSquaredDerivativeOfLossFunction / static_cast<double>(numSegmentsTotal);
	double denominator = pow( sumOfSecondDerivativeOfLossFunction / static_cast<double>(numSegmentsTotal), 2 );
	if (denominator < 1.0e-10) {
		denominator = 1.0e-10;
	}
	if (sumfOfBrxNorm < 1.0e-10) {
		sumfOfBrxNorm = 1.0e-10;
	}
	if (sumfOfBryNorm < 1.0e-10) {
		sumfOfBryNorm = 1.0e-10;
	}

	double** istfErr = new double* [numOfOutputAndInputVariables];
	for (int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar) {
		istfErr[iVar] = new double[2];
		istfErr[iVar][0] = variance[iVar] * numerator / denominator / sumfOfBrxNorm;
		istfErr[iVar][1] = variance[iVar] * numerator / denominator / sumfOfBryNorm;
	}

	const std::complex<double> czero = std::complex<double>(0.0, 0.0);
	assert(numOfInputVariables == 2);
	const int in0 = ptrControl->getChannelIndex(CommonParameters::INPUT, 0);
	const int in1 = ptrControl->getChannelIndex(CommonParameters::INPUT, 1);
	const std::complex<double> Txx = respOrg[in0][0];
	const std::complex<double> Txy = respOrg[in0][1];
	const std::complex<double> Tyx = respOrg[in1][0];
	const std::complex<double> Tyy = respOrg[in1][1];
	std::complex<double> TMatrix[2][2] = { Txx, Txy, Tyx, Tyy };
#ifdef _DEBUG_WRITE
	std::cout << "[";
	for (int row = 0; row < numOfInputVariables; ++row) {
		for (int col = 0; col < numOfReferenceVariables; ++col) {
			std::cout << TMatrix[row][col].real() << "+" << TMatrix[row][col].imag() << "im ";
		}
		if (row + 1 < numOfInputVariables) {
			std::cout << ";";
		}
	}
	std::cout << "]" << std::endl;
#endif	
	const std::complex<double> det = TMatrix[0][0] * TMatrix[1][1] - TMatrix[1][0] * TMatrix[0][1];
	const std::complex<double> TInvMatrix[2][2] = { TMatrix[1][1] / det, -TMatrix[0][1] / det, -TMatrix[1][0] / det, TMatrix[0][0] / det };
	const std::complex<double> TInvTMatrix[2][2] = { TInvMatrix[0][0], TInvMatrix[1][0], TInvMatrix[0][1], TInvMatrix[1][1] };
#ifdef _DEBUG_WRITE
	std::cout << "[";
	for (int row = 0; row < numOfReferenceVariables; ++row) {
		for (int col = 0; col < numOfInputVariables; ++col) {
			std::cout << TInvTMatrix[row][col].real() << "+" << TInvTMatrix[row][col].imag() << "im ";
		}
		if (row + 1 < numOfReferenceVariables) {
			std::cout << ";";
		}
	}
	std::cout << "]" << std::endl;
#endif	
	std::complex<double>** TInvTTInvMatrix = new std::complex<double>*[numOfInputVariables * numOfReferenceVariables];
	for (int irow = 0; irow < numOfInputVariables * numOfReferenceVariables; ++irow) {
		TInvTTInvMatrix[irow] = new std::complex<double>[numOfInputVariables * numOfReferenceVariables];
		for (int icol = 0; icol < numOfInputVariables * numOfReferenceVariables; ++icol) {
			TInvTTInvMatrix[irow][icol] = czero;// Zero clear
		}
	}
	for (int irow = 0; irow < numOfReferenceVariables; ++irow) {
		for (int icol = 0; icol < numOfInputVariables; ++icol) {
			const std::complex<double> factor = TInvTMatrix[irow][icol];
			for (int irow2 = 0; irow2 < numOfReferenceVariables; ++irow2) {
				for (int icol2 = 0; icol2 < numOfInputVariables; ++icol2) {
					const int irowOut = irow2 + irow * numOfReferenceVariables;
					const int icolOut = icol2 + icol * numOfInputVariables;
					TInvTTInvMatrix[irowOut][icolOut] = factor * TInvMatrix[irow2][icol2];
				}
			}
		}
	}
#ifdef _DEBUG_WRITE
	std::cout << "[";
	for (int row = 0; row < numOfInputVariables * numOfReferenceVariables; ++row) {
		for (int col = 0; col < numOfInputVariables * numOfReferenceVariables; ++col) {
			std::cout << TInvTTInvMatrix[row][col].real() << "+" << TInvTTInvMatrix[row][col].imag() << "im ";
		}
		if (row + 1 < numOfInputVariables * numOfReferenceVariables) {
			std::cout << ";";
		}
	}
	std::cout << "]" << std::endl;
#endif

	for (int iOut = 0; iOut < numOfOutputVariables; ++iOut) {
		const int out = ptrControl->getChannelIndex(CommonParameters::OUTPUT, iOut);
		std::complex<double> UMatrix[2] = { respOrg[out][0], respOrg[out][1] };
#ifdef _DEBUG_WRITE
		std::cout << "[";
		for (int col = 0; col < numOfReferenceVariables; ++col) {
			std::cout << UMatrix[col].real() << "+" << UMatrix[col].imag() << "im ";
		}
		std::cout << "]" << std::endl;
#endif
		std::complex<double>** IUMatrix = new std::complex<double>*[numOfInputVariables];
		for (int irow = 0; irow < numOfInputVariables; ++irow) {
			IUMatrix[irow] = new std::complex<double>[numOfInputVariables * numOfReferenceVariables];
			for (int icol = 0; icol < numOfInputVariables * numOfReferenceVariables; ++icol) {
				IUMatrix[irow][icol] = czero;// Zero clear
			}
		}
		for (int irow = 0; irow < numOfInputVariables; ++irow) {
			for (int icol = 0; icol < numOfInputVariables; ++icol) {
				const double factor = irow == icol ? 1.0 : 0.0;
				for (int icol2 = 0; icol2 < numOfReferenceVariables; ++icol2) {
					const int irowOut = irow;
					const int icolOut = icol2 + icol * numOfReferenceVariables;
					IUMatrix[irowOut][icolOut] = factor * UMatrix[icol2];
				}
			}
		}
#ifdef _DEBUG_WRITE
		std::cout << "[";
		for (int row = 0; row < numOfInputVariables; ++row) {
			for (int col = 0; col < numOfInputVariables * numOfReferenceVariables; ++col) {
				std::cout << IUMatrix[row][col].real() << "+" << IUMatrix[row][col].imag() << "im ";
			}
			if (row + 1 < numOfInputVariables) {
				std::cout << ";";
			}
		}
		std::cout << "]" << std::endl;
#endif
		std::complex<double>** IUTInvTTInvMatrix = new std::complex<double>*[numOfInputVariables];
		for (int irow = 0; irow < numOfInputVariables; ++irow) {
			IUTInvTTInvMatrix[irow] = new std::complex<double>[numOfInputVariables * numOfReferenceVariables];
			for (int icol = 0; icol < numOfInputVariables * numOfReferenceVariables; ++icol) {
				IUTInvTTInvMatrix[irow][icol] = czero;// Zero clear
			}
		}
		for (int irow = 0; irow < numOfInputVariables; ++irow) {
			for (int icol = 0; icol < numOfInputVariables * numOfReferenceVariables; ++icol) {
				std::complex<double> value = czero;
				for (int i = 0; i < numOfInputVariables * numOfReferenceVariables; ++i) {
					value += IUMatrix[irow][i] * TInvTTInvMatrix[i][icol];
				}
				IUTInvTTInvMatrix[irow][icol] = value;
			}
		}
		for (int irow = 0; irow < numOfInputVariables; ++irow) {
			delete[] IUMatrix[irow];
		}
		delete[] IUMatrix;
#ifdef _DEBUG_WRITE
		std::cout << "[";
		for (int row = 0; row < numOfInputVariables; ++row) {
			for (int col = 0; col < numOfInputVariables * numOfReferenceVariables; ++col) {
				std::cout << IUTInvTTInvMatrix[row][col].real() << "+" << IUTInvTTInvMatrix[row][col].imag() << "im ";
			}
			if (row + 1 < numOfInputVariables) {
				std::cout << ";";
			}
		}
		std::cout << "]" << std::endl;
#endif
		std::vector< std::complex<double> > delta[2];
		const std::complex<double> imaginaryUnit = std::complex<double>(0.0, 1.0);
		for (int i = 0; i < numOfInputVariables; ++i) {
			for (int row = 0; row < numOfInputVariables; ++row) {
				delta[i].push_back(TInvTMatrix[i][row] * istfErr[out][row]);
				delta[i].push_back(TInvTMatrix[i][row] * istfErr[out][row] * imaginaryUnit);
			}
			for (int row = 0; row < numOfInputVariables; ++row) {
				const int in = ptrControl->getChannelIndex(CommonParameters::INPUT, row);
				for (int col = 0; col < numOfReferenceVariables; ++col) {
					const int index = numOfInputVariables * col + row;
					delta[i].push_back( - IUTInvTTInvMatrix[i][index] * istfErr[in][col]);
					delta[i].push_back( - IUTInvTTInvMatrix[i][index] * istfErr[in][col] * imaginaryUnit);
				}
			}
		}
#ifdef _DEBUG_WRITE
		const std::complex<double> Ux = respOrg[out][0];
		const std::complex<double> Uy = respOrg[out][1];
		const std::complex<double> det = Txx * Tyy - Txy * Tyx;
		const std::complex<double> resp0 = (Ux * Tyy - Uy * Tyx) / det;
		const std::complex<double> resp1 = (Uy * Txx - Ux * Txy) / det;
		for (int i = 0; i < 3; ++i) {
			for (int j = 0; j < numOfReferenceVariables; ++j) {
				for (int isign = 0; isign < 2; ++isign){
					std::complex<double>** respMod = new std::complex<double>*[numOfOutputAndInputVariables];
					for (int i2 = 0; i2 < numOfOutputAndInputVariables; ++i2) {
						respMod[i2] = new std::complex<double>[numOfReferenceVariables];
						for (int j2 = 0; j2 < numOfReferenceVariables; ++j2) {
							respMod[i2][j2] = respOrg[i2][j2];
						}
					}
					const double deltaRatio = 0.001;
					double deltaValue(0.0);
					int index = i;
					if (i == 1) {
						index = in0;
					}
					else if (i == 2) {
						index = in1;
					}
					if (isign == 0) {
						deltaValue = respOrg[index][j].real() * deltaRatio;
						respMod[index][j] += deltaValue;
					}
					else {
						deltaValue = respOrg[index][j].imag() * deltaRatio;
						respMod[index][j] += imaginaryUnit * deltaValue;
					}
					// Slighty modify
					const std::complex<double> UxMod = respMod[out][0];
					const std::complex<double> UyMod = respMod[out][1];
					const std::complex<double> TxxMod = respMod[in0][0];
					const std::complex<double> TxyMod = respMod[in0][1];
					const std::complex<double> TyxMod = respMod[in1][0];
					const std::complex<double> TyyMod = respMod[in1][1];
					const std::complex<double> detMod = TxxMod * TyyMod - TxyMod * TyxMod;
					const std::complex<double> resp0Mod = (UxMod * TyyMod - UyMod * TyxMod) / detMod;
					const std::complex<double> resp1Mod = (UyMod * TxxMod - UxMod * TxyMod) / detMod;
					// Calculate derivative
					std::complex<double> derivative0(0.0, 0.0);
					std::complex<double> derivative1(0.0, 0.0);
					derivative0 = (resp0Mod - resp0) / deltaValue;
					derivative1 = (resp1Mod - resp1) / deltaValue;
					std::cout << i << " " << j << " " << isign
						<< " " << derivative0.real() << " " << derivative0.imag() 
						<< " " << derivative1.real() << " " << derivative1.imag() << std::endl;
					for (int i2 = 0; i2 < numOfOutputAndInputVariables; ++i2) {
						delete[] respMod[i2];
					}
					delete[] respMod;
				}
			}
		}

#endif
		for (int irow = 0; irow < numOfInputVariables; ++irow) {
			delete[] IUTInvTTInvMatrix[irow];
		}
		delete[] IUTInvTTInvMatrix;
		for (int i = 0; i < numOfInputVariables; ++i) {
			double squareSumReal(0.0);
			double squareSumImag(0.0);
			for (std::vector< std::complex<double> >::const_iterator itr = delta[i].begin(); itr != delta[i].end(); ++itr) {
				squareSumReal += pow(itr->real(), 2);
				squareSumImag += pow(itr->imag(), 2);
			}
			respErr[out][i] = sqrt(squareSumReal);
		}
	}

	// Release memory
	for (int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar) {
		delete[] complexResidualsOrg[iVar];
	}
	delete[] complexResidualsOrg;
	delete[] MDOrg;
	delete[] variance;
	for (int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar) {
		delete[] istfErr[iVar];
	}
	delete[] istfErr;
	for (int irow = 0; irow < numOfInputVariables * numOfReferenceVariables; ++irow) {
		delete[] TInvTTInvMatrix[irow];
	}
	delete[] TInvTTInvMatrix;


}

// Write residuals
void AnalysisMultivariateRegression::writeResiduals( const std::string& fileName, const int numSegmentsTotal, 
	const int numOfOutputAndInputVariables, const std::vector< std::pair<std::string, std::string> >& times,
	std::complex<double>** complexResiduals, const double* const MD, const double* const weights ) const{

	std::ofstream ofs;
	ofs.open( fileName.c_str(), std::ios::out );
	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	if( ofs.fail() ){
		ptrOutputFiles->writeErrorMessage( "File open error : " + fileName );
	}
	ofs << "index";
	ofs << ",start_time,end_time";
	for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
		ofs << ",residual_real_var_" << iVar;
		ofs << ",residual_imag_var_" << iVar;
	}
	ofs << ",MD";
	ofs << ",weight";
	ofs << std::endl;
	int index(0);
	for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg, ++index ){
		ofs << index;
		const std::string timeStart = times[iSeg].first;
		const std::string timeEnd = times[iSeg].second;
		ofs << "," << timeStart << "," << timeEnd;
		for( int iVar = 0; iVar < numOfOutputAndInputVariables; ++iVar ){
			ofs << "," << std::setprecision(10) << std::scientific << complexResiduals[iVar][iSeg].real();
			ofs << "," << std::setprecision(10) << std::scientific << complexResiduals[iVar][iSeg].imag();
		}
		ofs << "," << std::setprecision(10) << std::scientific << MD[iSeg];
		ofs << "," << std::setprecision(10) << std::scientific << weights[iSeg];
		ofs << std::endl;
	}
	ofs.close();

}
