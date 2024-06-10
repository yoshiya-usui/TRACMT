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
#include "RobustPrewhitening.h"
#include "OutputFiles.h"
#include "Control.h"
#include "Util.h"
#include "UtilRobust.h"
#include "DoubleDenseSquareSymmetricPositiveDefiniteMatrix.h"

#include <assert.h>
#include <iostream>
#include <iomanip>

// Default constructer
RobustPrewhitening::RobustPrewhitening()
{
}

// Destructer
RobustPrewhitening::~RobustPrewhitening()
{
}

// Return the instance of the class
RobustPrewhitening* RobustPrewhitening::getInstance(){
   	static RobustPrewhitening instance;// The only instance
  	return &instance;
}

// Perform robust prewhitening
void RobustPrewhitening::robustPrewhitening( std::vector<CommonParameters::DataFileSet>& dataFileSets, std::vector<double>* coeffsAROutput ) const{

	const Control* const ptrControl = Control::getInstance();
	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	const Control::ParamsForPrewhitening params = ptrControl->getParamsForPrewhitening();
	const int maxDegreesOfAR = params.maxDegreeOfARModel;
	const Control::ParamsForRobustFilter paramsForFilter = ptrControl->getParamsForRobustFilter();

	int numDataPointsAll(0);
	int iSection(0);
	for( std::vector<CommonParameters::DataFileSet>::const_iterator itr = dataFileSets.begin(); itr != dataFileSets.end(); ++itr, ++iSection ){
		numDataPointsAll += itr->numDataPoints;
	}
	const int numSections = iSection;

	// Calculate candidates
	const int numCandidates = params.numCandidatesOfPartialAutocorrelationFunction;
	if( numCandidates % 2 == 0 ){
		ptrOutputFiles->writeErrorMessage("Numbef of candidates of partial autocorrelation function should be an odd number");
	}
	double* candidates = new double[numCandidates];
	for( int iCan = 0; iCan < numCandidates; ++iCan ){
		const int index = iCan % 2 == 1 ? iCan / 2 + 1 : -(iCan / 2);
		candidates[iCan] = static_cast<double>(index) / static_cast<double>(numCandidates - 1) * 2.0;
#ifdef _DEBUG_WRITE
		std::cout << candidates[iCan] << std::endl;
#endif
	}

	const int numChannels = ptrControl->getNumberOfChannels();
	for( int iChan = 0; iChan < numChannels; ++iChan ){
		ptrOutputFiles->writeLogMessage("Perform prewhitening for channel " + Util::toString(iChan));
		double* coeffsAR= new double[maxDegreesOfAR];
		for( int iDeg = 0; iDeg < maxDegreesOfAR; ++iDeg ){
			coeffsAR[iDeg] = 0.0;
		}
		// Copy data to temporal arrays
		const int numOfDataSets = static_cast<int>(dataFileSets.size());
		int* numOfData = new int[numOfDataSets];;
		double** dataOrg = new double*[numOfDataSets];
		double** dataMod = new double*[numOfDataSets];
		int iSet(0);
		for( std::vector<CommonParameters::DataFileSet>::const_iterator itr = dataFileSets.begin(); itr != dataFileSets.end(); ++itr, ++iSet ){
			const int numData = itr->numDataPoints;
			numOfData[iSet] = numData;
			dataOrg[iSet] = new double[numData];
			dataMod[iSet] = new double[numData];
			for( int iData = 0; iData < numData; ++iData ){
				dataOrg[iSet][iData] = itr->dataFile[iChan].data[iData];
				dataMod[iSet][iData] = dataOrg[iSet][iData]; 
			}
		}
		DoubleDenseSquareSymmetricMatrix autoCovarianceMatrixMaxDegreesOfAR;
		if(paramsForFilter.applyRobustFilter){
			// Calculate auto-covariance matrix
			ptrOutputFiles->writeCvgAndLogMessage("Calculate auto-covariance matrix of channel " + Util::toString(iChan));
			calculateRobustAutoCovarianceMatrix( maxDegreesOfAR, numOfDataSets, numOfData, dataOrg, autoCovarianceMatrixMaxDegreesOfAR );
		}
		assert( iSet == numOfDataSets );
		double minAIC(1.0e20);
		int iDegAROfMinAIC(-1);
		double scaleOfMinAIC(-1.0);
		double* coeffsAROfMinAIC = new double[maxDegreesOfAR];
		for( int iDegAR = 1; iDegAR <= maxDegreesOfAR; ++iDegAR ){
			const int numDataForSEstimator = numDataPointsAll - iDegAR * numSections;
			// Calculate forward MMSE prediction residual and backward MMSE prediction residual
			double* forwardMMSEResidual = new double[numDataForSEstimator];
			double* backwardMMSEResidual = new double[numDataForSEstimator];
			int icount(0);
			for( int iSet = 0; iSet < numOfDataSets; ++iSet ){
				const int numData = numOfData[iSet];
				for( int iData = iDegAR; iData < numData; ++iData ){
					forwardMMSEResidual[icount] = dataMod[iSet][iData];
					for( int i = 0; i < iDegAR - 1; ++i ){
						const int index = iData - 1 - i;
						forwardMMSEResidual[icount] -= coeffsAR[i] * dataMod[iSet][index];
					}
					backwardMMSEResidual[icount] = dataMod[iSet][iData - iDegAR];
					for( int i = 0; i < iDegAR - 1; ++i ){
						const int index = iData - iDegAR + 1 + i;
						backwardMMSEResidual[icount] -= coeffsAR[i] * dataMod[iSet][index];
					}
					++icount;
				}
			}
			assert(icount == numDataForSEstimator);
			// Calculate AR coefficients
			double scale(0.0);
			double partialAutocorrelationFunction(0.0);
			double AIC(0.0);
			if( params.typeOfEstimator == Control::USE_S_ESTIMATOR_FOR_PREWHITENING ){
				UtilRobust::computeSEstimatorForUnivariateLinearRegression( numDataForSEstimator, numCandidates,
					forwardMMSEResidual, backwardMMSEResidual, candidates, -1.0, 1.0, scale, partialAutocorrelationFunction, AIC );
				ptrOutputFiles->writeLogMessage("Degree of AR model: " + Util::toString(iDegAR) + ", Sigma: " + Util::toString(scale)
					+ ", AICS: " + Util::toString(AIC));
			}else if( params.typeOfEstimator == Control::USE_LEAST_SQUARE_ESTIMATOR_FOR_PREWHITENING ){
				UtilRobust::computeLSEstimatorForUnivariateLinearRegression( numDataForSEstimator, forwardMMSEResidual, backwardMMSEResidual,
					-1.0, 1.0, scale, partialAutocorrelationFunction );
				AIC = 2.0 * static_cast<double>(numDataForSEstimator) * log(scale) + 2.0 * ( static_cast<double>(iDegAR) + 1.0 ) 
					+ static_cast<double>(numDataForSEstimator) + static_cast<double>(numDataForSEstimator) * log(2.0 * CommonParameters::PI);
				ptrOutputFiles->writeLogMessage("Degree of AR model: " + Util::toString(iDegAR) + ", Sigma: " + Util::toString(scale)
					+ ", AIC: " + Util::toString(AIC));
			}else{
				ptrOutputFiles->writeErrorMessage("Wrong type of estimator for prewhitening: " + Util::toString(params.typeOfEstimator));
			}
			delete [] forwardMMSEResidual;
			delete [] backwardMMSEResidual;
			coeffsAR[iDegAR - 1] = partialAutocorrelationFunction;
			if( iDegAR > 1 ){
				double* coeffsARPre = new double[iDegAR - 1];
				for( int i = 0; i < iDegAR - 1; ++i ){
					coeffsARPre[i] = coeffsAR[i];
				}
				for( int i = 0; i < iDegAR - 1; ++i ){
					const int index = iDegAR - 2 - i;
					coeffsAR[i] -= partialAutocorrelationFunction * coeffsARPre[index];
				}
#ifdef _DEBUG_WRITE
				for( int i = 0; i < iDegAR - 1; ++i ){
					std::cout << "AR old new: " << coeffsARPre[i] << " " << coeffsAR[i] << std::endl;
				}
#endif
				delete [] coeffsARPre;
			}
#ifdef _DEBUG_WRITE
			for( int i = 0; i < iDegAR; ++i ){
				std::cout << "AR coeffs: " << coeffsAR[i] << std::endl;
			}
#endif
			if( AIC < minAIC ){
				minAIC = AIC;
				iDegAROfMinAIC = iDegAR;
				scaleOfMinAIC = scale;
				for( int i = 0; i < iDegAR; ++i ){
					coeffsAROfMinAIC[i] = coeffsAR[i];
				}
				for( int i = iDegAR; i < maxDegreesOfAR; ++i ){
					coeffsAROfMinAIC[i] = 0.0;
				}
			}
			if(paramsForFilter.applyRobustFilter){
				// Calculate robust-filtered value
				double* autoCovariance = new double[iDegAR * iDegAR];
				// Copy matrix
				icount = 0;
				for( int col = 0; col < iDegAR; ++col ){
					for( int row = 0; row < iDegAR; ++row ){
						// Column major
						autoCovariance[icount] = autoCovarianceMatrixMaxDegreesOfAR.getValue(row, col);
						++icount;
					}
				}
#ifdef _DEBUG_WRITE
				Util::debugWriteRealMatrix(iDegAR, iDegAR, autoCovariance);
#endif
				int iSet(0);
				for( std::vector<CommonParameters::DataFileSet>::const_iterator itr = dataFileSets.begin(); itr != dataFileSets.end(); ++itr, ++iSet ){
					const int numData = itr->numDataPoints;
					calculateRobustFilteredValue( iChan, numData, iDegAR, coeffsAR, scale, dataOrg[iSet], autoCovariance, dataMod[iSet] );
				}
				delete [] autoCovariance;
			}else{
				int iSet = 0;
				for( std::vector<CommonParameters::DataFileSet>::const_iterator itr = dataFileSets.begin(); itr != dataFileSets.end(); ++itr, ++iSet ){
					const int numData = itr->numDataPoints;
					for( int iData = 0; iData < numData; ++iData ){
						dataMod[iSet][iData] = dataOrg[iSet][iData];
					}
				}
			}
		}
		delete [] coeffsAR;
		if( params.typeOfEstimator == Control::USE_S_ESTIMATOR_FOR_PREWHITENING ){
			ptrOutputFiles->writeLogMessage("The AR model of " + Util::toString(iDegAROfMinAIC) + " degress gives the minimum AICS (" + Util::toString(minAIC) + ")");
		}else if( params.typeOfEstimator == Control::USE_LEAST_SQUARE_ESTIMATOR_FOR_PREWHITENING ){
			ptrOutputFiles->writeLogMessage("The AR model of " + Util::toString(iDegAROfMinAIC) + " degress gives the minimum AIC (" + Util::toString(minAIC) + ")");
		}
		for( int i = 0; i < iDegAROfMinAIC; ++i ){
			coeffsAROutput[iChan].push_back(coeffsAROfMinAIC[i]);
		}
		std::ostringstream msg;
		msg << "AR coefficients: ";
		for( int i = 0; i < iDegAROfMinAIC; ++i ){
			msg << coeffsAROfMinAIC[i] << " ";
		}
		ptrOutputFiles->writeLogMessage(msg.str());
		if( paramsForFilter.applyRobustFilter && paramsForFilter.replaceTimeSeriesWithFilteredValues ){
			ptrOutputFiles->writeLogMessage("Replace observed time-series with filtered values");
			double* autoCovariance = new double[iDegAROfMinAIC * iDegAROfMinAIC];
			// Copy matrix
			int icount = 0;
			for( int col = 0; col < iDegAROfMinAIC; ++col ){
				for( int row = 0; row < iDegAROfMinAIC; ++row ){
					// Column major
					autoCovariance[icount] = autoCovarianceMatrixMaxDegreesOfAR.getValue(row, col);
					++icount;
				}
			}
#ifdef _DEBUG_WRITE
			Util::debugWriteRealMatrix(iDegAROfMinAIC, iDegAROfMinAIC, autoCovariance);
#endif
			// Calculate robust-filtered value
			int iSet = 0;
			for( std::vector<CommonParameters::DataFileSet>::const_iterator itr = dataFileSets.begin(); itr != dataFileSets.end(); ++itr, ++iSet ){
				const int numData = itr->numDataPoints;
				// Calculate robust-filtered value
				calculateRobustFilteredValue( iChan, numData, iDegAROfMinAIC, coeffsAROfMinAIC, scaleOfMinAIC, dataOrg[iSet], autoCovariance, dataMod[iSet] );
				// Output robust-filtered time serieas
				if( ptrControl->doesOutputTimeSeriesToCsv() ){
					std::ostringstream oss;
					oss << "time_series_robust_filtered_sect_" << iSet << "_chan_" << iChan << ".csv"; 
					std::ofstream ofs;
					ofs.open( oss.str().c_str(), std::ios::out );
					if( ofs.fail() ){
						ptrOutputFiles->writeLogMessage("File open error !! : " + oss.str());
					}
					for( int iData = 0; iData < numData; ++iData ){
						ofs << std::setprecision(12) << std::scientific << dataMod[iSet][iData] << std::endl;
					}
					ofs.close();
				}
			}
			delete [] autoCovariance;
		}else{
			int iSet = 0;
			for( std::vector<CommonParameters::DataFileSet>::const_iterator itr = dataFileSets.begin(); itr != dataFileSets.end(); ++itr, ++iSet ){
				const int numData = itr->numDataPoints;
				for( int iData = 0; iData < numData; ++iData ){
					dataMod[iSet][iData] = dataOrg[iSet][iData];
				}
			}
		}
		// Prewhitening
		iSet = 0;
		for( std::vector<CommonParameters::DataFileSet>::iterator itr = dataFileSets.begin(); itr != dataFileSets.end(); ++itr, ++iSet ){
			const int numData = itr->numDataPoints;
			// Replace observed time-series with filtered values
			for( int iData = 0; iData < iDegAROfMinAIC; ++iData ){
				itr->dataFile[iChan].data[iData] = 0.0;
			}
			for( int iData = iDegAROfMinAIC; iData < numData; ++iData ){
				double work = dataMod[iSet][iData];
				for( int iAR = 0; iAR < iDegAROfMinAIC; ++iAR ){
					const int index = iData - iAR - 1;
					work -= coeffsAROfMinAIC[iAR] * dataMod[iSet][index];
				}
				itr->dataFile[iChan].data[iData] = work;
			}
		}
		// Release memory
		delete [] coeffsAROfMinAIC;
		delete [] numOfData;
		for( int iSet = 0; iSet < numOfDataSets; ++iSet ){
			delete [] dataOrg[iSet];
			delete [] dataMod[iSet];
		}
		delete [] dataOrg;
		delete [] dataMod;
	}

	delete [] candidates;

}

// Perform prewhitening using user-defined AR coefficients
void RobustPrewhitening::prewhiteningUsingUserDefinedARCoeffs( std::vector<CommonParameters::DataFileSet>& dataFileSets, std::vector<double>* coeffsAROutput ) const{

	const Control* const ptrControl = Control::getInstance();
	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();

	const std::string fileName = "AR_coefficients.dat";
	std::ifstream ifs( fileName.c_str(), std::ios::in );
	if( ifs.fail() ){
		ptrOutputFiles->writeErrorMessage("File open error : " + fileName);
	}

	const int numChannels = ptrControl->getNumberOfChannels();
	for( int iChan = 0; iChan < numChannels; ++iChan ){
		ptrOutputFiles->writeLogMessage("Perform prewhitening for channel " + Util::toString(iChan));
		// Copy data to temporal arrays
		const int numOfDataSets = static_cast<int>(dataFileSets.size());
		double** dataOrg = new double*[numOfDataSets];
		int iSet(0);
		for( std::vector<CommonParameters::DataFileSet>::const_iterator itr = dataFileSets.begin(); itr != dataFileSets.end(); ++itr, ++iSet ){
			const int numData = itr->numDataPoints;
			dataOrg[iSet] = new double[numData];
			for( int iData = 0; iData < numData; ++iData ){
				dataOrg[iSet][iData] = itr->dataFile[iChan].data[iData];
			}
		}
		int numCoeffs(0);
		ifs >> numCoeffs;
		for( int i = 0; i < numCoeffs; ++i ){
			double dbuf(0.0);
			ifs >> dbuf;
			coeffsAROutput[iChan].push_back(dbuf);
		}
		std::ostringstream msg;
		msg << "AR coefficients: ";
		for( int i = 0; i < numCoeffs; ++i ){
			msg << coeffsAROutput[iChan][i] << " ";
		}
		ptrOutputFiles->writeLogMessage(msg.str());
		// Prewhitening
		iSet = 0;
		for( std::vector<CommonParameters::DataFileSet>::iterator itr = dataFileSets.begin(); itr != dataFileSets.end(); ++itr, ++iSet ){
			const int numData = itr->numDataPoints;
			// Replace observed time-series with filtered values
			for( int iData = 0; iData < numCoeffs; ++iData ){
				itr->dataFile[iChan].data[iData] = 0.0;
			}
			for( int iData = numCoeffs; iData < numData; ++iData ){
				double work = dataOrg[iSet][iData];
				for( int iAR = 0; iAR < numCoeffs; ++iAR ){
					const int index = iData - iAR - 1;
					work -= coeffsAROutput[iChan][iAR] * dataOrg[iSet][index];
				}
				itr->dataFile[iChan].data[iData] = work;
			}
		}
		// Release memory
		for( int iSet = 0; iSet < numOfDataSets; ++iSet ){
			delete [] dataOrg[iSet];
		}
		delete [] dataOrg;
	}

	ifs.close();

}

// Calculate robust-filtered value
void RobustPrewhitening::calculateRobustFilteredValue( const int iChan, const int numOfData, const int degreesOfAR, const double* const coeffsOfAR,
	double sigma, const double* const yOrg, const double* const autoCovariance, double* yMod ) const{

	const Control* const ptrControl = Control::getInstance();
	const Control::ParamsForRobustFilter params = ptrControl->getParamsForRobustFilter();

	const int maxNumOfConsecutiveReplacements = params.thresholds[iChan].maxNumOfConsecutiveReplacements;
	const double paramA = params.thresholds[iChan].firstThresholdFactor;
	const double paramB = params.thresholds[iChan].secondThresholdFactor;
	const bool doesGackToPreviousData = maxNumOfConsecutiveReplacements > 0;

	// Initial values of state vector and covariance matrix
	const double median = Util::calculateMedian( degreesOfAR, yOrg );
	double* xHat = new double[degreesOfAR];
	double* xHatPre = new double[degreesOfAR];
	for( int iDeg = 0; iDeg < degreesOfAR; ++iDeg ){
		xHat[iDeg] = median;
		xHatPre[iDeg] = xHat[iDeg];
	}
	double* covP = new double[degreesOfAR * degreesOfAR];
	double* covPPre = new double[degreesOfAR * degreesOfAR];
	for( int i = 0; i < degreesOfAR * degreesOfAR; ++i ){
		covP[i] = autoCovariance[i];
		covPPre[i] = covP[i];
	}

	for( int iDeg = 0; iDeg < degreesOfAR; ++iDeg ){
		yMod[iDeg] = yOrg[iDeg];
	}	
	double* xTilde = new double[degreesOfAR];
	double* covM = new double[degreesOfAR*degreesOfAR];
	double* workVector = new double[degreesOfAR];
	int numOfConsecutiveReplacements(0);
	int iDataPre = degreesOfAR - 1;
	bool residualAssumedToBeZero(false);
	for( int iData = degreesOfAR; iData < numOfData; ++iData ){
		// x_tilde = phi * x_hat
		xTilde[0] = 0.0;
		for( int iDeg = 0; iDeg < degreesOfAR; ++iDeg ){
			xTilde[0] += coeffsOfAR[iDeg] * xHat[iDeg];
		}
		for( int iDeg = 1; iDeg < degreesOfAR; ++iDeg ){
			xTilde[iDeg] = xHat[iDeg-1];
		}
#ifdef _DEBUG_WRITE
		std::cout << "[";
		for( int iDeg = 0; iDeg < degreesOfAR; ++iDeg ){
			std::cout << coeffsOfAR[iDeg] << " ";
		}
		std::cout << "; ";
		for( int row = 1; row < degreesOfAR; ++row ){
			for( int col = 0; col < degreesOfAR; ++col ){
				const double value = row - 1 == col ? 1.0 : 0.0;
				std::cout << value << " ";
			}
			std::cout << "; ";
		}
		std::cout << "]" << std::endl;
		for( int iDeg = 0; iDeg < degreesOfAR; ++iDeg ){
			std::cout << "iDeg xHat: " << iDeg << " " << xHat[iDeg] << std::endl;
		}
		for( int iDeg = 0; iDeg < degreesOfAR; ++iDeg ){
			std::cout << "iDeg xTilde: " << iDeg << " " << xTilde[iDeg] << std::endl;
		}
		Util::debugWriteRealMatrix(degreesOfAR, degreesOfAR, covP);
#endif
		for( int col = 0; col < degreesOfAR; ++col ){
			double work(0.0);
			for( int row = 0; row < degreesOfAR; ++row ){
				// P matrix is column major
				const int index = col * degreesOfAR + row;
				work += coeffsOfAR[row] * covP[index];
			}
			workVector[col] = work;
		}
#ifdef _DEBUG_WRITE
		for( int iDeg = 0; iDeg < degreesOfAR; ++iDeg ){
			std::cout << "iDeg workVector: " << iDeg << " " << workVector[iDeg] << std::endl;
		}
#endif
		// Zero clear
		for( int col = 0; col < degreesOfAR; ++col ){
			for( int row = 0; row < degreesOfAR; ++row ){
				const int index = col * degreesOfAR + row;
				covM[index] = 0.0;// Zero clear
			}
		}
#ifdef _DEBUG_WRITE
		Util::debugWriteRealMatrix(degreesOfAR, degreesOfAR, covM);
#endif
		// The (1.1) component of M_t+1 matrix
		covM[0] = 0.0;// Zero clear
		for( int i = 0; i < degreesOfAR; ++i ){
			covM[0] += workVector[i] * coeffsOfAR[i];
		}
		if( covM[0] < 0.0 ){
			covM[0] = 0.0; // Zero clear
			// Add (1.1) component of Q matrix
			covM[0] += pow(sigma, 2);
		}else{
			// Add (1.1) component of Q matrix
			covM[0] += pow(sigma, 2);
#ifdef _DEBUG_WRITE
			Util::debugWriteRealMatrix(degreesOfAR, degreesOfAR, covM);
#endif
			// The first column of M_t+1 matrix except (1.1) component
			for( int row = 1; row < degreesOfAR; ++row ){
				covM[row] = workVector[row - 1];
			}
#ifdef _DEBUG_WRITE
			Util::debugWriteRealMatrix(degreesOfAR, degreesOfAR, covM);
#endif
			// The first row of M_t+1 matrix except (1.1) component
			for( int col = 1; col < degreesOfAR; ++col ){
				const int index = col * degreesOfAR;
				covM[index] = workVector[col - 1];
			}
#ifdef _DEBUG_WRITE
			Util::debugWriteRealMatrix(degreesOfAR, degreesOfAR, covM);
#endif
			// The remaing components of M_t+1 matrix
			for( int col = 1; col < degreesOfAR; ++col ){
				for( int row = 1; row < degreesOfAR; ++row ){
					const int indexM = col * degreesOfAR + row;
					const int indexP = (col - 1) * degreesOfAR + (row - 1);
					covM[indexM] = covP[indexP];
				}
			}
		}
#ifdef _DEBUG_WRITE
		Util::debugWriteRealMatrix(degreesOfAR, degreesOfAR, covM);
#endif
		const double mt = std::max( covM[0], 1.0e-20 );
		const double normalizedResidual = ( yOrg[iData] - xTilde[0] ) / sqrt(mt);
		const double weightForStateVector = calculateWeightForStateVectorOfRobustFilter( normalizedResidual, paramA, paramB, residualAssumedToBeZero );
		for( int i = 0; i < degreesOfAR; ++i ){
			// M_t matrix is column major
			xHat[i] = xTilde[i] + covM[i] / sqrt(mt) * weightForStateVector;
		}
		yMod[iData] = xHat[0];
		const double weightForCovariance = calculateWeightForCovarianceMatrixOfRobustFilter( normalizedResidual, paramA, paramB, false );
		// P_t = M_t - chi * m_t * m_t^T  / m_t
		// First, copy M_t matrx to P_t matrx 
		for( int col = 0; col < degreesOfAR; ++col ){
			for( int row = 0; row < degreesOfAR; ++row ){
				// Column major
				const int index = col * degreesOfAR + row;
				covP[index] = covM[index];
			}
		}
#ifdef _DEBUG_WRITE
		Util::debugWriteRealMatrix(degreesOfAR, degreesOfAR, covP);
#endif
		for( int col = 0; col < degreesOfAR; ++col ){
			for( int row = 0; row < degreesOfAR; ++row ){
				// Column major
				const double work = covM[col] * covM[row] / mt * weightForCovariance;
				const int index = col * degreesOfAR + row;
				covP[index] -= work;
			}
		}
#ifdef _DEBUG_WRITE
		Util::debugWriteRealMatrix(degreesOfAR, degreesOfAR, covP);
#endif
		if( Util::calculateDeterminantOfMatrix( degreesOfAR, covP ) <= 0.0 ){
			for( int col = 0; col < degreesOfAR; ++col ){
				for( int row = 0; row < degreesOfAR; ++row ){
					// Column major
					const int index = col * degreesOfAR + row;
					covP[index] = 0.0;
				}
			}
		}
		// Check the number of consecutive replacements
		if( !residualAssumedToBeZero && fabs(yMod[iData] - yOrg[iData]) / std::max(fabs(yOrg[iData]), CommonParameters::EPS) > 1.0e-3 ){
#ifdef _DEBUG_WRITE
			std::cout << "org mod: " << yOrg[iData] << " " << yMod[iData] << std::endl;
#endif
			++numOfConsecutiveReplacements;
		}else{
			numOfConsecutiveReplacements = 0;
			iDataPre = iData;
			for( int iDeg = 0; iDeg < degreesOfAR; ++iDeg ){
				xHatPre[iDeg] = xHat[iDeg];
			}
			for( int i = 0; i < degreesOfAR * degreesOfAR; ++i ){
				covPPre[i] = covP[i];
			}
		}
		residualAssumedToBeZero = false;
		if( doesGackToPreviousData && numOfConsecutiveReplacements >= maxNumOfConsecutiveReplacements ){
			// Restore state vector and covariance matrix
			for( int iDeg = 0; iDeg < degreesOfAR; ++iDeg ){
				xHat[iDeg] = xHatPre[iDeg];
			}
			for( int i = 0; i < degreesOfAR * degreesOfAR; ++i ){
				covP[i] = covPPre[i];
			}
			numOfConsecutiveReplacements = 0;
			residualAssumedToBeZero = true;
			iData = iDataPre;
			// Go back to previous data
		}
	}

	delete [] xHat;
	delete [] xHatPre;
	delete [] xTilde;
	delete [] covP;
	delete [] covPPre;
	delete [] covM;
	delete [] workVector;

}

// Calculate robust auto-covariance matrix
void RobustPrewhitening::calculateRobustAutoCovarianceMatrix( const int degreesOfAR, const int numOfDataSets, 
	const int* const numOfData, double** data, DoubleDenseSquareSymmetricMatrix& covarianceMatrix ) const{

	const Control* const ptrControl = Control::getInstance();
	const Control::ParamsForRobustAutoCovariance params = ptrControl->getParamsForRobustAutoCovariance();

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();

	int numOfDataAll = 0; 
	for( int iSet = 0; iSet < numOfDataSets; ++iSet ){
		numOfDataAll += numOfData[iSet];
	}

	const int numOfSamples = numOfDataAll - degreesOfAR * numOfDataSets;
	covarianceMatrix.setDegreeOfEquation(degreesOfAR + 1);

	// Add to matrix
	for( int iSet = 0; iSet < numOfDataSets; ++iSet ){
		const int numData = numOfData[iSet];
		for( int iData = degreesOfAR; iData < numData; ++iData ){
			for( int iDiff = 0; iDiff < degreesOfAR + 1; ++iDiff ){
				const double value = data[iSet][iData] * data[iSet][iData - iDiff];
				for( int iCol = 0; iCol < degreesOfAR + 1; ++iCol ){
					const int iRow = iCol + iDiff;
					if( iRow < degreesOfAR + 1 ){
						covarianceMatrix.addValue(iRow, iCol, value);
					}
				}
			}
		}
	}
#ifdef _DEBUG_WRITE
	covarianceMatrix.debugWriteMatrix();
#endif
	// Add small number to the diagonals for making the matrix positive-definite
	for( int iCol = 0; iCol < degreesOfAR + 1; ++iCol ){
		const double value = covarianceMatrix.getValue(iCol, iCol) * params.percentageOfSmallNumberAddedToDiagonals / 100.0;
		covarianceMatrix.addValue(iCol, iCol, value);
	}
	covarianceMatrix.multiplyScalarToMatrix( 1.0/static_cast<double>(numOfSamples) );
#ifdef _DEBUG_WRITE
	covarianceMatrix.debugWriteMatrix();
#endif

	const double determinant = covarianceMatrix.calculateDeterminant();
	// Force det(C) to 1 
	covarianceMatrix.multiplyScalarToMatrix( pow(1.0/fabs(determinant), 1.0/static_cast<double>(degreesOfAR+1)) );

#ifdef _DEBUG_WRITE
	std::cout << "det: " << covarianceMatrix.calculateDeterminant() << std::endl;
	covarianceMatrix.debugWriteMatrix();
#endif

	double paramB(0.0);
	double paramC(0.0);
	RobustWeightTukeysBiweights::calculateParams(degreesOfAR + 1, numOfSamples, paramB, paramC);

	double* MD = new double[numOfSamples];
	calculateMD(degreesOfAR, numOfDataSets, numOfData, data, covarianceMatrix, MD);
	// Initial estimation of scale
	double scale = Util::calculateMedian(numOfSamples, MD);
	scale = RobustWeightTukeysBiweights::calculateRobustScale(numOfSamples, MD, scale, paramB, paramC);
#ifdef _DEBUG_WRITE
	std::cout << "scale: " << scale << std::endl;
#endif
	ptrOutputFiles->writeCvgMessage("Number of data = " + Util::toString(numOfSamples) +
		", Initial scale = " + Util::toString(scale) + ", Determinant = " + Util::toString(determinant));

	double scalePre = scale;
	for( int iter = 0; iter < params.maxNumOfIterations; ++iter ){
		covarianceMatrix.zeroClearMatrix();
#ifdef _DEBUG_WRITE
		covarianceMatrix.debugWriteMatrix();
#endif
		// Add to matrix
		double sumOfWeights(0.0);
		double denominator(0.0);
		int icount(0);
		for( int iSet = 0; iSet < numOfDataSets; ++iSet ){
			const int numData = numOfData[iSet];
			for( int iData = degreesOfAR; iData < numData; ++iData ){
				const double weight = RobustWeightTukeysBiweights::calculateWeights(MD[icount] / scale, paramC);
				for( int iDiff = 0; iDiff < degreesOfAR + 1; ++iDiff ){
					const double value = data[iSet][iData] * data[iSet][iData - iDiff] * weight;
					for( int iCol = 0; iCol < degreesOfAR + 1; ++iCol ){
						const int iRow = iCol + iDiff;
						if( iRow < degreesOfAR + 1 ){
							covarianceMatrix.addValue(iRow, iCol, value);
						}
					}
				}
				sumOfWeights += weight;
				denominator += RobustWeightTukeysBiweights::calculateTermInDenominatorOfRobustCovariance(MD[icount] / scale, paramC);
				++icount;
			}
		}
		// Add small number to the diagonals for making the matrix positive-definite
		for( int iCol = 0; iCol < degreesOfAR + 1; ++iCol ){
			const double value = covarianceMatrix.getValue(iCol, iCol) * params.percentageOfSmallNumberAddedToDiagonals / 100.0;
			covarianceMatrix.addValue(iCol, iCol, value);
		}
		if( denominator < CommonParameters::EPS ){
			ptrOutputFiles->writeWarningMessage("Denominator of robust covariance is too small (" +  Util::toString(denominator) + ")");
		}
		const double factor = static_cast<double>(degreesOfAR + 1) / std::max(denominator, CommonParameters::EPS);
		//covarianceMatrix.multiplyScalarToMatrix( 1.0/sumOfWeights );
		covarianceMatrix.multiplyScalarToMatrix(factor);
#ifdef _DEBUG_WRITE
		covarianceMatrix.debugWriteMatrix();
#endif
		const double determinant = covarianceMatrix.calculateDeterminant();
		// Force det(C) to 1
		covarianceMatrix.multiplyScalarToMatrix( pow(1.0/fabs(determinant), 1.0/static_cast<double>(degreesOfAR+1)) );
#ifdef _DEBUG_WRITE
		std::cout << "det: " << covarianceMatrix.calculateDeterminant() << std::endl;
		covarianceMatrix.debugWriteMatrix();
#endif
		calculateMD(degreesOfAR, numOfDataSets, numOfData, data, covarianceMatrix, MD);
		scale = RobustWeightTukeysBiweights::calculateRobustScale(numOfSamples, MD, scale, paramB, paramC);
		ptrOutputFiles->writeCvgMessage("Iteration number = " + Util::toString(iter) + ", Sum of weights = " + Util::toString(sumOfWeights) + 
			", Scale = " + Util::toString(scale) + ", Determinant = " + Util::toString(determinant));
		if( fabs(scale-scalePre)/fabs(scalePre) < params.convergenceCriteria ){
			break;
		}
		scalePre = scale;
	}

	covarianceMatrix.multiplyScalarToMatrix( pow(scale,2) );
#ifdef _DEBUG_WRITE
	covarianceMatrix.debugWriteMatrix();
#endif

	delete [] MD;

}

// Calculate Mahalanobis distances
// @note covariance matrix is factorized in this function
void RobustPrewhitening::calculateMD( const int degreesOfAR, const int numOfDataSets, const int* const numOfData, double** data,
	DoubleDenseSquareSymmetricMatrix& covarianceMatrix, double* MD ) const{

	double* residualVector = new double[degreesOfAR + 1];
	double* workVector = new double[degreesOfAR + 1];
	covarianceMatrix.factorizeMatrix();
#ifdef _DEBUG_WRITE
	covarianceMatrix.debugWriteMatrix();
#endif

	int icount(0);
	for( int iSet = 0; iSet < numOfDataSets; ++iSet ){
		const int numData = numOfData[iSet];
		for( int iData = degreesOfAR; iData < numData; ++iData ){
			for( int iDiff = 0; iDiff < degreesOfAR + 1; ++iDiff ){
				residualVector[iDiff] = data[iSet][iData - iDiff];
			}
			covarianceMatrix.solveLinearEquation(residualVector, workVector);
			double mahalanobisDistances(0.0);
			for( int iDiff = 0; iDiff < degreesOfAR + 1; ++iDiff ){
				mahalanobisDistances += residualVector[iDiff] * workVector[iDiff];
			}
			MD[icount] = sqrt(fabs(mahalanobisDistances));
//#ifdef _DEBUG_WRITE
//			for( int iDiff = 0; iDiff < degreesOfAR + 1; ++iDiff ){
//				std::cout << "residualVector:" << residualVector[iDiff] << std::endl;
//			}
//			for( int iDiff = 0; iDiff < degreesOfAR + 1; ++iDiff ){
//				std::cout << "workVector:" << workVector[iDiff] << std::endl;
//			}
//			std::cout << "mahalanobisDistances:" << mahalanobisDistances << std::endl;
//#endif
			++icount;
		}
	}

	delete [] residualVector;
	delete [] workVector;

}

// Calculate weight for state vector in robust filter
double RobustPrewhitening::calculateWeightForStateVectorOfRobustFilter( const double val, const double paramA, const double paramB, const bool residualAssumedToBeZero ) const{

	assert( paramA < paramB );

	if( residualAssumedToBeZero || fabs(val) <= paramA ){
		return val;
	}
	else if( val > paramA && val <= paramB ){
		return paramA * (paramB - val) / (paramB - paramA);
	}
	else if( val >= -paramB && val < -paramA ){
		return - paramA * (paramB + val) / (paramB - paramA);
	}
	else{
		return 0.0;
	}

}

// Calculate weight for covariance matrix in robust filter
double RobustPrewhitening::calculateWeightForCovarianceMatrixOfRobustFilter( const double val, const double paramA, const double paramB, const bool residualAssumedToBeZero  ) const{

	assert( paramA < paramB );

	if( residualAssumedToBeZero || fabs(val) <= paramA ){
		return 1.0;
	}
	else if( val > paramA && val <= paramB ){
		return paramA * (paramB / val - 1.0) / (paramB - paramA);
	}
	else if( val >= -paramB && val < -paramA ){
		return - paramA * (paramB / val + 1.0) / (paramB - paramA);
	}
	else{
		return 0.0;
	}

}
