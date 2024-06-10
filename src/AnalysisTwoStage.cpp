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
#include "AnalysisTwoStage.h"
#include "Control.h"
#include "OutputFiles.h"

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
AnalysisTwoStage::AnalysisTwoStage()
{
}

// Destructer
AnalysisTwoStage::~AnalysisTwoStage()
{
}

// Calculate response functions
void AnalysisTwoStage::calculateResponseFunctions( const int iSegLen, const int freqDegree, const double timeLength, const double freq, 
	const int numSegmentsTotal, std::complex<double>** ftval, const std::vector< std::pair<std::string, std::string> >& times, 
	std::ofstream& ofsResp, std::ofstream& ofsRhoaPhs ){

	const Control* const ptrControl = Control::getInstance();
	const int numOutputVariables = ptrControl->getNumOutputVariables();

	std::complex<double>* resp0 = new std::complex<double>[numOutputVariables];
	std::complex<double>* resp1 = new std::complex<double>[numOutputVariables];
	double** hatDiagonals = new double* [numOutputVariables];
	for( int iOut = 0; iOut < numOutputVariables; ++iOut ){
		hatDiagonals[iOut] = new double [numSegmentsTotal];
	}
	// Estimate response functions
	calculateResponseFunctionsAux( iSegLen, freqDegree, timeLength, freq, numSegmentsTotal, ftval, times, ofsResp, ofsRhoaPhs, false, resp0, resp1, hatDiagonals );
	
	if( ptrControl->getErrorEstimationMethod() == Control::SUBSET_DELETION_JACKKNIFE ){
		subsetDeletionJackknife( iSegLen, freqDegree, timeLength, freq, numSegmentsTotal, ftval, times, ofsResp, ofsRhoaPhs, resp0, resp1, hatDiagonals );
	}
	else if( ptrControl->getErrorEstimationMethod() == Control::STRICT_BOOTSTRAP ){
		strictBootstrap( iSegLen, freqDegree, timeLength, freq, numSegmentsTotal, ftval, times, ofsResp, ofsRhoaPhs, resp0, resp1 );
	}

	delete [] resp0;
	delete [] resp1;
	for( int iOut = 0; iOut < numOutputVariables; ++iOut ){
		delete [] hatDiagonals[iOut];
	}
	delete [] hatDiagonals;

}

void AnalysisTwoStage::calculateResponseFunctionsAux( const int iSegLen,const int freqDegree, const double timeLength, const double freq,
	const int numSegmentsTotal, std::complex<double>** ftval, const std::vector< std::pair<std::string, std::string> >& times, 
	std::ofstream& ofsResp, std::ofstream& ofsRhoaPhs, const bool forJackknife, 
	std::complex<double>* respOut0, std::complex<double>* respOut1, double** hatDiagsOut ) const{

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();

	ptrOutputFiles->writeLogMessage("Calculate response functions by two stage procedure");
	ptrOutputFiles->writeCvgMessage("================================================================================");
	ptrOutputFiles->writeCvgMessage("Now Frequency(Hz): " + Util::toString(freq) + ", Period(s): " + Util::toString(1.0/freq));
	ptrOutputFiles->writeCvgMessage("================================================================================");
	const Control* const ptrControl = Control::getInstance();

	// First stage
	ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
	ptrOutputFiles->writeCvgAndLogMessage("Start the first stage processing");
	ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
	const int numInputVariables = ptrControl->getNumInputVariables();
	const int numRemoteReferenceVariables = ptrControl->getNumRemoteReferenceVariables();
	const bool outputResidual = !forJackknife && ptrControl->getOutputLevel() >= 2;
	// Parameters used in the second stage
	double* weights1stStage = new double[numSegmentsTotal];
	for( int iSeg = 0 ; iSeg < numSegmentsTotal; ++iSeg){
		weights1stStage[iSeg] = 1.0;
	}
	std::complex<double>** dataSyn = new std::complex<double>*[numInputVariables];
	for( int iInVar = 0; iInVar < numInputVariables; ++iInVar ){
		dataSyn[iInVar] = new std::complex<double>[numSegmentsTotal];
	}
	for( int iInVar = 0; iInVar < numInputVariables; ++iInVar ){
		assert( numRemoteReferenceVariables == 2 );
		ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
		ptrOutputFiles->writeCvgAndLogMessage("Calculate response functions for input variable " + Util::toString(iInVar));
		ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
		// Calculate response functions by ordinary least square method
		double coherence(0.0);
		std::complex<double> resp0(0.0, 0.0);
		std::complex<double> resp1(0.0, 0.0);
		double* weights = new double[numSegmentsTotal];
		double sumOfWeights(0.0);
		for( int iSeg = 0 ; iSeg < numSegmentsTotal; ++iSeg){
			weights[iSeg] = 1.0;
			sumOfWeights += weights[iSeg];
		}
#ifdef _DEBUG_WRITE
		for( int iSeg = 0 ; iSeg < numSegmentsTotal; ++iSeg){
			std::cout << "iSeg weights: " << iSeg << " " << weights[iSeg] << std::endl;
		}
#endif
		if( sumOfWeights < CommonParameters::EPS ){
			ptrOutputFiles->writeErrorMessage("Sum of weights is too small: " + Util::toString(sumOfWeights));
		}
		std::complex<double>* residuals = new std::complex<double>[numSegmentsTotal];
		ptrOutputFiles->writeCvgAndLogMessage("Calculate response functions by the ordinary least square method");
		calculateResponseFunctionsByWLSForFirstStage( iInVar, numSegmentsTotal, ftval, weights, residuals, resp0, resp1, coherence );
		std::vector<std::string> titles;
		std::vector<double>* outputValues = new std::vector<double>[numSegmentsTotal];
		if( outputResidual ){
			titles.push_back("OLS");
			for( int iSeg = 0 ; iSeg < numSegmentsTotal; ++iSeg){
				outputValues[iSeg].push_back(residuals[iSeg].real());
				outputValues[iSeg].push_back(residuals[iSeg].imag());
				outputValues[iSeg].push_back(weights[iSeg]);
			}
		}
		double absMaxResidual(0.0);
		for( int iSeg = 0 ; iSeg < numSegmentsTotal; ++iSeg){
			const double absResidual = std::abs(residuals[iSeg]);
			if( absResidual > absMaxResidual ){
				absMaxResidual = absResidual;
			}
		}
		bool skipRobustMethod(false);
		if( absMaxResidual < CommonParameters::EPS ){
			// Such a case where the remote reference field is equivalent to the input field
			ptrOutputFiles->writeCvgMessage("Robust method is not performed because residuals is nearly zero");
			skipRobustMethod = true;
		}
		const int rr0 = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, 0 );
		const int rr1 = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, 1 );
		if(!skipRobustMethod){
			double* hatDiagonals = new double[numSegmentsTotal];
			const double maxHatDiagOrg = calculateDiagonalComponentsOfHatMatrix( numSegmentsTotal, rr0, rr1, ftval, weights, hatDiagonals );
			double sumOfWeights(0.0);
			for( int iSeg = 0 ; iSeg < numSegmentsTotal; ++iSeg){
				sumOfWeights += weights[iSeg];
			}
			const double maxYOrg = maxHatDiagOrg * sumOfWeights / 2;
			double* leverageWeights = new double[numSegmentsTotal];
			for( int iSeg = 0 ; iSeg < numSegmentsTotal; ++iSeg){
				// The leverage weights are initialized with original weights
				leverageWeights[iSeg] = 1.0;
			}
			// Outer loop by gradual downweighting of leverage points
			const Control::ParamsForTreatmentOfHatMatrix paramsLeverageWeights = ptrControl->getParamsForTreatmentOfHatMatrix();
			const double chiFinal = paramsLeverageWeights.threshold;
			int numOuterIteration = 1;
			double baseNumber = 2;
			if( paramsLeverageWeights.applyLeverageWeights && maxYOrg > chiFinal ){
				assert( maxYOrg > 0.0 );
				numOuterIteration += static_cast<int>( log(maxYOrg/chiFinal)/log(baseNumber) ) + 1;	
			}
			if( numOuterIteration > paramsLeverageWeights.maxNumberOfOuterIteration ){
				baseNumber = pow( maxYOrg/chiFinal, 1.0/(static_cast<double>(paramsLeverageWeights.maxNumberOfOuterIteration) - 1.0) );
				numOuterIteration = paramsLeverageWeights.maxNumberOfOuterIteration;
			}
			double scale = RobustWeight::calculateScaleByMADN(numSegmentsTotal, residuals);
			// Calculate hat matris diagonals
			ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
			ptrOutputFiles->writeCvgMessage("Maximum number of the outer iteration: " + Util::toString(numOuterIteration));
			for( int outerIter = 0; outerIter < numOuterIteration; ++outerIter ){
				const double chi = std::max( maxYOrg / pow(baseNumber, outerIter), chiFinal );
				double sumOfWeights(0.0);
				for( int iSeg = 0 ; iSeg < numSegmentsTotal; ++iSeg){
					sumOfWeights += weights[iSeg];
				}
				if( sumOfWeights < CommonParameters::EPS ){
					const std::string msg = "Iteration is finished because the sum of weights (" +  Util::toString(sumOfWeights) + ") is too small";
					ptrOutputFiles->writeCvgMessage(msg);
					ptrOutputFiles->writeWarningMessage(msg);
					break;
				}
				const double hatDiagExpected = 2.0 / sumOfWeights;
				const double maxHatDiag = calculateDiagonalComponentsOfHatMatrix( numSegmentsTotal, rr0, rr1, ftval, weights, hatDiagonals );
				const double maxY = maxHatDiag / hatDiagExpected;
				ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
				ptrOutputFiles->writeCvgMessage("Outer iteration " + Util::toString(outerIter));
				ptrOutputFiles->writeCvgMessage("Maximum value of statistic y: " + Util::toString(maxY));
				ptrOutputFiles->writeCvgMessage("Target chi value: " + Util::toString(chi));
				ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
				if( paramsLeverageWeights.applyLeverageWeights ){
					// Calculate leverage weights
					for( int iSeg = 0 ; iSeg < numSegmentsTotal; ++iSeg){
						// The leverage weights are initialized with unity
						const double term0 = exp( exp( - chi * chi ) );
						const double y = hatDiagonals[iSeg] / hatDiagExpected;
						const double term1 = exp( - exp( chi * (y - chi) ) );
						leverageWeights[iSeg] *= term0 * term1;
#ifdef _DEBUG_WRITE
						std::cout << "iSeg hatDiagonals y term0*term1 leverageWeights: " <<  
							iSeg << " " << hatDiagonals[iSeg] << " "  << y << " " << term0*term1 << " " << leverageWeights[iSeg] << std::endl;
#endif
					}
				}
				if( outputResidual ){
					titles.push_back("outer_"+Util::toString(outerIter));
					for( int iSeg = 0 ; iSeg < numSegmentsTotal; ++iSeg){
						outputValues[iSeg].push_back(hatDiagonals[iSeg]/hatDiagExpected);
						outputValues[iSeg].push_back(leverageWeights[iSeg]);
					}
				}
				// Calculate response functions by linear regression using the first M-estimator
				calculateResponseFunctionsByIRWLSForFirstStage(0, iInVar, numSegmentsTotal, false, scale, ftval, leverageWeights, weights, residuals,
					resp0, resp1, coherence, titles, outputValues);
				// Calculate response functions by linear regression using the second M-estimator
				calculateResponseFunctionsByIRWLSForFirstStage(1, iInVar, numSegmentsTotal, true, scale, ftval, leverageWeights, weights, residuals,
					resp0, resp1, coherence, titles, outputValues);
			}
			delete [] hatDiagonals;
			delete [] leverageWeights;
			if( ptrControl->getOutputLevel() > 0 ){
				// Output spectral density functions to cvg file
				const int inp = ptrControl->getChannelIndex( CommonParameters::INPUT, iInVar );
				outputSpectralDensityFunctionsToCvgFile( numSegmentsTotal, timeLength, ftval[inp], ftval[rr0], ftval[rr1], weights );
			}
		}
		// Output responses at the first stage
		std::ostringstream msg;
		msg << "Responses at the 1st stage: ";
		msg << "("  << std::setw(12) << std::setprecision(4) << std::scientific << resp0.real() << "," 
					<< std::setw(12) << std::setprecision(4) << std::scientific << resp0.imag() << "), ";
		msg << "("  << std::setw(12) << std::setprecision(4) << std::scientific << resp1.real() << "," 
					<< std::setw(12) << std::setprecision(4) << std::scientific << resp1.imag() << ")";
		ptrOutputFiles->writeCvgMessage(msg.str());
		// Copy parameters used for the second stage
		for( int iSeg = 0 ; iSeg < numSegmentsTotal; ++iSeg){
			if( weights[iSeg] < weights1stStage[iSeg] ){
				weights1stStage[iSeg] = weights[iSeg];
			}
			dataSyn[iInVar][iSeg] = resp0 * ftval[rr0][iSeg] + resp1 * ftval[rr1][iSeg];
#ifdef _DEBUG_WRITE
			const int inp = ptrControl->getChannelIndex( CommonParameters::INPUT, iInVar );
			std::cout << iSeg  << " " << ftval[rr0][iSeg] << " " << ftval[rr1][iSeg] << " " << resp0 << " " << resp1 << " " << dataSyn[iInVar][iSeg] << std::endl;
#endif
		}
		if( outputResidual ){
			std::ostringstream oss;
			oss << "segm" << iSegLen << "_index" << freqDegree << "_output" << iInVar << "_residuals_first_stage.csv";
			writeResiduals( oss.str(), numSegmentsTotal, times, titles, outputValues );
		}
		// Release memory
		delete [] residuals;
		delete [] outputValues;
		delete [] weights;
	}
	if( !forJackknife && ptrControl->doesOutputFreqDomainDataToCsv() ){
		outputSyntheticInputData( iSegLen, freqDegree, numSegmentsTotal, dataSyn );
	}
	// Second stage
	ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
	ptrOutputFiles->writeCvgAndLogMessage("Start the second stage processing");
	ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
	if( !forJackknife ){
		ofsResp << std::setprecision(10) << std::scientific << freq;
		ofsResp << "," << std::setprecision(10) << std::scientific << 1.0/freq;	
		if( ptrControl->doesOutputApparentResistivityAndPhase() ){
			ofsRhoaPhs << std::setprecision(10) << std::scientific << freq;
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << 1.0/freq;
		}
	}
	const int numOutputVariables = ptrControl->getNumOutputVariables();
	double** weightsFinal = new double*[numOutputVariables];
	for( int iOutVar = 0; iOutVar < numOutputVariables; ++iOutVar ){
		weightsFinal[iOutVar] = NULL;// Initialize
	}
	for( int iOutVar = 0; iOutVar < numOutputVariables; ++iOutVar ){
		ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
		ptrOutputFiles->writeCvgAndLogMessage("Calculate response functions for output variable " + Util::toString(iOutVar));
		ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
		// Calculate response functions by ordinary least square method
		double coherence(0.0);
		std::complex<double> resp0(0.0, 0.0);
		std::complex<double> resp1(0.0, 0.0);
		double* weights = new double[numSegmentsTotal];
		double sumOfWeights(0.0);
		for( int iSeg = 0 ; iSeg < numSegmentsTotal; ++iSeg){
			weights[iSeg] = weights1stStage[iSeg];
			sumOfWeights += weights[iSeg];
		}
		if( sumOfWeights < CommonParameters::EPS ){
			ptrOutputFiles->writeErrorMessage("Sum of weights is too small: " + Util::toString(sumOfWeights));
		}
		std::complex<double>* residuals = new std::complex<double>[numSegmentsTotal];
		ptrOutputFiles->writeCvgAndLogMessage("Calculate response functions by the ordinary least square method");
		calculateResponseFunctionsByWLSForSecondStage( iOutVar, numSegmentsTotal, ftval, dataSyn, weights, residuals, resp0, resp1, coherence );
		std::vector<std::string> titles;
		std::vector<double>* outputValues = new std::vector<double>[numSegmentsTotal];
		if( outputResidual ){
			titles.push_back("OLS");
			for( int iSeg = 0 ; iSeg < numSegmentsTotal; ++iSeg){
				outputValues[iSeg].push_back(residuals[iSeg].real());
				outputValues[iSeg].push_back(residuals[iSeg].imag());
				outputValues[iSeg].push_back(weights[iSeg]);
			}
		}
		double absMaxResidual(0.0);
		for( int iSeg = 0 ; iSeg < numSegmentsTotal; ++iSeg){
			const double absResidual = std::abs(residuals[iSeg]);
			if( absResidual > absMaxResidual ){
				absMaxResidual = absResidual;
			}
		}
		bool skipRobustMethod(false);
		if( absMaxResidual < CommonParameters::EPS ){
			// Such a case where the remote reference field is equivalent to the input field
			ptrOutputFiles->writeCvgMessage("Robust method is not performed because residuals is nearly zero");
			skipRobustMethod = true;
		}
		if(!skipRobustMethod){
			// Calculate hat matris diagonals
			double* hatDiagonals = new double[numSegmentsTotal];
			assert( numInputVariables == 2 );
			const double maxHatDiagOrg = calculateDiagonalComponentsOfHatMatrix( numSegmentsTotal, 0, 1, dataSyn, weights, hatDiagonals );
			double sumOfWeights(0.0);
			for( int iSeg = 0 ; iSeg < numSegmentsTotal; ++iSeg){
				sumOfWeights += weights[iSeg];
			}
			const double maxYOrg = maxHatDiagOrg * sumOfWeights / 2;
			double* leverageWeights = new double[numSegmentsTotal];
			for( int iSeg = 0 ; iSeg < numSegmentsTotal; ++iSeg){
				leverageWeights[iSeg] = weights1stStage[iSeg];
#ifdef _DEBUG_WRITE
				std::cout << "iSeg leverageWeights: " << iSeg << " " << leverageWeights[iSeg] << std::endl;
#endif
			}
			// Outer loop by gradual downweighting of leverage points
			const Control::ParamsForTreatmentOfHatMatrix paramsLeverageWeights = ptrControl->getParamsForTreatmentOfHatMatrix();
			const double chiFinal = paramsLeverageWeights.threshold;
			int numOuterIteration = 1;
			double baseNumber = 2;
			if( paramsLeverageWeights.applyLeverageWeights && maxYOrg > chiFinal ){
				assert( maxYOrg > 0.0 );
				numOuterIteration += static_cast<int>( log(maxYOrg/chiFinal)/log(baseNumber) ) + 1;	
			}
			if( numOuterIteration > paramsLeverageWeights.maxNumberOfOuterIteration ){
				baseNumber = pow( maxYOrg/chiFinal, 1.0/(static_cast<double>(paramsLeverageWeights.maxNumberOfOuterIteration) - 1.0) );
				numOuterIteration = paramsLeverageWeights.maxNumberOfOuterIteration;
			}
			double scale = RobustWeight::calculateScaleByMADN(numSegmentsTotal, residuals);
			ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
			ptrOutputFiles->writeCvgMessage("Maximum number of the outer iteration: " + Util::toString(numOuterIteration));
			for( int outerIter = 0; outerIter < numOuterIteration; ++outerIter ){
				const double chi = std::max( maxYOrg / pow(baseNumber, outerIter), chiFinal );
				double sumOfWeights(0.0);
				for( int iSeg = 0 ; iSeg < numSegmentsTotal; ++iSeg){
					sumOfWeights += weights[iSeg];
				}
				if( sumOfWeights < CommonParameters::EPS ){
					const std::string msg = "Iteration is finished because the sum of weights (" +  Util::toString(sumOfWeights) + ") is too small";
					ptrOutputFiles->writeCvgMessage(msg);
					ptrOutputFiles->writeWarningMessage(msg);
					break;
				}
				const double hatDiagExpected = 2.0 / sumOfWeights;
				const double maxHatDiag = calculateDiagonalComponentsOfHatMatrix( numSegmentsTotal, 0, 1, dataSyn, weights, hatDiagonals );
				const double maxY = maxHatDiag / hatDiagExpected;
				ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
				ptrOutputFiles->writeCvgMessage("Outer iteration " + Util::toString(outerIter));
				ptrOutputFiles->writeCvgMessage("Maximum value of statistic y: " + Util::toString(maxY));
				ptrOutputFiles->writeCvgMessage("Target chi value: " + Util::toString(chi));
				ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
#ifdef _DEBUG_WRITE
				std::cout << "hatDiagExpected maxHatDiag maxY chi: " << hatDiagExpected << " " << maxHatDiag << " " <<maxY << " " << chi << std::endl;
#endif
				if( paramsLeverageWeights.applyLeverageWeights ){
					// Calculate leverage weights
					for( int iSeg = 0 ; iSeg < numSegmentsTotal; ++iSeg){
						const double term0 = exp( exp( - chi * chi ) );
						const double y = hatDiagonals[iSeg] / hatDiagExpected;
						const double term1 = exp( - exp( chi * (y - chi) ) );
						leverageWeights[iSeg] *= term0 * term1;
#ifdef _DEBUG_WRITE
						std::cout << "iSeg y term0*term1 leverageWeights: " << iSeg << " " << term0 * term1 << " " << leverageWeights[iSeg] << std::endl;
#endif
					}
				}
				if( outputResidual ){
					titles.push_back("outer_"+Util::toString(outerIter));
					for( int iSeg = 0 ; iSeg < numSegmentsTotal; ++iSeg){
						outputValues[iSeg].push_back(hatDiagonals[iSeg]/hatDiagExpected);
						outputValues[iSeg].push_back(leverageWeights[iSeg]);
					}
				}
				// Calculate response functions by linear regression using the first M-estimator
				calculateResponseFunctionsByIRWLSForSecondStage(0, iOutVar, numSegmentsTotal, false, scale, ftval, dataSyn, 
					leverageWeights, weights, residuals, resp0, resp1, coherence, titles, outputValues);
				// Calculate response functions by linear regression using the second M-estimator
				calculateResponseFunctionsByIRWLSForSecondStage(1, iOutVar, numSegmentsTotal, true, scale, ftval, dataSyn,
					leverageWeights, weights, residuals, resp0, resp1, coherence, titles, outputValues);
			}		
			delete [] hatDiagonals;
			delete [] leverageWeights;
			if( ptrControl->getOutputLevel() > 0 ){
				// Output spectral density functions to cvg file
				const int out = ptrControl->getChannelIndex( CommonParameters::OUTPUT, iOutVar );
				outputSpectralDensityFunctionsToCvgFile( numSegmentsTotal, timeLength, ftval[out], dataSyn[0], dataSyn[1], weights );
			}
		}
		respOut0[iOutVar] = resp0;
		respOut1[iOutVar] = resp1;
		if( !forJackknife ){
			const double maxHatDiag = calculateDiagonalComponentsOfHatMatrix( numSegmentsTotal, 0, 1, dataSyn, weights, hatDiagsOut[iOutVar] );
			// Output results
			ofsResp    << "," << std::setprecision(10) << std::scientific << resp0.real();
			ofsResp    << "," << std::setprecision(10) << std::scientific << resp0.imag();
			ofsResp    << "," << std::setprecision(10) << std::scientific << resp1.real();
			ofsResp    << "," << std::setprecision(10) << std::scientific << resp1.imag();
			ofsResp    << "," << std::setprecision(10) << std::scientific << coherence;
			if( ptrControl->doesOutputApparentResistivityAndPhase() ){
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcApparentResistivity(freq, resp0);
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcPhase(resp0);
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcApparentResistivity(freq, resp1);
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcPhase(resp1);
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << coherence;
			}
			weightsFinal[iOutVar] = new double[numSegmentsTotal];
			memcpy(weightsFinal[iOutVar], weights, sizeof(double)*numSegmentsTotal);
		}
		if( outputResidual ){
			std::ostringstream oss;
			oss << "segm" << iSegLen << "_index" << freqDegree << "_output" << iOutVar << "_residuals_second_stage.csv";
			writeResiduals( oss.str(), numSegmentsTotal, times, titles, outputValues );
		}
		// Release memory
		delete [] residuals;
		delete [] outputValues;
		delete [] weights;
	}
	// Release memory
	delete [] weights1stStage;

	const int typeOfErrorEstimationMethod = ptrControl->getErrorEstimationMethod();
	switch(typeOfErrorEstimationMethod){
	case Control::FIXED_WEIGHTS_JACKKNIFE:
		// Fixed-weights jackknife
		fixedWeightsJackknife( freq, numSegmentsTotal, weightsFinal, ftval, dataSyn, ofsResp, ofsRhoaPhs, respOut0, respOut1, hatDiagsOut );
		break;
	case Control::FIXED_WEIGHTS_BOOTSTRAP:
		// Fixed-weights bootstrap
		fixedWeightsBootstrap( freq, numSegmentsTotal, weightsFinal, ftval, dataSyn, ofsResp, ofsRhoaPhs, respOut0, respOut1 );
		break;
	case Control::SUBSET_DELETION_JACKKNIFE:
		// Go through
	case Control::STRICT_BOOTSTRAP:
		break;
	default:
		ptrOutputFiles->writeErrorMessage("Unsupported error estimation method : " + Util::toString(typeOfErrorEstimationMethod));
		break;
	}

	// Release memory
	for( int iInVar = 0; iInVar < numInputVariables; ++iInVar ){
		delete [] dataSyn[iInVar];
	}
	delete [] dataSyn;
	for( int iOutVar = 0; iOutVar < numOutputVariables; ++iOutVar ){
		if( weightsFinal[iOutVar] != NULL ){
			delete [] weightsFinal[iOutVar];
		}
	}
	delete [] weightsFinal;

}

// Calculate response functions by iteratively reweighted least squares for the first stage
void AnalysisTwoStage::calculateResponseFunctionsByIRWLSForFirstStage(const int iRobustWeight, const int inputVariableIndex,
	const int numSegments, const bool fixScale, double& scale, std::complex<double>** data, const double* const leverageWeights,
	double* weights, std::complex<double>* residuals, std::complex<double>& resp0, std::complex<double>& resp1,
	double& coherence, std::vector<std::string>& titles, std::vector<double>* outputValues) const {

	const Control* const ptrControl = Control::getInstance();
	const int inp = ptrControl->getChannelIndex(CommonParameters::INPUT, inputVariableIndex);
	const int numRemoteReferenceVariables = ptrControl->getNumRemoteReferenceVariables();
	assert(numRemoteReferenceVariables >= 2);
	const int rr0 = ptrControl->getChannelIndex(CommonParameters::REMOTE_REFERENCE, 0);
	const int rr1 = ptrControl->getChannelIndex(CommonParameters::REMOTE_REFERENCE, 1);
	const bool priorityOnFirst = (inputVariableIndex % 2 == 0) ? true : false;// Priority on diagonals

	calculateResponseFunctionsByIRWLS(iRobustWeight, data[inp], data[rr0], data[rr1], numSegments, fixScale, scale, leverageWeights,
		weights, residuals, resp0, resp1, coherence, titles, outputValues, priorityOnFirst);

}

// Calculate response functions by iteratively reweighted least squares for the second stage
void AnalysisTwoStage::calculateResponseFunctionsByIRWLSForSecondStage(const int iRobustWeight, const int outputVariableIndex,
	const int numSegments, const bool fixScale, double& scale, std::complex<double>** data, std::complex<double>** dataSyn, const double* const leverageWeights,
	double* weights, std::complex<double>* residuals, std::complex<double>& resp0, std::complex<double>& resp1, double& coherence,
	std::vector<std::string>& titles, std::vector<double>* outputValues) const {

	const Control* const ptrControl = Control::getInstance();
	const int out = ptrControl->getChannelIndex(CommonParameters::OUTPUT, outputVariableIndex);
	const bool priorityOnFirst = (outputVariableIndex % 2 == 0) ? false : true;// Priority on off-diagonals

	calculateResponseFunctionsByIRWLS(iRobustWeight, data[out], dataSyn[0], dataSyn[1], numSegments, fixScale, scale, leverageWeights,
		weights, residuals, resp0, resp1, coherence, titles, outputValues, priorityOnFirst);

}

// Calculate response function by the weighted least square method for the first stage
void AnalysisTwoStage::calculateResponseFunctionsByWLSForFirstStage( const int inputVariableIndex, const int numSegments, std::complex<double>** data, 
	 double* weights, std::complex<double>* residuals, std::complex<double>& resp0, std::complex<double>& resp1, double& coherence ) const{

	const Control* const ptrControl = Control::getInstance();
	const int inp = ptrControl->getChannelIndex( CommonParameters::INPUT, inputVariableIndex );
	const int numRemoteReferenceVariables = ptrControl->getNumRemoteReferenceVariables();
	assert( numRemoteReferenceVariables >= 2 );
	const int rr0 = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, 0 );
	const int rr1 = ptrControl->getChannelIndex( CommonParameters::REMOTE_REFERENCE, 1 );
	const bool priorityOnFirst = ( inputVariableIndex % 2 == 0 ) ? true : false;// Priority on diagonals

	calculateResponseFunctionByWLS( data[inp], data[rr0], data[rr1], numSegments, weights, residuals, resp0, resp1, coherence, priorityOnFirst );

}

// Calculate response function by the weighted least square method for the second stage
void AnalysisTwoStage::calculateResponseFunctionsByWLSForSecondStage( const int outputVariableIndex, const int numSegments, std::complex<double>** data, 
	 std::complex<double>** dataSyn, double* weights, std::complex<double>* residuals, std::complex<double>& resp0, std::complex<double>& resp1, double& coherence ) const{

	const Control* const ptrControl = Control::getInstance();
	const int out = ptrControl->getChannelIndex( CommonParameters::OUTPUT, outputVariableIndex );
	const bool priorityOnFirst = ( outputVariableIndex % 2 == 0 ) ? false : true;// Priority on off-diagonals

	calculateResponseFunctionByWLS( data[out], dataSyn[0], dataSyn[1], numSegments, weights, residuals, resp0, resp1, coherence, priorityOnFirst );

}

// Perform fixed-weights jackknife
void AnalysisTwoStage::fixedWeightsJackknife( const double freq, const int numSegmentsTotal, double** weightsOrg, std::complex<double>** ftval,
	std::complex<double>** dataSyn, std::ofstream& ofsResp, std::ofstream& ofsRhoaPhs, 
	const std::complex<double>* const resp0, const std::complex<double>* const resp1, double** hatDiagonals ) const{

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeLogMessage("Fixed-weights jackknife is performed to estimate errors");

	const Control* const ptrControl = Control::getInstance();
	const int numOutputVariables = ptrControl->getNumOutputVariables();

	double** weights = new double*[numOutputVariables];
	for( int iOut = 0; iOut < numOutputVariables; ++iOut ){
		weights[iOut] = new double[numSegmentsTotal];
		memcpy(weights[iOut], weightsOrg[iOut], sizeof(double)*numSegmentsTotal);
	}
#ifdef _DEBUG_WRITE
	for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
		for( int iOut = 0; iOut < numOutputVariables; ++iOut ){
			std::cout << weightsOrg[iOut][iSeg] << " ";
		}
		std::cout << std::endl;
		for( int iOut = 0; iOut < numOutputVariables; ++iOut ){
			std::cout << weights[iOut][iSeg] << " ";
		}
		std::cout << std::endl;
	}
#endif
		
	std::complex<double>** pseudoResp0 = new std::complex<double>*[numSegmentsTotal];
	std::complex<double>** pseudoResp1 = new std::complex<double>*[numSegmentsTotal];
	for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
		pseudoResp0[iSeg] = new std::complex<double>[numOutputVariables];
		pseudoResp1[iSeg] = new std::complex<double>[numOutputVariables];
		for( int iOut = 0; iOut < numOutputVariables; ++iOut ){
			weights[iOut][iSeg] = 0.0;// Replace
			const int index = ptrControl->getChannelIndex( CommonParameters::OUTPUT, iOut );
			const bool priorityOnFirst = ( iOut % 2 == 0 ) ? false : true;// Priority on off-diagonals
			std::complex<double> temp0(0.0,0.0);
			std::complex<double> temp1(0.0,0.0);
			calculateResponseFunctionByWLSAux( numSegmentsTotal, ftval[index], dataSyn[0], dataSyn[1], weights[iOut], temp0, temp1, priorityOnFirst );
			const double hatMatrixDiagonal = hatDiagonals[iOut][iSeg];
			double factor = static_cast<double>(numSegmentsTotal) * (1.0 - hatMatrixDiagonal);
			if( hatMatrixDiagonal > 1.0 ){
				factor= 0.0;
			}else if( hatMatrixDiagonal < 0.0 ){
				factor = static_cast<double>(numSegmentsTotal);
			}
			pseudoResp0[iSeg][iOut] = resp0[iOut] + factor * (resp0[iOut] - temp0);
			pseudoResp1[iSeg][iOut] = resp1[iOut] + factor * (resp1[iOut] - temp1);
			weights[iOut][iSeg] = weightsOrg[iOut][iSeg];// Restore
		}
	}
	// Release memory
	for( int iOut = 0; iOut < numOutputVariables; ++iOut ){
		delete [] weights[iOut];
	}
	delete [] weights;

	// Calculate & output error bars
	for( int iOut = 0; iOut < numOutputVariables; ++iOut ){
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
			const double dResp0 = sqrt(variance0);
			const double dResp1 = sqrt(variance1);
			ofsResp    << "," << std::setprecision(10) << std::scientific << dResp0;		
			ofsResp    << "," << std::setprecision(10) << std::scientific << dResp1;
			if( ptrControl->doesOutputApparentResistivityAndPhase() ){
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcApparentResistivityError(freq, resp0[iOut], dResp0);
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcPhaseError(resp0[iOut], dResp0);
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcApparentResistivityError(freq, resp1[iOut], dResp1);
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcPhaseError(resp1[iOut], dResp1);
			}
		}else{
			ofsResp    << "," << std::setprecision(10) << std::scientific << 1.0e10;
			ofsResp    << "," << std::setprecision(10) << std::scientific << 1.0e10;
			if( ptrControl->doesOutputApparentResistivityAndPhase() ){
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << 1.0e10;
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << 180.0;
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << 1.0e10;
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << 180.0;
			}
		}
	}
	ofsResp << std::endl;
	ofsResp.flush();
	if( ptrControl->doesOutputApparentResistivityAndPhase() ){
		ofsRhoaPhs << std::endl;
		ofsRhoaPhs.flush();
	}

	// Release memory
	for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
		delete [] pseudoResp0[iSeg];
		delete [] pseudoResp1[iSeg];
	}
	delete [] pseudoResp0;
	delete [] pseudoResp1;

}

// Perform fixed-weights bootstrap
void AnalysisTwoStage::fixedWeightsBootstrap( const double freq, const int numSegmentsTotal, double** weightsOrg, std::complex<double>** ftval,
	std::complex<double>** dataSyn, std::ofstream& ofsResp, std::ofstream& ofsRhoaPhs,
	const std::complex<double>* const resp0, const std::complex<double>* const resp1 ) const{

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeLogMessage("Fixed-weights bootstrap is performed to estimate errors");

	const Control* const ptrControl = Control::getInstance();
	const int numOutputVariables = ptrControl->getNumOutputVariables();

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
		resp0Sample[iSample] = new std::complex<double>[numOutputVariables];
		resp1Sample[iSample] = new std::complex<double>[numOutputVariables];
		for( int iOut = 0; iOut < numOutputVariables; ++iOut ){
			const int index = ptrControl->getChannelIndex( CommonParameters::OUTPUT, iOut );
			const bool priorityOnFirst = ( iOut % 2 == 0 ) ? false : true;// Priority on off-diagonals
			std::complex<double> temp0(0.0,0.0);
			std::complex<double> temp1(0.0,0.0);
			calculateResponseFunctionByWLSForBootstrap( numSegmentsTotal, segmentIndexes, ftval[index], dataSyn[0], dataSyn[1], weightsOrg[iOut],
				temp0, temp1, priorityOnFirst );
			resp0Sample[iSample][iOut] = temp0;
			resp1Sample[iSample][iOut] = temp1;
		}
	}
	delete [] segmentIndexes;

	// Calculate error of response functions
	assert(numOfSamples > 2);
	for( int iOut = 0; iOut < numOutputVariables; ++iOut ){
		// Calculate average
		std::complex<double> average0 = std::complex<double>(0.0, 0.0);
		std::complex<double> average1 = std::complex<double>(0.0, 0.0);
		for( int iSample = 0; iSample < numOfSamples; ++iSample ){
			average0 += resp0Sample[iSample][iOut];
			average1 += resp1Sample[iSample][iOut];
		}
		average0 /= static_cast<double>(numOfSamples);
		average1 /= static_cast<double>(numOfSamples);
		// Calculate variance
		double variance0(0.0);
		double variance1(0.0);
		for( int iSample = 0; iSample < numOfSamples; ++iSample ){
			variance0 += std::norm(resp0Sample[iSample][iOut] - average0);
			variance1 += std::norm(resp1Sample[iSample][iOut] - average1);
		}
		variance0 /= static_cast<double>(2 * numOfSamples - 4);
		variance1 /= static_cast<double>(2 * numOfSamples - 4);
		const double dResp0 = sqrt(variance0);
		const double dResp1 = sqrt(variance1);
		ofsResp    << "," << std::setprecision(10) << std::scientific << dResp0;	
		ofsResp    << "," << std::setprecision(10) << std::scientific << dResp1;
		if( ptrControl->doesOutputApparentResistivityAndPhase() ){
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcApparentResistivityError(freq, resp0[iOut], dResp0);
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcPhaseError(resp0[iOut], dResp0);
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcApparentResistivityError(freq, resp1[iOut], dResp1);
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcPhaseError(resp1[iOut], dResp1);
		}
	}
	ofsResp << std::endl;
	ofsResp.flush();
	if( ptrControl->doesOutputApparentResistivityAndPhase() ){
		ofsRhoaPhs << std::endl;
		ofsRhoaPhs.flush();
	}
	// Release memory
	for( int iSample = 0; iSample < numOfSamples; ++iSample ){
		delete [] resp0Sample[iSample];
		delete [] resp1Sample[iSample];
	}
	delete [] resp0Sample;
	delete [] resp1Sample;

}

// Output synthetic input data
void AnalysisTwoStage::outputSyntheticInputData( const int iSegLen, const int freqDegree, const int numSegmentsTotal, std::complex<double>** dataSyn ) const{

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeLogMessage("Output synthetic input data");

	const Control* const ptrControl = Control::getInstance();
	for( int iInVar = 0; iInVar < ptrControl->getNumInputVariables(); ++iInVar ){
		const int iChan = ptrControl->getChannelIndex( CommonParameters::INPUT, iInVar );
		std::ostringstream fileName;
		fileName << "syn_inp_data_segm" << iSegLen << "_index" << freqDegree << "_chan_" << iChan << ".csv";
		std::ofstream ofs;
		ofs.open( fileName.str().c_str(), std::ios::out );
		if( ofs.fail() ){
			ptrOutputFiles->writeLogMessage("File open error !! : " + fileName.str());
		}
		ofs << "real,imaginary" << std::endl;
		for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
			ofs << std::setprecision(12) << std::scientific << dataSyn[iInVar][iSeg].real() << ",";
			ofs << std::setprecision(12) << std::scientific << dataSyn[iInVar][iSeg].imag() << std::endl;
		}
		ofs.close();
	}

}

// Estimate error by strict bootstrap
void AnalysisTwoStage::strictBootstrap( const int iSegLen,const int freqDegree, const double timeLength, const double freq,
	const int numSegmentsTotal, std::complex<double>** ftval, const std::vector< std::pair<std::string, std::string> >& times,
	std::ofstream& ofsResp, std::ofstream& ofsRhoaPhs,
	const std::complex<double>* const resp0, const std::complex<double>* const resp1 ) const{

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeLogMessage("Strict bootstrap is performed to estimate errors");

	const Control* const ptrControl = Control::getInstance();
	const int numOutputVariables = ptrControl->getNumOutputVariables();

	// Copy Fourier transformed values
	const int numChannels = ptrControl->getNumberOfChannels();
	std::complex<double>** ftvalForBootstrap = new std::complex<double>*[numChannels];
	for( int iChan = 0; iChan < numChannels; ++iChan ){
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
				ftvalForBootstrap[iChan][icount] = ftval[iChan][iSeg];
			}
		}
		resp0Sample[iSample] = new std::complex<double>[numOutputVariables];
		resp1Sample[iSample] = new std::complex<double>[numOutputVariables];
		calculateResponseFunctionsAux( iSegLen, freqDegree, timeLength, freq, numSegmentsTotal, ftvalForBootstrap, times, 
			ofsResp, ofsRhoaPhs, true, resp0Sample[iSample], resp1Sample[iSample], NULL );
	}
	ptrOutputFiles->restartToWriteCvgMessage();
	ptrOutputFiles->restartToWriteLogMessage();
	ptrOutputFiles->restartToWriteWarningMessage();

	// Calculate error of response functions
	assert(numOfSamples > 2);
	for( int iOut = 0; iOut < numOutputVariables; ++iOut ){
		// Calculate average
		std::complex<double> average0 = std::complex<double>(0.0, 0.0);
		std::complex<double> average1 = std::complex<double>(0.0, 0.0);
		for( int iSample = 0; iSample < numOfSamples; ++iSample ){
			average0 += resp0Sample[iSample][iOut];
			average1 += resp1Sample[iSample][iOut];
		}
		average0 /= static_cast<double>(numOfSamples);
		average1 /= static_cast<double>(numOfSamples);
		// Calculate variance
		double variance0(0.0);
		double variance1(0.0);
		for( int iSample = 0; iSample < numOfSamples; ++iSample ){
			variance0 += std::norm(resp0Sample[iSample][iOut] - average0);
			variance1 += std::norm(resp1Sample[iSample][iOut] - average1);
		}
		variance0 /= static_cast<double>(2 * numOfSamples - 4);
		variance1 /= static_cast<double>(2 * numOfSamples - 4);
		const double dResp0 = sqrt(variance0);
		const double dResp1 = sqrt(variance1);
		ofsResp    << "," << std::setprecision(10) << std::scientific << dResp0;		
		ofsResp    << "," << std::setprecision(10) << std::scientific << dResp1;
		if( ptrControl->doesOutputApparentResistivityAndPhase() ){
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcApparentResistivityError(freq, resp0[iOut], dResp0);
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcPhaseError(resp0[iOut], dResp0);
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcApparentResistivityError(freq, resp1[iOut], dResp1);
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcPhaseError(resp1[iOut], dResp1);
		}
	}
	ofsResp << std::endl;
	ofsResp.flush();
	if( ptrControl->doesOutputApparentResistivityAndPhase() ){
		ofsRhoaPhs << std::endl;
		ofsRhoaPhs.flush();
	}

	// Release memory
	for( int iChan = 0; iChan < numChannels; ++iChan ){
		delete [] ftvalForBootstrap[iChan];
	}
	delete [] ftvalForBootstrap;
	for( int iSample = 0; iSample < numOfSamples; ++iSample ){
		delete [] resp0Sample[iSample];
		delete [] resp1Sample[iSample];
	}
	delete [] resp0Sample;
	delete [] resp1Sample;
	delete [] segmentIndexes;

}

// Perform subset deletion jackknife
void AnalysisTwoStage::subsetDeletionJackknife( const int iSegLen,const int freqDegree, const double timeLength, const double freq,
	const int numSegmentsTotal, std::complex<double>** ftval, const std::vector< std::pair<std::string, std::string> >& times,
	std::ofstream& ofsResp, std::ofstream& ofsRhoaPhs, 
	const std::complex<double>* const resp0, const std::complex<double>* const resp1, double** hatDiagonals ) const{

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
	ptrOutputFiles->stopToWriteCvgMessage();
	ptrOutputFiles->stopToWriteLogMessage();
	ptrOutputFiles->stopToWriteWarningMessage();

	const int numOutputVariables = ptrControl->getNumOutputVariables();

#ifdef _DEBUG_WRITE
	for( int iSeg = 0;iSeg < numSegmentsTotal; ++iSeg){
		for( int iOut = 0; iOut < numOutputVariables; ++iOut ){
			std::cout << std::setw(15) << hatDiagonals[iOut][iSeg];
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	for( int iChan = 0; iChan < ptrControl->getNumberOfChannels(); ++iChan ){
		for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
			std::cout << "iChan iSeg val: " << iChan << " " << iSeg << " " << ftval[iChan][iSeg] << std::endl;
		}
	}
#endif

	// Copy Fourier transformed values
	const int numChannels = ptrControl->getNumberOfChannels();
	std::complex<double>** ftvalForJackknife = new std::complex<double>*[numChannels];
	for( int iChan = 0; iChan < numChannels; ++iChan ){
		ftvalForJackknife[iChan] = new std::complex<double>[numSegmentsTotal - numOmittedData];
	}
	std::complex<double>** pseudoResp0 = new std::complex<double>*[numOfSubsets];
	std::complex<double>** pseudoResp1 = new std::complex<double>*[numOfSubsets];
	for( int iSubset = 0 ; iSubset < numOfSubsets; ++iSubset ){
		double* averageHatDiags = new double[numOutputVariables];
		for( int iOut = 0; iOut < numOutputVariables; ++iOut ){
			// Zero clear
			averageHatDiags[iOut] = 0.0;
		}
		const int iSegOmitStart = iSubset * numOmittedData;
		const int iSegOmitEnd = iSegOmitStart + numOmittedData;
		assert(iSegOmitEnd <= numSegmentsTotal);
		for( int iSeg = iSegOmitStart; iSeg < iSegOmitEnd; ++iSeg ){
			for( int iOut = 0; iOut < numOutputVariables; ++iOut ){
				averageHatDiags[iOut] += hatDiagonals[iOut][iSeg];
			}
		}
		for( int iOut = 0; iOut < numOutputVariables; ++iOut ){
			averageHatDiags[iOut] /= static_cast<double>(numOmittedData);
		}
#ifdef _DEBUG_WRITE
		for( int iOut = 0; iOut < numOutputVariables; ++iOut ){
			std::cout << std::setw(20) << averageHatDiags[iOut];
		}
		std::cout << std::endl;
#endif
		for( int iChan = 0; iChan < numChannels; ++iChan ){
			int icount(0);
			for( int iSeg = 0; iSeg < iSegOmitStart; ++iSeg, ++icount ){
				ftvalForJackknife[iChan][icount] = ftval[iChan][iSeg];
			}
			for( int iSeg = iSegOmitEnd; iSeg < numSegmentsTotal; ++iSeg, ++icount ){
				ftvalForJackknife[iChan][icount] = ftval[iChan][iSeg];
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
		std::complex<double>* temp0 = new std::complex<double>[numOutputVariables];
		std::complex<double>* temp1 = new std::complex<double>[numOutputVariables];
		calculateResponseFunctionsAux( iSegLen, freqDegree, timeLength, freq, numSegmentsTotal - numOmittedData,
			ftvalForJackknife, times, ofsResp, ofsRhoaPhs, true, temp0, temp1, NULL );
		pseudoResp0[iSubset] = new std::complex<double>[numOutputVariables];
		pseudoResp1[iSubset] = new std::complex<double>[numOutputVariables];
		for( int iOut = 0; iOut < numOutputVariables; ++iOut ){
			double factor = static_cast<double>(numOfSubsets) * (1.0 - averageHatDiags[iOut]);
			if( averageHatDiags[iOut] > 1.0 ){
				factor= 0.0;
			}else if( averageHatDiags[iOut] < 0.0 ){
				factor = static_cast<double>(numOfSubsets);
			}
			pseudoResp0[iSubset][iOut] = resp0[iOut] + factor * (resp0[iOut] - temp0[iOut]);
			pseudoResp1[iSubset][iOut] = resp1[iOut] + factor * (resp1[iOut] - temp1[iOut]);
		}
		delete [] averageHatDiags;
		delete [] temp0;
		delete [] temp1;
	}
	ptrOutputFiles->restartToWriteCvgMessage();
	ptrOutputFiles->restartToWriteLogMessage();
	ptrOutputFiles->restartToWriteWarningMessage();
	for( int iChan = 0; iChan < numChannels; ++iChan ){
		delete [] ftvalForJackknife[iChan];
	}
	delete [] ftvalForJackknife;

	// Calculate & output error bars
	for( int iOut = 0; iOut < numOutputVariables; ++iOut ){
		if( numOfSubsets > 2 ){
			std::complex<double> avgResp0(0.0, 0.0);
			std::complex<double> avgResp1(0.0, 0.0);
			for( int iSubset = 0 ; iSubset < numOfSubsets; ++iSubset ){
				avgResp0 += pseudoResp0[iSubset][iOut];
				avgResp1 += pseudoResp1[iSubset][iOut];
			}
			const double factor = 1.0 / static_cast<double>(numOfSubsets);
			avgResp0 *= factor;
			avgResp1 *= factor;
			double variance0(0.0);
			double variance1(0.0);
			for( int iSubset = 0 ; iSubset < numOfSubsets; ++iSubset ){
				variance0 += std::norm( pseudoResp0[iSubset][iOut] - avgResp0 );
				variance1 += std::norm( pseudoResp1[iSubset][iOut] - avgResp1 );
			}
			const double factor2 = factor / static_cast<double>(2 * numOfSubsets - 4);
			variance0 *= factor2;
			variance1 *= factor2;
			const double dResp0 = sqrt(variance0);
			const double dResp1 = sqrt(variance1);
			ofsResp    << "," << std::setprecision(10) << std::scientific << dResp0;
			ofsResp    << "," << std::setprecision(10) << std::scientific << dResp1;
			if( ptrControl->doesOutputApparentResistivityAndPhase() ){
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcApparentResistivityError(freq, resp0[iOut], dResp0);
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcPhaseError(resp0[iOut], dResp0);
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcApparentResistivityError(freq, resp1[iOut], dResp1);
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcPhaseError(resp1[iOut], dResp1);
			}
		}else{
			ofsResp    << "," << std::setprecision(10) << std::scientific << 1.0e10;
			ofsResp    << "," << std::setprecision(10) << std::scientific << 1.0e10;
			if( ptrControl->doesOutputApparentResistivityAndPhase() ){
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << 1.0e10;
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << 180.0;
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << 1.0e10;
				ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << 180.0;
			}
		}
	}
	ofsResp << std::endl;
	ofsResp.flush();
	if( ptrControl->doesOutputApparentResistivityAndPhase() ){
		ofsRhoaPhs << std::endl;
		ofsRhoaPhs.flush();
	}

	for( int iSubset = 0 ; iSubset < numOfSubsets; ++iSubset ){
		delete [] pseudoResp0[iSubset];
		delete [] pseudoResp1[iSubset];
	}
	delete [] pseudoResp0;
	delete [] pseudoResp1;

}

// Write residuals
void AnalysisTwoStage::writeResiduals( const std::string& fileName, const int numSegmentsTotal,
	const std::vector< std::pair<std::string, std::string> >& times, const std::vector<std::string>& titles, const std::vector<double>* outputValues ) const{

	std::ofstream ofs;
	ofs.open( fileName.c_str(), std::ios::out );
	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	if( ofs.fail() ){
		ptrOutputFiles->writeErrorMessage( "File open error : " + fileName );
	}
	ofs << "index";
	ofs << ",start_time,end_time";
	for( std::vector<std::string>::const_iterator itr = titles.begin(); itr != titles.end(); ++itr ){
		if( itr->substr(0,6) == "outer_" ){
			ofs << ",statistic_y_" << *itr;
			ofs << ",leverage_weight_" << *itr;
		}else{
			ofs << ",residual_real_" << *itr;
			ofs << ",residual_imag_" << *itr;
			ofs << ",weight_" << *itr;
		}
	}
	ofs << std::endl;

	int index(0);
	for( int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg ){
		ofs << index;
		const std::string timeStart = times[iSeg].first;
		const std::string timeEnd = times[iSeg].second;
		ofs << "," << timeStart << "," << timeEnd;
		for( std::vector<double>::const_iterator itr = outputValues[index].begin(); itr != outputValues[index].end(); ++itr ){
			ofs << "," << std::setprecision(10) << std::scientific << *itr;
		}
		ofs << std::endl;
		++index;
	}
	ofs.close();

}
