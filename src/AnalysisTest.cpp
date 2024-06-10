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
#include "AnalysisTest.h"
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
AnalysisTest::AnalysisTest()
{
}

// Destructer
AnalysisTest::~AnalysisTest()
{
}

// Calculate response functions
void AnalysisTest::calculateResponseFunctions(const int iSegLen, const int freqDegree, const double timeLength, const double freq,
	const int numSegmentsTotal, std::complex<double>** ftval, const std::vector< std::pair<std::string, std::string> >& times,
	std::ofstream& ofsResp, std::ofstream& ofsRhoaPhs) {

	// Estimate response functions
	const Control* const ptrControl = Control::getInstance();
	const int numOutputVariables = ptrControl->getNumOutputVariables();
	const int numInputVariables = ptrControl->getNumInputVariables();
	assert(numInputVariables == 2);
	std::complex<double>* istf0 = new std::complex<double>[numOutputVariables + numInputVariables];
	std::complex<double>* istf1 = new std::complex<double>[numOutputVariables + numInputVariables];
	double coherenceMin(1.0);
	int icount(0);
	for (int iOut = 0; iOut < numOutputVariables; ++iOut, ++icount) {
		const int index = ptrControl->getChannelIndex(CommonParameters::OUTPUT, iOut);
		double coherence(0.0);
		calculateResponseFunctionsAux(iSegLen, freqDegree, timeLength, freq, numSegmentsTotal, ftval, times, ofsResp, ofsRhoaPhs, false, index, istf0[icount], istf1[icount], coherence);
		if (coherence < coherenceMin) {
			coherenceMin = coherence;
		}
	}
	for (int iInp = 0; iInp < numInputVariables; ++iInp, ++icount ) {
		const int index = ptrControl->getChannelIndex(CommonParameters::INPUT, iInp);
		double coherence(0.0);
		calculateResponseFunctionsAux(iSegLen, freqDegree, timeLength, freq, numSegmentsTotal, ftval, times, ofsResp, ofsRhoaPhs, false, index, istf0[icount], istf1[icount], coherence);
		if (coherence < coherenceMin) {
			coherenceMin = coherence;
		}
	}

	// Output response functions
	const std::complex<double> Txx = istf0[numOutputVariables];
	const std::complex<double> Txy = istf1[numOutputVariables];
	const std::complex<double> Tyx = istf0[numOutputVariables + 1];
	const std::complex<double> Tyy = istf1[numOutputVariables + 1];
	const std::complex<double> det = Txx * Tyy - Txy * Tyx;
	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	if (std::abs(det) < CommonParameters::EPS) {
		ptrOutputFiles->writeErrorMessage("Determinant is too small: " + Util::toString(std::abs(det)));
	}
	ofsResp << std::setprecision(10) << std::scientific << freq;
	ofsResp << "," << std::setprecision(10) << std::scientific << 1.0 / freq;
	if (ptrControl->doesOutputApparentResistivityAndPhase()) {
		ofsRhoaPhs << std::setprecision(10) << std::scientific << freq;
		ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << 1.0 / freq;
	}
	std::complex<double>* respOut0 = new std::complex<double>[numOutputVariables];
	std::complex<double>* respOut1 = new std::complex<double>[numOutputVariables];
	for (int iOut = 0; iOut < numOutputVariables; ++iOut) {
		const int index = ptrControl->getChannelIndex(CommonParameters::OUTPUT, iOut);
		const std::complex<double> U_x = istf0[index];
		const std::complex<double> U_y = istf1[index];
		respOut0[iOut] = (U_x * Tyy - U_y * Tyx) / det;
		respOut1[iOut] = (U_y * Txx - U_x * Txy) / det;
		ofsResp << "," << std::setprecision(10) << std::scientific << respOut0[iOut].real();
		ofsResp << "," << std::setprecision(10) << std::scientific << respOut0[iOut].imag();
		ofsResp << "," << std::setprecision(10) << std::scientific << respOut1[iOut].real();
		ofsResp << "," << std::setprecision(10) << std::scientific << respOut1[iOut].imag();
		ofsResp << "," << std::setprecision(10) << std::scientific << coherenceMin;
		if (ptrControl->doesOutputApparentResistivityAndPhase()) {
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcApparentResistivity(freq, respOut0[iOut]);
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcPhase(respOut0[iOut]);
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcApparentResistivity(freq, respOut1[iOut]);
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcPhase(respOut1[iOut]);
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << coherenceMin;
		}
	}

	if (ptrControl->getErrorEstimationMethod() == Control::STRICT_BOOTSTRAP) {
		strictBootstrap(iSegLen, freqDegree, timeLength, freq, numSegmentsTotal, ftval, times, ofsResp, ofsRhoaPhs, respOut0, respOut1);
	}
	else {
		ptrOutputFiles->writeErrorMessage("Unsupported error estimation method : " + Util::toString(ptrControl->getErrorEstimationMethod()));
	}

	delete[] istf0;
	delete[] istf1;
	delete[] respOut0;
	delete[] respOut1;

}

void AnalysisTest::calculateResponseFunctionsAux(const int iSegLen, const int freqDegree, const double timeLength, const double freq,
	const int numSegmentsTotal, std::complex<double>** ftval, const std::vector< std::pair<std::string, std::string> >& times,
	std::ofstream& ofsResp, std::ofstream& ofsRhoaPhs, const bool forJackknife, const int out, std::complex<double>& respOut0, std::complex<double>& respOut1,
	double& coherenceOut ) const {

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();

	ptrOutputFiles->writeLogMessage("Calculate response functions by ordinary remote reference");
	ptrOutputFiles->writeCvgMessage("================================================================================");
	ptrOutputFiles->writeCvgMessage("Now Frequency(Hz): " + Util::toString(freq) + ", Period(s): " + Util::toString(1.0 / freq));
	ptrOutputFiles->writeCvgMessage("================================================================================");
	const Control* const ptrControl = Control::getInstance();

	const int numRemoteReferenceVariables = ptrControl->getNumRemoteReferenceVariables();
	assert(numRemoteReferenceVariables >= 2);
	const int rr0 = ptrControl->getChannelIndex(CommonParameters::REMOTE_REFERENCE, 0);
	const int rr1 = ptrControl->getChannelIndex(CommonParameters::REMOTE_REFERENCE, 1);

	double* unitWeights = new double[numSegmentsTotal];
	for (int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg) {
		unitWeights[iSeg] = 1.0;
	}

	ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
	ptrOutputFiles->writeCvgAndLogMessage("Calculate response functions for output variable " + Util::toString(out));
	ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
	// Calculate response functions by the ordinary least square method
	double coherence(0.0);
	std::complex<double> resp0(0.0, 0.0);
	std::complex<double> resp1(0.0, 0.0);
	double* weights = new double[numSegmentsTotal];
	double sumOfWeights(0.0);
	for (int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg) {
		weights[iSeg] = 1.0;
		sumOfWeights += weights[iSeg];
	}
	if (sumOfWeights < CommonParameters::EPS) {
		ptrOutputFiles->writeErrorMessage("Sum of weights is too small: " + Util::toString(sumOfWeights));
	}
	std::complex<double>* residuals = new std::complex<double>[numSegmentsTotal];
	ptrOutputFiles->writeCvgAndLogMessage("Calculate response functions by the ordinary least square method");
	calculateResponseFunctionByWLSRemoteReference(ftval[out], ftval[rr0], ftval[rr1], ftval[rr0], ftval[rr1],
		numSegmentsTotal, weights, residuals, resp0, resp1, coherence);
	coherenceOut = coherence;
	std::vector<std::string> titles;
	std::vector<double>* outputValues = new std::vector<double>[numSegmentsTotal];
	const bool outputResidual = !forJackknife && ptrControl->getOutputLevel() >= 2;
	if (outputResidual) {
		titles.push_back("OLS");
		for (int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg) {
			outputValues[iSeg].push_back(residuals[iSeg].real());
			outputValues[iSeg].push_back(residuals[iSeg].imag());
			outputValues[iSeg].push_back(weights[iSeg]);
		}
	}
	double absMaxResidual(0.0);
	for (int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg) {
		const double absResidual = std::abs(residuals[iSeg]);
		if (absResidual > absMaxResidual) {
			absMaxResidual = absResidual;
		}
	}
	if (absMaxResidual < CommonParameters::EPS) {
		// Such a case where the remote reference field is equivalent to the input field
		ptrOutputFiles->writeCvgMessage("Robust method is not performed because residuals is nearly zero");
	}else{
		double scale = RobustWeight::calculateScaleByMADN(numSegmentsTotal, residuals);
		// Calculate response functions by regression using the first M-estimator
		calculateResponseFunctionsByIRWLSRemoteReference(0, ftval[out], ftval[rr0], ftval[rr1], ftval[rr0], ftval[rr1],
			numSegmentsTotal, false, scale, unitWeights, weights, residuals, resp0, resp1, coherence, titles, outputValues);
		// Calculate response functions by regression using the second M-estimator
		calculateResponseFunctionsByIRWLSRemoteReference(1, ftval[out], ftval[rr0], ftval[rr1], ftval[rr0], ftval[rr1],
			numSegmentsTotal, true, scale, unitWeights, weights, residuals, resp0, resp1, coherence, titles, outputValues);
		coherenceOut = coherence;
	}
	if (!forJackknife && ptrControl->getOutputLevel() > 0) {
		// Output spectral density functions to cvg file
		outputSpectralDensityFunctionsToCvgFile(numSegmentsTotal, timeLength, ftval[out], ftval[rr0], ftval[rr1], weights);
	}
	respOut0 = resp0;
	respOut1 = resp1;

	if (outputResidual) {
		std::ostringstream oss;
		oss << "segm" << iSegLen << "_index" << freqDegree << "_output" << out << "_residuals.csv";
		writeResiduals(oss.str(), numSegmentsTotal, times, titles, outputValues);
	}

	// Release memory
	delete[] residuals;
	delete[] outputValues;
	delete[] weights;
	delete[] unitWeights;

}

// Estimate error by strict bootstrap
void AnalysisTest::strictBootstrap(const int iSegLen, const int freqDegree, const double timeLength, const double freq,
	const int numSegmentsTotal, std::complex<double>** ftval, const std::vector< std::pair<std::string, std::string> >& times,
	std::ofstream& ofsResp, std::ofstream& ofsRhoaPhs, const std::complex<double>* const resp0, const std::complex<double>* const resp1) const {

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeLogMessage("Strict bootstrap is performed to estimate errors");

	const Control* const ptrControl = Control::getInstance();
	const int numOutputVariables = ptrControl->getNumOutputVariables();
	const int numInputVariables = ptrControl->getNumInputVariables();
	assert(numInputVariables == 2);

	// Copy Fourier transformed values
	const int numChannels = ptrControl->getNumberOfChannels();
	std::complex<double>** ftvalForBootstrap = new std::complex<double>*[numChannels];
	for (int iChan = 0; iChan < numChannels; ++iChan) {
		ftvalForBootstrap[iChan] = new std::complex<double>[numSegmentsTotal];
	}
#ifdef _DEBUG_WRITE
	for (int iChan = 0; iChan < numChannels; ++iChan) {
		for (int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg) {
			std::cout << iChan << " " << iSeg << " " << ftval[iChan][iSeg] << std::endl;
		}
	}
#endif

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
#ifdef _DEBUG_WRITE
		for (int iSeg = 0; iSeg < numSegmentsTotal; ++iSeg) {
			std::cout << iSeg << " " << segmentIndexes[iSeg] << std::endl;
		}
#endif
		// Copy data
		for (int iChan = 0; iChan < numChannels; ++iChan) {
			for (int icount = 0; icount < numSegmentsTotal; ++icount) {
				const int iSeg = segmentIndexes[icount];
				ftvalForBootstrap[iChan][icount] = ftval[iChan][iSeg];
			}
		}
		std::complex<double>* istf0 = new std::complex<double>[numOutputVariables + numInputVariables];
		std::complex<double>* istf1 = new std::complex<double>[numOutputVariables + numInputVariables];
		int icount(0);
		for (int iOut = 0; iOut < numOutputVariables; ++iOut, ++icount) {
			const int index = ptrControl->getChannelIndex(CommonParameters::OUTPUT, iOut);
			double coherence(0.0);
			calculateResponseFunctionsAux(iSegLen, freqDegree, timeLength, freq, numSegmentsTotal, ftvalForBootstrap, times,
				ofsResp, ofsRhoaPhs, true, index, istf0[icount], istf1[icount], coherence);
		}
		for (int iInp = 0; iInp < numInputVariables; ++iInp, ++icount) {
			const int index = ptrControl->getChannelIndex(CommonParameters::INPUT, iInp);
			double coherence(0.0);
			calculateResponseFunctionsAux(iSegLen, freqDegree, timeLength, freq, numSegmentsTotal, ftvalForBootstrap, times,
				ofsResp, ofsRhoaPhs, true, index, istf0[icount], istf1[icount], coherence);
		}
#ifdef _DEBUG_WRITE
		for (int i = 0; i < numOutputVariables + numInputVariables; ++i) {
			std::cout << i << " " << iSample << " " << istf0[i] << " " << istf1[i] << std::endl;
		}
#endif
		resp0Sample[iSample] = new std::complex<double>[numOutputVariables];
		resp1Sample[iSample] = new std::complex<double>[numOutputVariables];
		const std::complex<double> Txx = istf0[numOutputVariables];
		const std::complex<double> Txy = istf1[numOutputVariables];
		const std::complex<double> Tyx = istf0[numOutputVariables + 1];
		const std::complex<double> Tyy = istf1[numOutputVariables + 1];
		const std::complex<double> det = Txx * Tyy - Txy * Tyx;
		OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
		if (std::abs(det) < CommonParameters::EPS) {
			ptrOutputFiles->writeErrorMessage("Determinant is too small: " + Util::toString(std::abs(det)));
		}
		for (int iOut = 0; iOut < numOutputVariables; ++iOut) {
			const int index = ptrControl->getChannelIndex(CommonParameters::OUTPUT, iOut);
			const std::complex<double> U_x = istf0[index];
			const std::complex<double> U_y = istf1[index];
			resp0Sample[iSample][iOut] = (U_x * Tyy - U_y * Tyx) / det;
			resp1Sample[iSample][iOut] = (U_y * Txx - U_x * Txy) / det;
		}
		delete [] istf0;
		delete [] istf1;
	}
	ptrOutputFiles->restartToWriteCvgMessage();
	ptrOutputFiles->restartToWriteLogMessage();
	ptrOutputFiles->restartToWriteWarningMessage();

#ifdef _DEBUG_WRITE
	for (int iOut = 0; iOut < numOutputVariables; ++iOut) {
		for (int iSample = 0; iSample < numOfSamples; ++iSample) {
			std::cout << iOut << " " << iSample << " " << resp0Sample[iSample][iOut] << " " << resp1Sample[iSample][iOut] << std::endl;
		}
	}
#endif

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

// Write residuals
void AnalysisTest::writeResiduals(const std::string& fileName, const int numSegmentsTotal,
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