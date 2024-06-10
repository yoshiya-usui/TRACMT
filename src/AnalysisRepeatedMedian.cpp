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
#include "AnalysisRepeatedMedian.h"
#include "Control.h"
#include "OutputFiles.h"

#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <iomanip>

#include "Util.h"


// Default constructer
AnalysisRepeatedMedian::AnalysisRepeatedMedian()
{
}

// Destructer
AnalysisRepeatedMedian::~AnalysisRepeatedMedian()
{
}

// Calculate response functions
void AnalysisRepeatedMedian::calculateResponseFunctions(const int iSegLen, const int freqDegree, const double timeLength, const double freq,
	const int numSegmentsTotal, std::complex<double>** ftval, const std::vector< std::pair<std::string, std::string> >& times,
	std::ofstream& ofsResp, std::ofstream& ofsRhoaPhs) {

	OutputFiles* ptrOutputFiles = OutputFiles::getInstance();
	ptrOutputFiles->writeLogMessage("Calculate response functions by repeated median estimator");
	ptrOutputFiles->writeCvgMessage("================================================================================");
	ptrOutputFiles->writeCvgMessage("Now Frequency(Hz): " + Util::toString(freq) + ", Period(s): " + Util::toString(1.0 / freq));
	ptrOutputFiles->writeCvgMessage("================================================================================");
	const Control* const ptrControl = Control::getInstance();
	const int in0 = ptrControl->getChannelIndex(CommonParameters::INPUT, 0);
	const int in1 = ptrControl->getChannelIndex(CommonParameters::INPUT, 1);
	const int numRemoteReferenceVariables = ptrControl->getNumRemoteReferenceVariables();
	assert(numRemoteReferenceVariables >= 2);
	const int rr0 = ptrControl->getChannelIndex(CommonParameters::REMOTE_REFERENCE, 0);
	const int rr1 = ptrControl->getChannelIndex(CommonParameters::REMOTE_REFERENCE, 1);
	const int numOutputVariables = ptrControl->getNumOutputVariables();

	ofsResp << std::setprecision(10) << std::scientific << freq;
	ofsResp << "," << std::setprecision(10) << std::scientific << 1.0 / freq;
	if (ptrControl->doesOutputApparentResistivityAndPhase()) {
		ofsRhoaPhs << std::setprecision(10) << std::scientific << freq;
		ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << 1.0 / freq;
	}

	std::complex<double>* respOut0 = new std::complex<double>[numOutputVariables];
	std::complex<double>* respOut1 = new std::complex<double>[numOutputVariables];
	double* dResp0 = new double[numOutputVariables];
	double* dResp1 = new double[numOutputVariables];
	for (int iOutVar = 0; iOutVar < numOutputVariables; ++iOutVar) {
		ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
		ptrOutputFiles->writeCvgAndLogMessage("Calculate response functions for output variable " + Util::toString(iOutVar));
		ptrOutputFiles->writeCvgMessage("--------------------------------------------------------------------------------");
		const double weights[2] = { 1.0, 1.0 };
		std::complex<double> ftvalOut[2];
		std::complex<double> ftvalIn0[2];
		std::complex<double> ftvalIn1[2];
		std::complex<double> ftvalRr0[2];
		std::complex<double> ftvalRr1[2];
		double* resp0RealMedian = new double[numSegmentsTotal];
		double* resp0ImagMedian = new double[numSegmentsTotal];
		double* resp1RealMedian = new double[numSegmentsTotal];
		double* resp1ImagMedian = new double[numSegmentsTotal];
		const int out = ptrControl->getChannelIndex(CommonParameters::OUTPUT, iOutVar);
		for (int iSeg0 = 0; iSeg0 < numSegmentsTotal; ++iSeg0) {
			double* resp0Real = new double[numSegmentsTotal - 1];
			double* resp0Imag = new double[numSegmentsTotal - 1];
			double* resp1Real = new double[numSegmentsTotal - 1];
			double* resp1Imag = new double[numSegmentsTotal - 1];
			int icount(0);
			ftvalOut[0] = ftval[out][iSeg0];
			ftvalIn0[0] = ftval[in0][iSeg0];
			ftvalIn1[0] = ftval[in1][iSeg0];
			ftvalRr0[0] = ftval[rr0][iSeg0];
			ftvalRr1[0] = ftval[rr1][iSeg0];
			for (int iSeg1 = 0; iSeg1 < numSegmentsTotal; ++iSeg1) {
				if (iSeg0 != iSeg1) {
					ftvalOut[1] = ftval[out][iSeg1];
					ftvalIn0[1] = ftval[in0][iSeg1];
					ftvalIn1[1] = ftval[in1][iSeg1];
					ftvalRr0[1] = ftval[rr0][iSeg1];
					ftvalRr1[1] = ftval[rr1][iSeg1];
					std::complex<double> resp0(0.0, 0.0);
					std::complex<double> resp1(0.0, 0.0);
					calculateResponseFunctionByWLSRemoteReferenceAux(2, ftvalOut, ftvalIn0, ftvalIn1, ftvalRr0, ftvalRr1,
						weights, resp0, resp1);
					resp0Real[icount] = resp0.real();
					resp0Imag[icount] = resp0.imag();
					resp1Real[icount] = resp1.real();
					resp1Imag[icount] = resp1.imag();
					++icount;
				}
			}
			resp0RealMedian[iSeg0] = Util::calculateMedian(numSegmentsTotal - 1, resp0Real);
			resp0ImagMedian[iSeg0] = Util::calculateMedian(numSegmentsTotal - 1, resp0Imag);
			resp1RealMedian[iSeg0] = Util::calculateMedian(numSegmentsTotal - 1, resp1Real);
			resp1ImagMedian[iSeg0] = Util::calculateMedian(numSegmentsTotal - 1, resp1Imag);
			delete[] resp0Real;
			delete[] resp0Imag;
			delete[] resp1Real;
			delete[] resp1Imag;
		}
		const double respOut0Real = Util::calculateMedian(numSegmentsTotal, resp0RealMedian);
		const double respOut0Imag = Util::calculateMedian(numSegmentsTotal, resp0ImagMedian);
		const double respOut1Real = Util::calculateMedian(numSegmentsTotal, resp1RealMedian);
		const double respOut1Imag = Util::calculateMedian(numSegmentsTotal, resp1ImagMedian);
		// Output results
		const double coherence = 1.0;
		ofsResp << "," << std::setprecision(10) << std::scientific << respOut0Real;
		ofsResp << "," << std::setprecision(10) << std::scientific << respOut0Imag;
		ofsResp << "," << std::setprecision(10) << std::scientific << respOut1Real;
		ofsResp << "," << std::setprecision(10) << std::scientific << respOut1Imag;
		ofsResp << "," << std::setprecision(10) << std::scientific << coherence;
		respOut0[iOutVar] = std::complex<double>(respOut0Real, respOut0Imag);
		respOut1[iOutVar] = std::complex<double>(respOut1Real, respOut1Imag);
		if (ptrControl->doesOutputApparentResistivityAndPhase()) {
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcApparentResistivity(freq, respOut0[iOutVar]);
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcPhase(respOut0[iOutVar]);
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcApparentResistivity(freq, respOut1[iOutVar]);
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcPhase(respOut1[iOutVar]);
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << coherence;
		}
		const double factor = 1.483 / sqrt(numSegmentsTotal);
		const double dResp0Real = factor * Util::calculateMAD(numSegmentsTotal, respOut0Real, resp0RealMedian);
		const double dResp0Imag = factor * Util::calculateMAD(numSegmentsTotal, respOut0Imag, resp0ImagMedian);
		const double dResp1Real = factor * Util::calculateMAD(numSegmentsTotal, respOut1Real, resp1RealMedian);
		const double dResp1Imag = factor * Util::calculateMAD(numSegmentsTotal, respOut1Imag, resp1ImagMedian);
		dResp0[iOutVar] = sqrt(0.5 * (pow(dResp0Real, 2) + pow(dResp0Imag, 2)));
		dResp1[iOutVar] = sqrt(0.5 * (pow(dResp1Real, 2) + pow(dResp1Imag, 2)));
		delete[] resp0RealMedian;
		delete[] resp0ImagMedian;
		delete[] resp1RealMedian;
		delete[] resp1ImagMedian;
	}
	for (int iOutVar = 0; iOutVar < numOutputVariables; ++iOutVar) {
		ofsResp << "," << std::setprecision(10) << std::scientific << dResp0[iOutVar];
		ofsResp << "," << std::setprecision(10) << std::scientific << dResp1[iOutVar];
		if (ptrControl->doesOutputApparentResistivityAndPhase()) {
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcApparentResistivityError(freq, respOut0[iOutVar], dResp0[iOutVar]);
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcPhaseError(respOut0[iOutVar], dResp0[iOutVar]);
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcApparentResistivityError(freq, respOut1[iOutVar], dResp1[iOutVar]);
			ofsRhoaPhs << "," << std::setprecision(10) << std::scientific << calcPhaseError(respOut1[iOutVar], dResp1[iOutVar]);
		}
	}
	ofsResp << std::endl;
	ofsResp.flush();
	if (ptrControl->doesOutputApparentResistivityAndPhase()) {
		ofsRhoaPhs << std::endl;
		ofsRhoaPhs.flush();
	}
	delete[] respOut0;
	delete[] respOut1;
	delete[] dResp0;
	delete[] dResp1;

}
