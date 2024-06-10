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
#ifndef DBLDEF_ANALYSIS_TEST
#define DBLDEF_ANALYSIS_TEST

#include "Analysis.h"

// Class of test
class AnalysisTest : public Analysis {

public:

	// Default constructer
	AnalysisTest();

	// Destructer
	virtual ~AnalysisTest();

private:

	// Calculate response functions
	virtual void calculateResponseFunctions( const int iSegLen,const int freqDegree, const double timeLength, const double freq,
		const int numSegmentsTotal, std::complex<double>** ftval, const std::vector< std::pair<std::string, std::string> >& times,
		std::ofstream& ofsResp, std::ofstream& ofsRhoaPhs );

	void parametricErrorEstimation() const;

	void calculateResponseFunctionsAux(const int iSegLen, const int freqDegree, const double timeLength, const double freq,
		const int numSegmentsTotal, std::complex<double>** ftval, const std::vector< std::pair<std::string, std::string> >& times,
		std::ofstream& ofsResp, std::ofstream& ofsRhoaPhs, const bool forJackknife, const int out, std::complex<double>& respOut0, std::complex<double>& respOut1,
		double& coherenceOut ) const;

	// Estimate errors by strict bootstrap
	void strictBootstrap(const int iSegLen, const int freqDegree, const double timeLength, const double freq,
		const int numSegmentsTotal, std::complex<double>** ftval, const std::vector< std::pair<std::string, std::string> >& times,
		std::ofstream& ofsResp, std::ofstream& ofsRhoaPhs, const std::complex<double>* const resp0, const std::complex<double>* const resp1 ) const;

	// Write residuals
	void writeResiduals(const std::string& fileName, const int numSegmentsTotal,
		const std::vector< std::pair<std::string, std::string> >& times, const std::vector<std::string>& titles, const std::vector<double>* outputValues) const;

};

#endif
