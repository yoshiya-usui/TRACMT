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
#ifndef DBLDEF_COMMON_PARAMETERS
#define DBLDEF_COMMON_PARAMETERS

#include <string>
#include <vector>

namespace CommonParameters{

enum DataType{
	OUTPUT = 0,
	INPUT,
	REMOTE_REFERENCE
};

struct DataFile{
	std::string fileName;
	std::string mth5GroupName;
	int numSkipData;
	double* data;
};

struct DataFileSet{
	int numDataPoints;
	std::vector<DataFile> dataFile;
};

// Circular constant
const static double PI = 3.14159265358979323846;

// Factor converting values from radians to degrees
const static double RAD2DEG = 180.0 / PI;

// Factor converting values from degrees to radians
const static double DEG2RAD = PI / 180.0;

const static double EPS = 1.0e-20;

static char programName[]="TRACMT";

static char version[] = "v2.3.0";

}

#endif
