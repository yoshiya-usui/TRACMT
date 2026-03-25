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
#include "MTH5Filter.h"
#include "OutputFiles.h"

// Default constructer
MTH5Filter::MTH5Filter() :
	m_sequenceNumber(0)
{
}

// Destructer
MTH5Filter::~MTH5Filter(){
}

// Copy constructer
MTH5Filter::MTH5Filter(const MTH5Filter& rhs) :
	m_name(rhs.m_name),
	m_type(rhs.m_type),
	m_unitsIn(rhs.m_unitsIn),
	m_unitsOut(rhs.m_unitsOut),
	m_sequenceNumber(rhs.m_sequenceNumber),
	m_data(rhs.m_data)
{
}

// Assignment operator
MTH5Filter& MTH5Filter::operator=(const MTH5Filter& rhs) {

	if (this == &rhs) {
		return *this;
	}
	m_name = rhs.m_name;
	m_type = rhs.m_type;
	m_unitsIn = rhs.m_unitsIn;
	m_unitsOut = rhs.m_unitsOut;
	m_sequenceNumber = rhs.m_sequenceNumber;
	m_data = rhs.m_data;

	return *this;

}

// Get filter name
std::string MTH5Filter::getName() const {
	return m_name;
}

// Get filter type (zpk, coefficient, time_delay, fap, fir)
std::string MTH5Filter::getType() const {
	return m_type;
}

// Input units
std::string MTH5Filter::getUnitsIn() const {
	return m_unitsIn;
}

// Output units
std::string MTH5Filter::getUnitsOut() const {
	return m_unitsOut;
}

// Get order in filter chain
int MTH5Filter::getSequenceNumber() const {
	return m_sequenceNumber;
}

// Set filter name
void MTH5Filter::setName(const std::string& name) {
	m_name = name;
}

// Set filter type (zpk, coefficient, time_delay, fap, fir)
void MTH5Filter::setType(const std::string& type) {
	m_type = type;
}

// Set input units
void MTH5Filter::setUnitsIn(const std::string& unitsIn) {
	m_unitsIn = unitsIn;
}

// Set output units
void MTH5Filter::setUnitsOut(const std::string& unitsOut) {
	m_unitsOut = unitsOut;
}

// Set order in filter chain
void MTH5Filter::setSequenceNumber(const int sequenceNumber) {
	m_sequenceNumber = sequenceNumber;
}

// Get frequency response functions using the requency response functions of filter
std::complex<double> MTH5Filter::getFrequencyResponse(const double freq) const {
	OutputFiles::getInstance()->writeWarningMessage("Filter is not properly defined");
	return 1.0;
}
