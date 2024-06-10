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
#include <iostream>
#include <sstream>
#include <iomanip>

#include "CommonParameters.h"
#include "OutputFiles.h"

// Return the the instance of the class
OutputFiles* OutputFiles::getInstance(){
   	static OutputFiles instance;// The only instance
  	return &instance;
}

// Restart to write cvg message
void OutputFiles::restartToWriteCvgMessage(){
	m_isWritingCvgMessageStopped = false;
}

// Restart to write log message
void OutputFiles::restartToWriteLogMessage(){
	m_isWritingLogMessageStopped = false;
}

// Restart to write warning essage
void OutputFiles::restartToWriteWarningMessage(){
	m_isWritingWarningMessageStopped = false;
}

// Set flag specifing whether messages are shown to console
void OutputFiles::setOutputToConsole(const bool outputToConsole) {
	m_outputToConsole = outputToConsole;
}

// Stop to write cvg message
void OutputFiles::stopToWriteCvgMessage(){
	m_isWritingCvgMessageStopped = true;
}

// Stop to write log message
void OutputFiles::stopToWriteLogMessage(){
	m_isWritingLogMessageStopped = true;
}

// Stop to write warning essage
void OutputFiles::stopToWriteWarningMessage(){
	m_isWritingWarningMessageStopped = true;
}

// Write to cvg message
void OutputFiles::writeCvgMessage( const std::string& msg ){
	if(m_isWritingCvgMessageStopped){
		return;
	}
	m_cvgFile << msg << std::endl;
}

// Write to cvg message
void OutputFiles::writeCvgAndLogMessage( const std::string& msg, const bool withElapsedTime ){
	writeLogMessage(msg, withElapsedTime);
	writeCvgMessage(msg);
}

// Write log message
void OutputFiles::writeLogMessage( const std::string& msg, const bool withElapsedTime ){
	if(m_isWritingLogMessageStopped){
		return;
	}
	if( withElapsedTime ){
		m_logFile << msg << " " << outputElapsedTime() << std::endl;
	}else{
		m_logFile << msg << std::endl;
	}
	if (m_outputToConsole) {
		if (withElapsedTime) {
			std::cout << msg << " " << outputElapsedTime() << std::endl;
		}
		else {
			std::cout << msg << std::endl;
		}
	}
}

// Write error message
void OutputFiles::writeErrorMessage( const std::string& msg ){
	m_logFile << "[Error] " << msg << std::endl;
	if (m_outputToConsole) {
		std::cout << "[Error] " << msg << std::endl;
	}
	exit(1);
}

// Write warning message
void OutputFiles::writeWarningMessage( const std::string& msg ){
	if(m_isWritingWarningMessageStopped){
		return;
	}
	m_logFile << "[Warning] " << msg << std::endl;
	if (m_outputToConsole) {
		std::cout << "[Warning] " << msg << std::endl;
	}
}

// Constructer
OutputFiles::OutputFiles():
	m_isWritingCvgMessageStopped(false),
	m_isWritingLogMessageStopped(false),
	m_isWritingWarningMessageStopped(false),
	m_outputToConsole(false),
	m_startTime(NULL)
{

	// Open cvg file
	std::ostringstream cvgFileName;
	cvgFileName << CommonParameters::programName <<".cvg";
	m_cvgFile.open( cvgFileName.str().c_str(), std::ios::out );
	if( m_cvgFile.fail() ){
		std::cerr << "File open error !! : " << cvgFileName.str() << std::endl;
		exit(1);
	}

	// Open log file
	std::ostringstream logFileName;
	logFileName << CommonParameters::programName <<".log";
	m_logFile.open( logFileName.str().c_str(), std::ios::out );
	if( m_logFile.fail() ){
		std::cerr << "File open error !! : " << logFileName.str() << std::endl;
		exit(1);
	}

	// Measure the start time
	time(&m_startTime);

	OutputFiles::m_logFile << "Start " << CommonParameters::programName << " Version " << CommonParameters::version << std::endl;

}

// Destructer
OutputFiles::~OutputFiles(){

	if( m_cvgFile.is_open() ){
		m_cvgFile.close();
	}

	if( m_logFile.is_open() ){
		m_logFile.close();
	}

}

// Copy constructer
OutputFiles::OutputFiles(const OutputFiles& rhs){
	std::cerr << "Error : Copy constructer of the class LogFile is not implemented." << std::endl;
	exit(1);
}

// Copy assignment operator
OutputFiles& OutputFiles::operator=(const OutputFiles& rhs){
	std::cerr << "Error : Copy assignment operator of the class LogFile is not implemented." << std::endl;
	exit(1);
}

// Calculate elapsed time
std::string OutputFiles::outputElapsedTime() const{

	time_t curTime(NULL);
	time(&curTime);

	std::ostringstream output;
	output << "( " << difftime(curTime, m_startTime) << " sec )";

	return output.str();

}
