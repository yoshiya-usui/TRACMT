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
#ifndef DBLDEF_OUTPUT_FILES
#define DBLDEF_OUTPUT_FILES

#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <ctime>

// Class of output files
class OutputFiles{

public:
	// Return the the instance of the class
	static OutputFiles* getInstance();

	// Restart to write cvg message
	void restartToWriteCvgMessage();

	// Restart to write log message
	void restartToWriteLogMessage();

	// Restart to write warning essage
	void restartToWriteWarningMessage();

	// Set flag specifing whether messages are shown to console
	void setOutputToConsole(const bool outputToConsole);

	// Stop to write cvg message
	void stopToWriteCvgMessage();

	// Stop to write log message
	void stopToWriteLogMessage();

	// Stop to write warning essage
	void stopToWriteWarningMessage();

	// Write to cvg message
	void writeCvgMessage( const std::string& msg );

	// Write to cvg message
	void writeCvgAndLogMessage( const std::string& msg, const bool withElapsedTime = true );

	// Write log message
	void writeLogMessage( const std::string& msg, const bool withElapsedTime = true );

	// Write error message
	void writeErrorMessage( const std::string& msg );

	// Write warning message
	void writeWarningMessage( const std::string& msg );

private:
	// Conrvergence file
	std::ofstream m_cvgFile;

	// Log file
	std::ofstream m_logFile;

	// Flag specifing whether writing cvg message is stopped
	bool m_isWritingCvgMessageStopped;

	// Flag specifing whether writing log message is stopped
	bool m_isWritingLogMessageStopped;

	// Flag specifing whether writing warning message is stopped
	bool m_isWritingWarningMessageStopped;

	// Flag specifing whether messages are shown to console
	bool m_outputToConsole;

	// The time the class instanced
	time_t m_startTime;

	// Constructer
	OutputFiles();

	// Destructer
	~OutputFiles();

	// Copy constructer
	OutputFiles(const OutputFiles& rhs);

	// Assignment operator
	OutputFiles& operator=(const OutputFiles& rhs);

	// Calculate and output elapsed time
	std::string outputElapsedTime() const;
		
};

#endif
