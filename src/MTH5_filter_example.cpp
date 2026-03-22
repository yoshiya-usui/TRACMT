// Example usage of MTH5 filter reading functionality
// This demonstrates how to read filters from an MTH5 file for channel calibration

#include "MTH5.h"
#include <iostream>
#include <vector>

int main() {
    MTH5* mth5 = MTH5::getInstance();
    
    // Example MTH5 file and channel path
    std::string fileName = "example_data.mth5";
    std::string channelPath = "/Survey/Stations/MT001/MT001a/Ex";
    
    // Method 1: Get all filters for a channel (recommended)
    // This is similar to Python's channel.channel_response property
    std::vector<FilterInfo> filters = mth5->getChannelFilters(fileName, channelPath);
    
    std::cout << "Channel: " << channelPath << std::endl;
    std::cout << "Number of filters: " << filters.size() << std::endl;
    std::cout << std::endl;
    
    // Display filter information
    for (size_t i = 0; i < filters.size(); i++) {
        std::cout << "Filter " << (i + 1) << ":" << std::endl;
        std::cout << "  Name: " << filters[i].name << std::endl;
        std::cout << "  Type: " << filters[i].type << std::endl;
        std::cout << "  Sequence: " << filters[i].sequence_number << std::endl;
        std::cout << "  Units In: " << filters[i].units_in << std::endl;
        std::cout << "  Units Out: " << filters[i].units_out << std::endl;
        std::cout << std::endl;
    }
    
    // Method 2: Read filter names separately (if you need just the names)
    std::vector<std::string> filterNames = mth5->readFilterNamesFromChannel(fileName, channelPath);
    
    std::cout << "Filter names from channel attribute:" << std::endl;
    for (const auto& name : filterNames) {
        std::cout << "  - " << name << std::endl;
    }
    std::cout << std::endl;
    
    // Method 3: Read a specific filter individually
    if (!filterNames.empty()) {
        std::string filtersPath = "/Survey/Filters";  // Adjust based on file version
        FilterInfo singleFilter = mth5->readFilter(fileName, filtersPath, 
                                                   filterNames[0], 1);
        
        std::cout << "Single filter details:" << std::endl;
        std::cout << "  Name: " << singleFilter.name << std::endl;
        std::cout << "  Type: " << singleFilter.type << std::endl;
    }
    
    /* 
     * Using filters for calibration:
     * 
     * 1. Read the channel time series data using readMTH5File()
     * 2. Get the filters using getChannelFilters()
     * 3. Apply filters in sequence based on sequence_number
     * 4. Filter types and their typical use:
     *    - zpk: Zeros-poles-gain transfer function (instrument response)
     *    - coefficient: FIR filter coefficients
     *    - time_delay: Time shift correction
     *    - fap: Frequency-amplitude-phase lookup table
     *    - fir: Finite impulse response filter
     * 
     * For calibration, filters are typically applied in the frequency domain:
     * 1. FFT the time series
     * 2. Apply each filter's transfer function
     * 3. Inverse FFT back to time domain
     */
    
    return 0;
}

/*
 * Compilation example (adjust paths as needed):
 * 
 * g++ -std=c++11 MTH5_filter_example.cpp MTH5.cpp Util.cpp OutputFiles.cpp \
 *     -I/path/to/hdf5/include -L/path/to/hdf5/lib -lhdf5 -o filter_example
 * 
 * 
 * Expected output:
 * 
 * Channel: /Survey/Stations/MT001/MT001a/Ex
 * Number of filters: 3
 * 
 * Filter 1:
 *   Name: dipole_104.0m
 *   Type: coefficient
 *   Sequence: 1
 *   Units In: volts
 *   Units Out: millivolts per kilometer
 * 
 * Filter 2:
 *   Name: v to counts (24v)
 *   Type: zpk
 *   Sequence: 2
 *   Units In: millivolts per kilometer
 *   Units Out: counts
 * 
 * Filter 3:
 *   Name: counts to v (24v)
 *   Type: zpk
 *   Sequence: 3
 *   Units In: counts
 *   Units Out: millivolts per kilometer
 */
