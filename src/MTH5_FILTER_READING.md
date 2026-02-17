# MTH5 Filter Reading Functionality

## Overview

The MTH5 class now supports reading filter information from MTH5 files for channel calibration. This functionality mirrors the Python implementation in `mth5.groups.channel_dataset.ChannelDataset.channel_response`.

## New Data Structure

### FilterInfo

A structure to hold filter information:

```cpp
struct FilterInfo {
    std::string name;           // Filter name
    std::string type;           // Filter type (zpk, coefficient, time_delay, fap, fir)
    std::string units_in;       // Input units
    std::string units_out;      // Output units
    int sequence_number;        // Order in filter chain
    std::vector<double> data;   // Generic data array (for future use)
};
```

## New Public Methods

### 1. `getChannelFilters()` - Recommended Method

Gets all filters for a channel, similar to Python's `channel.channel_response` property.

```cpp
std::vector<FilterInfo> getChannelFilters(
    const std::string& fileName,     // MTH5 file path
    const std::string& channelPath   // Full path to channel dataset
) const;
```

**Example:**
```cpp
MTH5* mth5 = MTH5::getInstance();
std::vector<FilterInfo> filters = mth5->getChannelFilters(
    "data.mth5", 
    "/Survey/Stations/MT001/MT001a/Ex"
);
```

### 2. `readFilterNamesFromChannel()`

Reads filter names from channel dataset attributes.

```cpp
std::vector<std::string> readFilterNamesFromChannel(
    const std::string& fileName,
    const std::string& channelPath
) const;
```

**Example:**
```cpp
std::vector<std::string> names = mth5->readFilterNamesFromChannel(
    "data.mth5",
    "/Survey/Stations/MT001/MT001a/Ex"
);
```

### 3. `readFilter()`

Reads a specific filter from the Filters group.

```cpp
FilterInfo readFilter(
    const std::string& fileName,
    const std::string& filterPath,     // Path to Filters group
    const std::string& filterName,
    int sequenceNumber
) const;
```

**Example:**
```cpp
FilterInfo filter = mth5->readFilter(
    "data.mth5",
    "/Survey/Filters",
    "dipole_104.0m",
    1
);
```

## MTH5 Filter Architecture

### File Structure

MTH5 files store filters in a hierarchical structure:

```
/Survey/
├── Stations/
│   └── {station_id}/
│       └── {run_id}/
│           └── {channel}          <- Channel dataset
│               └── @filter.name   <- Attribute with filter names
└── Filters/                        <- Filters group
    ├── zpk/                        <- Zero-pole-gain filters
    │   └── {filter_name}
    ├── coefficient/                <- Coefficient filters
    │   └── {filter_name}
    ├── time_delay/                 <- Time delay filters
    │   └── {filter_name}
    ├── fap/                        <- Frequency-amplitude-phase filters
    │   └── {filter_name}
    └── fir/                        <- FIR filters
        └── {filter_name}
```

Or for v0.2.0 files:

```
/Experiment/
└── Surveys/
    └── {survey_id}/
        ├── Stations/
        │   └── {station_id}/
        │       └── {run_id}/
        │           └── {channel}
        └── Filters/
            ├── zpk/
            ├── coefficient/
            ├── time_delay/
            ├── fap/
            └── fir/
```

### Filter Types

1. **zpk** - Zeros, Poles, and Gain
   - Represents transfer functions
   - Common for instrument responses
   - Contains zeros, poles, and gain values

2. **coefficient** - FIR Coefficients
   - Finite Impulse Response filter
   - Array of filter coefficients
   - Used for simple filtering operations

3. **time_delay** - Time Delay
   - Represents time shift
   - Single delay value in seconds

4. **fap** - Frequency-Amplitude-Phase
   - Lookup table format
   - Arrays of frequency, amplitude, and phase

5. **fir** - Finite Impulse Response
   - Similar to coefficient type
   - More detailed metadata

## Usage Workflow

### Basic Usage

```cpp
#include "MTH5.h"
#include <iostream>

int main() {
    MTH5* mth5 = MTH5::getInstance();
    
    // Get all filters for a channel
    std::vector<FilterInfo> filters = mth5->getChannelFilters(
        "my_data.mth5",
        "/Survey/Stations/MT001/MT001a/Ex"
    );
    
    // Process filters
    for (const auto& filter : filters) {
        std::cout << "Filter: " << filter.name << std::endl;
        std::cout << "  Type: " << filter.type << std::endl;
        std::cout << "  Sequence: " << filter.sequence_number << std::endl;
        std::cout << "  " << filter.units_in << " -> " 
                  << filter.units_out << std::endl;
    }
    
    return 0;
}
```

### Calibration Workflow

```cpp
// 1. Read time series data
double* data = new double[npts];
mth5->readMTH5File("data.mth5", channelPath, 0, npts, data);

// 2. Get filters
std::vector<FilterInfo> filters = mth5->getChannelFilters("data.mth5", channelPath);

// 3. Apply filters in sequence based on sequence_number
// Typically done in frequency domain:
//    - FFT time series
//    - Apply each filter's transfer function
//    - Inverse FFT
for (const auto& filter : filters) {
    if (filter.type == "zpk") {
        // Apply zeros-poles-gain transfer function
    } else if (filter.type == "coefficient") {
        // Apply FIR coefficients
    } else if (filter.type == "time_delay") {
        // Apply time shift
    }
    // ... additional filter types
}

delete[] data;
```

## Implementation Details

### Filter Name Normalization

Filter names are normalized to match the Python implementation:
- Forward slashes (`/`) replaced with `" per "`
- Converted to lowercase

Example: `"V/m"` becomes `"v per m"`

### Attribute Reading

The implementation reads the `filter.name` attribute from channel datasets. This is stored as an array of strings in HDF5.

For newer MTH5 file formats that use structured filter objects, the implementation can be extended to parse those structures.

### Error Handling

The methods use the existing `OutputFiles` class for logging:
- **Error messages**: Critical failures (file not found, dataset not accessible)
- **Warning messages**: Non-critical issues (filter not found)
- **Log messages**: Informational (filters found, reading progress)

## Comparison with Python

The C++ implementation parallels the Python version:

**Python (`channel_dataset.py`):**
```python
@property
def channel_response(self) -> ChannelResponse:
    filters_group = FiltersGroup(
        self.hdf5_dataset.parent.parent.parent.parent["Filters"]
    )
    f_list = []
    for count, name in enumerate(self.metadata.filter.name, 1):
        name = name.replace("/", " per ").lower()
        filt_obj = filters_group.to_filter_object(name)
        filt_obj.sequence_number = count
        f_list.append(filt_obj)
    return ChannelResponse(filters_list=f_list)
```

**C++ (MTH5.cpp):**
```cpp
std::vector<FilterInfo> MTH5::getChannelFilters(
    const std::string& fileName, 
    const std::string& channelPath
) const {
    std::vector<FilterInfo> filters;
    
    // Read filter names from channel
    std::vector<std::string> filterNames = 
        readFilterNamesFromChannel(fileName, channelPath);
    
    // Get Filters group path
    std::string filtersPath = getFiltersGroupPath(channelPath);
    
    // Read each filter
    for (size_t i = 0; i < filterNames.size(); i++) {
        FilterInfo filter = readFilter(fileName, filtersPath, 
                                       filterNames[i], i + 1);
        filters.push_back(filter);
    }
    
    return filters;
}
```

## Future Enhancements

### Reading Filter Data

The current implementation reads filter metadata. To fully support calibration, extend `readFilter()` to read filter-specific datasets:

```cpp
// For ZPK filters
if (ftype == "zpk") {
    // Read zeros, poles, gain datasets
    H5LTread_dataset_double(filter_id, "zeros", ...);
    H5LTread_dataset_double(filter_id, "poles", ...);
    // Store in FilterInfo.data or create ZPKFilter struct
}

// For coefficient filters  
if (ftype == "coefficient") {
    // Read coefficients dataset
    H5LTread_dataset_double(filter_id, "coefficients", ...);
}
```

### Structured Filter Types

Create specialized filter structures:

```cpp
struct ZPKFilter : public FilterInfo {
    std::vector<std::complex<double>> zeros;
    std::vector<std::complex<double>> poles;
    double gain;
};

struct CoefficientFilter : public FilterInfo {
    std::vector<double> coefficients;
};
```

### Filter Application

Implement filter application methods:

```cpp
void applyFilters(
    double* data, 
    int npts, 
    const std::vector<FilterInfo>& filters
);
```

## See Also

- **Python Reference**: `mth5/groups/channel_dataset.py` lines 343-410
- **Filter Groups**: `mth5/groups/filters.py`
- **Example Usage**: `MTH5_filter_example.cpp`
