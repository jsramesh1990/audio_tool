Simple C++ audio tool for tuning (pitch shift via resampling) and filtering (biquad filters including low-shelf, high-pass, peaking EQ). Processes WAV files using libsndfile and writes output WAV. Not real-time — offline processing for clarity and portability.

Requirements
-----------
- C++17
- CMake
- libsndfile (for WAV I/O)

On Debian/Ubuntu:
  sudo apt update && sudo apt install build-essential cmake libsndfile1-dev

Build
-----
mkdir build && cd build
cmake ..
make
