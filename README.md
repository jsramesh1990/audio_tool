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

If successful, you’ll get an executable:
./audio_tool

Usage
-----
./audio_tool -i input.wav -o output.wav [--pitch <semitones>] [--lowshelf <freq> <gain_db> <Q>] [--peak <freq> <gain_db> <Q>] [--hp <freq> <Q>]

Examples
--------
# Raise pitch by 3 semitones and add bass boost at 100 Hz +6 dB
./audio_tool -i in.wav -o out.wav --pitch 3 --lowshelf 100 6 0.7

# Apply a peaking EQ at 1000 Hz -3 dB
./audio_tool -i in.wav -o out.wav --peak 1000 -3 1.0
