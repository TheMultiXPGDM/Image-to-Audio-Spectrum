# Image-to-Audio-Spectrum
Convert any image into a WAV audio file containing the image on its spectrum

## Building
In order to compile this project, you'll need:
- CMake
- C++ Compiler (MSVC, MinGW)

To build this project, you only have to run build.bat

## Usage
--min-freq (Minimum Frequency: Hz)
The frequency at which the image will start drawing.

--max-freq (Maximum Frequency: Hz)
The frequency at which the image will finish drawing.

--duration (Duration: seconds)
The duration of the resulting audio file.

--sample-rate (Sample Rate: Hz)
The sample rate of the resulting audio file.

--gamma (Gamma: float)
Change the contrast of the image.

--window_size (Window Size: int (power of 2))
Size of the FFT window.
Bigger window will increase resolution in the frequency domain (vertical) but decrease resolution in the time domain (horizontal)

--window_overlap (Window Overlap: int (power of 2))
How many windows will overlap.
More overlaps will result in more high quality spectrum (less spectral leakage) but it will also take more time to process.

--input (Input Image: '.png', '.jpg', '.bmp')
Path to the image file (remember to add the extension)

--output (Output Wave File: '.wav')
Path to the output file (remember to add the extension)
