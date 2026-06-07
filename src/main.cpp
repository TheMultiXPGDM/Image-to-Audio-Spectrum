#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <cstring>

#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>
#define DR_WAV_IMPLEMENTATION
#include <dr_wav.h>

#ifndef M_PIf
    #define M_PIf 3.14159265358979323846f
#endif

int log2int(int n) {
    if (n == 0) return 0;
    // check if __builtin_clz is available
#if defined(__GNUC__) || defined(__clang__)
    return 31 - __builtin_clz(n);
#elif defined(_MSC_VER)
    unsigned long index;
    _BitScanReverse(&index, n);
    return index;
#else
    int log = 0;
    if (n >= 1 << 16) { n >>= 16; log |= 16; }
    if (n >= 1 <<  8) { n >>=  8; log |=  8; }
    if (n >= 1 <<  4) { n >>=  4; log |=  4; }
    if (n >= 1 <<  2) { n >>=  2; log |=  2; }
    if (n >= 1 <<  1) { log |=  1; }
    return log;
#endif
}

// asumes that both a and b are powers of 2
int fastDivPow2(int a, int b) {
    if (b == 0) {
        printf("Error: division by zero\n");
        return 0;
    }
    if (a == b) return 1;
    return a >> log2int(b);
}

int exp2int(int n) {
    return 1 << n;
}

int forceExp2int(int n) {
    if (n == 0) return 1;
    // check if n isn't already a power of 2
    if (n & (n - 1)) {
        printf("Warning: %d is not a power of 2, rounding up to the next power of 2\n", n);
        return exp2int(log2int(n) + 1);
    }
    return n;
}

struct RingBuffer {
    float *buffer; // pointer to the buffer
    int capacity = 0; // size of the buffer
    int read = 0; // start of the valid data
    int write = 0; // next position to write to

    RingBuffer(int capacity) : capacity(capacity) {
        buffer = new float[capacity];
        for (int i = 0; i < capacity; i++) {
            buffer[i] = 0.0f;
        }
    }

    ~RingBuffer() {
        delete[] buffer;
    }

    void push(const float value) {
        buffer[write] = value;
        write = (write + 1) % capacity;
    }

    float pop() {
        float value = buffer[read];
        read = (read + 1) % capacity;
        return value;
    }


    void set (const int index, const float value) {
        buffer[(write + index) % capacity] = value;
    }

    void add (const int index, const float value) {
        buffer[(write + index) % capacity] += value;
    }

    float peek(const int index) const {
        return buffer[(read + index) % capacity];
    }

    void jump_read(const int offset) {
        read = (read + offset) % capacity;
    }

    void jump_write(const int offset) {
        write = (write + offset) % capacity;
    }
};

int reverseBits(int n, int bits) {
    int reversed = 0;
    for (int i = 0; i < bits; i++) {
        reversed <<= 1;
        reversed |= (n & 1);
        n >>= 1;
    }
    return reversed;
}

struct FFT {
    int size = 0; // must be a power of 2
    int bits = 0; // log2(size)
    int *rev = nullptr; // bit-reversed index
    float *real = nullptr;
    float *imag = nullptr;
    float *twiddle_real = nullptr; // precomputed twiddle factors
    float *twiddle_imag = nullptr;

    FFT(int size) : size(size) {
        if ((size & (size - 1)) != 0) {
            std::cout << "Error: size must be a power of 2\n";
            return;
        }

        real = new float[size];
        imag = new float[size];

        twiddle_real = new float[size / 2 + 1];
        twiddle_imag = new float[size / 2 + 1];

        rev = new int[size];
        bits = log2int(size);

        for (int i = 0; i < size; i++)
            rev[i] = reverseBits(i, bits);

        for (int i = 0; i <= size / 2; i++) {
            float theta = -2.0f * M_PIf * i / size;
            twiddle_real[i] = cosf(theta);
            twiddle_imag[i] = sinf(theta);
        }
    }

    ~FFT() {
        delete[] real;
        delete[] imag;
        delete[] twiddle_real;
        delete[] twiddle_imag;
        delete[] rev;
    }

    void set(int index, float real, float imag) {
        this->real[index] = real;
        this->imag[index] = imag;
    }

    void process(bool inverse = false) {
        for (int i = 0; i < size; i++) {
            if (i < rev[i]) {
                std::swap(real[i], real[rev[i]]);
                std::swap(imag[i], imag[rev[i]]);
            }
        }

        for (int s = 1; s <= bits; ++s) {
            int m = 1 << s;
            int m2 = m >> 1;
            float w_r = 1.0f;
            float w_i = 0.0f;
            float w_m_r = twiddle_real[size / m];
            float w_m_i = inverse ? -twiddle_imag[size / m] : twiddle_imag[size / m];

            for (int j = 0; j < m2; ++j) {
                for (int k = j; k < size; k += m) {
                    int t = k + m2;
                    float u_r = real[k];
                    float u_i = imag[k];
                    float v_r = real[t] * w_r - imag[t] * w_i;
                    float v_i = real[t] * w_i + imag[t] * w_r;
                    real[k] = u_r + v_r;
                    imag[k] = u_i + v_i;
                    real[t] = u_r - v_r;
                    imag[t] = u_i - v_i;
                }
                float temp_w_r = w_r * w_m_r - w_i * w_m_i;
                w_i = w_r * w_m_i + w_i * w_m_r;
                w_r = temp_w_r;
            }
        }

        if (inverse) {
            float inv_size = 1.0f / size;
            for (int i = 0; i < size; i++) {
                real[i] *= inv_size;
                imag[i] *= inv_size;
            }
        }
    }
};

struct whiteNoise {
    static constexpr float _scale = 2.0f / 0xffffffff;
    
    unsigned int g_x1 = 0x67452301;
    unsigned int g_x2 = 0xefcdab89;

    float operator()() {
        g_x1 ^= g_x2;
        float res = g_x2 * _scale - 1.0f;
        g_x2 += g_x1;
        return res;
    }
};

int main(int argc, char *argv[]) {
    int window_size = 512;
    int window_overlap = 16;
    int sample_rate = 44100;
    float min_freq = 0.0f;
    float max_freq = (float) sample_rate * 0.5f; // Nyquist
    double duration = 3.0;
    float gamma = 2.2f;
    char *input_path = "\0";
    char *output_path = "\0";

    const char* help_text = "Usage: \n--min-freq (Minimum Frequency: Hz)\n--max-freq (Maximum Frequency: Hz)\n--duration (Duration: seconds)\n--sample_rate (Sample Rate: Hz)\n--gamma (Gamma: float)\n--window_size (Window Size: int (power of 2))\n--window_overlap (Window Overlap: int (power of 2))\n--input (Input Image: '.png', '.jpg', '.bmp')\n--output (Output Wave File: '.wav')\n";

    for (int arg = 1; arg < argc; arg++) {
        if (!strcmp(argv[arg], "--min-freq")) {
            min_freq = (float) atof(argv[++arg]);
        } else if (!strcmp(argv[arg], "--max-freq")) {
            max_freq = (float) atof(argv[++arg]);
        } else if (!strcmp(argv[arg], "--duration")) {
            duration = atof(argv[++arg]);
        } else if (!strcmp(argv[arg], "--sample-rate")) {
            sample_rate = (int) atoi(argv[++arg]);
        } else if (!strcmp(argv[arg], "--gamma")) {
            gamma = (float) atof(argv[++arg]);
        } else if (!strcmp(argv[arg], "--input")) {
            input_path = argv[++arg];
        } else if (!strcmp(argv[arg], "--output")) {
            output_path = argv[++arg];
        } else if (!strcmp(argv[arg], "--window-size")) {
            window_size = forceExp2int((int) atoi(argv[++arg]));
        } else if (!strcmp(argv[arg], "--window-overlap")) {
            window_overlap = forceExp2int((int) atoi(argv[++arg]));
        } else if (!strcmp(argv[arg], "--help")) {
            std::cout << help_text << std::endl;
            return 0;
        }
    }

    // check if input or output paths are empty
    if (input_path[0] == '\0') {
        std::cout << "Missing input path\n" << help_text << std::endl;
        return 1;
    }
    if (output_path[0] == '\0') {
        std::cout << "Missing output path\n" << help_text << std::endl;
        return 1;
    }

    // check if min_freq and max_freq are valid

    // flip if min_freq > max_freq
    if (min_freq > max_freq) {
        float temp = min_freq;
        min_freq = max_freq;
        max_freq = temp;
        std::cout << "Swapping min_freq and max_freq\n";
    }

    // throw error if min_freq and max_freq are equal (division by zero)
    if (std::abs(min_freq - max_freq) < 1e-6) {
        std::cout << "min_freq and max_freq cannot be equal\n";
        return 1;
    }

    int num_samples = static_cast<int>(duration * sample_rate);
    int hop_size = fastDivPow2(window_size, window_overlap);

    // load image
    stbi_set_flip_vertically_on_load(true);
    int width, height, channels;
    unsigned char *image = stbi_load(input_path, &width, &height, &channels, 0);
    if (!image) {
        std::cout << "Failed to load image\n";
        return 1;
    }
    
    std::cout << "Image size: " << width << "x" << height << "x" << channels << "\n";
    std::cout << "Window size: " << window_size << "\n";
    std::cout << "Window overlap: " << window_overlap << "\n";
    std::cout << "Number of samples: " << num_samples << "\n";
    std::cout << "Hop size: " << hop_size << "\n";
    std::cout << "Min freq: " << min_freq << "\n";
    std::cout << "Max freq: " << max_freq << "\n";
    std::cout << "Duration: " << duration << "\n";
    std::cout << "Sample rate: " << sample_rate << "\n";
    std::cout << "Gamma: " << gamma << "\n";

    std::cout << "\nProcessing...\n";

    float inv_scale = 1.0f / (float) (window_overlap >> 1);

    // fetch luminance of a pixel from RGB with float coordinates
    // take into account alpha channel if present (multiply in such case)
    // use bilinear interpolation
    // orientation is: +x = right, +y = up
    // range is: [0, 1] for each axis
    auto sample_image = [&](float x, float y) -> float {
        // clamp coordinates
        x = std::max(std::min(x, 1.0f), 0.0f);
        y = std::max(std::min(y, 1.0f), 0.0f);

        // map to pixel coordinates
        x *= width - 1;
        y *= height - 1;

        // integer pixels
        int i0 = (int) std::floor(x);
        int j0 = (int) std::floor(y);
        int i1 = i0 + 1;
        int j1 = j0 + 1;

        // clamp to valid coordinates
        i0 = std::max(std::min(i0, width - 1), 0);
        j0 = std::max(std::min(j0, height - 1), 0);
        i1 = std::max(std::min(i1, width - 1), 0);
        j1 = std::max(std::min(j1, height - 1), 0);

        // fractional parts
        float dx = x - i0;
        float dy = y - j0;

        // helper function to get luminance
        auto get_luminance = [&](int ix, int iy) -> float {
            int idx = (iy * width + ix) * channels;
            float r = image[idx + 0] / 255.0f;
            float g = image[idx + 1] / 255.0f;
            float b = image[idx + 2] / 255.0f;
            float lum = 0.2126f * r + 0.7152f * g + 0.0722f * b;
            if (channels == 4)
                lum *= image[idx + 3] / 255.0f;
            return lum;
        };

        // sample pixels
        float v00 = get_luminance(i0, j0);
        float v01 = get_luminance(i0, j1);
        float v10 = get_luminance(i1, j0);
        float v11 = get_luminance(i1, j1);

        // interpolate
        float v0 = v00 * (1.0f - dx) + v10 * dx;
        float v1 = v01 * (1.0f - dx) + v11 * dx;
        float v = v0 * (1.0f - dy) + v1 * dy;

        return pow(v, gamma);
    };

    // precompute hann window
    float *hann = new float[window_size];
    for (int i = 0; i < window_size; i++) {
        hann[i] = 0.5f * (1.0f - cosf(2.0f * M_PIf * i / (float) (window_size - 1)));
    }

    // bins shaping (from 0 to nyquist)
    int synth_size = window_size / 2 + 1;
    float *synth = new float[synth_size];

    // streamed buffer
    RingBuffer buff_in(window_size << 1);
    // accumulated buffer
    RingBuffer buff_acc(window_size << 1);
    
    FFT fft(window_size);
    whiteNoise noise;

    auto signal = noise;

    // fill the buffer with first window of noise
    for (int i = 0; i < window_size; i++) {
        buff_in.push(signal());
    }

    float *output_arr = new float[hop_size];

    // open output file
    /*std::ofstream output_file(output_path, std::ios::binary);
    if (!output_file) {
        printf("Error: could not open output file\n");
        return 1;
    }*/
    drwav wav;
    drwav_data_format format;
    format.container = drwav_container_riff; // standard WAV container
    format.format = DR_WAVE_FORMAT_IEEE_FLOAT; // 32-bit float
    format.channels = 1;
    format.sampleRate = sample_rate;
    format.bitsPerSample = 32;
    if (!drwav_init_file_write(&wav, output_path, &format, nullptr)) {
        printf("Error: could not open output file\n");
        return 1;
    }


    // process each window of samples
    for (int cursor = 0; cursor < num_samples; cursor += hop_size) {
        // read the next window of samples, apply the hann window beforehand
        for (int i = 0; i < window_size; i++) {
            fft.set(i, buff_in.peek(i) * hann[i], 0.0f);
        }

        // apply FFT
        fft.process();
        
        // set bin shaping curve
        for (int i = 0; i < synth_size; i++) {
            // x -> cursor / total samples
            float x = (float) cursor / (float) num_samples;
            // y -> (freq - min_freq) / (max_freq - min_freq)
            float freq = i * sample_rate / (float) window_size;
            float y = (freq - min_freq) / (max_freq - min_freq);
            
            // sample the shaping curve
            synth[i] = sample_image(x, y);

        }

        // scale the first half of the spectrum
        for (int i = 0; i < synth_size; i++) {
            fft.real[i] *= synth[i];
            fft.imag[i] *= synth[i];
        }

        // scale the redundant second half of the spectrum (conjugate symmetric)
        for (int i = synth_size; i < window_size; i++) {
            fft.real[i] *= synth[window_size - i];
            fft.imag[i] *= synth[window_size - i];
        }
        
        // apply inverse FFT
        fft.process(true);

        // accumulate the output of each window
        for (int i = 0; i < window_size; i++) {
            buff_acc.add(i, fft.real[i]);
        }

        // write the current hop of samples to the output file
        int samples_to_write = std::min(hop_size, num_samples - cursor);
        for (int i = 0; i < samples_to_write; i++)
            output_arr[i] = buff_acc.peek(i) * inv_scale;
        drwav_uint64 samples_written = drwav_write_pcm_frames(&wav, samples_to_write, output_arr);

        // clear oldest contributions
        for (int i = 0; i < hop_size; i++) {
            buff_acc.set(i, 0.0f);
        }

        // move the cursor for the next hop
        buff_acc.jump_read(hop_size);
        buff_acc.jump_write(hop_size);

        // move the cursor for the next hop
        buff_in.jump_read(hop_size); 

        // generate the next hop of noise samples and push them into the input buffer
        for (int i = 0; i < hop_size; i++) {
            buff_in.push(signal());
        }
    }

    // close the output file
    // output_file.close();
    drwav_uninit(&wav);

    delete[] hann;
    delete[] synth;

    std::cout << "Done\n";
    
    return 0;
}
