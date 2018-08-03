/*
    Copyright (c) 2018 Christopher A. Taylor.  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.
    * Neither the name of Solinas64 nor the names of its contributors may be
      used to endorse or promote products derived from this software without
      specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
    LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    POSSIBILITY OF SUCH DAMAGE.
*/

#include "../solinas64.h"
#include "gf256.h"

#define SOLINAS64_ENABLE_GF256_COMPARE

/**
    Solinas64 Benchmarks

    The goal of the benchmarks is to determine how fast Solinas prime field
    arithmetic is for the purpose of implementing erasure codes in software.
*/

#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
using namespace std;


#ifdef _WIN32
    #ifndef NOMINMAX
        #define NOMINMAX
    #endif
    #include <windows.h>
#elif __MACH__
    #include <mach/mach_time.h>
    #include <mach/mach.h>
    #include <mach/clock.h>

    extern mach_port_t clock_port;
#else
    #include <time.h>
    #include <sys/time.h>
#endif


//------------------------------------------------------------------------------
// Timing

#ifdef _WIN32
// Precomputed frequency inverse
static double PerfFrequencyInverseUsec = 0.;
static double PerfFrequencyInverseMsec = 0.;

static void InitPerfFrequencyInverse()
{
    LARGE_INTEGER freq = {};
    if (!::QueryPerformanceFrequency(&freq) || freq.QuadPart == 0)
        return;
    const double invFreq = 1. / (double)freq.QuadPart;
    PerfFrequencyInverseUsec = 1000000. * invFreq;
    PerfFrequencyInverseMsec = 1000. * invFreq;
}
#elif __MACH__
static bool m_clock_serv_init = false;
static clock_serv_t m_clock_serv = 0;

static void InitClockServ()
{
    m_clock_serv_init = true;
    host_get_clock_service(mach_host_self(), SYSTEM_CLOCK, &m_clock_serv);
}
#endif // _WIN32

uint64_t GetTimeUsec()
{
#ifdef _WIN32
    LARGE_INTEGER timeStamp = {};
    if (!::QueryPerformanceCounter(&timeStamp))
        return 0;
    if (PerfFrequencyInverseUsec == 0.)
        InitPerfFrequencyInverse();
    return (uint64_t)(PerfFrequencyInverseUsec * timeStamp.QuadPart);
#elif __MACH__
    if (!m_clock_serv_init)
        InitClockServ();

    mach_timespec_t tv;
    clock_get_time(m_clock_serv, &tv);

    return 1000000 * tv.tv_sec + tv.tv_nsec / 1000;
#else
    struct timeval tv;
    gettimeofday(&tv, nullptr);
    return 1000000 * tv.tv_sec + tv.tv_usec;
#endif
}


//------------------------------------------------------------------------------
// Fp61 Erasure Code Encoder

/**
    Encode()

    This function implements the encoder for an erasure code.
    It accepts a set of equal-sized data packets and outputs one recovery packet
    that can repair one lost original packet.

    The recovery packet must be GetRecoveryBytes() in size.

    Returns the number of bytes written.
*/
unsigned Encode(
    const std::vector<std::vector<uint8_t>>& originals,
    unsigned N,
    unsigned bytes,
    uint64_t seed,
    uint8_t* workspace,
    uint8_t* recovery,
    unsigned maxRecoveryBytes)
{
    // Set up row seed
    uint64_t seedMix = solinas64::HashU64(seed);

    // Unroll first column
    uint64_t coeff0 = solinas64::HashToNonzeroFp(seedMix + 0);
    unsigned recoveryBytes = solinas64::MultiplyRegion(
        &originals[0][0],
        bytes,
        coeff0,
        workspace,
        recovery);

    // Pad with zeros in case others overflow more
    memset(recovery + recoveryBytes, 0, maxRecoveryBytes - recoveryBytes);

    // For each remaining column:
    for (unsigned i = 1; i < N; ++i)
    {
        uint64_t coeff_i = solinas64::HashToNonzeroFp(seedMix + i);

        unsigned written = solinas64::MultiplyAddRegion(
            &originals[i][0],
            bytes,
            coeff_i,
            workspace,
            recovery);

        if (recoveryBytes < written) {
            recoveryBytes = written;
        }
    }

    return recoveryBytes;
}

void EncodeGF256(
    const std::vector<std::vector<uint8_t>>& originals,
    unsigned N,
    unsigned bytes,
    uint64_t seed,
    uint8_t* recovery)
{
    uint64_t seedMix = solinas64::HashU64(seed);

    uint8_t coeff = (uint8_t)solinas64::HashToNonzeroFp(seedMix + 0);
    if (coeff == 0) {
        coeff = 1;
    }

    gf256_mul_mem(recovery, &originals[0][0], coeff, bytes);

    for (unsigned i = 1; i < N; ++i)
    {
        coeff = (uint8_t)solinas64::HashToNonzeroFp(seedMix + 0);
        if (coeff == 0) {
            coeff = 1;
        }

        gf256_muladd_mem(recovery, coeff, &originals[i][0], bytes);
    }
}


//------------------------------------------------------------------------------
// Benchmarks

static const unsigned kFileSizes[] = {
    10, 100, 1000, 10000, 100000
};
static const unsigned kFileSizesCount = static_cast<unsigned>(sizeof(kFileSizes) / sizeof(kFileSizes[0]));

static const unsigned kFileN[] = {
    2, 4, 8, 16, 32, 64, 128, 256, 512
};
static const unsigned kFileNCount = static_cast<unsigned>(sizeof(kFileN) / sizeof(kFileN[0]));

static const unsigned kTrials = 1000;

void RunBenchmarks()
{
    solinas64::Random prng;
    prng.Seed(0);

    std::vector<std::vector<uint8_t>> original_data;
    std::vector<uint8_t> recovery_data;
    std::vector<uint8_t> workspace_data;

    for (unsigned i = 0; i < kFileSizesCount; ++i)
    {
        unsigned fileSizeBytes = kFileSizes[i];

        cout << "Testing file size = " << fileSizeBytes << " bytes" << endl;

        for (unsigned j = 0; j < kFileNCount; ++j)
        {
            unsigned N = kFileN[j];

            cout << "N = " << N << " : ";

            uint64_t sizeSum = 0, timeSum = 0;
            uint64_t timeSum_gf256 = 0;

            for (unsigned k = 0; k < kTrials; ++k)
            {
                /*
                    File pieces: f0, f1, f3, f4, ...
                    Coefficients: m0, m1, m2, m3, ...

                    R = m0 * f0 + m1 * f1 + m2 * f2 + ...

                    R = sum(m_i * f_i) (mod 2^61-1)

                    To compute the recovery packet R we process the calculations
                    for the first word from all of the file pieces to produce a
                    single word of output.  This is a matrix-vector product
                    between file data f_i (treated as Fp words) and randomly
                    chosen generator matrix coefficients m_i.

                    Lazy reduction can be used to simplify the add steps.

                    Then we continue to the next word for all the file pieces,
                    producing the next word of output.

                    It is possible to interleave the calculations for output
                    words, and for input words to achieve higher throughput.

                    The number of words for each file piece can vary slightly
                    based on the data (if the data bytes do not fit evenly into
                    the Fp words, we have to add extra bits to resolve
                    ambiguities).

                    The result is a set of 61-bit Fp words serialized to bytes,
                    that is about 8 bytes more than the original file sizes.

                    The erasure code decoder (not implemented) would be able
                    to take these recovery packets and fix lost data.
                    The decoder performance would be fairly similar to the
                    encoder performance for this type of erasure code, since
                    the runtime is dominated by this matrix-vector product.
                */

                original_data.resize(N);
                for (unsigned s = 0; s < N; ++s)
                {
                    // Add 8 bytes padding to simplify tester
                    original_data[s].resize(fileSizeBytes + 8);

                    // Fill the data with random bytes
                    for (unsigned r = 0; r < i; r += 8)
                    {
                        uint64_t w;
                        if (prng.Next() % 100 <= 3) {
                            w = ~(uint64_t)0;
                        }
                        else {
                            w = prng.Next();
                        }
                        solinas64::WriteU64_LE(&original_data[s][r], w);
                    }
                }

                const unsigned maxRecoveryBytes = solinas64::AppDataReader::GetMaxOutputBytes(fileSizeBytes);
                const unsigned workspaceBytes = solinas64::AppDataReader::GetWorkspaceBytes(fileSizeBytes);
                recovery_data.resize(maxRecoveryBytes);
                workspace_data.resize(workspaceBytes);

                {
                    uint64_t t0 = GetTimeUsec();

                    unsigned recoveryBytes = Encode(
                        original_data,
                        N,
                        fileSizeBytes,
                        k,
                        &workspace_data[0],
                        &recovery_data[0],
                        maxRecoveryBytes);

                    uint64_t t1 = GetTimeUsec();

                    sizeSum += recoveryBytes;
                    timeSum += t1 - t0;
                }

#ifdef SOLINAS64_ENABLE_GF256_COMPARE
                {
                    uint64_t t0 = GetTimeUsec();

                    EncodeGF256(original_data, N, fileSizeBytes, k, &recovery_data[0]);

                    uint64_t t1 = GetTimeUsec();

                    timeSum_gf256 += t1 - t0;
                }
#endif // SOLINAS64_ENABLE_GF256_COMPARE
            }

#ifdef SOLINAS64_ENABLE_GF256_COMPARE
            cout << " gf256_MBPS=" << (uint64_t)fileSizeBytes * N * kTrials / timeSum_gf256;
#endif // SOLINAS64_ENABLE_GF256_COMPARE
            cout << " Solinas64_MBPS=" << (uint64_t)fileSizeBytes * N * kTrials / timeSum;
            cout << " Solinas64_OutputBytes=" << sizeSum / (float)kTrials;
            cout << endl;
        }
    }
}


//------------------------------------------------------------------------------
// Entrypoint

int main()
{
    cout << "Benchmarks for Fp61 erasure codes.  Before running the benchmarks please run the tests to make sure everything's working on your PC.  It's going to run quite a bit faster with 64-bit builds because it takes advantage of the speed of 64-bit multiplications." << endl;
    cout << endl;

    gf256_init();

    RunBenchmarks();

    cout << endl;
    return 0;
}
