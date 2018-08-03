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

#ifndef CAT_SOLINAS64_H
#define CAT_SOLINAS64_H

#include <stdint.h>

/** \mainpage
    Solinas64 : Finite field arithmetic modulo Solinas prime p = 2^64-2^32+1 in C++

    The Solinas64 software takes advantage of the commonly available fast
    64x64->128 multiplier to accelerate finite (base) field arithmetic.
    So it runs a lot faster when built into a 64-bit executable.
*/

// Define this to avoid any unaligned memory accesses while reading data.
// This is useful as a quick-fix for mobile applications.
// A preferred solution is to ensure that the data provided is aligned.
// Another reason to do this is if the platform is big-endian.
//#define SOLINAS64_SAFE_MEMORY_ACCESSES


//------------------------------------------------------------------------------
// Portability Macros

// Compiler-specific force inline keyword
#ifdef _MSC_VER
# define SOLINAS64_FORCE_INLINE inline __forceinline
#else
# define SOLINAS64_FORCE_INLINE inline __attribute__((always_inline))
#endif


//------------------------------------------------------------------------------
// Portable 64x64->128 Multiply
// CAT_MUL128: r{hi,lo} = x * y

// Returns low part of product, and high part is set in r_hi
SOLINAS64_FORCE_INLINE uint64_t Emulate64x64to128(
    uint64_t& r_hi,
    const uint64_t x,
    const uint64_t y)
{
    // Form temporary 32-bit words
    const uint32_t x0 = static_cast<uint32_t>(x);
    const uint32_t x1 = static_cast<uint32_t>(x >> 32);
    const uint32_t y0 = static_cast<uint32_t>(y);
    const uint32_t y1 = static_cast<uint32_t>(y >> 32);

    // Calculate 32x32->64 bit products
    const uint64_t p11 = static_cast<uint64_t>(x1) * y1;
    const uint64_t p01 = static_cast<uint64_t>(x0) * y1;
    const uint64_t p10 = static_cast<uint64_t>(x1) * y0;
    const uint64_t p00 = static_cast<uint64_t>(x0) * y0;

    /*
        This is implementing schoolbook multiplication:

                x1 x0
        X       y1 y0
        -------------
                   00  LOW PART
        -------------
                00
             10 10     MIDDLE PART
        +       01
        -------------
             01 
        + 11 11        HIGH PART
        -------------
    */

    // 64-bit product + two 32-bit values
    const uint64_t middle = p10
        + static_cast<uint32_t>(p00 >> 32)
        + static_cast<uint32_t>(p01);

    /*
        Proof that 64-bit products can accumulate two more 32-bit values
        without overflowing:

        Max 32-bit value is 2^32 - 1.
        PSum = (2^32-1) * (2^32-1) + (2^32-1) + (2^32-1)
             = 2^64 - 2^32 - 2^32 + 1 + 2^32 - 1 + 2^32 - 1
             = 2^64 - 1
        Therefore it cannot overflow regardless of input.
    */

    // 64-bit product + two 32-bit values
    r_hi = p11
        + static_cast<uint32_t>(middle >> 32)
        + static_cast<uint32_t>(p01 >> 32);

    // Add LOW PART and lower half of MIDDLE PART
    return (middle << 32) | static_cast<uint32_t>(p00);
}

#if defined(_MSC_VER) && defined(_WIN64)
// Visual Studio 64-bit

# include <intrin.h>
# pragma intrinsic(_umul128)
# define CAT_MUL128(r_hi, r_lo, x, y) \
    r_lo = _umul128(x, y, &(r_hi));

#elif defined(__SIZEOF_INT128__)
// Compiler supporting 128-bit values (GCC/Clang)

# define CAT_MUL128(r_hi, r_lo, x, y)                   \
    {                                                   \
        unsigned __int128 w = (unsigned __int128)x * y; \
        r_lo = (uint64_t)w;                             \
        r_hi = (uint64_t)(w >> 64);                     \
    }

#else
// Emulate 64x64->128-bit multiply with 64x64->64 operations

# define CAT_MUL128(r_hi, r_lo, x, y) \
    r_lo = Emulate64x64to128(r_hi, x, y);

#endif // End CAT_MUL128


namespace solinas64 {


//------------------------------------------------------------------------------
// Constants

// p = 2^64 - 2^32 + 1
static const int64_t kPrimeSubC = ((uint64_t)1 << 32) - 1;
static const uint64_t kPrime = 0 - kPrimeSubC;


//------------------------------------------------------------------------------
// API

SOLINAS64_FORCE_INLINE bool adc(uint64_t& x, uint64_t y)
{
    if ((x += y) >= y) {
        return false;
    }
    return true;
}

SOLINAS64_FORCE_INLINE bool sbb(uint64_t& x, uint64_t y)
{
    const uint64_t x0 = x;
    x = x0 - y;
    if (x0 >= y) {
        return false;
    }
    return true;
}

/**
    r = solinas64::Add(x, y)

    Returns sum x + y (mod p).
*/
SOLINAS64_FORCE_INLINE uint64_t Add(uint64_t x, uint64_t y)
{
    if (adc(x, y)) {
        if (adc(x, kPrimeSubC)) {
            adc(x, kPrimeSubC);
        }
    }
    return x;
}

/**
    r = solinas64::Subtract(x, y)

    Returns difference x - y (mod p).
*/
SOLINAS64_FORCE_INLINE uint64_t Subtract(uint64_t x, uint64_t y)
{
    if (sbb(x, y)) {
        if (sbb(x, kPrimeSubC)) {
            sbb(x, kPrimeSubC);
        }
    }
    return x;
}

/**
    r = solinas64::Multiply(x, y)

    r = x * y (mod p)
*/
SOLINAS64_FORCE_INLINE uint64_t Multiply(uint64_t x, uint64_t y)
{
    uint64_t p_lo, p_hi;
    CAT_MUL128(p_hi, p_lo, x, y);

    /*
      2^64 - 2^32 + 1

       3  2  1  0
            -3 -2
         +3 +2

             1  0
            +2 -2
               -3
    */

    uint32_t a2 = static_cast<uint32_t>(p_hi);
    uint32_t a3 = static_cast<uint32_t>(p_hi >> 32);

    uint64_t t = (static_cast<uint64_t>(a2) << 32) - a2;

    if (adc(p_lo, t)) {
        adc(p_lo, kPrimeSubC);
    }
    sbb(p_lo, a3);
    return p_lo;
}

/**
    r = solinas64::Inverse(x)

    r = x^-1 (mod p)
    The input value x can be any 64-bit value.

    This operation is kind of heavy so it should be avoided where possible.

    This operation is not constant-time.
    A constant-time version can be implemented using Euler's totient method and
    a straight line similar to
    https://github.com/catid/snowshoe/blob/master/src/fp.inc#L545

    Returns the multiplicative inverse of x modulo p.
    0 < result < p

    If the inverse does not exist, it returns 0.
*/
uint64_t Inverse(uint64_t x);


//------------------------------------------------------------------------------
// Memory Reading

/// Read 8 bytes in little-endian byte order
SOLINAS64_FORCE_INLINE uint64_t ReadU64_LE(const uint8_t* data)
{
#ifdef SOLINAS64_SAFE_MEMORY_ACCESSES
    return ((uint64_t)data[7] << 56) | ((uint64_t)data[6] << 48) | ((uint64_t)data[5] << 40) |
        ((uint64_t)data[4] << 32) | ((uint64_t)data[3] << 24) | ((uint64_t)data[2] << 16) |
        ((uint64_t)data[1] << 8) | data[0];
#else
    const uint64_t* wordPtr = reinterpret_cast<const uint64_t*>(data);
    return *wordPtr;
#endif
}

/// Read 4 bytes in little-endian byte order
SOLINAS64_FORCE_INLINE uint32_t ReadU32_LE(const uint8_t* data)
{
#ifdef SOLINAS64_SAFE_MEMORY_ACCESSES
    return ((uint32_t)data[3] << 24) | ((uint32_t)data[2] << 16) | ((uint32_t)data[1] << 8) | data[0];
#else
    const uint32_t* wordPtr = reinterpret_cast<const uint32_t*>(data);
    return *wordPtr;
#endif
}

/// Read between 0..8 bytes in little-endian byte order
/// Returns 0 for any other value for `bytes`
uint64_t ReadBytes_LE(const uint8_t* data, unsigned bytes);

/// Values with all these bits set are ambiguous and need an extra bit to disambiguate the high bit.
static const uint64_t kAmbiguityMask = 0x7fffffff00000000ULL;

/// Returns true if the word provided needs an extra bit to represent it.
SOLINAS64_FORCE_INLINE bool IsU64Ambiguous(uint64_t u64_word)
{
    return (u64_word & kAmbiguityMask) == kAmbiguityMask;
}

/// Mask to clear the high bit
static const uint64_t kHighBitMask = 0x7fffffffffffffffULL;


//------------------------------------------------------------------------------
// Memory Writing

/// Write 4 bytes in little-endian byte order
SOLINAS64_FORCE_INLINE void WriteU32_LE(uint8_t* data, uint32_t value)
{
#ifdef SOLINAS64_SAFE_MEMORY_ACCESSES
    data[3] = (uint8_t)(value >> 24);
    data[2] = (uint8_t)(value >> 16);
    data[1] = (uint8_t)(value >> 8);
    data[0] = (uint8_t)value;
#else
    uint32_t* wordPtr = reinterpret_cast<uint32_t*>(data);
    *wordPtr = value;
#endif
}

/// Write 8 bytes in little-endian byte order
SOLINAS64_FORCE_INLINE void WriteU64_LE(uint8_t* data, uint64_t value)
{
#ifdef SOLINAS64_SAFE_MEMORY_ACCESSES
    data[7] = (uint8_t)(value >> 56);
    data[6] = (uint8_t)(value >> 48);
    data[5] = (uint8_t)(value >> 40);
    data[4] = (uint8_t)(value >> 32);
    data[3] = (uint8_t)(value >> 24);
    data[2] = (uint8_t)(value >> 16);
    data[1] = (uint8_t)(value >> 8);
    data[0] = (uint8_t)value;
#else
    uint64_t* wordPtr = reinterpret_cast<uint64_t*>(data);
    *wordPtr = value;
#endif
}

/// Write between 0..8 bytes in little-endian byte order
void WriteBytes_LE(uint8_t* data, unsigned bytes, uint64_t value);


//------------------------------------------------------------------------------
// Reading Data

/**
    AppDataReader

    This is used to read byte data into 64-bit field words.
    The complication is that 64-bit words >= p are not in the field, so to pack
    them in efficiently ReadNext8Bytes() emits an extra bit when the top half of
    the word is all ones (except the high bit).  The extra bits must be
    processed as additional data once all the original data is processed.
*/
struct AppDataReader
{
    uint8_t* Data;
    uint8_t* DataWritePtr;
    uint64_t Workspace;
    int Available;


    /// Returns the number of extra temporary workspace bytes needed.
    static SOLINAS64_FORCE_INLINE unsigned GetWorkspaceBytes(unsigned bytes)
    {
        const unsigned inputBits = bytes * 8;

        // All words may be expanded by one bit, hence the (bits / 64) factor,
        // but only full words can be too large, so we round down.
        const unsigned maxExtraBits = inputBits / 64;

        // Round up to the number of words needed.
        const unsigned words = (maxExtraBits + 63) / 64;

        // Return the number of bytes needed in the Data pointer.
        return words * 8;
    }

    /// Returns the number of bytes overall that will be produced.
    /// This includes the original data converted to words plus the extra words
    /// that are generated by this class to handle input overflows.
    static SOLINAS64_FORCE_INLINE unsigned GetMaxOutputBytes(unsigned bytes)
    {
        const unsigned originalWords = (bytes + 7) / 8;
        return GetWorkspaceBytes(bytes) + originalWords * 8;
    }

    /// Provide the workspace data buffer of size GetWorkspaceBytes().
    SOLINAS64_FORCE_INLINE void SetupWorkspace(uint8_t* workspace)
    {
        Data = workspace;
        DataWritePtr = workspace;
        Workspace = 0;
        Available = 0;
    }

    /// Fits a word from file/packet data into a Fp element and emits an extra
    /// bit if needed.  This should be called from the first word of data until
    /// the last word of data.
    SOLINAS64_FORCE_INLINE uint64_t ReadNext8Bytes(const uint8_t* data)
    {
        // Read the word of data
        const uint64_t* wordPtr = reinterpret_cast<const uint64_t*>(data);
        uint64_t word = *wordPtr;

        // If word needs a bit emitted:
        if (IsU64Ambiguous(word))
        {
            // Emit words that have the high bit cleared so they are all in Fp:

            // If we ran out of space:
            if (Available >= 63)
            {
                WriteU64_LE(DataWritePtr, Workspace);

                DataWritePtr += 8;
                Workspace = word >> 63;
                Available = 1;
            }
            else
            {
                Workspace |= (word >> 63) << Available;
                ++Available;
            }

            // Clear high bit either way
            word &= kHighBitMask;
        }

        return word;
    }

    /// Read the final few bytes of data.
    /// Precondition: Bytes > 0.
    SOLINAS64_FORCE_INLINE uint64_t ReadFinalBytes(const uint8_t* data, unsigned bytes)
    {
        // No need for reduction because at least high 8 bits are 0.
        return ReadBytes_LE(data, bytes);
    }

    /// Flush any remaining bits.
    /// Returns the number of words starting at Data, up to DataWritePtr.
    /// The application can either loop based on the return value or choose
    /// to loop until Data exceeds DataWritePtr.
    /// Each word can be read via ReadU64_LE().
    SOLINAS64_FORCE_INLINE unsigned FlushAndGetWordCount()
    {
        // If there are some available bytes:
        if (Available != 0)
        {
            // Write final workspace
            WriteU64_LE(DataWritePtr, Workspace);
            DataWritePtr += 8;
        }

        // Calculate number of words written overall
        const uintptr_t writtenWords = static_cast<uintptr_t>(DataWritePtr - Data) / 8;
        return static_cast<unsigned>(writtenWords);
    }
};


//------------------------------------------------------------------------------
// Random Numbers

#define CAT_ROL64(x, bits) ( ((uint64_t)(x) << (bits)) | ((uint64_t)(x) >> (64 - (bits))) )

/**
    Random

    Xoshiro256+ based pseudo-random number generator (PRNG) for 64-bit output.

    Call Seed() to provide a 64-bit generator seed.
    Call Next() to produce a random 64-bit number.
*/
struct Random
{
    uint64_t State[4];


    /// Seed the generator.  (This is somewhat expensive.)
    void Seed(uint64_t x);

    /// Get the next 64-bit random number.
    /// The low 3 bits are slightly weak according to the authors.
    // From http://xoshiro.di.unimi.it/xoshiro256plus.c
    // Written in 2018 by David Blackman and Sebastiano Vigna (vigna@acm.org)
    SOLINAS64_FORCE_INLINE uint64_t Next()
    {
        uint64_t s0 = State[0], s1 = State[1], s2 = State[2], s3 = State[3];

        const uint64_t result = s0 + s3;

        const uint64_t t = s1 << 17;
        s2 ^= s0;
        s3 ^= s1;
        s1 ^= s2;
        s0 ^= s3;
        s2 ^= t;
        s3 = CAT_ROL64(s3, 45);

        State[0] = s0, State[1] = s1, State[2] = s2, State[3] = s3;

        return result;
    }
};

/// Hash a 64-bit value to another 64-bit value
uint64_t HashU64(uint64_t x);

/// Hash a seed into a value from 1..p-1
SOLINAS64_FORCE_INLINE uint64_t HashToNonzeroFp(uint64_t word)
{
    // Run a simple mixer based on HashU64()
    word += 0x9e3779b97f4a7c15;
    word = (word ^ (word >> 30)) * 0xbf58476d1ce4e5b9;

    // Take the top 61 bits
    word >>= 3;

    // Eliminate values = p
    word -= (word + 1) >> 61;

    // Eliminate values = 0
    word += (word - 1) >> 63;

    return word;
}


//------------------------------------------------------------------------------
// Bulk Operations

/**
    MultiplyRegion()

    output[] = data[] * coeff

    Preconditions:
        0 <= coeff < p.
        data != null, output != null, workspace != null, bytes > 0

    Note: This expands the input data up to solinas64::GetMaxOutputBytes() bytes.

    Returns the number of bytes written.
*/
unsigned MultiplyRegion(
    const uint8_t* data,    ///< Input data
    unsigned bytes,         ///< Number of input data bytes
    uint64_t coeff,         ///< Coefficient to multiply the data by
    uint8_t* workspace,     ///< Size calculated by solinas64::GetWorkspaceBytes()
    uint8_t* output);       ///< Size calculated by solinas64::GetMaxOutputBytes()

/**
    MultiplyAddRegion()

    output[] = output[] + data[] * coeff

    Preconditions:
        0 <= coeff < p.
        data != null, output != null, workspace != null, bytes > 0

    Note: This expands the input data up to solinas64::GetMaxOutputBytes() bytes.

    Returns the number of bytes written.
*/
unsigned MultiplyAddRegion(
    const uint8_t* data,    ///< Input data
    unsigned bytes,         ///< Number of input data bytes
    uint64_t coeff,         ///< Coefficient to multiply the data by
    uint8_t* workspace,     ///< Size calculated by solinas64::GetWorkspaceBytes()
    uint8_t* output);       ///< Size calculated by solinas64::GetMaxOutputBytes()


} // namespace solinas64


#endif // CAT_SOLINAS64_H
