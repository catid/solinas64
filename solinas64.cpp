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

#include "solinas64.h"

#include <string.h>

namespace solinas64 {


// This is an unrolled implementation of Knuth's unsigned version of the eGCD,
// specialized for the prime.  It handles any input.
uint64_t Inverse(uint64_t u)
{
    uint64_t u1, u3, v1, v3, qt;

    qt = u / kPrime;
    u3 = u % kPrime;
    u1 = 1;

    if (u3 == 0) {
        return 0; // No inverse
    }

    qt = kPrime / u3;
    v3 = kPrime % u3;
    v1 = qt;

    for (;;)
    {
        if (v3 == 0) {
            return u3 == 1 ? u1 : 0;
        }

        qt = u3 / v3;
        u3 %= v3;
        u1 += qt * v1;

        if (u3 == 0) {
            return v3 == 1 ? kPrime - v1 : 0;
        }

        qt = v3 / u3;
        v3 %= u3;
        v1 += qt * u1;
    }
}


//------------------------------------------------------------------------------
// Memory Reading

uint64_t ReadBytes_LE(const uint8_t* data, unsigned bytes)
{
    switch (bytes)
    {
    case 8: return ReadU64_LE(data);
    case 7: return ((uint64_t)data[6] << 48) | ((uint64_t)data[5] << 40) | ((uint64_t)data[4] << 32) | ReadU32_LE(data);
    case 6: return ((uint64_t)data[5] << 40) | ((uint64_t)data[4] << 32) | ReadU32_LE(data);
    case 5: return ((uint64_t)data[4] << 32) | ReadU32_LE(data);
    case 4: return ReadU32_LE(data);
    case 3: return ((uint32_t)data[2] << 16) | ((uint32_t)data[1] << 8) | data[0];
    case 2: return ((uint32_t)data[1] << 8) | data[0];
    case 1: return data[0];
    default: break;
    }
    return 0;
}


//------------------------------------------------------------------------------
// Memory Writing

void WriteBytes_LE(uint8_t* data, unsigned bytes, uint64_t value)
{
    switch (bytes)
    {
    case 8: WriteU64_LE(data, value);
        return;
    case 7: data[6] = (uint8_t)(value >> 48);
    case 6: data[5] = (uint8_t)(value >> 40);
    case 5: data[4] = (uint8_t)(value >> 32);
    case 4: WriteU32_LE(data, static_cast<uint32_t>(value));
        return;
    case 3: data[2] = (uint8_t)(value >> 16);
    case 2: data[1] = (uint8_t)(value >> 8);
    case 1: data[0] = (uint8_t)value;
    default: break;
    }
}


//------------------------------------------------------------------------------
// Random

// From http://xoshiro.di.unimi.it/splitmix64.c
// Written in 2015 by Sebastiano Vigna (vigna@acm.org)
uint64_t HashU64(uint64_t x)
{
    x += 0x9e3779b97f4a7c15;
    uint64_t z = x;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
    z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
    return z ^ (z >> 31);
}

void Random::Seed(uint64_t x)
{
    // Fill initial state as recommended by authors
    uint64_t h = HashU64(x);
    State[0] = h;
    h = HashU64(h);
    State[1] = h;
    h = HashU64(h);
    State[2] = h;
    h = HashU64(h);
    State[3] = h;
}


//------------------------------------------------------------------------------
// Bulk Operations

unsigned MultiplyRegion(
    const uint8_t* data,
    unsigned bytes,
    uint64_t coeff,
    uint8_t* workspace,
    uint8_t* output)
{
    const unsigned minimumOutputBytes = (bytes + 7) & ~7u;

    // Special fast cases
    if (coeff <= 1)
    {
        if (coeff == 0) {
            memset(output, 0, minimumOutputBytes);
        }
        else
        {
            memcpy(output, data, bytes);
            memset(output + bytes, 0, minimumOutputBytes - bytes);
        }
        return minimumOutputBytes;
    }

    AppDataReader reader;
    reader.SetupWorkspace(workspace);

    while (bytes >= 32)
    {
        bytes -= 32;

        uint64_t x0 = Multiply(coeff, reader.ReadNext8Bytes(data));
        uint64_t x1 = Multiply(coeff, reader.ReadNext8Bytes(data + 8));
        uint64_t x2 = Multiply(coeff, reader.ReadNext8Bytes(data + 16));
        uint64_t x3 = Multiply(coeff, reader.ReadNext8Bytes(data + 24));

        data += 32;

        WriteU64_LE(output, x0);
        WriteU64_LE(output + 8, x1);
        WriteU64_LE(output + 16, x2);
        WriteU64_LE(output + 24, x3);

        output += 32;
    }

    while (bytes >= 8)
    {
        bytes -= 8;

        uint64_t x0 = Multiply(coeff, reader.ReadNext8Bytes(data));
        data += 8;

        WriteU64_LE(output, x0);
        output += 8;
    }

    if (bytes > 0)
    {
        uint64_t x0 = Multiply(coeff, reader.ReadFinalBytes(data, bytes));
        WriteU64_LE(output, x0);
        output += 8;
    }

    // Finalize the overflow bits
    const unsigned extraWordBytes = reader.FlushAndGetWordCount() * 8;
    const uint8_t* readPtr = reader.Data;

    // Also work on the overflow bits
    for (unsigned i = 0; i < extraWordBytes; i += 8)
    {
        WriteU64_LE(
            output + i,
            Multiply(
                coeff,
                ReadU64_LE(readPtr + i)));
    }

    return minimumOutputBytes + extraWordBytes;
}

unsigned MultiplyAddRegion(
    const uint8_t* data,
    unsigned bytes,
    uint64_t coeff,
    uint8_t* workspace,
    uint8_t* output)
{
    const unsigned minimumOutputBytes = (bytes + 7) & ~7u;

    // Special fast case
    if (coeff == 0)
    {
        // TODO: Add a special case for coeff = 1
        return minimumOutputBytes;
    }

    AppDataReader reader;
    reader.SetupWorkspace(workspace);

    /**** This loop takes over 95% of the execution time. ****/
    while (bytes >= 32)
    {
        bytes -= 32;

        WriteU64_LE(output, Add(Multiply(coeff, reader.ReadNext8Bytes(data)), ReadU64_LE(output)));
        WriteU64_LE(output + 8, Add(Multiply(coeff, reader.ReadNext8Bytes(data + 8)), ReadU64_LE(output + 8)));
        WriteU64_LE(output + 16, Add(Multiply(coeff, reader.ReadNext8Bytes(data + 16)), ReadU64_LE(output + 16)));
        WriteU64_LE(output + 24, Add(Multiply(coeff, reader.ReadNext8Bytes(data + 24)), ReadU64_LE(output + 24)));

        data += 32;
        output += 32;
    }

    while (bytes >= 8)
    {
        bytes -= 8;

        uint64_t x0 = Add(Multiply(coeff, reader.ReadNext8Bytes(data)), ReadU64_LE(output));
        data += 8;

        WriteU64_LE(output, x0);
        output += 8;
    }

    if (bytes > 0)
    {
        uint64_t x0 = Add(Multiply(coeff, reader.ReadFinalBytes(data, bytes)), ReadU64_LE(output));

        WriteU64_LE(output, x0);
        output += 8;
    }

    // Finalize the overflow bits
    const unsigned extraWordBytes = reader.FlushAndGetWordCount() * 8;
    const uint8_t* readPtr = reader.Data;

    // Also work on the overflow bits
    for (unsigned i = 0; i < extraWordBytes; i += 8)
    {
        const uint64_t x = ReadU64_LE(output + i);
        WriteU64_LE(
            output + i,
            Add(
                Multiply(
                    coeff,
                    ReadU64_LE(readPtr + i)),
                x));
    }

    return minimumOutputBytes + extraWordBytes;
}


} // namespace solinas64
