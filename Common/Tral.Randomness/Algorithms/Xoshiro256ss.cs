//--------------------------------------------------------------------------------------------------
// PROJECT      : TRAL
// COPYRIGHT    : Andy Thomas / https://kuiper.zone
// LICENSE      : Apache License, Version 2.0
//--------------------------------------------------------------------------------------------------

using System.Runtime.CompilerServices;
using Tral.Randomness.Utility;

namespace Tral.Randomness.Algorithms;

/// <summary>
/// Xoshiro256** is a fast 64-bit all purposes (non-cryptographic) pseudo random number
/// generator with 256 bits of internal state. It supports fast jumping with <see cref="Jump"/>
/// equivalent to 2^192 calls of <see cref="Next"/>; it can be used to generate 2^64 starting
/// points, from each of which each jump will generate 2^64 non-overlapping subsequences for
/// parallel distributed computations. On construction, the initial state is populated with a
/// fixed arbitrary value.
/// </summary>
/// <remarks>
/// The xoshiro256** is 4-dimensionally equidistributed, whereas xoshiro256++ is 3-dimensionally
/// equidistributed. However, xoshiro256++ has higher linear complexity. This implementation is
/// derived from code put in the public domain by Sebastiano Vigna. Refer
/// to: http://xoroshiro.di.unimi.it
/// </remarks>
public class Xoshiro256ss : IJumpableAlgorithm
{
    // Long jump variant
    private static readonly ulong[] JumpConstants = { 0x76E15D3EFEFDCBBFUL, 0xC5004E441C522FB3UL, 0x77710069854EE241UL, 0x39109BB02ACBE635UL };

    private ulong _s0;
    private ulong _s1;
    private ulong _s2;
    private ulong _s3;

    /// <summary>
    /// Default constructor.
    /// </summary>
    public Xoshiro256ss()
    {
        Initialize();
    }

    private Xoshiro256ss(Xoshiro256ss other)
    {
        _s0 = other._s0;
        _s1 = other._s1;
        _s2 = other._s2;
        _s3 = other._s3;
    }

    /// <summary>
    /// Implements <see cref="IRandomAlgorithm.AlgorithmName"/>.
    /// </summary>
    public string AlgorithmName => "xoshiro256** 1.0";

    /// <summary>
    /// Implements <see cref="IRandomAlgorithm.MaxNext"/>.
    /// </summary>
    public ulong MaxNext => ulong.MaxValue;

    /// <summary>
    /// Implements <see cref="ISeedableAlgorithm.SeedLength"/>.
    /// </summary>
    public int SeedLength => 32;

    /// <summary>
    /// Implements <see cref="IRandomAlgorithm.Next"/>.
    /// </summary>
    public ulong Next()
    {
        ulong rslt = Rotl(_s1 * 5, 7) * 9;
        ulong t = _s1 << 17;

        _s2 ^= _s0;
        _s3 ^= _s1;
        _s1 ^= _s2;
        _s0 ^= _s3;

        _s2 ^= t;
        _s3 = Rotl(_s3, 45);

        return rslt;
    }

    /// <summary>
    /// Implements <see cref="ISeedableAlgorithm.SetSeed(byte[])"/>.
    /// </summary>
    public void SetSeed(byte[] seed)
    {
        SetSeed(SeedHelper.ToUInt64Array(seed, SeedLength));
    }

    /// <summary>
    /// Sets the internal state from an array of integers. The array must be of length 4 or greater.
    /// </summary>
    public void SetSeed(ulong[] seed)
    {
        _s0 = seed[0];
        _s1 = seed[1];
        _s2 = seed[2];
        _s3 = seed[3];

        switch (_s0)
        {
            case 0 when _s1 == 0 && _s2 == 0 && _s3 == 0:
                // Cannot all be 0
                Initialize();
                break;
        }
    }

    /// <summary>
    /// Implements <see cref="ISeedableAlgorithm.Clone"/>.
    /// </summary>
    public ISeedableAlgorithm Clone()
    {
        return new Xoshiro256ss(this);
    }

    /// <summary>
    /// Implements <see cref="IJumpableAlgorithm.Jump"/>.
    /// </summary>
    public void Jump()
    {
        ulong s0 = 0;
        ulong s1 = 0;
        ulong s2 = 0;
        ulong s3 = 0;

        foreach (ulong t in JumpConstants)
        {
            for (int b = 0; b < 64; ++b)
            {
                if ((t & 0x01UL << b) != 0)
                {
                    s0 ^= _s0;
                    s1 ^= _s1;
                    s2 ^= _s2;
                    s3 ^= _s3;
                }

                Next();
            }
        }

        _s0 = s0;
        _s1 = s1;
        _s2 = s2;
        _s3 = s3;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static ulong Rotl(ulong x, int k)
    {
        return (x << k) | (x >> (64 - k));
    }

    private void Initialize()
    {
        // Internal state can be any value except 0 everywhere.
        // We choose an arbitrary starting state because there is no published example.
        _s0 = 0x5A17EB1D35400E31UL;
        _s1 = 0x51F6183DF8203CDBUL;
        _s2 = 0x40E396FFFA72E9A8UL;
        _s3 = 0x7D6942DB29A87560UL;
    }
}