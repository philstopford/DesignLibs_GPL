//--------------------------------------------------------------------------------------------------
// PROJECT      : TRAL
// COPYRIGHT    : Andy Thomas / https://kuiper.zone
// LICENSE      : Apache License, Version 2.0
//--------------------------------------------------------------------------------------------------

using System;
using Tral.Randomness.Utility;

namespace Tral.Randomness.Algorithms;

/// <summary>
/// MT19937-64 is a 64-bit variant of MT19937-32. On 64-bit architectures, this variant will
/// typically be faster than <see cref="MersenneTwister32"/>. It is not suited to cryptographic
/// purposes. On construction, the initial state is seeded with the integer value 5489 according
/// to the MT seeding routine.
/// </summary>
/// <remarks>
/// Mersenne Twister (MT) was originally developed by Makoto Matsumoto and Takuji Nishimura in 1996/1997.
/// </remarks>
public class MersenneTwister64 : ISeedableAlgorithm
{
    private const uint DefaultSeed = 5489U;
    private const int StateN = 312;

    private const int BitWidth = 64;
    private const int BitWidthM2 = 62;

    private const ulong InitSeedA = 19650218UL;
    private const ulong InitSeedB = 3935559000370003845UL;
    private const ulong InitSeedC = 2862933555777941757UL;

    private const ulong UMask = 0xFFFFFFFF80000000UL;
    private const ulong LMask = 0x000000007FFFFFFFUL;
    private const ulong Matrix = 0xB5026F5AA96619E9UL;
    private const int ShiftU = 29;
    private const int ShiftS = 17;
    private const ulong MaskB = 0x71D67FFFEDA60000UL;
    private const int ShiftT = 37;
    private const ulong MaskC = 0xFFF7EEE000000000UL;
    private const ulong MaskD = 0x5555555555555555U;
    private const int ShiftL = 43;
    private const ulong InitF = 0x5851F42D4C957F2DUL;
    private const int Offset = 156;

    private int _idx;
    private readonly ulong[] _mt = new ulong[StateN];

    /// <summary>
    /// Default constructor. The initial state is seeded with the integer value 5489 according
    /// to the MT seeding routine.
    /// </summary>
    public MersenneTwister64()
    {
        SetSeed(DefaultSeed);
    }

    private MersenneTwister64(MersenneTwister64 other)
    {
        _idx = other._idx;
        Buffer.BlockCopy(other._mt, 0, _mt, 0, SeedLength);
    }

    /// <summary>
    /// Implements <see cref="IRandomAlgorithm.AlgorithmName"/>.
    /// </summary>
    public string AlgorithmName => "MT19937-64";

    /// <summary>
    /// Implements <see cref="IRandomAlgorithm.MaxNext"/>.
    /// </summary>
    public ulong MaxNext => ulong.MaxValue;

    /// <summary>
    /// Implements <see cref="ISeedableAlgorithm.SeedLength"/>.
    /// </summary>
    public int SeedLength => 2496;

    /// <summary>
    /// Implements <see cref="IRandomAlgorithm.Next"/>.
    /// </summary>
    public ulong Next()
    {
        switch (_idx)
        {
            case StateN:
                Twist();
                _idx = 0;
                break;
        }

        ulong y = _mt[_idx++];

        y ^= (y >> ShiftU) & MaskD;
        y ^= (y << ShiftS) & MaskB;
        y ^= (y << ShiftT) & MaskC;
        return y ^ (y >> ShiftL);
    }

    /// <summary>
    /// Implements <see cref="ISeedableAlgorithm.SetSeed(byte[])"/>.
    /// </summary>
    public void SetSeed(byte[] seed)
    {
        SetSeed(SeedHelper.ToUInt64Array(seed));
    }

    /// <summary>
    /// Sets the internal state from an array of integers. The array can be any length, but
    /// should be at least 1 or more.
    /// </summary>
    public void SetSeed(ulong[] seed)
    {
        // Pre-initialize state
        // This will call CopyCacheFrom()
        SetSeed(InitSeedA);

        int i = 1;
        uint j = 0;
        int jlen = seed.Length;

        int klen = jlen;
        klen = klen switch
        {
            < StateN => StateN,
            _ => klen
        };

        // Use as local
        ulong[] sloc = _mt;

        for (int k = 0; k < klen; ++k)
        {
            ulong a = sloc[i - 1];
            sloc[i] = (sloc[i] ^ ((a ^ (a >> BitWidthM2)) * InitSeedB)) + seed[j] + j;

            if (++i == StateN)
            {
                sloc[0] = sloc[StateN - 1];
                i = 1;
            }

            if (++j == jlen)
            {
                j = 0;
            }
        }

        for (int k = 1; k < StateN; ++k)
        {
            ulong a = sloc[i - 1];
            sloc[i] = (sloc[i] ^ ((a ^ (a >> BitWidthM2)) * InitSeedC)) - (uint)i;

            if (++i != StateN)
            {
                continue;
            }

            sloc[0] = sloc[StateN - 1];
            i = 1;
        }

        sloc[0] = 0x01UL << (BitWidth - 1);
    }

    /// <summary>
    /// Sets the internal state using a single integer value according to the MT seeding routine.
    /// </summary>
    public void SetSeed(ulong seed)
    {
        ulong[] sloc = _mt;

        sloc[0] = seed;

        for (int n = 1; n < StateN; ++n)
        {
            ulong a = sloc[n - 1];
            sloc[n] = InitF * (a ^ (a >> BitWidthM2)) + (ulong)n;
        }

        _idx = StateN;
    }

    /// <summary>
    /// Implements <see cref="ISeedableAlgorithm.Clone"/>.
    /// </summary>
    public ISeedableAlgorithm Clone()
    {
        return new MersenneTwister64(this);
    }

    private void Twist()
    {
        int nrot = 1;
        int mrot = Offset;

        ulong[] sloc = _mt;

        for (int n = 0; n < StateN; ++n)
        {
            ulong x = (sloc[n] & UMask) | (sloc[nrot] & LMask);
            ulong xa = x >> 1;

            switch (x & 0x01UL)
            {
                case 1:
                    xa ^= Matrix;
                    break;
            }

            sloc[n] = sloc[mrot] ^ xa;

            if (++nrot == StateN)
            {
                nrot = 0;
            }

            if (++mrot == StateN)
            {
                mrot = 0;
            }
        }

    }
}