//--------------------------------------------------------------------------------------------------
// PROJECT      : TRAL
// COPYRIGHT    : Andy Thomas / https://kuiper.zone
// LICENSE      : Apache License, Version 2.0
//--------------------------------------------------------------------------------------------------

using System;
using Tral.Randomness.Utility;

namespace Tral.Randomness.Algorithms;

/// <summary>
/// WELL512 is a fast 32-bit all purposes pseudo random number generator with 512 bits of
/// internal state. It is not suited to cryptographic purposes. On construction, the initial
/// state is seeded with arbitrary (hard-coded) values.
/// </summary>
/// <remarks>
/// This implementation is derived from code put in the public domain by Chris Lomont. Refer
/// to: http://www.lomont.org/Math/Papers/2008/Lomont_PRNG_2008.pdf
/// </remarks>
public class Well512 : ISeedableAlgorithm
{
    private int _index;
    private readonly uint[] _s = new uint[16];

    /// <summary>
    /// Default constructor.
    /// </summary>
    public Well512()
    {
        // We choose an arbitrary starting state.
        _s[0] = 0x94C24F02U;
        _s[1] = 0x2AD46D70U;
        _s[2] = 0xA41E5A66U;
        _s[3] = 0x5C8EA936U;
        _s[4] = 0xE1023BC2U;
        _s[5] = 0xF926F3B2U;
        _s[6] = 0x018AF11FU;
        _s[7] = 0x45495ADBU;
        _s[8] = 0xB1C4540CU;
        _s[9] = 0x531D2701U;
        _s[10] = 0x86FED955U;
        _s[11] = 0x4A47440FU;
        _s[12] = 0x5B8B06F1U;
        _s[13] = 0xD0C8F8CBU;
        _s[14] = 0x950110E9U;
        _s[15] = 0x2C797BA8U;
    }

    private Well512(Well512 other)
    {
        _index = other._index;
        Buffer.BlockCopy(other._s, 0, _s, 0, SeedLength);
    }

    /// <summary>
    /// Implements <see cref="IRandomAlgorithm.AlgorithmName"/>.
    /// </summary>
    public string AlgorithmName => "WELL512";

    /// <summary>
    /// Implements <see cref="IRandomAlgorithm.MaxNext"/>.
    /// </summary>
    public ulong MaxNext => uint.MaxValue;

    /// <summary>
    /// Implements <see cref="ISeedableAlgorithm.SeedLength"/>.
    /// </summary>
    public int SeedLength => 64;

    /// <summary>
    /// Implements <see cref="IRandomAlgorithm.Next"/>.
    /// </summary>
    public ulong Next()
    {
        uint a = _s[_index];
        uint c = _s[(_index + 13) & 0x0F];
        uint b = a ^ c ^ (a << 16) ^ (c << 15);

        c = _s[(_index + 9) & 15];
        c ^= c >> 11;
        a = _s[_index] = b ^ c;

        uint d = a ^ ((a << 5) & 0xDA442D24U);

        _index = (_index + 15) & 0x0F;

        a = _s[_index];
        _s[_index] = a ^ b ^ d ^ (a << 2) ^ (b << 18) ^ (c << 28);

        return _s[_index];
    }

    /// <summary>
    /// Implements <see cref="ISeedableAlgorithm.SetSeed(byte[])"/>.
    /// </summary>
    public void SetSeed(byte[] seed)
    {
        SetSeed(SeedHelper.ToUInt32Array(seed, SeedLength));
    }

    /// <summary>
    /// Sets the internal state from an array of 32-bit integers. The array must be of length 16 or greater.
    /// </summary>
    public void SetSeed(uint[] seed)
    {
        _index = 0;
        Buffer.BlockCopy(seed, 0, _s, 0, SeedLength);
    }

    /// <summary>
    /// Implements <see cref="ISeedableAlgorithm.Clone"/>.
    /// </summary>
    public ISeedableAlgorithm Clone()
    {
        return new Well512(this);
    }

}