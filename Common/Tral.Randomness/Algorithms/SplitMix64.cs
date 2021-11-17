//--------------------------------------------------------------------------------------------------
// PROJECT      : TRAL
// COPYRIGHT    : Andy Thomas / https://kuiper.zone
// LICENSE      : Apache License, Version 2.0
//--------------------------------------------------------------------------------------------------

using System.Runtime.CompilerServices;

namespace Tral.Randomness.Algorithms;

/// <summary>
/// SplitMix64 is a fast all purposes (non-cryptographic) pseudo random number generator, but
/// with only 64 bits of internal state. The implementation is a fixed increment version of Java
/// 8's SplittableRandom generator. It is passes "Big Crush" and is useful where only 64 bits of
/// state is required. On construction, the initial state is populated with a fixed arbitrary value.
/// </summary>
/// <remarks>
/// This implementation is derived from code put in the public domain by Sebastiano Vigna. Refer
/// to: http://xoroshiro.di.unimi.it
/// </remarks>
public class SplitMix64 : ISeedableAlgorithm
{
    private ulong _s = 0x675187BC468FEF54UL;

    /// <summary>
    /// Implements <see cref="IRandomAlgorithm.AlgorithmName"/>.
    /// </summary>
    public string AlgorithmName => "SplitMix64";

    /// <summary>
    /// Implements <see cref="IRandomAlgorithm.MaxNext"/>.
    /// </summary>
    public ulong MaxNext => ulong.MaxValue;

    /// <summary>
    /// Implements <see cref="ISeedableAlgorithm.SeedLength"/>.
    /// </summary>
    public int SeedLength => 8;

    /// <summary>
    /// Implements <see cref="IRandomAlgorithm.Next"/>.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public ulong Next()
    {
        ulong z = _s += 0x9E3779B97F4A7C15UL;
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9UL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBUL;
        return z ^ (z >> 31);
    }

    /// <summary>
    /// Implements <see cref="ISeedableAlgorithm.SetSeed(byte[])"/>.
    /// </summary>
    public void SetSeed(byte[] seed)
    {
        ulong s = seed[0];

        s = (s << 8) | seed[1];
        s = (s << 8) | seed[2];
        s = (s << 8) | seed[3];
        s = (s << 8) | seed[4];
        s = (s << 8) | seed[5];
        s = (s << 8) | seed[6];

        _s = (s << 8) | seed[7];
    }

    /// <summary>
    /// Overload with 64-bit integer.
    /// </summary>
    public void SetSeed(ulong seed)
    {
        _s = seed;
    }

    /// <summary>
    /// Implements <see cref="ISeedableAlgorithm.Clone"/>.
    /// </summary>
    public ISeedableAlgorithm Clone()
    {
        SplitMix64 clone = new()
        {
            _s = _s
        };
        return clone;
    }

}