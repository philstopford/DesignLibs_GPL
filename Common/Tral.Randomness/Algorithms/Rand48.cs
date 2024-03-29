﻿//--------------------------------------------------------------------------------------------------
// PROJECT      : TRAL
// COPYRIGHT    : Andy Thomas / https://kuiper.zone
// LICENSE      : Apache License, Version 2.0
//--------------------------------------------------------------------------------------------------

using System.Runtime.CompilerServices;

namespace Tral.Randomness.Algorithms;

/// <summary>
/// The RAND48 linear congruential algorithm (lrand48). This is a widely used but legacy
/// generator with 48 bits of internal state. It is not suitable to cryptographic applications.
/// On construction, the initial state is populated with a fixed arbitrary value.
/// </summary>
public class Rand48 : ISeedableAlgorithm
{
    private ulong _s;

    /// <summary>
    /// Constructor.
    /// </summary>
    public Rand48()
    {
        SetSeed(0xD0C3F8CBU);
    }

    /// <summary>
    /// Implements <see cref="IRandomAlgorithm.AlgorithmName"/>.
    /// </summary>
    public string AlgorithmName => "RAND48";

    /// <summary>
    /// Implements <see cref="IRandomAlgorithm.MaxNext"/>.
    /// </summary>
    public ulong MaxNext => int.MaxValue;

    /// <summary>
    /// Implements <see cref="ISeedableAlgorithm.SeedLength"/>.
    /// </summary>
    public int SeedLength => 6;

    /// <summary>
    /// Implements <see cref="IRandomAlgorithm.Next"/>.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public ulong Next()
    {
        // Modulus of 2^48
        _s = (25214903917UL * _s + 11UL) & 0x0000FFFFFFFFFFFFUL;
        return _s >> 17;
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

        SetSeed((s << 8) | seed[5]);
    }

    /// <summary>
    /// Sets the internal state from the seed value.
    /// </summary>
    public void SetSeed(ulong seed)
    {
        _s = (seed << 16) | 0x0000330EU;
    }

    /// <summary>
    /// Implements <see cref="ISeedableAlgorithm.Clone"/>.
    /// </summary>
    public ISeedableAlgorithm Clone()
    {
        Rand48 clone = new()
        {
            _s = _s
        };
        return clone;
    }

}