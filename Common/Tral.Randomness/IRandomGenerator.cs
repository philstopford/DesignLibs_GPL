//--------------------------------------------------------------------------------------------------
// PROJECT      : TRAL
// COPYRIGHT    : Andy Thomas / https://kuiper.zone
// LICENSE      : Apache License, Version 2.0
//--------------------------------------------------------------------------------------------------

using System;
using Tral.Randomness.Algorithms;

namespace Tral.Randomness;

/// <summary>
/// Random number generator interface. This extends <see cref="IRandomRoutines"/> with seed and jump routines.
/// </summary>
public interface IRandomGenerator : IRandomRoutines
{
    /// <summary>
    /// Gets whether the generator is seedable. The value is true if the underlying algorithm
    /// implements <see cref="ISeedableAlgorithm"/>. If false, methods which require a seedable
    /// algorithm throw <see cref="NotSupportedException"/>.
    /// </summary>
    bool IsSeedable { get; }

    /// <summary>
    /// Gets whether the generator is jumpable. The value is true if the underlying algorithm
    /// implements <see cref="IJumpableAlgorithm"/>. If false, methods which require a jumpable
    /// algorithm throw <see cref="NotSupportedException"/>.
    /// </summary>
    bool IsJumpable { get; }

    /// <summary>
    /// Gets the seed length in bytes. The value is to be immutable.
    /// </summary>
    int SeedLength { get; }

    /// <summary>
    /// Creates a clone of the generator with the same internal state. The <see cref="IsSeedable"/> property
    /// must be true in order for the instance to be clonable.
    /// </summary>
    /// <exception cref="NotSupportedException"><see cref="IsSeedable"/> is false</exception>
    IRandomGenerator Clone();

    /// <summary>
    /// Advances the internal state "n" times. Where n is large, this call may be expensive,
    /// however, it can be used where <see cref="IsJumpable"/> is false. Note that it also
    /// resets internal cache and, therefore, may not be directly equivalent to calling <see
    /// cref="IRandomRoutines.Next"/> the specified number of times. A value of 0 or less
    /// resets internal cache but does not advance the underlying generator. See also <see cref="Jump"/>.
    /// </summary>
    void Discard(int n);

    /// <summary>
    /// Performs a "fast jump", i.e. an efficient operation equivalent to calling <see
    /// cref="Discard(int)"/> with a very large "n" value. Note that this call resets internal
    /// cache and, therefore, may not be directly equivalent to calling <see
    /// cref="IRandomRoutines.Next"/> an equivalent number of times.
    /// </summary>
    /// <exception cref="NotSupportedException"><see cref="IsJumpable"/> is false</exception>
    void Jump();

    /// <summary>
    /// Equivalent to NewInstance(true).
    /// </summary>
    IRandomGenerator NewInstance();

    /// <summary>
    /// Creates a new instance of the generator using the same underlying <see
    /// cref="IRandomAlgorithm"/> implementation. The new instance is randomized if "randomize"
    /// is true, otherwise it will be in its default state. Unlike <see cref="Clone"/>, internal
    /// state is not copied.
    /// </summary>
    IRandomGenerator NewInstance(bool randomize);

    /// <summary>
    /// Randomizes the internal state.
    /// </summary>
    /// <exception cref="NotSupportedException"><see cref="IsSeedable"/> is false</exception>
    void Randomize();

    /// <summary>
    /// Seeds the internal state of the generator using a single "seed" integer value. The <see
    /// cref="SetSeed(byte[])"/> is called with a byte sequence derived from the supplied
    /// integer value. The algorithm for extending "seed" is not defined (i.e. implementation
    /// specific subject to future change). Use <see cref="SetSeed(byte[])"/> where it is
    /// required to seed the algorithm in a well defined way.
    /// </summary>
    /// <exception cref="NotSupportedException"><see cref="IsSeedable"/> is false</exception>
    void SetSeed(ulong seed);

    /// <summary>
    /// Signed overload of <see cref="SetSeed(ulong)"/>.
    /// </summary>
    /// <exception cref="NotSupportedException"><see cref="IsSeedable"/> is false</exception>
    void SetSeed(long seed);

    /// <summary>
    /// Sets the internal seeds of the generator. The "seed" array must be of length <see
    /// cref="ISeedableAlgorithm.SeedLength"/> or greater (whether additional bytes are ignored
    /// or used is algorithm dependent). The seed is supplied directly to the underlying <see
    /// cref="ISeedableAlgorithm"/> instance. The same byte sequence shall yield the same
    /// internal state irrespective of endian.
    /// </summary>
    /// <exception cref="NotSupportedException"><see cref="IsSeedable"/> is false</exception>
    /// <exception cref="ArgumentException">Array length too short</exception>
    void SetSeed(byte[] seed);
}