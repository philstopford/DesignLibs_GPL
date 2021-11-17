//--------------------------------------------------------------------------------------------------
// PROJECT      : TRAL
// COPYRIGHT    : Andy Thomas / https://kuiper.zone
// LICENSE      : Apache License, Version 2.0
//--------------------------------------------------------------------------------------------------

namespace Tral.Randomness.Algorithms;

/// <summary>
/// Interface for an underlying random number generation algorithm. Unless the implementation
/// is to utilize hardware generation, it is typically expected that classes will implement <see
/// cref="ISeedableAlgorithm"/> or, where applicable, <see cref="IJumpableAlgorithm"/>. Both of these
/// inherit this as a base interface.
/// </summary>
public interface IRandomAlgorithm
{
    /// <summary>
    /// Gets a short name for the underlying algorithm. The value is to be immutable and should not be null or empty.
    /// </summary>
    string AlgorithmName { get; }

    /// <summary>
    /// Gets the maximum inclusive value of <see cref="Next"/>. It may be any value greater than
    /// 0, however, a value of <see cref="ulong.MaxValue"/> will generally yield the best
    /// performance when used with the <see cref="RandomGenerator{TAlgo}"/> class. The value is
    /// to be immutable.
    /// </summary>
    ulong MaxNext { get; }

    /// <summary>
    /// Returns a random unsigned integer value uniformly distributed in the range [0, <see
    /// cref="MaxNext"/>]. Results are those of the native generation routine which, for
    /// seedable generators, are expected to match test vectors from a known starting state.
    /// </summary>
    ulong Next();
}