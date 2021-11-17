//--------------------------------------------------------------------------------------------------
// PROJECT      : TRAL
// COPYRIGHT    : Andy Thomas / https://kuiper.zone
// LICENSE      : Apache License, Version 2.0
//--------------------------------------------------------------------------------------------------

namespace Tral.Randomness.Algorithms;

/// <summary>
/// Interface for a jumpable random number generation algorithm. Its <see cref="Jump"/> method
/// is to provide an efficient means to advance internal state by a large number of equivalent
/// calls to <see cref="IRandomAlgorithm.Next"/>. It can be used to ensure that generators,
/// initially seeded with the same value, do not produce the same or overlapping results in
/// parallel scenarios where repeatability is required. Jumpable generators must be seedable
/// and, therefore, the interface inherits <see cref="ISeedableAlgorithm"/>.
/// </summary>
public interface IJumpableAlgorithm : ISeedableAlgorithm
{
    /// <summary>
    /// Performs a "fast jump", i.e. an efficient operation equivalent to calling <see
    /// cref="IRandomAlgorithm.Next"/> a very large number of times.
    /// </summary>
    void Jump();
}