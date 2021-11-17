//--------------------------------------------------------------------------------------------------
// PROJECT      : TRAL
// COPYRIGHT    : Andy Thomas / https://kuiper.zone
// LICENSE      : Apache License, Version 2.0
//--------------------------------------------------------------------------------------------------

using Tral.Randomness.Utility;

namespace Tral.Randomness.Algorithms;

/// <summary>
/// Interface for a seedable (pseudo) random number algorithm.
/// </summary>
public interface ISeedableAlgorithm : IRandomAlgorithm
{
    /// <summary>
    /// Gets the seed length in bytes. The value is to be immutable.
    /// </summary>
    int SeedLength { get; }

    /// <summary>
    /// Sets the internal seeds of the generator. Where used with <see cref="RandomGenerator"/>,
    /// the supplied array will always contain at least <see cref="SeedLength"/> bytes. Whether
    /// the algorithm accepts shorter lengths or discards bytes of longer arrays may be
    /// implementation dependent. The array data is to be copied (not referenced assigned).
    /// See <see cref="SeedHelper"/> for conversions.
    /// </summary>
    void SetSeed(byte[] seed);

    /// <summary>
    /// Creates a clone of the generator with the same internal state.
    /// </summary>
    ISeedableAlgorithm Clone();
}