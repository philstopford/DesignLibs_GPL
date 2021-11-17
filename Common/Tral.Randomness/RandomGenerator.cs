//--------------------------------------------------------------------------------------------------
// PROJECT      : TRAL
// COPYRIGHT    : Andy Thomas / https://kuiper.zone
// LICENSE      : Apache License, Version 2.0
//--------------------------------------------------------------------------------------------------

using System;
using Tral.Randomness.Algorithms;

namespace Tral.Randomness;

/// <summary>
/// A pseudo random number generator utilizing a general purpose (non-cryptographic) jumpable
/// algorithm. It should be considered be that the algorithm used is undefined by the class
/// specification and may be subject to future change. Use <see cref="RandomGenerator{TAlgo}"/>
/// where a specific underlying algorithm is required. See also <see cref="IRandomGenerator"/>
/// and <see cref="IRandomRoutines"/> for information.
/// </summary>
public sealed class RandomGenerator : RandomGenerator<Xoshiro256pp>
{
    [ThreadStatic]
    private static IRandomRoutines? st_global;

    /// <summary>
    /// Constructor. <see cref="IRandomGenerator.Randomize"/> is called if "randomize" is true, otherwise
    /// the generator will be in its default state.
    /// </summary>
    public RandomGenerator(bool randomize = true)
        : base(randomize)
    {
    }

    /// <summary>
    /// Gets global (singleton) random routines for convenience. The result is a thread-local
    /// instance and is randomized on its construction.
    /// </summary>
    public static IRandomRoutines Global
    {
        get
        {
            IRandomRoutines? g = st_global;

            if (g != null)
            {
                return g;
            }

            g = new RandomGenerator();
            st_global = g;

            return g;
        }
    }
}