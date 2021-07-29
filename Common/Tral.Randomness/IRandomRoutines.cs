//--------------------------------------------------------------------------------------------------
// PROJECT      : TRAL
// COPYRIGHT    : Andy Thomas / https://kuiper.zone
// LICENSE      : Apache License, Version 2.0
//--------------------------------------------------------------------------------------------------

using System;
using System.Collections.Generic;

namespace Tral.Randomness
{
    /// <summary>
    /// Interface providing random number utility routines. This interface does not allowing seeding
    /// or randomization, but see <see cref="IRandomGenerator"/>.
    /// </summary>
    public interface IRandomRoutines
    {
        /// <summary>
        /// Constant string of lowercase alpha characters for use with <see cref="NextString(int, string)"/>.
        /// </summary>
        const string AlphaLower = "abcdefghijklmnopqrstuvwxyz";

        /// <summary>
        /// Constant string of mixed case alpha characters for use with <see cref="NextString(int, string)"/>.
        /// </summary>
        const string AlphaMixed = AlphaLower + AlphaUpper;

        /// <summary>
        /// Constant string of lowercase alpha numeric characters for use with <see cref="NextString(int, string)"/>.
        /// </summary>
        const string AlphaNumericLower = AlphaLower + Numeric;

        /// <summary>
        /// Constant string of mixed case alpha numeric characters for use with <see cref="NextString(int, string)"/>.
        /// </summary>
        const string AlphaNumericMixed = AlphaMixed + Numeric;

        /// <summary>
        /// Constant string of uppercase alpha numeric characters for use with <see cref="NextString(int, string)"/>.
        /// </summary>
        const string AlphaNumericUpper = AlphaUpper + Numeric;

        /// <summary>
        /// Constant string of uppercase alpha characters for use with <see cref="NextString(int, string)"/>.
        /// </summary>
        const string AlphaUpper = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

        /// <summary>
        /// Constant string of numeric characters for use with <see cref="NextString(int, string)"/>.
        /// </summary>
        const string Numeric = "0123456789";

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

        /// <summary>
        /// Generates an unsigned random integer value in the range [0, 2^32).
        /// </summary>
        uint Next32();

        /// <summary>
        /// Generates an unsigned random integer value in the range [0, 2^64).
        /// </summary>
        ulong Next64();

        /// <summary>
        /// Generates a byte array of "count" length and populates it with random values.
        /// </summary>
        byte[] NextBytes(int count);

        /// <summary>
        /// Populates the array with random byte values. The array should not be null.
        /// </summary>
        void NextBytes(byte[] dest);

        /// <summary>
        /// Populates the array with random byte values, starting from "offset". Returns the number
        /// bytes assigned. The array should not be null.
        /// </summary>
        /// <exception cref="ArgumentOutOfRangeException">Array offset out of range</exception>
        int NextBytes(byte[] dest, int offset);

        /// <summary>
        /// Populates the array with "count" random values, starting from the given "offset". If
        /// count is -1, values are written from offset to the end of the array. Returns the number
        /// bytes assigned. The array should not be null.
        /// </summary>
        /// <exception cref="ArgumentException">Array count exceeds length of array</exception>
        /// <exception cref="ArgumentOutOfRangeException">Array offset out of range</exception>
        int NextBytes(byte[] dest, int offset, int count);

        /// <summary>
        /// Generates a random double value uniformly distributed in the range [0, 1.0).
        /// </summary>
        double NextDouble();

        /// <summary>
        /// Generates a random double value uniformly distributed in the range [min, max). Either or
        /// both "min" or "max" may be negative provided "max" is greater than "min".
        /// </summary>
        /// <exception cref="ArgumentOutOfRangeException">Invalid range</exception>
        double NextDouble(double min, double max);

        /// <summary>
        /// Returns true or false, with equal probability.
        /// </summary>
        bool NextFlip();

        /// <summary>
        /// Generates a random integer value in the range [0, max). The range is not restricted by
        /// <see cref="MaxNext"/>.
        /// </summary>
        /// <exception cref="ArgumentOutOfRangeException">Invalid range</exception>
        int NextInt(int max);

        /// <summary>
        /// Generates a random integer value in the range [min, max). Either or both "min" or "max"
        /// may be negative provided "max" is greater than "min". The range is not restricted by
        /// <see cref="MaxNext"/>.
        /// </summary>
        /// <exception cref="ArgumentOutOfRangeException">Invalid range</exception>
        int NextInt(int min, int max);

        /// <summary>
        /// Generates a random integer value in the range [min, max). The result will exclude the
        /// "exclude" value. The range must greater than 1 in magnitude.
        /// </summary>
        /// <exception cref="ArgumentOutOfRangeException">Invalid range</exception>
        int NextIntExclude(int min, int max, int exclude);

        /// <summary>
        /// A <see cref="long"/> variant of <see cref="NextInt(int)"/>.
        /// </summary>
        long NextLong(long max);

        /// <summary>
        /// A <see cref="long"/> variant of <see cref="NextInt(int, int)"/>.
        /// </summary>
        long NextLong(long min, long max);

        /// <summary>
        /// Generates a random integer value in the range [min, max). The result will exclude the
        /// "exclude" value. The range must greater than 1 in magnitude.
        /// </summary>
        /// <exception cref="ArgumentOutOfRangeException">Invalid range</exception>
        long NextLongExclude(long min, long max, long exclude);

        /// <summary>
        /// Generates a random double value uniformly distributed in the range (0, 1.0).
        /// </summary>
        double NextOpenDouble();

        /// <summary>
        /// Returns a normally distributed (Gaussian) random double value, with a mean of 0 and standard
        /// deviation of +1.0.
        /// </summary>
        double NextStdNormal();

        /// <summary>
        /// Generates a string of "count" length with characters randomly drawn from "alphabet". The
        /// result may contain repeated characters. A range of default alphabet constants are defined
        /// such as <see cref="AlphaUpper"/>.
        /// </summary>
        /// <exception cref="ArgumentOutOfRangeException">Count less than 0</exception>
        string NextString(int count, string alphabet);

        /// <summary>
        /// Shuffles the items in the array.
        /// </summary>
        void Shuffle<T>(T[] items);

        /// <summary>
        /// Shuffles the items in the list.
        /// </summary>
        void Shuffle<T>(IList<T> items);

        /// <summary>
        /// Shuffles "count" items in the array, starting from the given offset. If count is -1,
        /// items are shuffled from "offset" to the end of the array. If the array is empty, the
        /// method throws unless "offset" is 0 and "count" is 0 or less. It does nothing where
        /// "count" is 0 or 1.
        /// </summary>
        /// <exception cref="ArgumentException">Count exceeds length</exception>
        /// <exception cref="ArgumentOutOfRangeException">Offset out of range</exception>
        void Shuffle<T>(T[] items, int offset, int count);

        /// <summary>
        /// Shuffles "count" items in the list, starting from the given offset. See <see
        /// cref="Shuffle{T}(T[], int, int)"/> for information.
        /// </summary>
        /// <exception cref="ArgumentException">Count exceeds length</exception>
        /// <exception cref="ArgumentOutOfRangeException">Offset out of range</exception>
        void Shuffle<T>(IList<T> items, int offset, int count);

        /// <summary>
        /// Shuffles the string characters and returns the result.
        /// </summary>
        string Shuffle(string s);
    }
}