//--------------------------------------------------------------------------------------------------
// PROJECT      : TRAL
// COPYRIGHT    : Andy Thomas / https://kuiper.zone
// LICENSE      : Apache License, Version 2.0
//--------------------------------------------------------------------------------------------------

using System;
using System.Runtime.CompilerServices;
using System.Security.Cryptography;
using Tral.Randomness.Algorithms;

namespace Tral.Randomness.Utility
{
    /// <summary>
    /// Seeding utility routines.
    /// </summary>
    public static class SeedHelper
    {
        private const int RunSplitMix = 8;

        /// <summary>
        /// Gets whether the machine is little endian.
        /// </summary>
        public static bool IsLittleEndian { get; } = BitConverter.IsLittleEndian;

        /// <summary>
        /// Generates "count" bytes of seed according using <see cref="RNGCryptoServiceProvider"/>.
        /// </summary>
        public static byte[] GenerateSeed(int count)
        {
            if (count != 0)
            {
                byte[] bytes = new byte[count];

                using (var c = new RNGCryptoServiceProvider())
                {
                    c.GetBytes(bytes);
                }

                return bytes;
            }

            return Array.Empty<byte>();
        }

        /// <summary>
        /// Extends the given integer seed value and returns and array of "count" bytes. The method
        /// should be considered implementation specific (i.e. may change between updates) and
        /// is unsuitable for cryptography.
        /// </summary>
        public static byte[] ExtendSeed(ulong seed, int count)
        {
            var bytes = new byte[count];

            // Extend the seed
            var smix = new SplitMix64();
            smix.SetSeed(seed);

            // Run it in a little to get rid of 0 state
            for (int n = 0; n < RunSplitMix; ++n)
            {
                smix.Next();
            }

            for (int n = 0; n < count; ++n)
            {
                bytes[n] = (byte)(smix.Next() >> 8);
            }

            return bytes;
        }

        /// <summary>
        /// Converts a maximum of "max" bytes from the array to a corresponding array of 32-bit
        /// unsigned integers. Note that it is assumed the byte data is big endian and results are
        /// transposed where <see cref="IsLittleEndian"/> is true. Therefore, the same byte values
        /// result in the same numeric integers irrespective of endian. The resulting length is the
        /// byte length divided by 4 (remaining bytes are discarded). If "max" is negative, all
        /// bytes are converted. If the byte length is less than 4, the result is an empty array.
        /// </summary>
        public static uint[] ToUInt32Array(byte[] bytes, int max = -1)
        {
            const int SizeOfInt = 4;

            if (max < 0)
            {
                max = bytes.Length;
            }

            int len = Math.Min(bytes.Length, max) / SizeOfInt;

            var rslt = new uint[len];
            Buffer.BlockCopy(bytes, 0, rslt, 0, len * SizeOfInt);

            if (IsLittleEndian)
            {
                for (int n = 0; n < len; ++n)
                {
                    rslt[n] = SwapEndian(rslt[n]);
                }
            }

            return rslt;
        }

        /// <summary>
        /// Converts a maximum of "max" bytes from the array to a corresponding array of 64-bit
        /// unsigned integers. Note that it is assumed the byte data is big endian and results are
        /// transposed where <see cref="IsLittleEndian"/> is true. Therefore, the same byte values
        /// result in the same numeric integers irrespective of endian. The resulting length is the
        /// byte length divided by 8 (remaining bytes are discarded). If "max" is negative, all
        /// bytes are converted. If the byte length is less than 8, the result is an empty array.
        /// </summary>
        public static ulong[] ToUInt64Array(byte[] bytes, int max = -1)
        {
            const int SizeOfLong = 8;

            if (max < 0)
            {
                max = bytes.Length;
            }

            int len = Math.Min(bytes.Length, max) / SizeOfLong;

            var rslt = new ulong[len];
            Buffer.BlockCopy(bytes, 0, rslt, 0, len * SizeOfLong);

            if (IsLittleEndian)
            {
                for (int n = 0; n < len; ++n)
                {
                    rslt[n] = SwapEndian(rslt[n]);
                }
            }

            return rslt;
        }

        /// <summary>
        /// Array bound check. Confirms validity of "offset" and "count" within the array length or
        /// throws appropriately. On success, the return is "count", or if count is negative, the
        /// number of bytes from "offset" to the end of the array. If the array is empty, the method
        /// throws unless "offset" is 0 and "count" is 0 or less.
        /// </summary>
        /// <exception cref="ArgumentException">Count exceeds length</exception>
        /// <exception cref="ArgumentOutOfRangeException">Offset out of range</exception>
        public static int AssertBounds(int length, int offset, int count)
        {
            if (offset >= length || offset < 0)
            {
                if (offset == 0 && count <= 0)
                {
                    // Legal exception for empty array
                    return 0;
                }

                throw new ArgumentOutOfRangeException("offset", "Offset out of range");
            }

            if (count < 0)
            {
                return length - offset;
            }

            if ((long)offset + count > length)
            {
                throw new ArgumentException("count", "Count exceeds length");
            }

            return count;
        }

        /// <summary>
        /// Exchanges the endian of "x" and returns the result.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static uint SwapEndian(uint x)
        {
            return (x & 0x000000FF) << 24 | (x & 0x0000FF00) << 8 |
                (x & 0x00FF0000) >> 8 | (x & 0xFF000000) >> 24;
        }

        /// <summary>
        /// Exchanges the endian of "x" and returns the result.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ulong SwapEndian(ulong x)
        {
            return (x & 0x00000000000000FFUL) << 56 | (x & 0x000000000000FF00UL) << 40 |
                (x & 0x0000000000FF0000UL) << 24 | (x & 0x00000000FF000000UL) << 8 |
                (x & 0x000000FF00000000UL) >> 8 | (x & 0x0000FF0000000000UL) >> 24 |
                (x & 0x00FF000000000000UL) >> 40 | (x & 0xFF00000000000000UL) >> 56;
        }

    }
}
