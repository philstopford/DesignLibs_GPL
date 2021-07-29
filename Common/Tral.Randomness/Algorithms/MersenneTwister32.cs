//--------------------------------------------------------------------------------------------------
// PROJECT      : TRAL
// COPYRIGHT    : Andy Thomas / https://kuiper.zone
// LICENSE      : Apache License, Version 2.0
//--------------------------------------------------------------------------------------------------

using System;
using Tral.Randomness.Utility;

namespace Tral.Randomness.Algorithms
{
    /// <summary>
    /// Mersenne Twister, or MT19937-32, is a widely used 32-bit pseudo random number generator. It
    /// is not suited to cryptographic purposes. On construction, the initial state is seeded with
    /// the integer value 5489 according to the MT seeding routine.
    /// </summary>
    /// <remarks>
    /// Mersenne Twister (MT) was originally developed by Makoto Matsumoto and Takuji Nishimura in 1996/1997.
    /// </remarks>
    public class MersenneTwister32 : ISeedableAlgorithm
    {
        private const uint DefaultSeed = 5489U;
        private const int StateN = 624;

        private const int BitWidth = 32;
        private const int BitWidthM2 = 30;

        private const uint InitSeedA = 19650218U;
        private const uint InitSeedB = 1664525U;
        private const uint InitSeedC = 1566083941U;

        private const uint UMask = 0x80000000U;
        private const uint LMask = 0x7FFFFFFFU;
        private const uint Matrix = 0x9908B0DFU;
        private const int ShiftU = 11;
        private const int ShiftS = 7;
        private const uint MaskB = 0x9D2C5680U;
        private const int ShiftT = 15;
        private const uint MaskC = 0xEFC60000U;
        private const int ShiftL = 18;
        private const uint InitF = 0x6C078965U;
        private const int Offset = 397;

        private int _idx;
        private readonly uint[] _mt = new uint[StateN];

        /// <summary>
        /// Default constructor. The initial state is seeded with the integer value 5489 according
        /// to the MT seeding routine.
        /// </summary>
        public MersenneTwister32()
        {
            SetSeed(DefaultSeed);
        }

        private MersenneTwister32(MersenneTwister32 other)
        {
            _idx = other._idx;
            Buffer.BlockCopy(other._mt, 0, _mt, 0, SeedLength);
        }

        /// <summary>
        /// Implements <see cref="IRandomAlgorithm.AlgorithmName"/>.
        /// </summary>
        public string AlgorithmName { get; } = "MT19937-32";

        /// <summary>
        /// Implements <see cref="IRandomAlgorithm.MaxNext"/>.
        /// </summary>
        public ulong MaxNext { get; } = uint.MaxValue;

        /// <summary>
        /// Implements <see cref="ISeedableAlgorithm.SeedLength"/>.
        /// </summary>
        public int SeedLength { get; } = 2496;

        /// <summary>
        /// Implements <see cref="IRandomAlgorithm.Next"/>.
        /// </summary>
        public ulong Next()
        {
            if (_idx == StateN)
            {
                Twist();
                _idx = 0;
            }

            uint y = _mt[_idx++];

            y ^= y >> ShiftU;
            y ^= (y << ShiftS) & MaskB;
            y ^= (y << ShiftT) & MaskC;
            return y ^ (y >> ShiftL);
        }

        /// <summary>
        /// Implements <see cref="ISeedableAlgorithm.SetSeed(byte[])"/>.
        /// </summary>
        public void SetSeed(byte[] seed)
        {
            SetSeed(SeedHelper.ToUInt32Array(seed));
        }

        /// <summary>
        /// Sets the internal state from an array of integers. The array can be any length, but
        /// should be at least 1 or more.
        /// </summary>
        public void SetSeed(uint[] seed)
        {
            // Pre-initialize state
            // This will call CopyCacheFrom()
            SetSeed(InitSeedA);

            int i = 1;
            uint j = 0;
            int jlen = seed.Length;

            int klen = jlen;
            if (klen < StateN) klen = StateN;

            // Use as local
            var sloc = _mt;

            for (int k = 0; k < klen; ++k)
            {
                uint a = sloc[i - 1];
                sloc[i] = (sloc[i] ^ ((a ^ (a >> BitWidthM2)) * InitSeedB)) + seed[j] + j;

                if (++i == StateN)
                {
                    sloc[0] = sloc[StateN - 1];
                    i = 1;
                }

                if (++j == jlen) j = 0;
            }

            for (int k = 1; k < StateN; ++k)
            {
                uint a = sloc[i - 1];
                sloc[i] = (sloc[i] ^ ((a ^ (a >> BitWidthM2)) * InitSeedC)) - (uint)i;

                if (++i == StateN)
                {
                    sloc[0] = sloc[StateN - 1];
                    i = 1;
                }
            }

            sloc[0] = 0x01U << (BitWidth - 1);
        }

        /// <summary>
        /// Sets the internal state using a single integer value according to the MT seeding routine.
        /// </summary>
        public void SetSeed(uint seed)
        {
            uint a;
            var sloc = _mt;

            sloc[0] = seed;

            for (int n = 1; n < StateN; ++n)
            {
                a = sloc[n - 1];
                sloc[n] = InitF * (a ^ (a >> BitWidthM2)) + (uint)n;
            }

            _idx = StateN;
        }

        /// <summary>
        /// Implements <see cref="ISeedableAlgorithm.Clone"/>.
        /// </summary>
        public ISeedableAlgorithm Clone()
        {
            return new MersenneTwister32(this);
        }

        private void Twist()
        {
            int nrot = 1;
            int mrot = Offset;

            var sloc = _mt;

            for (int n = 0; n < StateN; ++n)
            {
                uint x = (sloc[n] & UMask) | (sloc[nrot] & LMask);
                uint xa = x >> 1;

                if ((x & 0x01U) == 1) xa ^= Matrix;

                sloc[n] = sloc[mrot] ^ xa;

                if (++nrot == StateN) nrot = 0;
                if (++mrot == StateN) mrot = 0;
            }

        }
    }
}