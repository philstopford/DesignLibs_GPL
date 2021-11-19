//--------------------------------------------------------------------------------------------------
// PROJECT      : TRAL
// COPYRIGHT    : Andy Thomas / https://kuiper.zone
// LICENSE      : Apache License, Version 2.0
//--------------------------------------------------------------------------------------------------

using System;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using System.Text;
using Tral.Randomness.Algorithms;
using Tral.Randomness.Utility;

namespace Tral.Randomness;

/// <summary>
/// Implements <see cref="IRandomGenerator"/> utilizing the underlying algorithm specified by
/// the generic type "TAlgo". TAlgo must implement <see cref="IRandomAlgorithm"/> and may
/// optionally implement <see cref="ISeedableAlgorithm"/> and <see cref="IJumpableAlgorithm"/>
/// to yield a seedable and jumpable generator respectively. Where "TAlgo" implements only <see
/// cref="IRandomAlgorithm"/>, the resulting instance will not be capable of being seeded. Where
/// it is necessary to pass an instance of this class without declaring the generic type, pass
/// it either as a type of <see cref="IRandomGenerator"/> or <see cref="IRandomRoutines"/>.
/// </summary>
public class RandomGenerator<TAlgo> : IRandomGenerator
    where TAlgo : class, IRandomAlgorithm, new()
{
    private readonly IRandomAlgorithm _algorithm;
    private readonly IJumpableAlgorithm? _jumpable;
    private readonly ISeedableAlgorithm? _seedable;
    private readonly int _shiftBits;
    private readonly ulong _shiftMax;

    // Improves performance in many scenarios
    private ulong _flipCache;
    private int _flipsRemain;
    private double _gaussCache;
    private bool _gaussFlag;
    private uint _next32Cache;
    private bool _next32Flag;

    /// <summary>
    /// Constructor. If "TAlgo" inherits <see cref="ISeedableAlgorithm"/>, <see
    /// cref="IRandomGenerator.Randomize"/> is called if "randomize" is true, otherwise the
    /// generator will be in its default state. The "randomize" value is ignored where TAlgo
    /// inherits only <see cref="IRandomAlgorithm"/>.
    /// </summary>
    /// <exception cref="InvalidOperationException">Invalid algorithm</exception>
    public RandomGenerator(bool randomize = true)
        : this(new TAlgo(), randomize)
    {
    }

    private RandomGenerator(IRandomAlgorithm algo, bool randomize)
    {
        _algorithm = algo;
        AlgorithmName = algo.AlgorithmName;
        MaxNext = algo.MaxNext;

        switch (MaxNext)
        {
            case < 1:
                throw new InvalidOperationException(nameof(MaxNext) + " cannot be 0 or less");
        }

        if (MaxNext != uint.MaxValue && MaxNext != ulong.MaxValue)
        {
            // Get number of bits to use for GetUnbiasedBits().
            // A NextMax of 2^31-1 gives a shift of 31.
            _shiftBits = BitSize(MaxNext);
            _shiftMax = ulong.MaxValue >> (64 - _shiftBits);

            if (MaxNext != _shiftMax)
            {
                _shiftBits -= 1;
                _shiftMax >>= 1;
            }
        }

        if (algo is not ISeedableAlgorithm seedable)
        {
            return;
        }

        IsSeedable = true;
        _seedable = seedable;
        SeedLength = seedable.SeedLength;

        switch (SeedLength)
        {
            case < 1:
                throw new InvalidOperationException(nameof(SeedLength) + " cannot be 0 or less");
        }

        switch (algo)
        {
            case IJumpableAlgorithm jumpable:
                IsJumpable = true;
                _jumpable = jumpable;
                break;
        }

        switch (randomize)
        {
            case true:
                Randomize();
                break;
        }
    }

    /// <summary>
    /// Implements <see cref="IRandomAlgorithm.AlgorithmName"/>.
    /// </summary>
    public string AlgorithmName { get; }

    /// <summary>
    /// Implements <see cref="IRandomGenerator.IsJumpable"/>. It is true if "TAlgo" inherits <see cref="IJumpableAlgorithm"/>.
    /// </summary>
    public bool IsJumpable { get; }

    /// <summary>
    /// Implements <see cref="IRandomGenerator.IsSeedable"/>. It is true if "TAlgo" inherits <see cref="ISeedableAlgorithm"/>.
    /// </summary>
    public bool IsSeedable { get; }

    /// <summary>
    /// Implements <see cref="IRandomAlgorithm.MaxNext"/>.
    /// </summary>
    public ulong MaxNext { get; }

    /// <summary>
    /// Implements <see cref="ISeedableAlgorithm.SeedLength"/>.
    /// </summary>
    public int SeedLength { get; }

    /// <summary>
    /// Implements <see cref="IRandomGenerator.Clone"/>.
    /// </summary>
    public IRandomGenerator Clone()
    {
        switch (_seedable)
        {
            case null:
                throw new NotSupportedException(nameof(ISeedableAlgorithm) + " operation not supported for " + AlgorithmName);
            default:
            {
                RandomGenerator<TAlgo> clone = new(_seedable.Clone(), false)
                {
                    // Clone state of this instance also
                    _next32Flag = _next32Flag,
                    _next32Cache = _next32Cache,
                    _flipCache = _flipCache,
                    _flipsRemain = _flipsRemain,
                    _gaussFlag = _gaussFlag,
                    _gaussCache = _gaussCache
                };

                return clone;
            }
        }
    }

    /// <summary>
    /// Implements <see cref="IRandomGenerator.Discard(int)"/>.
    /// </summary>
    public void Discard(int n)
    {
        for (int i = 0; i < n; ++i)
        {
            _algorithm.Next();
        }

        ZeroCache();
    }

    /// <summary>
    /// Implements <see cref="IRandomGenerator.NewInstance()"/>.
    /// </summary>
    public IRandomGenerator NewInstance()
    {
        return new RandomGenerator<TAlgo>();
    }

    /// <summary>
    /// Implements <see cref="IRandomGenerator.NewInstance(bool)"/>.
    /// </summary>
    public IRandomGenerator NewInstance(bool randomize)
    {
        return new RandomGenerator<TAlgo>(randomize);
    }

    /// <summary>
    /// Implements <see cref="IJumpableAlgorithm.Jump"/>.
    /// </summary>
    public void Jump()
    {
        switch (_jumpable)
        {
            case null:
                throw new NotSupportedException(nameof(IJumpableAlgorithm) + " operation not supported for " + AlgorithmName);
            default:
                _jumpable.Jump();
                ZeroCache();
                break;
        }
    }

    /// <summary>
    /// Implements <see cref="IRandomAlgorithm.Next"/>.
    /// </summary>
    public ulong Next()
    {
        return _algorithm.Next();
    }

    /// <summary>
    /// Implements <see cref="IRandomRoutines.Next32"/>.
    /// </summary>
    public uint Next32()
    {
        switch (_next32Flag)
        {
            case true:
                // Caching improves the fastest
                // 64-bit generators by around 20%
                _next32Flag = false;
                return _next32Cache;
            default:
                switch (MaxNext)
                {
                    case ulong.MaxValue:
                    {
                        ulong r64 = _algorithm.Next();
                        _next32Flag = true;
                        _next32Cache = (uint)(r64 >> 32);
                        return (uint)r64;
                    }
                    case uint.MaxValue:
                        return (uint)_algorithm.Next();
                    default:
                        // Build unbiased 64-bit for non 32/64
                        // bit types. Relatively slow.
                        return (uint)GetUnbiasedBits(32);
                }
        }
    }

    /// <summary>
    /// Implements <see cref="IRandomRoutines.Next64"/>.
    /// </summary>
    public ulong Next64()
    {
        return MaxNext switch
        {
            ulong.MaxValue => _algorithm.Next(),
            uint.MaxValue => (_algorithm.Next() << 32) | _algorithm.Next(),
            _ => GetUnbiasedBits(64)
        };
    }

    /// <summary>
    /// Implements <see cref="IRandomRoutines.NextBytes(int)"/>.
    /// </summary>
    public byte[] NextBytes(int count)
    {
        byte[] rslt = new byte[count];
        NextBytes(rslt, 0, -1);
        return rslt;
    }

    /// <summary>
    /// Implements <see cref="IRandomRoutines.NextBytes(byte[])"/>.
    /// </summary>
    public void NextBytes(byte[] dest)
    {
        NextBytes(dest, 0, -1);
    }

    /// <summary>
    /// Implements <see cref="IRandomRoutines.NextBytes(byte[], int)"/>.
    /// </summary>
    public int NextBytes(byte[] dest, int offset)
    {
        return NextBytes(dest, offset, -1);
    }

    /// <summary>
    /// Implements <see cref="IRandomRoutines.NextBytes(byte[], int, int)"/>.
    /// </summary>
    public int NextBytes(byte[] dest, int offset, int count)
    {
        // Buffer size in ulongs (i.e. multiples of 8 bytes).
        const int MaxBufferSize = 512;

        // This will throw on bounds error
        count = SeedHelper.AssertBounds(dest.Length, offset, count);

        int i8Count = count;
        int i64Count = i8Count / 8;
        if (i8Count % 8 != 0)
        {
            ++i64Count;
        }

        int bufSize = Math.Min(i64Count, MaxBufferSize);
        ulong[] buf = new ulong[bufSize];

        int byteSize = bufSize * 8;
        if (byteSize > i8Count)
        {
            byteSize = i8Count;
        }

        while (true)
        {
            int bufIdx = 0;
            do { buf[bufIdx] = Next64(); } while (++bufIdx < bufSize);

            Buffer.BlockCopy(buf, 0, dest, offset, byteSize);

            if ((i8Count -= byteSize) == 0)
            {
                return count;
            }

            offset += byteSize;

            if ((i64Count -= bufSize) >= bufSize)
            {
                continue;
            }

            // Short last block
            bufSize = i64Count;
            byteSize = i8Count;
        }
    }

    /// <summary>
    /// Implements <see cref="IRandomRoutines.NextDouble()"/>.
    /// </summary>
    public double NextDouble()
    {
        // Multiplier equals 0x1.0p-53, or as 64 bits: 0x3CA0000000000000
        return (Next64() >> 11) * 1.1102230246251565E-16;
    }

    /// <summary>
    /// Implements <see cref="IRandomRoutines.NextDouble(double, double)"/>.
    /// </summary>
    public double NextDouble(double min, double max)
    {
        if (max > min)
        {
            return min + NextDouble() * (max - min);
        }

        throw new ArgumentOutOfRangeException(nameof(max), "Invalid range");
    }

    /// <summary>
    /// Implements <see cref="IRandomRoutines.NextFlip"/>.
    /// </summary>
    public bool NextFlip()
    {
        switch (_flipsRemain)
        {
            case 0:
                _flipCache = Next32();
                _flipsRemain = 32;
                break;
        }

        bool rslt = (_flipCache & 0x01UL) == 1U;
        _flipCache >>= 1;
        _flipsRemain -= 1;
        return rslt;
    }

    /// <summary>
    /// Implements <see cref="IRandomRoutines.NextInt(int)"/>.
    /// </summary>
    public int NextInt(int max)
    {
        return max switch
        {
            > 0 => (int) (((ulong) max * Next32()) >> 32),
            _ => throw new ArgumentOutOfRangeException(nameof(max))
        };
    }

    /// <summary>
    /// Implements <see cref="IRandomRoutines.NextInt(int, int)"/>.
    /// </summary>
    public int NextInt(int min, int max)
    {
        if (max > min)
        {
            return min + (int)(((ulong)Next32() * (uint)(max - min)) >> 32);
        }

        throw new ArgumentOutOfRangeException(nameof(max), "Invalid range");
    }

    /// <summary>
    /// Implements <see cref="IRandomRoutines.NextIntExclude(int, int, int)"/>.
    /// </summary>
    public int NextIntExclude(int min, int max, int exclude)
    {
        if (max <= min + 1)
        {
            throw new ArgumentOutOfRangeException(nameof(max), "Invalid range");
        }

        int rslt = min + (int)(((ulong)Next32() * (uint)(max - min)) >> 32);

        while (rslt == exclude)
        {
            rslt = min + (int)(((ulong)Next32() * (uint)(max - min)) >> 32);
        }

        return rslt;

    }

    /// <summary>
    /// Implements <see cref="IRandomRoutines.NextLong(long)"/>.
    /// </summary>
    public long NextLong(long max)
    {
        return max switch
        {
            > 0 => (long) NextLongInternal((ulong) max),
            _ => throw new ArgumentOutOfRangeException(nameof(max))
        };
    }

    /// <summary>
    /// Implements <see cref="IRandomRoutines.NextLong(long, long)"/>.
    /// </summary>
    public long NextLong(long min, long max)
    {
        if (max > min)
        {
            return min + (long)NextLongInternal((ulong)(max - min));
        }

        throw new ArgumentException("Invalid range");
    }

    /// <summary>
    /// Implements <see cref="IRandomRoutines.NextLongExclude(long, long, long)"/>.
    /// </summary>
    public long NextLongExclude(long min, long max, long exclude)
    {
        if (max <= min + 1)
        {
            throw new ArgumentOutOfRangeException(nameof(max), "Invalid range");
        }

        long rslt = min + (long)NextLongInternal((ulong)(max - min));

        while (rslt == exclude)
        {
            rslt = min + (long)NextLongInternal((ulong)(max - min));
        }

        return rslt;

    }

    /// <summary>
    /// Implements <see cref="IRandomRoutines.NextOpenDouble"/>.
    /// </summary>
    public double NextOpenDouble()
    {
        double rslt = NextDouble();

        while (rslt == 0)
        {
            rslt = NextDouble();
        }

        return rslt;
    }

    /// <summary>
    /// Implements <see cref="IRandomRoutines.NextStdNormal"/>.
    /// </summary>
    public double NextStdNormal()
    {
        switch (_gaussFlag)
        {
            case true:
                // Two values are generated, with
                // one being cached for the next call.
                _gaussFlag = false;
                return _gaussCache;
        }

        // Polar form of Box-Muller.
        // Ref: http://www.design.caltech.edu/erik/Misc/Gaussian.html
        double r0, r1, w;

        do
        {
            r0 = 2.0 * NextDouble() - 1.0;
            r1 = 2.0 * NextDouble() - 1.0;
            w = r0 * r0 + r1 * r1;

        } while (w >= 1 || w <= double.Epsilon);

        w = Math.Sqrt(-2.0 * Math.Log(w) / w);
        _gaussCache = r0 * w;
        _gaussFlag = true;

        return r1 * w;
    }

    /// <summary>
    /// Implements <see cref="IRandomRoutines.NextString(int, string)"/>.
    /// </summary>
    public string NextString(int count, string alphabet)
    {
        StringBuilder sb = new(count);

        for (int n = 0; n < count; ++n)
        {
            sb.Append(alphabet[NextInt(0, alphabet.Length)]);
        }

        return sb.ToString();
    }

    /// <summary>
    /// Implements <see cref="IRandomGenerator.Randomize()"/>.
    /// </summary>
    public void Randomize()
    {
        switch (_seedable)
        {
            case null:
                throw new NotSupportedException(nameof(ISeedableAlgorithm) + " operation not supported for " + AlgorithmName);
            default:
                SetSeed(SeedHelper.GenerateSeed(SeedLength));
                break;
        }
    }

    /// <summary>
    /// Implements <see cref="IRandomGenerator.SetSeed(ulong)"/>.
    /// </summary>
    public void SetSeed(ulong seed)
    {
        switch (_seedable)
        {
            case null:
                throw new NotSupportedException(nameof(ISeedableAlgorithm) + " operation not supported for " + AlgorithmName);
            default:
                _seedable.SetSeed(SeedHelper.ExtendSeed(seed, SeedLength));
                ZeroCache();
                break;
        }
    }

    /// <summary>
    /// Implements <see cref="IRandomGenerator.SetSeed(long)"/>.
    /// </summary>
    public void SetSeed(long seed)
    {
        SetSeed((ulong)seed);
    }

    /// <summary>
    /// Implements <see cref="ISeedableAlgorithm.SetSeed(byte[])"/>. The "seed" array must be of
    /// length <see cref="ISeedableAlgorithm.SeedLength"/> or greater (whether additional bytes
    /// are ignored or used is algorithm dependent). The seed is supplied directly to the
    /// underlying <see cref="ISeedableAlgorithm"/> instance. The same byte sequence shall yield
    /// the same internal state irrespective of endian.
    /// </summary>
    /// <exception cref="NotSupportedException"><see cref="IsSeedable"/> is false</exception>
    /// <exception cref="ArgumentException">Array length too short</exception>
    public void SetSeed(byte[] seed)
    {
        switch (_seedable)
        {
            case null:
                throw new NotSupportedException(nameof(ISeedableAlgorithm) + " operation not supported for " + AlgorithmName);
        }

        if (seed.Length < SeedLength)
        {
            throw new ArgumentException("Array length too short");
        }

        _seedable.SetSeed(seed);
        ZeroCache();
    }

    /// <summary>
    /// Implements <see cref="IRandomRoutines.Shuffle{T}(T[])"/>.
    /// </summary>
    public void Shuffle<T>(T[] items)
    {
        Shuffle(items, 0, -1);
    }

    /// <summary>
    /// Implements <see cref="IRandomRoutines.Shuffle{T}(IList{T})"/>.
    /// </summary>
    public void Shuffle<T>(IList<T> items)
    {
        Shuffle(items, 0, -1);
    }

    /// <summary>
    /// Implements <see cref="IRandomRoutines.Shuffle{T}(T[], int, int)"/>.
    /// </summary>
    public void Shuffle<T>(T[] items, int offset, int count)
    {
        // Fisher-Yates algorithm
        count = SeedHelper.AssertBounds(items.Length, offset, count);

        int end = count + offset;

        while (count > 0)
        {
            int j = offset + NextInt(count--);

            T temp = items[j];
            items[j] = items[--end];
            items[end] = temp;
        }
    }

    /// <summary>
    /// Implements <see cref="IRandomRoutines.Shuffle{T}(T[], int, int)"/>.
    /// </summary>
    public void Shuffle<T>(IList<T> items, int offset, int count)
    {
        // Fisher-Yates algorithm
        count = SeedHelper.AssertBounds(items.Count, offset, count);

        int end = count + offset;

        while (count > 0)
        {
            int j = offset + NextInt(count--);

            T temp = items[j];
            items[j] = items[--end];
            items[end] = temp;
        }
    }

    /// <summary>
    /// Implements <see cref="IRandomRoutines.Shuffle(string)"/>.
    /// </summary>
    public string Shuffle(string s)
    {
        switch (s.Length)
        {
            case <= 1:
                return s;
        }

        char[] arr = s.ToCharArray();
        Shuffle(arr);
        return new string(arr);

    }

    /// <summary>
    /// Zero (reset) internal cache values belonging to <see cref="RandomGenerator"/>. It does
    /// not modify the state of the internal algorithm.
    /// </summary>
    protected void ZeroCache()
    {
        _next32Flag = false;
        _next32Cache = 0;
        _flipCache = 0;
        _flipsRemain = 0;
        _gaussFlag = false;
        _gaussCache = 0;
    }

    private static int BitSize(ulong x)
    {
        // Calculate the number of bits needed
        // for x. I.e, for 0x70, it would be 7.

        switch (x)
        {
            // Shortcuts to common values.
            case int.MaxValue:
                return 31;
            case long.MaxValue:
                return 63;
        }

        switch (x)
        {
            case 0:
                return 0;
        }

        int bits = 64;
        while (x >> --bits == 0)
        {
        }

        return bits + 1;

    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private ulong NextLongInternal(ulong max)
    {
        switch (max)
        {
            case <= uint.MaxValue:
                // Typically faster
                return (max * Next32()) >> 32;
        }

        ulong rxh = max >> 32;
        ulong rxl = max & uint.MaxValue;

        ulong rav = Next64();
        ulong ravh = rav >> 32;
        ulong ravl = (uint)rav;

        return ((rxl * ravh) >> 32) + ((rxh * ravl) >> 32) + rxh * ravh;
    }

    private ulong GetUnbiasedBits(int bits)
    {
        // Build unbiased integer.
        ulong rslt = 0;
        int sb = _shiftBits;
        ulong sm = _shiftMax;

        do
        {
            ulong x;
            do { x = Next(); } while (x > sm);

            rslt <<= sb;
            rslt |= x;

        } while ((bits -= sb) > 0);

        return rslt;
    }

}