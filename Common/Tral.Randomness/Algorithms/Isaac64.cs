//--------------------------------------------------------------------------------------------------
// PROJECT      : TRAL
// COPYRIGHT    : Andy Thomas / https://kuiper.zone
// LICENSE      : Apache License, Version 2.0
//--------------------------------------------------------------------------------------------------

using System;
using System.Runtime.CompilerServices;
using Tral.Randomness.Utility;

namespace Tral.Randomness.Algorithms;

/// <summary>
/// ISAAC-64 is an all purpose cryptographic pseudo random number generator with good
/// performance. It is 64-bit variant of ISAAC-32, but gives different output. On construction,
/// the initial state is populated with a fixed arbitrary value.
/// </summary>
/// <remarks>
/// The ISAAC algorithm was created by Bob Jenkins and is in the public domain. For more
/// information, see: http://burtleburtle.net/bob/rand/isaacafa.html
/// </remarks>
public class Isaac64 : ISeedableAlgorithm
{
    private const int StateSize = 256;
    private const uint StateMask = 0xFF;
    private const ulong GoldenRatio = 0x9E3779B97F4A7C13UL;

    private int _stateIdx;
    private ulong _aa, _bb, _cc;
    private readonly ulong[] _state = new ulong[StateSize];
    private readonly ulong[] _ready = new ulong[StateSize];

    /// <summary>
    /// Implements <see cref="IRandomAlgorithm.AlgorithmName"/>.
    /// </summary>
    public string AlgorithmName => "ISAAC-64";

    /// <summary>
    /// Implements <see cref="IRandomAlgorithm.MaxNext"/>.
    /// </summary>
    public ulong MaxNext => ulong.MaxValue;

    /// <summary>
    /// Implements <see cref="ISeedableAlgorithm.SeedLength"/>.
    /// </summary>
    public int SeedLength => 2048;

    /// <summary>
    /// Default constructor.
    /// </summary>
    public Isaac64()
    {
        Initialize();
    }

    /// <summary>
    /// Copy constructor.
    /// </summary>
    private Isaac64(Isaac64 other)
    {
        _stateIdx = other._stateIdx;

        _aa = other._aa;
        _bb = other._bb;
        _cc = other._cc;

        Buffer.BlockCopy(other._state, 0, _state, 0, SeedLength);
        Buffer.BlockCopy(other._ready, 0, _ready, 0, SeedLength);
    }

    /// <summary>
    /// Implements <see cref="IRandomAlgorithm.Next"/>.
    /// </summary>
    public ulong Next()
    {
        if (_stateIdx != 0)
        {
            return _ready[--_stateIdx];
        }

        Reload();
        _stateIdx = StateSize;

        return _ready[--_stateIdx];
    }

    /// <summary>
    /// Implements <see cref="ISeedableAlgorithm.SetSeed(byte[])"/>.
    /// </summary>
    public void SetSeed(byte[] seed)
    {
        Buffer.BlockCopy(SeedHelper.ToUInt64Array(seed, SeedLength), 0, _ready, 0, SeedLength);
        Initialize();
    }

    /// <summary>
    /// Implements <see cref="ISeedableAlgorithm.Clone"/>.
    /// </summary>
    public ISeedableAlgorithm Clone()
    {
        return new Isaac64(this);
    }

    private void Initialize()
    {
        _stateIdx = 0;

        _aa = 0;
        _bb = 0;
        _cc = 0;

        ulong a = GoldenRatio;
        ulong b = a, c = a, d = a, e = a, f = a, g = a, h = a;

        // Scramble
        for (int n = 0; n < 4; ++n)
        {
            Mixer(ref a, ref b, ref c, ref d, ref e, ref f, ref g, ref h);
        }

        for (int n = 0; n < StateSize; n += 8)
        {
            // Populated with seed
            a += _ready[n];
            b += _ready[n + 1];
            c += _ready[n + 2];
            d += _ready[n + 3];
            e += _ready[n + 4];
            f += _ready[n + 5];
            g += _ready[n + 6];
            h += _ready[n + 7];

            Mixer(ref a, ref b, ref c, ref d, ref e, ref f, ref g, ref h);

            // Initialize state
            _state[n] = a;
            _state[n + 1] = b;
            _state[n + 2] = c;
            _state[n + 3] = d;
            _state[n + 4] = e;
            _state[n + 5] = f;
            _state[n + 6] = g;
            _state[n + 7] = h;
        }

        // Second pass
        for (int n = 0; n < StateSize; n += 8)
        {
            a += _state[n];
            b += _state[n + 1];
            c += _state[n + 2];
            d += _state[n + 3];
            e += _state[n + 4];
            f += _state[n + 5];
            g += _state[n + 6];
            h += _state[n + 7];

            Mixer(ref a, ref b, ref c, ref d, ref e, ref f, ref g, ref h);

            _state[n] = a;
            _state[n + 1] = b;
            _state[n + 2] = c;
            _state[n + 3] = d;
            _state[n + 4] = e;
            _state[n + 5] = f;
            _state[n + 6] = g;
            _state[n + 7] = h;
        }
    }

    private void Reload()
    {
        ulong a = _aa;
        ulong b = _bb + ++_cc;

        int nrot = StateSize / 2 - 1;

        // Local copies
        ulong[] sloc = _state;
        ulong[] rloc = _ready;

        for (int n = 0; n < StateSize; ++n)
        {
            ulong x = sloc[n];
            int idx = n & 0x03;

            switch (idx)
            {
                case 0:
                    a = ~(a ^ (a << 21));
                    break;
                case 1:
                    a ^= a >> 5;
                    break;
                case 2:
                    a ^= a << 12;
                    break;
                default:
                    a ^= a >> 33;
                    break;
            }

            a += sloc[++nrot & StateMask];

            ulong y;
            sloc[n] = y = sloc[(x >> 3) & StateMask] + a + b;
            rloc[n] = b = sloc[(y >> 11) & StateMask] + x;
        }

        _aa = a;
        _bb = b;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void Mixer(ref ulong a, ref ulong b, ref ulong c, ref ulong d,
        ref ulong e, ref ulong f, ref ulong g, ref ulong h)
    {
        a -= e; f ^= h >> 9; h += a;
        b -= f; g ^= a << 9; a += b;
        c -= g; h ^= b >> 23; b += c;
        d -= h; a ^= c << 15; c += d;
        e -= a; b ^= d >> 14; d += e;
        f -= b; c ^= e << 20; e += f;
        g -= c; d ^= f >> 17; f += g;
        h -= d; e ^= g << 14; g += h;
    }

}