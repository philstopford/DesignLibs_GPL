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
/// ISAAC-32 is an all purpose cryptographic pseudo random number generator with good performance.
/// On construction, the initial state is populated with a fixed arbitrary value.
/// </summary>
/// <remarks>
/// The ISAAC algorithm was created by Bob Jenkins and is in the public domain. For more
/// information, see: http://burtleburtle.net/bob/rand/isaacafa.html
/// </remarks>
public class Isaac32 : ISeedableAlgorithm
{
    private const int StateSize = 256;
    private const uint StateMask = 0xFF;
    private const uint GoldenRatio = 0x9E3779B9U;

    private int _stateIdx;
    private uint _aa, _bb, _cc;
    private readonly uint[] _state = new uint[StateSize];
    private readonly uint[] _ready = new uint[StateSize];

    /// <summary>
    /// Implements <see cref="IRandomAlgorithm.AlgorithmName"/>.
    /// </summary>
    public string AlgorithmName => "ISAAC-32";

    /// <summary>
    /// Implements <see cref="IRandomAlgorithm.MaxNext"/>.
    /// </summary>
    public ulong MaxNext => uint.MaxValue;

    /// <summary>
    /// Implements <see cref="ISeedableAlgorithm.SeedLength"/>.
    /// </summary>
    public int SeedLength => 1024;

    /// <summary>
    /// Default constructor.
    /// </summary>
    public Isaac32()
    {
        Initialize();
    }

    /// <summary>
    /// Copy constructor.
    /// </summary>
    private Isaac32(Isaac32 other)
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
        Buffer.BlockCopy(SeedHelper.ToUInt32Array(seed, SeedLength), 0, _ready, 0, SeedLength);
        Initialize();
    }

    /// <summary>
    /// Implements <see cref="ISeedableAlgorithm.Clone"/>.
    /// </summary>
    public ISeedableAlgorithm Clone()
    {
        return new Isaac32(this);
    }

    private void Initialize()
    {
        _stateIdx = 0;

        _aa = 0;
        _bb = 0;
        _cc = 0;

        uint a = GoldenRatio;
        uint b = a, c = a, d = a, e = a, f = a, g = a, h = a;

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
        uint a = _aa;
        uint b = _bb + ++_cc;

        int nrot = StateSize / 2 - 1;

        // Local copies
        uint[] sloc = _state;
        uint[] rloc = _ready;

        for (int n = 0; n < StateSize; ++n)
        {
            uint x = sloc[n];
            int idx = n & 0x03;

            a ^= idx switch
            {
                0 => a << 13,
                1 => a >> 6,
                2 => a << 2,
                _ => a >> 16
            };

            a += sloc[++nrot & StateMask];

            uint y;
            sloc[n] = y = sloc[(x >> 2) & StateMask] + a + b;
            rloc[n] = b = sloc[(y >> 10) & StateMask] + x;
        }

        _aa = a;
        _bb = b;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void  Mixer(ref uint a, ref uint b, ref uint c, ref uint d,
        ref uint e, ref uint f, ref uint g, ref uint h)
    {
        a ^= b << 11; d += a; b += c;
        b ^= c >> 2; e += b; c += d;
        c ^= d << 8; f += c; d += e;
        d ^= e >> 16; g += d; e += f;
        e ^= f << 10; h += e; f += g;
        f ^= g >> 4; a += f; g += h;
        g ^= h << 8; b += g; h += a;
        h ^= a >> 9; c += h; a += b;
    }

}