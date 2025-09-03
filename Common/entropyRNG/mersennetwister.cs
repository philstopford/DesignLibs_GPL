#nullable enable
using MersenneTwisterRNG;
using System;
using System.Runtime.CompilerServices;
using System.Security.Cryptography;
using System.Threading;

namespace entropyRNG;

public static class MersenneTwister_RNG
{
    /*
     * This class is interesting. Originally, the intent was to have per-thread RNGs, but it became apparent that threads that instantiated an RNG
     * would get the same random number distribution when the RNGs were initialized at the same system time.
     * To avoid this, earlier systems made a common RNG and the threads would query from that RNG.
     * However, this was not thread-safe, such that the calls to the RNG would start returning 0 and the application enters a spiral of death as the RNG was 
     * continuously, and unsuccessfully, polled for a non-zero value.
     * 
     * Locking the RNG was one option, so that only one thread could query at a time, but this caused severe performance issues.
     * 
     * So, to address this, I've gone back to a per-thread RNG (referenced in jobSettings()) and then use the RNGCryptoServiceProvider to provide a 
     * 'seed' value for a null RNG entity. This avoids some severe performance issues if the RNGCryptoServiceProvider is used for all random numbers.
     * 
     * Ref : http://blogs.msdn.com/b/pfxteam/archive/2009/02/19/9434171.aspx
     */
    private static readonly RandomNumberGenerator _global = RandomNumberGenerator.Create();
    private static readonly ThreadLocal<byte[]> _threadLocalBuffer = new(() => new byte[4]);

    // Performance optimization: Cache mathematical constants
    private const double TWO_PI = 2.0 * Math.PI;
    private const double MIN_RANDOM_VALUE = 1E-15;
    private const double ONE_THIRD = 1.0 / 3.0;

    [ThreadStatic] private static MersenneTwister? _local;
    
    // Performance optimization: Box-Muller cached value for better efficiency
    [ThreadStatic] private static double _cachedGaussian;
    [ThreadStatic] private static bool _hasValidCachedGaussian;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double[] random_gauss3()
    {
        double[] myReturn = random_gauss();
        // Performance optimization: Use pre-calculated constant
        return [myReturn[0] * ONE_THIRD, myReturn[1] * ONE_THIRD];
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double[] random_gauss()
    {
        // Performance optimization: Use cached Box-Muller value if available
        if (_hasValidCachedGaussian)
        {
            _hasValidCachedGaussian = false;
            return [_cachedGaussian, GenerateGaussianPair()[1]];
        }

        var pair = GenerateGaussianPair();
        _cachedGaussian = pair[1];
        _hasValidCachedGaussian = true;
        return [pair[0], pair[1]];
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double[] GenerateGaussianPair()
    {
        MersenneTwister random = GetThreadLocalRandom();

        // Box-Muller transform - optimized version
        // We aren't allowed 0, so we reject any values approaching zero.
        double U1, U2;
        
        // Performance optimization: Use do-while instead of while for better branch prediction
        do
        {
            U1 = random.NextDouble();
        } while (U1 < MIN_RANDOM_VALUE);
        
        do
        {
            U2 = random.NextDouble();
        } while (U2 < MIN_RANDOM_VALUE);

        // Performance optimization: Pre-calculate common sub-expressions
        double logU2 = Math.Log(U2);
        double sqrtNeg2LogU2 = Math.Sqrt(-2.0 * logU2);
        double twoPiU1 = TWO_PI * U1;
        
        // PAs are 3-sigma, so this needs to be divided by 3 to give single sigma value when used
        double A1 = sqrtNeg2LogU2 * Math.Cos(twoPiU1);
        double A2 = sqrtNeg2LogU2 * Math.Sin(twoPiU1);
        return [A1, A2];
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double nextdouble(int seed)
    {
        var buffer = _threadLocalBuffer.Value!;
        _global.GetBytes(buffer);
        var random = new MersenneTwister(seed);
        return random.NextDouble();
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double nextdouble()
    {
        MersenneTwister random = GetThreadLocalRandom();
        return random.NextDouble();
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int nextint()
    {
        MersenneTwister random = GetThreadLocalRandom();
        return nextint(random);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static int nextint(MersenneTwister random)
    {
        return random.Next();
    }
        
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int nextint(int min, int max, int seed)
    {
        var buffer = _threadLocalBuffer.Value!;
        _global.GetBytes(buffer);
        var random = new MersenneTwister(seed);
        return nextint(random, min, max);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int nextint(int min, int max)
    {
        MersenneTwister random = GetThreadLocalRandom();
        return nextint(random, min, max);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static int nextint(MersenneTwister random, int min, int max)
    {
        return random.Next(min, max);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static MersenneTwister GetThreadLocalRandom()
    {
        if (_local == null)
        {
            var buffer = _threadLocalBuffer.Value!;
            _global.GetBytes(buffer);
            _local = new MersenneTwister(BitConverter.ToInt32(buffer, 0));
        }
        return _local;
    }
}