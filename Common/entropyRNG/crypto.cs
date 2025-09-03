using System;
using System.Runtime.CompilerServices;
using System.Security.Cryptography;
using System.Threading;

namespace entropyRNG;

public static class Crypto_RNG
{
    private const double maxInt = int.MaxValue; // max int value.
    
    // Performance optimization: Cache mathematical constants
    private const double TWO_PI = 2.0 * Math.PI;
    private const double MIN_RANDOM_VALUE = 1E-15;
    private const double ONE_THIRD = 1.0 / 3.0;
    
    // Thread-local storage for optimal performance in multi-threaded scenarios
    private static readonly ThreadLocal<RandomNumberGenerator> _threadLocalRng = 
        new(() => RandomNumberGenerator.Create());
    
    private static readonly ThreadLocal<byte[]> _threadLocalBuffer = 
        new(() => new byte[sizeof(long)]); // Use long buffer for efficiency

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
        var random = _threadLocalRng.Value!;
        var tmp = _threadLocalBuffer.Value!;

        // Box-Muller transform - optimized version
        // We aren't allowed 0, so we reject any values approaching zero.
        double U1, U2;
        
        // Performance optimization: Use do-while for better branch prediction
        do
        {
            random.GetBytes(tmp.AsSpan(0, sizeof(int)));
            U1 = Math.Abs(BitConverter.ToInt32(tmp, 0)) / maxInt;
        } while (U1 < MIN_RANDOM_VALUE);
        
        do
        {
            random.GetBytes(tmp.AsSpan(0, sizeof(int)));
            U2 = Math.Abs(BitConverter.ToInt32(tmp, 0)) / maxInt;
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
    public static double nextdouble()
    {
        var random = _threadLocalRng.Value!;
        var tmp = _threadLocalBuffer.Value!;

        random.GetBytes(tmp.AsSpan(0, sizeof(int)));
        return Math.Abs(BitConverter.ToInt32(tmp, 0)) / maxInt;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int nextint()
    {
        var random = _threadLocalRng.Value!;
        var tmp = _threadLocalBuffer.Value!;

        random.GetBytes(tmp.AsSpan(0, sizeof(int)));
        return BitConverter.ToInt32(tmp, 0);
    }
}