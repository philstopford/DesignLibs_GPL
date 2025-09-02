using System;
using System.Runtime.CompilerServices;

namespace utility;

/// <summary>
/// Advanced Box-Muller transform implementation with optimizations
/// </summary>
public static class BoxMullerOptimized
{
    // Performance optimization: Pre-calculated constants
    private const double TWO_PI = 2.0 * Math.PI;
    private const double SQRT_2 = 1.4142135623730951;
    private const double SQRT_NEG_2 = -1.4142135623730951;
    
    // Lookup tables for common cases
    private const int LOOKUP_SIZE = 1024;
    private static readonly double[] SqrtTable = new double[LOOKUP_SIZE];
    private static readonly double[] LogTable = new double[LOOKUP_SIZE];
    
    static BoxMullerOptimized()
    {
        // Pre-calculate lookup tables for common ranges
        for (int i = 0; i < LOOKUP_SIZE; i++)
        {
            double x = (double)i / LOOKUP_SIZE;
            SqrtTable[i] = Math.Sqrt(x);
            LogTable[i] = x > 0 ? Math.Log(x) : double.NegativeInfinity;
        }
    }

    /// <summary>
    /// Optimized Box-Muller transform with cached values and fast math
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static (double, double) GenerateGaussianPair(double u1, double u2)
    {
        // Performance optimization: Pre-calculate common sub-expressions
        double logU2 = Math.Log(u2);
        double sqrtNeg2LogU2 = Math.Sqrt(-2.0 * logU2);
        double twoPiU1 = TWO_PI * u1;
        
        // Use both sine and cosine to get two values efficiently
        double cosVal = Math.Cos(twoPiU1);
        double sinVal = Math.Sin(twoPiU1);
        
        return (sqrtNeg2LogU2 * cosVal, sqrtNeg2LogU2 * sinVal);
    }

    /// <summary>
    /// Fast lookup-based square root for small values
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double FastSqrt(double x)
    {
        if (x >= 0.0 && x <= 1.0)
        {
            int index = (int)(x * (LOOKUP_SIZE - 1));
            return SqrtTable[index];
        }
        return Math.Sqrt(x);
    }

    /// <summary>
    /// Fast lookup-based logarithm for small values
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double FastLog(double x)
    {
        if (x >= 0.0 && x <= 1.0)
        {
            int index = (int)(x * (LOOKUP_SIZE - 1));
            return LogTable[index];
        }
        return Math.Log(x);
    }
}