using System;
using System.Buffers;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;

namespace DesignLibs.Performance;

/// <summary>
/// Collection of performance optimization utilities and patterns
/// for use across the DesignLibs codebase.
/// </summary>
public static class PerformanceOptimizations
{
    /// <summary>
    /// Pool for double arrays to reduce GC pressure
    /// </summary>
    private static readonly ArrayPool<double> DoubleArrayPool = ArrayPool<double>.Shared;
    
    /// <summary>
    /// Pool for int arrays to reduce GC pressure
    /// </summary>
    private static readonly ArrayPool<int> IntArrayPool = ArrayPool<int>.Shared;

    /// <summary>
    /// Rent a double array from the pool. Remember to return it!
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double[] RentDoubleArray(int minimumLength)
    {
        return DoubleArrayPool.Rent(minimumLength);
    }

    /// <summary>
    /// Return a double array to the pool
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static void ReturnDoubleArray(double[] array, bool clearArray = false)
    {
        DoubleArrayPool.Return(array, clearArray);
    }

    /// <summary>
    /// Rent an int array from the pool. Remember to return it!
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int[] RentIntArray(int minimumLength)
    {
        return IntArrayPool.Rent(minimumLength);
    }

    /// <summary>
    /// Return an int array to the pool
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static void ReturnIntArray(int[] array, bool clearArray = false)
    {
        IntArrayPool.Return(array, clearArray);
    }

    /// <summary>
    /// Fast copy using Span&lt;T&gt; for better performance than Array.Copy
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static void FastCopy<T>(ReadOnlySpan<T> source, Span<T> destination)
    {
        source.CopyTo(destination);
    }

    /// <summary>
    /// Optimized distance calculation squared (avoids expensive sqrt)
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double DistanceSquared(double x1, double y1, double x2, double y2)
    {
        double dx = x2 - x1;
        double dy = y2 - y1;
        return dx * dx + dy * dy;
    }

    /// <summary>
    /// Optimized distance calculation squared for 3D points
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double DistanceSquared3D(double x1, double y1, double z1, double x2, double y2, double z2)
    {
        double dx = x2 - x1;
        double dy = y2 - y1;
        double dz = z2 - z1;
        return dx * dx + dy * dy + dz * dz;
    }

    /// <summary>
    /// Fast integer power for small exponents
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double FastPow(double value, int exponent)
    {
        return exponent switch
        {
            0 => 1.0,
            1 => value,
            2 => value * value,
            3 => value * value * value,
            4 => value * value * value * value,
            _ => Math.Pow(value, exponent)
        };
    }

    /// <summary>
    /// Optimized lerp that avoids multiplication when t is 0 or 1
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double FastLerp(double t, double a, double b)
    {
        return t switch
        {
            0.0 => a,
            1.0 => b,
            _ => a + t * (b - a)
        };
    }

    /// <summary>
    /// Optimized smoothstep function for interpolation
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double SmoothStep(double t)
    {
        // Optimized: t * t * (3 - 2 * t)
        return t * t * (3.0 - 2.0 * t);
    }

    /// <summary>
    /// Check if a value is within epsilon of zero (optimized comparison)
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static bool IsNearZero(double value, double epsilon = 1e-10)
    {
        return Math.Abs(value) < epsilon;
    }
}

/// <summary>
/// Performance-focused disposable wrapper for pooled arrays
/// </summary>
public readonly ref struct PooledArray<T>
{
    private readonly ArrayPool<T> _pool;
    private readonly T[] _array;
    private readonly bool _clearOnReturn;

    public PooledArray(ArrayPool<T> pool, int minimumLength, bool clearOnReturn = false)
    {
        _pool = pool;
        _array = pool.Rent(minimumLength);
        _clearOnReturn = clearOnReturn;
    }

    public Span<T> Span => _array.AsSpan();
    public T[] Array => _array;

    public void Dispose()
    {
        _pool.Return(_array, _clearOnReturn);
    }
}