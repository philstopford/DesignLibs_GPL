using System;
using System.Runtime.CompilerServices;
using System.Security.Cryptography;
using System.Threading;

namespace entropyRNG;

public static class RNG
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

    [ThreadStatic] private static Random? _local;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double[] random_gauss()
    {
        Random random = GetThreadLocalRandom();

        // Box-Muller transform
        // We aren't allowed 0, so we reject any values approaching zero.
        double U1 = random.NextDouble();
        while (U1 < 1E-15)
        {
            U1 = random.NextDouble();
        }
        double U2 = random.NextDouble();
        while (U2 < 1E-15)
        {
            U2 = random.NextDouble();
        }
        // PAs are 3-sigma, so this needs to be divided by 3 to give single sigma value when used
        double A1 = Math.Sqrt(-2 * Math.Log(U2)) * Math.Cos(2 * Math.PI * U1);
        double A2 = Math.Sqrt(-2 * Math.Log(U1)) * Math.Sin(2 * Math.PI * U2);
        return [A1, A2];
    }
    
    // This is our Gaussian RNG
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double[] random_gauss3()
    {
        double[] myReturn = random_gauss();
        return [myReturn[0] / 3.0, myReturn[1] / 3.0]; // Use double precision
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double nextdouble(int seed)
    {
        var buffer = _threadLocalBuffer.Value!;
        _global.GetBytes(buffer);
        var random = new Random(seed);
        return random.NextDouble();
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double nextdouble()
    {
        Random random = GetThreadLocalRandom();
        return random.NextDouble();
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int nextint()
    {
        Random random = GetThreadLocalRandom();
        return nextint(random);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static int nextint(Random random)
    {
        return random.Next();
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int nextint(int seed)
    {
        var buffer = _threadLocalBuffer.Value!;
        _global.GetBytes(buffer);
        var random = new Random(seed);
        return random.Next();
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int nextint(int min, int max, int seed)
    {
        var buffer = _threadLocalBuffer.Value!;
        _global.GetBytes(buffer);
        var random = new Random(seed);
        return nextint(random, min, max);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int nextint(int min, int max)
    {
        Random random = GetThreadLocalRandom();
        return nextint(random, min, max);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static int nextint(Random random, int min, int max)
    {
        return random.Next(min, max);
    }

    // This is a slightly different version of our Gaussian RNG
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double random_gauss2()
    {
        Random random = GetThreadLocalRandom();
        // We aren't allowed 0, so we reject any values approaching zero.
        double U1 = random.NextDouble();
        while (U1 < 1E-15)
        {
            U1 = random.NextDouble();
        }
        double U2 = random.NextDouble();
        while (U2 < 1E-15)
        {
            U2 = random.NextDouble();
        }
        // PAs are 3-sigma, so this needs to be divided by 3 to give single sigma value when used
        double A2 = Math.Sqrt(-2 * Math.Log(U1)) * Math.Sin(2 * Math.PI * U2) / 3;
        return A2;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static Random GetThreadLocalRandom()
    {
        if (_local == null)
        {
            var buffer = _threadLocalBuffer.Value!;
            _global.GetBytes(buffer);
            _local = new Random(BitConverter.ToInt32(buffer, 0));
        }
        return _local;
    }
}