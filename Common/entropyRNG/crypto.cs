using System;
using System.Runtime.CompilerServices;
using System.Security.Cryptography;
using System.Threading;

namespace entropyRNG;

public static class Crypto_RNG
{
    private const double maxInt = int.MaxValue; // max int value.
    
    // Thread-local storage for optimal performance in multi-threaded scenarios
    private static readonly ThreadLocal<RandomNumberGenerator> _threadLocalRng = 
        new(() => RandomNumberGenerator.Create());
    
    private static readonly ThreadLocal<byte[]> _threadLocalBuffer = 
        new(() => new byte[sizeof(long)]); // Use long buffer for efficiency

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double[] random_gauss3()
    {
        double[] myReturn = random_gauss();
        return [myReturn[0] / 3.0, myReturn[1] / 3.0]; // Use double precision
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double[] random_gauss()
    {
        var random = _threadLocalRng.Value!;
        var tmp = _threadLocalBuffer.Value!;

        // Box-Muller transform
        // We aren't allowed 0, so we reject any values approaching zero.
        random.GetBytes(tmp.AsSpan(0, sizeof(int)));
        double U1 = Math.Abs(BitConverter.ToInt32(tmp, 0)) / maxInt;
        while (U1 < 1E-15)
        {
            random.GetBytes(tmp.AsSpan(0, sizeof(int)));
            U1 = Math.Abs(BitConverter.ToInt32(tmp, 0)) / maxInt;
        }
        random.GetBytes(tmp.AsSpan(0, sizeof(int)));
        double U2 = Math.Abs(BitConverter.ToInt32(tmp, 0)) / maxInt;
        while (U2 < 1E-15)
        {
            random.GetBytes(tmp.AsSpan(0, sizeof(int)));
            U2 = Math.Abs(BitConverter.ToInt32(tmp, 0)) / maxInt;
        }
        // PAs are 3-sigma, so this needs to be divided by 3 to give single sigma value when used
        double A1 = Math.Sqrt(-2 * Math.Log(U2)) * Math.Cos(2 * Math.PI * U1);
        double A2 = Math.Sqrt(-2 * Math.Log(U1)) * Math.Sin(2 * Math.PI * U2);
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