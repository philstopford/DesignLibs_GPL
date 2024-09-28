using MersenneTwisterRNG;
using System;
using System.Security.Cryptography;

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

    [ThreadStatic] private static MersenneTwister _local;

    public static double[] random_gauss3()
    {
        double[] myReturn = random_gauss();
        return [myReturn[0] / 3.0f, myReturn[1] / 3.0f];
    }

    public static double[] random_gauss()
    {
        MersenneTwister random = _local;

        switch (random)
        {
            case null:
            {
                byte[] buffer = new byte[4];
                _global.GetBytes(buffer);
                _local = random = new MersenneTwister(BitConverter.ToInt32(buffer, 0));
                break;
            }
        }

        // Box-Muller transform
        // We aren't allowed 0, so we reject any values approaching zero.
        var U1 = random.NextDouble();
        while (U1 < 1E-15)
        {
            U1 = random.NextDouble();
        }
        var U2 = random.NextDouble();
        while (U2 < 1E-15)
        {
            U2 = random.NextDouble();
        }
        // PAs are 3-sigma, so this needs to be divided by 3 to give single sigma value when used
        double A1 = Math.Sqrt(-2 * Math.Log(U2, Math.E)) * Math.Cos(2 * Math.PI * U1);
        double A2 = Math.Sqrt(-2 * Math.Log(U1, Math.E)) * Math.Sin(2 * Math.PI * U2);
        double[] myReturn = [A1, A2];
        return myReturn;
    }

    public static double nextdouble(int seed)
    {
        MersenneTwister random = _local;

        switch (random)
        {
            case null:
            {
                byte[] buffer = new byte[4];
                _global.GetBytes(buffer);
                _local = random = new MersenneTwister(seed);
                break;
            }
        }

        return random.NextDouble();
    }

    public static double nextdouble()
    {
        MersenneTwister random = _local;

        switch (random)
        {
            case null:
            {
                byte[] buffer = new byte[4];
                _global.GetBytes(buffer);
                _local = random = new MersenneTwister(BitConverter.ToInt32(buffer, 0));
                break;
            }
        }

        return random.NextDouble();
    }

    public static int nextint()
    {
        MersenneTwister random = _local;

        switch (random)
        {
            case null:
            {
                byte[] buffer = new byte[4];
                _global.GetBytes(buffer);
                _local = random = new MersenneTwister(BitConverter.ToInt32(buffer, 0));
                break;
            }
        }

        return nextint(random);
    }

    private static int nextint(MersenneTwister random)
    {
        return random.Next();
    }
        
    public static int nextint(int min, int max, int seed)
    {
        MersenneTwister random = _local;

        switch (random)
        {
            case null:
            {
                byte[] buffer = new byte[4];
                _global.GetBytes(buffer);
                _local = random = new MersenneTwister(seed);
                break;
            }
        }

        return nextint(random, min, max);
    }

    public static int nextint(int min, int max)
    {
        MersenneTwister random = _local;

        switch (random)
        {
            case null:
            {
                byte[] buffer = new byte[4];
                _global.GetBytes(buffer);
                _local = random = new MersenneTwister(BitConverter.ToInt32(buffer, 0));
                break;
            }
        }

        return nextint(random, min, max);
    }

    private static int nextint(MersenneTwister random, int min, int max)
    {
        return random.Next(min, max);
    }
}