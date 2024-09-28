using System;
using System.Security.Cryptography;

namespace entropyRNG;

public static class Crypto_RNG
{
    private const double maxInt = int.MaxValue; // max int value.

    public static double[] random_gauss3()
    {
        double[] myReturn = random_gauss();
        return [myReturn[0] / 3.0f, myReturn[1] / 3.0f];
    }

    public static double[] random_gauss()
    {
        RandomNumberGenerator random = RandomNumberGenerator.Create();

        byte[] tmp = new byte[sizeof(int)];

        // Box-Muller transform
        // We aren't allowed 0, so we reject any values approaching zero.
        random.GetBytes(tmp);
        var U1 = BitConverter.ToInt32(tmp, 0) / maxInt;
        while (U1 < 1E-15)
        {
            random.GetBytes(tmp);
            U1 = BitConverter.ToInt32(tmp, 0) / maxInt;
        }
        random.GetBytes(tmp);
        var U2 = BitConverter.ToInt32(tmp, 0) / maxInt;
        while (U2 < 1E-15)
        {
            random.GetBytes(tmp);
            U2 = BitConverter.ToInt32(tmp, 0) / maxInt;
        }
        // PAs are 3-sigma, so this needs to be divided by 3 to give single sigma value when used
        double A1 = Math.Sqrt(-2 * Math.Log(U2, Math.E)) * Math.Cos(2 * Math.PI * U1);
        double A2 = Math.Sqrt(-2 * Math.Log(U1, Math.E)) * Math.Sin(2 * Math.PI * U2);
        double[] myReturn = [A1, A2];
        return myReturn;
    }

    public static double nextdouble()
    {
        RandomNumberGenerator random = RandomNumberGenerator.Create();

        byte[] tmp = new byte[sizeof(int)];

        random.GetBytes(tmp);
        return BitConverter.ToInt32(tmp, 0) / maxInt;
    }

    public static int nextint()
    {
        RandomNumberGenerator random = RandomNumberGenerator.Create();

        byte[] tmp = new byte[sizeof(int)];

        random.GetBytes(tmp);
        return BitConverter.ToInt32(tmp, 0);
    }
}