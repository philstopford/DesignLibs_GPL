using System;
using System.Linq;

namespace Noise;

/// <summary>
/// This class generates several noise "octaves", frequencies, and amplitudes, then combines them to make the final Simplex Noise.
/// </summary>
public class SimplexNoise
{
    private SimplexNoiseOctave[] octaves;
    private double[] frequencys;
    private double[] amplitudes;

    public SimplexNoise(int seed) : this(30, 0.7, seed)
    {
        // Defaults from original repository were 850, 0.7, seed.
        // 850 produces a very rough surface, so dialled it down.
    }

    public SimplexNoise(int largestFeature, double persistence, int seed)
    {
        // Receives a number (eg 128) and calculates what power of 2 it is (eg 2^7)
        int numberOfOctaves = (int)Math.Ceiling(Math.Log10(largestFeature) / Math.Log10(2));

        octaves = new SimplexNoiseOctave[numberOfOctaves];
        frequencys = new double[numberOfOctaves];
        amplitudes = new double[numberOfOctaves];

        Random rand = new(seed);

        for (int i = 0; i < numberOfOctaves; i++)
        {
            octaves[i] = new SimplexNoiseOctave(rand.Next());

            frequencys[i] = Math.Pow(2, i);
            amplitudes[i] = Math.Pow(persistence, octaves.Length - i);
        }
    }

    /// <summary>
    /// Combines the separate octaves, frequencies, and amplitudes into a single Simplex Noise value for the specified position.
    /// </summary>
    /// <param name="x">The x position in the noise.</param>
    /// <param name="y">The y position in the noise.</param>
    /// <returns></returns>
    public double GetNoise(int x, int y)
    {
        return octaves.Select((t, i) => t.Noise(x / frequencys[i], y / frequencys[i]) * amplitudes[i]).Sum();
    }

    public double GetNoise(double x, double y)
    {
        return octaves.Select((t, i) => t.Noise(x / frequencys[i], y / frequencys[i]) * amplitudes[i]).Sum();
    }
}

/// <summary>
/// A speed-improved simplex noise algorithm for 2D in C#.  Based on example code by Stefan Gustavson (stegu@itn.liu.se).
/// Optimizations by Peter Eastman (peastman@drizzle.stanford.edu).  Better rank ordering method by Stefan Gustavson in 2012.
/// 
/// This code was placed in the public domain by its original author, Stefan Gustavson. You may use it as you see fit, but
/// attribution is appreciated.
/// </summary>
public class SimplexNoiseOctave
{
    public const int RandomSeed = 0;        // Determines if a random seed should be used
    private const int NumberOfSwaps = 400;

    /// <summary>
    /// A private inner class to for gradient computations.
    /// </summary>
    private class Grad
    {
        public double x { get; private set; }
        public double y { get; private set; }
        public double z { get; private set; }
        public double w { get; private set; }

        public Grad(double x, double y, double z)
        {
            this.x = x;
            this.y = y;
            this.z = z;
        }

        public Grad(double x, double y, double z, double w)
        {
            this.x = x;
            this.y = y;
            this.z = z;
            this.w = w;
        }
    }

    private static Grad[] grad3 = { new(1,1,0),new(-1,1,0),new(1,-1,0),new(-1,-1,0),
        new(1,0,1),new(-1,0,1),new(1,0,-1),new(-1,0,-1),
        new(0,1,1),new(0,-1,1),new(0,1,-1),new(0,-1,-1) };

    // Contains all the numbers between 0 and 255, these are put in a random order depending upon the seed
    private static short[] p_supply = { 151,160,137,91,90,15,131,13,201,95,96,53,194,233,7,225,140,
        36,103,30,69,142,8,99,37,240,21,10,23,190, 6,148,247,120,234,
        75,0,26,197,62,94,252,219,203,117,35,11,32,57,177,33,88,237,
        149,56,87,174,20,125,136,171,168,68,175,74,165,71,134,139,
        48,27,166,77,146,158,231,83,111,229,122,60,211,133,230,220,
        105,92,41,55,46,245,40,244,102,143,54,65,25,63,161,1,216,80,
        73,209,76,132,187,208, 89,18,169,200,196,135,130,116,188,159,
        86,164,100,109,198,173,186,3,64,52,217,226,250,124,123,5,
        202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,
        182,189,28,42,223,183,170,213,119,248,152,2,44,154,163,70,221,
        153,101,155,167,43,172,9,129,22,39,253,19,98,108,110,79,113,
        224,232,178,185,112,104,218,246,97,228,251,34,242,193,238,210,
        144,12,191,179,162,241,81,51,145,235,249,14,239,107,49,192,
        214,31,181,199,106,157,184,84,204,176,115,121,50,45,127,4,150,
        254,138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,
        215,61,156,180 };

    private short[] p = new short[p_supply.Length];

    // To remove the need for index wrapping, double the permutation table length
    private short[] perm = new short[512];
    private short[] permMod12 = new short[512];

    public SimplexNoiseOctave(int seed)
    {
        p_supply.CopyTo(p, 0);

        seed = seed switch
        {
            RandomSeed => new Random().Next(),
            _ => seed
        };

        // The random for the swaps
        Random rand = new(seed);

        // The seed determines the swaps that occur between the default order and the order we're actually going to use
        for (int i = 0; i < NumberOfSwaps; i++)
        {
            int swapFrom = rand.Next(p.Length);
            int swapTo = rand.Next(p.Length);

            (p[swapFrom], p[swapTo]) = (p[swapTo], p[swapFrom]);
        }


        for (int i = 0; i < 512; i++)
        {
            perm[i] = p[i & 255];
            permMod12[i] = (short)(perm[i] % 12);
        }
    }

    // Skewing and unskewing factors for 2 dimensions
    private static double F2 = 0.5 * (Math.Sqrt(3.0) - 1.0);
    private static double G2 = (3.0 - Math.Sqrt(3.0)) / 6.0;

    /// <summary>
    /// Fast floor function - this method is a *lot* faster than using (int)Math.floor(x) - at least that's the case for the original Java
    /// code, not sure about .NET. :)
    /// </summary>
    /// <param name="x"></param>
    /// <returns></returns>
    private static int FastFloor(double x)
    {
        int xi = (int)x;
        return x < xi ? xi - 1 : xi;
    }

    /// <summary>
    ///  Computes the dot product for the gradient and specified position.
    /// </summary>
    /// <param name="g"></param>
    /// <param name="x"></param>
    /// <param name="y"></param>
    /// <returns></returns>
    private static double dot(Grad g, double x, double y)
    {
        return g.x * x + g.y * y;
    }

    /// <summary>
    /// Computes a 2D Simplex Noise value for the specified location in the noise.
    /// </summary>
    /// <param name="xin">The x position in the noise.</param>
    /// <param name="yin">The y position in the noise.</param>
    /// <returns></returns>
    public double Noise(double xin, double yin)
    {
        double n0, n1, n2; // Noise contributions from the three corners

        // Skew the input space to determine which simplex cell we're in
        double s = (xin + yin) * F2; // Hairy factor for 2D
        int i = FastFloor(xin + s);
        int j = FastFloor(yin + s);
        double t = (i + j) * G2;
        double X0 = i - t; // Unskew the cell origin back to (x,y) space
        double Y0 = j - t;
        double x0 = xin - X0; // The x,y distances from the cell origin
        double y0 = yin - Y0;

        // For the 2D case, the simplex shape is an equilateral triangle.
        // Determine which simplex we are in.
        int i1, j1; // Offsets for second (middle) corner of simplex in (i,j) coords
        if (x0 > y0)
        {
            // lower triangle, XY order: (0,0)->(1,0)->(1,1)
            i1 = 1;
            j1 = 0;
        }
        else
        {
            // upper triangle, YX order: (0,0)->(0,1)->(1,1)
            i1 = 0;
            j1 = 1;
        }

        // A step of (1,0) in (i,j) means a step of (1-c,-c) in (x,y), and
        // a step of (0,1) in (i,j) means a step of (-c,1-c) in (x,y), where
        // c = (3-sqrt(3))/6
        double x1 = x0 - i1 + G2; // Offsets for middle corner in (x,y) unskewed coords
        double y1 = y0 - j1 + G2;
        double x2 = x0 - 1.0 + 2.0 * G2; // Offsets for last corner in (x,y) unskewed coords
        double y2 = y0 - 1.0 + 2.0 * G2;

        // Work out the hashed gradient indices of the three simplex corners
        int ii = i & 255;
        int jj = j & 255;
        int gi0 = permMod12[ii + perm[jj]];
        int gi1 = permMod12[ii + i1 + perm[jj + j1]];
        int gi2 = permMod12[ii + 1 + perm[jj + 1]];

        // Calculate the contribution from the three corners
        double t0 = 0.5 - x0 * x0 - y0 * y0;
        switch (t0)
        {
            case < 0:
                n0 = 0.0;
                break;
            default:
                t0 *= t0;
                n0 = t0 * t0 * dot(grad3[gi0], x0, y0);  // (x,y) of grad3 used for 2D gradient
                break;
        }

        double t1 = 0.5 - x1 * x1 - y1 * y1;
        switch (t1)
        {
            case < 0:
                n1 = 0.0;
                break;
            default:
                t1 *= t1;
                n1 = t1 * t1 * dot(grad3[gi1], x1, y1);
                break;
        }

        double t2 = 0.5 - x2 * x2 - y2 * y2;
        switch (t2)
        {
            case < 0:
                n2 = 0.0;
                break;
            default:
                t2 *= t2;
                n2 = t2 * t2 * dot(grad3[gi2], x2, y2);
                break;
        }

        // Add contributions from each corner to get the final noise value.
        // The result is scaled to return values in the interval [-1,1].
        return 70.0 * (n0 + n1 + n2);
    }
}

public class SimplexNoiseOld
{
    private const double Div6 = 1.0 / 6.0;
    private const double Div3 = 1.0 / 3.0;

    /* Unit vectors for gradients to points on cube,equal distances apart (ie vector from center to the middle of each side */
    private static readonly int[][] Grad3 = new[]
    {
        new[] { 1, 1, 0 }, new[] { -1, 1, 0 }, new[] { 1, -1, 0 }, new[] { -1, -1, 0 },
        new[] { 1, 0, 1 }, new[] { -1, 0, 1 }, new[] { 1, 0, -1 }, new[] { -1, 0, -1 },
        new[] { 0, 1, 1 }, new[] { 0, -1, 1 }, new[] { 0, 1, -1 }, new[] { 0, -1, -1 }
    };

    //0..255, randomized
    private static readonly int[] P =
    {
        151, 160, 137, 91, 90, 15, 131, 13, 201, 95, 96, 53, 194, 233, 7, 225, 140, 36, 103,
        30, 69, 142, 8, 99, 37, 240, 21, 10, 23, 190, 6, 148, 247, 120, 234, 75, 0, 26, 197,
        62, 94, 252, 219, 203, 117, 35, 11, 32, 57, 177, 33, 88, 237, 149, 56, 87, 174, 20,
        125, 136, 171, 168, 68, 175, 74, 165, 71, 134, 139, 48, 27, 166, 77, 146, 158, 231,
        83, 111, 229, 122, 60, 211, 133, 230, 220, 105, 92, 41, 55, 46, 245, 40, 244, 102,
        143, 54, 65, 25, 63, 161, 1, 216, 80, 73, 209, 76, 132, 187, 208, 89, 18, 169, 200,
        196, 135, 130, 116, 188, 159, 86, 164, 100, 109, 198, 173, 186, 3, 64, 52, 217, 226,
        250, 124, 123, 5, 202, 38, 147, 118, 126, 255, 82, 85, 212, 207, 206, 59, 227, 47, 16,
        58, 17, 182, 189, 28, 42, 223, 183, 170, 213, 119, 248, 152, 2, 44, 154, 163, 70,
        221, 153, 101, 155, 167, 43, 172, 9, 129, 22, 39, 253, 19, 98, 108, 110, 79, 113, 224,
        232, 178, 185, 112, 104, 218, 246, 97, 228, 251, 34, 242, 193, 238, 210, 144, 12,
        191, 179, 162, 241, 81, 51, 145, 235, 249, 14, 239, 107, 49, 192, 214, 31, 181, 199,
        106, 157, 184, 84, 204, 176, 115, 121, 50, 45, 127, 4, 150, 254, 138, 236, 205, 93,
        222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180
    };

    private static readonly int[] Perm = new int[512];

    private static readonly double Sqrt3 = Math.Sqrt(3);
    private static readonly double F2 = 0.5 * (Sqrt3 - 1.0);
    private static readonly double G2 = (3.0 - Sqrt3) * Div6;

    static SimplexNoiseOld()
    {
        for (int i = 0; i < 512; i++)
        {
            Perm[i] = P[i & 255];
        }

        Singleton = new SimplexNoiseOld();
    }

    public static SimplexNoiseOld Singleton { get; private set; }

    public double Noise01(double x, double y)
    {
        return (Noise(x, y) + 1) * 0.5;
    }

    public double MultiNoise(int octaves, double x, double y)
    {
        double value = 0.0;
        float mul = 1;
        for (int i = 0; i < octaves; i++)
        {
            value += Noise((x + 10) * mul, (y + 15) * mul) / mul;

            mul *= 2;
        }
        return value;
    }

    public double MultiNoise01(int octaves, double x, double y)
    {
        return (MultiNoise(octaves, x, y) + 1.0) * 0.5;
    }

    public double RidgedMulti(int octaves, double x, double y)
    {
        double value = 0.0;
        double mul = 1;
        for (int i = 0; i < octaves; i++)
        {
            double added = Noise(x * mul, y * mul) / mul;
            value += Math.Abs(added);

            mul *= 2.18387276;
        }
        return value;
    }

    public double Noise(double xin, double yin)
    {
        double n0, n1, n2; // Noise contributions from the three corners

        // Skew the input space to a square to determine which simplex cell we're in

        double s = (xin + yin) * F2; // Hairy factor for 2D
        int i = FastFloor(xin + s);
        int j = FastFloor(yin + s);

        double t = (i + j) * G2;
        double x0p = i - t; // Unskew the cell origin back to (x,y) space
        double y0p = j - t;
        double x0 = xin - x0p; // The x,y distances from the cell origin
        double y0 = yin - y0p;

        // For the 2D case, the simplex shape is an equilateral triangle.
        // Determine which simplex we are in.
        int i1, j1; // Offsets for second (middle) corner of simplex in (i,j) coords
        if (x0 > y0)
        {
            // lower triangle, XY order: (0,0)->(1,0)->(1,1)
            i1 = 1;
            j1 = 0;
        }
        else
        {
            i1 = 0;
            j1 = 1;
        } // upper triangle, YX order: (0,0)->(0,1)->(1,1)

        // A step of (1,0) in (i,j) means a step of (1-c,-c) in (x,y), and
        // a step of (0,1) in (i,j) means a step of (-c,1-c) in (x,y), where
        // c = (3-sqrt(3))/6

        double x1 = x0 - i1 + G2; // Offsets for middle corner in (x,y) unskewed coords
        double y1 = y0 - j1 + G2;
        double x2 = x0 - 1.0 + 2.0 * G2; // Offsets for last corner in (x,y) unskewed coords
        double y2 = y0 - 1.0 + 2.0 * G2;

        // Work out the hashed gradient indices of the three simplex corners
        int ii = i & 255;
        int jj = j & 255;
        int gi0 = Perm[ii + Perm[jj]] % 12;
        int gi1 = Perm[ii + i1 + Perm[jj + j1]] % 12;
        int gi2 = Perm[ii + 1 + Perm[jj + 1]] % 12;

        // Calculate the contribution from the three corners
        double t0 = 0.5 - x0 * x0 - y0 * y0;
        switch (t0)
        {
            case < 0:
                n0 = 0.0;
                break;
            default:
                t0 *= t0;
                n0 = t0 * t0 * Dot(Grad3[gi0], x0, y0); // (x,y) of grad3 used for 2D gradient
                break;
        }
        double t1 = 0.5 - x1 * x1 - y1 * y1;
        switch (t1)
        {
            case < 0:
                n1 = 0.0;
                break;
            default:
                t1 *= t1;
                n1 = t1 * t1 * Dot(Grad3[gi1], x1, y1);
                break;
        }
        double t2 = 0.5 - x2 * x2 - y2 * y2;
        switch (t2)
        {
            case < 0:
                n2 = 0.0;
                break;
            default:
                t2 *= t2;
                n2 = t2 * t2 * Dot(Grad3[gi2], x2, y2);
                break;
        }
        // Add contributions from each corner to get the final noise value.
        // The result is scaled to return values in the interval [-1,1].
        return 70.0 * (n0 + n1 + n2);
    }

    public double Multi01(int octaves, double x, double y, double z)
    {
        return (Multi(octaves, x, y, z) + 1) * 0.5;
    }

    public double Multi(int octaves, double x, double y, double z)
    {
        double value = 0.0;
        double mul = 1;
        for (int i = 0; i < octaves; i++)
        {
            double added = Noise(x * mul, y * mul, z * mul) / mul;
            value += added;
            mul *= 2;
        }
        return value;
    }

    public double Noise01(double x, double y, double z)
    {
        // Noise  is in the range -1 to +1
        double val = Noise(x, y, z);
        return (val + 1) * 0.5;
    }

    public double RidgedMulti(int octaves, double x, double y, double z)
    {
        double value = 0.0;
        double mul = 1;
        for (int i = 0; i < octaves; i++)
        {
            double added = Noise(x * mul, y * mul, z * mul) / mul;
            value += Math.Abs(added);
            mul *= 2;
        }
        return value;
    }

    public double Noise(double xin, double yin, double zin)
    {
        double n0, n1, n2, n3; // Noise contributions from the four corners
        // Skew the input space to determine which simplex cell we're in

        double s = (xin + yin + zin) * Div3; // Very nice and simple skew factor for 3D
        int i = FastFloor(xin + s);
        int j = FastFloor(yin + s);
        int k = FastFloor(zin + s);

        double t = (i + j + k) * Div6;
        double ax0 = i - t; // Unskew the cell origin back to (x,y,z) space
        double ay0 = j - t;
        double az0 = k - t;
        double x0 = xin - ax0; // The x,y,z distances from the cell origin
        double y0 = yin - ay0;
        double z0 = zin - az0;
        // For the 3D case, the simplex shape is a slightly irregular tetrahedron.
        // Determine which simplex we are in.
        int i1, j1, k1; // Offsets for second corner of simplex in (i,j,k) coords
        int i2, j2, k2; // Offsets for third corner of simplex in (i,j,k) coords
        if (x0 >= y0)
        {
            if (y0 >= z0)
            {
                i1 = 1;
                j1 = 0;
                k1 = 0;
                i2 = 1;
                j2 = 1;
                k2 = 0;
            }
            else if (x0 >= z0)
            {
                i1 = 1;
                j1 = 0;
                k1 = 0;
                i2 = 1;
                j2 = 0;
                k2 = 1;
            }
            else
            {
                i1 = 0;
                j1 = 0;
                k1 = 1;
                i2 = 1;
                j2 = 0;
                k2 = 1;
            }
        }
        else
        {
            // x0<y0
            if (y0 < z0)
            {
                i1 = 0;
                j1 = 0;
                k1 = 1;
                i2 = 0;
                j2 = 1;
                k2 = 1;
            }
            else if (x0 < z0)
            {
                i1 = 0;
                j1 = 1;
                k1 = 0;
                i2 = 0;
                j2 = 1;
                k2 = 1;
            }
            else
            {
                i1 = 0;
                j1 = 1;
                k1 = 0;
                i2 = 1;
                j2 = 1;
                k2 = 0;
            }
        }
        // A step of (1,0,0) in (i,j,k) means a step of (1-c,-c,-c) in (x,y,z),
        // a step of (0,1,0) in (i,j,k) means a step of (-c,1-c,-c) in (x,y,z), and
        // a step of (0,0,1) in (i,j,k) means a step of (-c,-c,1-c) in (x,y,z), where
        // c = 1/6.
        double x1 = x0 - i1 + Div6; // Offsets for second corner in (x,y,z) coords
        double y1 = y0 - j1 + Div6;
        double z1 = z0 - k1 + Div6;
        double x2 = x0 - i2 + 2.0 * Div6; // Offsets for third corner in (x,y,z) coords
        double y2 = y0 - j2 + 2.0 * Div6;
        double z2 = z0 - k2 + 2.0 * Div6;
        double x3 = x0 - 1.0 + 3.0 * Div6; // Offsets for last corner in (x,y,z) coords
        double y3 = y0 - 1.0 + 3.0 * Div6;
        double z3 = z0 - 1.0 + 3.0 * Div6;
        // Work out the hashed gradient indices of the four simplex corners
        int ii = i & 255;
        int jj = j & 255;
        int kk = k & 255;
        int gi0 = Perm[ii + Perm[jj + Perm[kk]]] % 12;
        int gi1 = Perm[ii + i1 + Perm[jj + j1 + Perm[kk + k1]]] % 12;
        int gi2 = Perm[ii + i2 + Perm[jj + j2 + Perm[kk + k2]]] % 12;
        int gi3 = Perm[ii + 1 + Perm[jj + 1 + Perm[kk + 1]]] % 12;
        // Calculate the contribution from the four corners
        double t0 = 0.6 - x0 * x0 - y0 * y0 - z0 * z0;
        switch (t0)
        {
            case < 0:
                n0 = 0.0;
                break;
            default:
                t0 *= t0;
                n0 = t0 * t0 * Dot(Grad3[gi0], x0, y0, z0);
                break;
        }
        double t1 = 0.6 - x1 * x1 - y1 * y1 - z1 * z1;
        switch (t1)
        {
            case < 0:
                n1 = 0.0;
                break;
            default:
                t1 *= t1;
                n1 = t1 * t1 * Dot(Grad3[gi1], x1, y1, z1);
                break;
        }
        double t2 = 0.6 - x2 * x2 - y2 * y2 - z2 * z2;
        switch (t2)
        {
            case < 0:
                n2 = 0.0;
                break;
            default:
                t2 *= t2;
                n2 = t2 * t2 * Dot(Grad3[gi2], x2, y2, z2);
                break;
        }
        double t3 = 0.6 - x3 * x3 - y3 * y3 - z3 * z3;
        switch (t3)
        {
            case < 0:
                n3 = 0.0;
                break;
            default:
                t3 *= t3;
                n3 = t3 * t3 * Dot(Grad3[gi3], x3, y3, z3);
                break;
        }
        // Add contributions from each corner to get the final noise value.
        // The result is scaled to stay just inside [-1,1]
        return 32.0 * (n0 + n1 + n2 + n3);
    }

    private static int FastFloor(double x)
    {
        // This method is a *lot* faster than using (int)Math.floor(x)
        // This comment is regards to C, 
        return x > 0 ? (int)x : (int)x - 1;
    }

    private static double Dot(int[] g, double x, double y)
    {
        return g[0] * x + g[1] * y;
    }

    private static double Dot(int[] g, double x, double y, double z)
    {
        return g[0] * x + g[1] * y + g[2] * z;
    }
}