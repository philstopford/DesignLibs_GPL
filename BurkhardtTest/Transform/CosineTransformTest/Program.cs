using System;
using Burkardt.Transform;
using Burkardt.Uniform;

namespace CosineTransformTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for COSINE_TRANSFORM_TEST.
        //
        //  Discussion:
        //
        //    COSINE_TRANSFORM_TEST tests the COSINE_TRANSFORM library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 August 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("COSINE_TRANSFORM_TEST");
        Console.WriteLine("  Test the COSINE_TRANSFORM library.");

        cosine_transform_test01();

        Console.WriteLine("");
        Console.WriteLine("COSINE_TRANSFORM_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void cosine_transform_test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COSINE_TRANSFORM_TEST01 applies the DCT and its inverse.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 August 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int n = 10;
        int seed;
        double[] r;
        double[] s;
        double[] t;

        seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("COSINE_TRANSFORM_TEST01:");
        Console.WriteLine("  COSINE_TRANSFORM_DATA does a cosine transform of data");
        Console.WriteLine("  defined by a vector.");
        Console.WriteLine("");
        Console.WriteLine("  Apply the transform, then its inverse.");
        Console.WriteLine("  Let R be a random N vector.");
        Console.WriteLine("  Let S be the transform of D.");
        Console.WriteLine("  Let T be the transform of E.");
        Console.WriteLine("  Then R and T will be equal.");

        r = UniformRNG.r8vec_uniform_01_new(n, ref seed);
        s = Cosine.cosine_transform_data(n, r);
        t = Cosine.cosine_transform_inverse(n, s);

        Console.WriteLine("");
        Console.WriteLine("     I      R(I)        S(I)        T(I)");
        Console.WriteLine("");

        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString().PadLeft(4)
                                   + "  " + r[i].ToString().PadLeft(10)
                                   + "  " + s[i].ToString().PadLeft(10)
                                   + "  " + t[i].ToString().PadLeft(10) + "");
        }
    }
}