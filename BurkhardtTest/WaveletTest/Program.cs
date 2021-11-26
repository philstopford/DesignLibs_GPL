using System;
using System.Globalization;
using Burkardt.DaubechiesWavelet;
using Burkardt.Types;
using Burkardt.Uniform;

namespace WaveletTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for WAVELET_TEST.
        //
        //  Discussion:
        //
        //    WAVELET_TEST tests the WAVELET library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 May 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("WAVELET_TEST");
            
        Console.WriteLine("  Test the WAVELET library.");

        test01();
        test02();
        test03();
        test04();
        test05();
        test06();
        test07();
        test08();
        test09();

        test10();
        Console.WriteLine("");
        Console.WriteLine("WAVELET_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests DAUB2_TRANSFORM and DAUB2_TRANSFORM_INVERSE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  DAUB2_TRANSFORM computes the DAUB2 transform of a vector.");
        Console.WriteLine("  DAUB2_TRANSFORM_INVERSE inverts it.");
        //
        //  Random data.
        //
        int n = 16;
        int seed = 123456789;

        double[] u = UniformRNG.r8vec_uniform_01_new(n, ref seed);

        double[] v = Daub2.daub2_transform(n, u);

        double[] w = Daub2.daub2_transform_inverse(n, v);

        Console.WriteLine("");
        Console.WriteLine("   i      U          D2(U)    D2inv(D2(U))");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        //
        //  Constant signal.
        //
        n = 8;
        u = new double[8];
        for (i = 0; i < n; i++)
        {
            u[i] = 1.0;
        }

        v = Daub2.daub2_transform(n, u);

        w = Daub2.daub2_transform_inverse(n, v);

        Console.WriteLine("");
        Console.WriteLine("   i      U          D2(U)    D2inv(D2(U))");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        //
        //  Linear signal.
        //
        n = 16;
        double a_first = 1.0;
        double a_last = n;
        u = typeMethods.r8vec_linspace_new(n, a_first, a_last);

        v = Daub2.daub2_transform(n, u);

        w = Daub2.daub2_transform_inverse(n, v);

        Console.WriteLine("");
        Console.WriteLine("   i      U          D2(U)    D2inv(D2(U))");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        //
        //  Quadratic data.
        //
        n = 8;
        u = new double[n];
        for (i = 0; i < n; i++)
        {
            u[i] = Math.Pow(i - 5, 2);
        }

        v = Daub2.daub2_transform(n, u);

        w = Daub2.daub2_transform_inverse(n, v);

        Console.WriteLine("");
        Console.WriteLine("   i      U          D2(U)    D2inv(D2(U))");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests DAUB4_TRANSFORM and DAUB4_TRANSFORM_INVERSE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  DAUB4_TRANSFORM computes the DAUB4 transform of a vector.");
        Console.WriteLine("  DAUB4_TRANSFORM_INVERSE inverts it.");
        //
        //  Random data.
        //
        int n = 16;
        int seed = 123456789;

        double[] u = UniformRNG.r8vec_uniform_01_new(n, ref seed);

        double[] v = Daub4.daub4_transform(n, u);

        double[] w = Daub4.daub4_transform_inverse(n, v);

        Console.WriteLine("");
        Console.WriteLine("   i      U          D4(U)    D4inv(D4(U))");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        //
        //  Constant signal.
        //
        n = 8;
        u = new double[8];
        for (i = 0; i < n; i++)
        {
            u[i] = 1.0;
        }

        v = Daub4.daub4_transform(n, u);

        w = Daub4.daub4_transform_inverse(n, v);

        Console.WriteLine("");
        Console.WriteLine("   i      U          D4(U)    D4inv(D4(U))");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        //
        //  Linear signal.
        //
        n = 16;
        double a_first = 1.0;
        double a_last = n;
        u = typeMethods.r8vec_linspace_new(n, a_first, a_last);

        v = Daub4.daub4_transform(n, u);

        w = Daub4.daub4_transform_inverse(n, v);

        Console.WriteLine("");
        Console.WriteLine("   i      U          D4(U)    D4inv(D4(U))");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        //
        //  Quadratic data.
        //
        n = 8;
        u = new double[n];
        for (i = 0; i < n; i++)
        {
            u[i] = Math.Pow(i - 5, 2);
        }

        v = Daub4.daub4_transform(n, u);

        w = Daub4.daub4_transform_inverse(n, v);

        Console.WriteLine("");
        Console.WriteLine("   i      U          D4(U)    D4inv(D4(U))");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests DAUB6_TRANSFORM and DAUB6_TRANSFORM_INVERSE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  DAUB6_TRANSFORM computes the DAUB6 transform of a vector.");
        Console.WriteLine("  DAUB6_TRANSFORM_INVERSE inverts it.");
        //
        //  Random data.
        //
        int n = 16;
        int seed = 123456789;

        double[] u = UniformRNG.r8vec_uniform_01_new(n, ref seed);

        double[] v = Daub6.daub6_transform(n, u);

        double[] w = Daub6.daub6_transform_inverse(n, v);

        Console.WriteLine("");
        Console.WriteLine("   i      U          D6(U)    D6inv(D6(U))");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        //
        //  Constant signal.
        //
        n = 8;
        u = new double[8];
        for (i = 0; i < n; i++)
        {
            u[i] = 1.0;
        }

        v = Daub6.daub6_transform(n, u);

        w = Daub6.daub6_transform_inverse(n, v);

        Console.WriteLine("");
        Console.WriteLine("   i      U          D6(U)    D6inv(D6(U))");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        //
        //  Linear signal.
        //
        n = 16;
        double a_first = 1.0;
        double a_last = n;
        u = typeMethods.r8vec_linspace_new(n, a_first, a_last);

        v = Daub6.daub6_transform(n, u);

        w = Daub6.daub6_transform_inverse(n, v);

        Console.WriteLine("");
        Console.WriteLine("   i      U          D6(U)    D6inv(D6(U))");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        //
        //  Quadratic data.
        //
        n = 8;
        u = new double[n];
        for (i = 0; i < n; i++)
        {
            u[i] = Math.Pow(i - 5, 2);
        }

        v = Daub6.daub6_transform(n, u);

        w = Daub6.daub6_transform_inverse(n, v);

        Console.WriteLine("");
        Console.WriteLine("   i      U          D6(U)    D6inv(D6(U))");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }
    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests DAUB8_TRANSFORM and DAUB8_TRANSFORM_INVERSE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  DAUB8_TRANSFORM computes the DAUB8 transform of a vector.");
        Console.WriteLine("  DAUB8_TRANSFORM_INVERSE inverts it.");
        //
        //  Random data.
        //
        int n = 16;
        int seed = 123456789;

        double[] u = UniformRNG.r8vec_uniform_01_new(n, ref seed);

        double[] v = Daub8.daub8_transform(n, u);

        double[] w = Daub8.daub8_transform_inverse(n, v);

        Console.WriteLine("");
        Console.WriteLine("   i      U          D8(U)    D8inv(D8(U))");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        //
        //  Constant signal.
        //
        n = 8;
        u = new double[8];
        for (i = 0; i < n; i++)
        {
            u[i] = 1.0;
        }

        v = Daub8.daub8_transform(n, u);

        w = Daub8.daub8_transform_inverse(n, v);

        Console.WriteLine("");
        Console.WriteLine("   i      U          D8(U)    D8inv(D8(U))");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        //
        //  Linear signal.
        //
        n = 16;
        double a_first = 1.0;
        double a_last = n;
        u = typeMethods.r8vec_linspace_new(n, a_first, a_last);

        v = Daub8.daub8_transform(n, u);

        w = Daub8.daub8_transform_inverse(n, v);

        Console.WriteLine("");
        Console.WriteLine("   i      U          D8(U)    D8inv(D8(U))");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        //
        //  Quadratic data.
        //
        n = 8;
        u = new double[n];
        for (i = 0; i < n; i++)
        {
            u[i] = Math.Pow(i - 5, 2);
        }

        v = Daub8.daub8_transform(n, u);

        w = Daub8.daub8_transform_inverse(n, v);

        Console.WriteLine("");
        Console.WriteLine("   i      U          D8(U)    D8inv(D8(U))");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }
    }

    private static void test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 tests DAUB10_TRANSFORM and DAUB10_TRANSFORM_INVERSE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 May 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  DAUB10_TRANSFORM computes the DAUB10 transform of a vector.");
        Console.WriteLine("  DAUB10_TRANSFORM_INVERSE inverts it.");
        //
        //  Random data.
        //
        int n = 16;
        int seed = 123456789;

        double[] u = UniformRNG.r8vec_uniform_01_new(n, ref seed);

        double[] v = Daub10.daub10_transform(n, u);

        double[] w = Daub10.daub10_transform_inverse(n, v);

        Console.WriteLine("");
        Console.WriteLine("   i      U          D10(U)    D10inv(D10(U))");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        //
        //  Constant signal.
        //
        n = 8;
        u = new double[8];
        for (i = 0; i < n; i++)
        {
            u[i] = 1.0;
        }

        v = Daub10.daub10_transform(n, u);

        w = Daub10.daub10_transform_inverse(n, v);

        Console.WriteLine("");
        Console.WriteLine("   i      U          D10(U)    D10inv(D10(U))");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        //
        //  Linear signal.
        //
        n = 16;
        double a_first = 1.0;
        double a_last = n;
        u = typeMethods.r8vec_linspace_new(n, a_first, a_last);

        v = Daub10.daub10_transform(n, u);

        w = Daub10.daub10_transform_inverse(n, v);

        Console.WriteLine("");
        Console.WriteLine("   i      U          D10(U)    D10inv(D10(U))");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        //
        //  Quadratic data.
        //
        n = 8;
        u = new double[n];
        for (i = 0; i < n; i++)
        {
            u[i] = Math.Pow(i - 5, 2);
        }

        v = Daub10.daub10_transform(n, u);

        w = Daub10.daub10_transform_inverse(n, v);

        Console.WriteLine("");
        Console.WriteLine("   i      U          D10(U)    D10inv(D10(U))");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }
    }

    private static void test06()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST06 tests DAUB12_TRANSFORM and DAUB12_TRANSFORM_INVERSE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 May 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("TEST06");
        Console.WriteLine("  DAUB12_TRANSFORM computes the DAUB12 transform of a vector.");
        Console.WriteLine("  DAUB12_TRANSFORM_INVERSE inverts it.");
        //
        //  Random data.
        //
        int n = 16;
        int seed = 123456789;

        double[] u = UniformRNG.r8vec_uniform_01_new(n, ref seed);

        double[] v = Daub12.daub12_transform(n, u);

        double[] w = Daub12.daub12_transform_inverse(n, v);

        Console.WriteLine("");
        Console.WriteLine("   i      U          D12(U)    D12inv(D12(U))");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        //
        //  Constant signal.
        //
        n = 8;
        u = new double[8];
        for (i = 0; i < n; i++)
        {
            u[i] = 1.0;
        }

        v = Daub12.daub12_transform(n, u);

        w = Daub12.daub12_transform_inverse(n, v);

        Console.WriteLine("");
        Console.WriteLine("   i      U          D12(U)    D12inv(D12(U))");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        //
        //  Linear signal.
        //
        n = 16;
        double a_first = 1.0;
        double a_last = n;
        u = typeMethods.r8vec_linspace_new(n, a_first, a_last);

        v = Daub12.daub12_transform(n, u);

        w = Daub12.daub12_transform_inverse(n, v);

        Console.WriteLine("");
        Console.WriteLine("   i      U          D12(U)    D12inv(D12(U))");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        //
        //  Quadratic data.
        //
        n = 8;
        u = new double[n];
        for (i = 0; i < n; i++)
        {
            u[i] = Math.Pow(i - 5, 2);
        }

        v = Daub12.daub12_transform(n, u);

        w = Daub12.daub12_transform_inverse(n, v);

        Console.WriteLine("");
        Console.WriteLine("   i      U          D12(U)    D12inv(D12(U))");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }
    }

    private static void test07()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST07 tests DAUB14_TRANSFORM and DAUB14_TRANSFORM_INVERSE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 May 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("TEST07");
        Console.WriteLine("  DAUB14_TRANSFORM computes the DAUB14 transform of a vector.");
        Console.WriteLine("  DAUB14_TRANSFORM_INVERSE inverts it.");
        //
        //  Random data.
        //
        int n = 16;
        int seed = 123456789;

        double[] u = UniformRNG.r8vec_uniform_01_new(n, ref seed);

        double[] v = Daub14.daub14_transform(n, u);

        double[] w = Daub14.daub14_transform_inverse(n, v);

        Console.WriteLine("");
        Console.WriteLine("   i      U          D14(U)    D14inv(D14(U))");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        //
        //  Constant signal.
        //
        n = 8;
        u = new double[8];
        for (i = 0; i < n; i++)
        {
            u[i] = 1.0;
        }

        v = Daub14.daub14_transform(n, u);

        w = Daub14.daub14_transform_inverse(n, v);

        Console.WriteLine("");
        Console.WriteLine("   i      U          D14(U)    D14inv(D14(U))");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        //
        //  Linear signal.
        //
        n = 16;
        double a_first = 1.0;
        double a_last = n;
        u = typeMethods.r8vec_linspace_new(n, a_first, a_last);

        v = Daub14.daub14_transform(n, u);

        w = Daub14.daub14_transform_inverse(n, v);

        Console.WriteLine("");
        Console.WriteLine("   i      U          D14(U)    D14inv(D14(U))");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        //
        //  Quadratic data.
        //
        n = 8;
        u = new double[n];
        for (i = 0; i < n; i++)
        {
            u[i] = Math.Pow(i - 5, 2);
        }

        v = Daub14.daub14_transform(n, u);

        w = Daub14.daub14_transform_inverse(n, v);

        Console.WriteLine("");
        Console.WriteLine("   i      U          D14(U)    D14inv(D14(U))");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }
    }

    private static void test08()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST08 tests DAUB16_TRANSFORM and DAUB16_TRANSFORM_INVERSE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 May 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("TEST08");
        Console.WriteLine("  DAUB16_TRANSFORM computes the DAUB16 transform of a vector.");
        Console.WriteLine("  DAUB16_TRANSFORM_INVERSE inverts it.");
        //
        //  Random data.
        //
        int n = 16;
        int seed = 123456789;

        double[] u = UniformRNG.r8vec_uniform_01_new(n, ref seed);

        double[] v = Daub16.daub16_transform(n, u);

        double[] w = Daub16.daub16_transform_inverse(n, v);

        Console.WriteLine("");
        Console.WriteLine("   i      U          D16(U)    D16inv(D16(U))");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        //
        //  Constant signal.
        //
        n = 8;
        u = new double[8];
        for (i = 0; i < n; i++)
        {
            u[i] = 1.0;
        }

        v = Daub16.daub16_transform(n, u);

        w = Daub16.daub16_transform_inverse(n, v);

        Console.WriteLine("");
        Console.WriteLine("   i      U          D16(U)    D16inv(D16(U))");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        //
        //  Linear signal.
        //
        n = 16;
        double a_first = 1.0;
        double a_last = n;
        u = typeMethods.r8vec_linspace_new(n, a_first, a_last);

        v = Daub16.daub16_transform(n, u);

        w = Daub16.daub16_transform_inverse(n, v);

        Console.WriteLine("");
        Console.WriteLine("   i      U          D16(U)    D16inv(D16(U))");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        //
        //  Quadratic data.
        //
        n = 8;
        u = new double[n];
        for (i = 0; i < n; i++)
        {
            u[i] = Math.Pow(i - 5, 2);
        }

        v = Daub16.daub16_transform(n, u);

        w = Daub16.daub16_transform_inverse(n, v);

        Console.WriteLine("");
        Console.WriteLine("   i      U          D16(U)    D16inv(D16(U))");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }
    }

    private static void test09()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST09 tests DAUB18_TRANSFORM and DAUB18_TRANSFORM_INVERSE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 May 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("TEST09");
        Console.WriteLine("  DAUB18_TRANSFORM computes the DAUB18 transform of a vector.");
        Console.WriteLine("  DAUB18_TRANSFORM_INVERSE inverts it.");
        //
        //  Random data.
        //
        int n = 16;
        int seed = 123456789;

        double[] u = UniformRNG.r8vec_uniform_01_new(n, ref seed);

        double[] v = Daub18.daub18_transform(n, u);

        double[] w = Daub18.daub18_transform_inverse(n, v);

        Console.WriteLine("");
        Console.WriteLine("   i      U          D18(U)    D18inv(D18(U))");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        //
        //  Constant signal.
        //
        n = 8;
        u = new double[8];
        for (i = 0; i < n; i++)
        {
            u[i] = 1.0;
        }

        v = Daub18.daub18_transform(n, u);

        w = Daub18.daub18_transform_inverse(n, v);

        Console.WriteLine("");
        Console.WriteLine("   i      U          D18(U)    D18inv(D18(U))");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        //
        //  Linear signal.
        //
        n = 16;
        double a_first = 1.0;
        double a_last = n;
        u = typeMethods.r8vec_linspace_new(n, a_first, a_last);

        v = Daub18.daub18_transform(n, u);

        w = Daub18.daub18_transform_inverse(n, v);

        Console.WriteLine("");
        Console.WriteLine("   i      U          D18(U)    D18inv(D18(U))");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        //
        //  Quadratic data.
        //
        n = 8;
        u = new double[n];
        for (i = 0; i < n; i++)
        {
            u[i] = Math.Pow(i - 5, 2);
        }

        v = Daub18.daub18_transform(n, u);

        w = Daub18.daub18_transform_inverse(n, v);

        Console.WriteLine("");
        Console.WriteLine("   i      U          D18(U)    D18inv(D18(U))");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }
    }

    private static void test10()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST10 tests DAUB20_TRANSFORM and DAUB20_TRANSFORM_INVERSE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 May 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("TEST10");
        Console.WriteLine("  DAUB20_TRANSFORM computes the DAUB20 transform of a vector.");
        Console.WriteLine("  DAUB20_TRANSFORM_INVERSE inverts it.");
        //
        //  Random data.
        //
        int n = 16;
        int seed = 123456789;

        double[] u = UniformRNG.r8vec_uniform_01_new(n, ref seed);

        double[] v = Daub20.daub20_transform(n, u);

        double[] w = Daub20.daub20_transform_inverse(n, v);

        Console.WriteLine("");
        Console.WriteLine("   i      U          D20(U)    D20inv(D20(U))");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        //
        //  Constant signal.
        //
        n = 8;
        u = new double[8];
        for (i = 0; i < n; i++)
        {
            u[i] = 1.0;
        }

        v = Daub20.daub20_transform(n, u);

        w = Daub20.daub20_transform_inverse(n, v);

        Console.WriteLine("");
        Console.WriteLine("   i      U          D20(U)    D20inv(D20(U))");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        //
        //  Linear signal.
        //
        n = 16;
        double a_first = 1.0;
        double a_last = n;
        u = typeMethods.r8vec_linspace_new(n, a_first, a_last);

        v = Daub20.daub20_transform(n, u);

        w = Daub20.daub20_transform_inverse(n, v);

        Console.WriteLine("");
        Console.WriteLine("   i      U          D20(U)    D20inv(D20(U))");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        //
        //  Quadratic data.
        //
        n = 8;
        u = new double[n];
        for (i = 0; i < n; i++)
        {
            u[i] = Math.Pow(i - 5, 2);
        }

        v = Daub20.daub20_transform(n, u);

        w = Daub20.daub20_transform_inverse(n, v);

        Console.WriteLine("");
        Console.WriteLine("   i      U          D20(U)    D20inv(D20(U))");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }
    }
}