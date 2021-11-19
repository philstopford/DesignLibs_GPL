using System;
using System.Numerics;
using Burkardt.FourierTransform;
using Burkardt.Types;
using Burkardt.Uniform;

namespace SlowFourierTransformTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SFTPACK_TEST.
        //
        //  Discussion:
        //
        //    SFTPACK_TEST tests the SFTPACK library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //   27 June 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("SFTPACK_TEST");

        Console.WriteLine("  Test the SFTPACK library.");

        c4mat_sft_test();
        c4vec_sft_test();

        c8mat_sft_test();
        c8vec_sft_test();

        r4vec_sft_test();

        r8vec_sct_test();
        r8vec_sft_test();
        r8vec_sht_test();
        r8vec_sqct_test();
        r8vec_sqst_test();
        r8vec_sst_test();
        r8vec_swt_test();

        Console.WriteLine("");
        Console.WriteLine("SFTPACK_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void c4mat_sft_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C4MAT_SFT_TEST tests C4MAT_SFTB and C4MAT_SFTF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 June 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n1 = 10;
        int n2 = 4;
        int seed;
        Complex[] x;
        Complex[] x2;
        Complex[] y;

        Console.WriteLine("");
        Console.WriteLine("C4MAT_SFT_TEST");
        Console.WriteLine("  C4MAT_SFTF computes the forward slow Fourier transform.");
        Console.WriteLine("  C4MAT_SFTB computes the backward slow Fourier transform.");
        Console.WriteLine("");
        Console.WriteLine("  The data has dimensions N1 = " + n1 + ", N2 = " + n2 + "");

        seed = 123456789;

        x = UniformRNG.c4mat_uniform_01_new(n1, n2, ref seed);

        typeMethods.c4mat_print_some(n1, n2, x, 1, 1, 10, 10, "  The original data:");
        //
        //  Compute the slow Fourier transform of the data.
        //
        y = Slow.c4mat_sftf(n1, n2, x);

        typeMethods.c4mat_print_some(n1, n2, y, 1, 1, 10, 10, "  The Fourier coefficients:");
        //
        //  Now try to retrieve the data from the coefficients.
        //
        x2 = Slow.c4mat_sftb(n1, n2, y);

        typeMethods.c4mat_print_some(n1, n2, x2, 1, 1, 10, 10, "  The retrieved data:");
    }

    private static void c4vec_sft_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C4VEC_SFT_TEST tests C4VEC_SFTB and C4VEC_SFTF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 June 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n = 36;
        int seed;
        Complex[] x;
        Complex[] x2;
        Complex[] y;

        Console.WriteLine("");
        Console.WriteLine("C4VEC_SFT_TEST");
        Console.WriteLine("  C4VEC_SFTF computes the forward slow Fourier transform.");
        Console.WriteLine("  C4VEC_SFTB computes the backward slow Fourier transform.");
        Console.WriteLine("");
        Console.WriteLine("  The number of data values, N = " + n + "");

        seed = 123456789;

        x = UniformRNG.c4vec_uniform_01_new(n, ref seed);

        typeMethods.c4vec_print_part(n, x, 10, "  The original data:");
        //
        //  Compute the slow Fourier transform of the data.
        //
        y = Slow.c4vec_sftf(n, x);

        typeMethods.c4vec_print_part(n, y, 10, "  The Fourier coefficients:");
        //
        //  Now try to retrieve the data from the coefficients.
        //
        x2 = Slow.c4vec_sftb(n, y);

        typeMethods.c4vec_print_part(n, x2, 10, "  The retrieved data:");
    }

    private static void c8mat_sft_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8MAT_SFT_TEST tests C8MAT_SFTB and C8MAT_SFTF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 June 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n1 = 10;
        int n2 = 4;
        int seed;
        Complex[] x;
        Complex[] x2;
        Complex[] y;

        Console.WriteLine("");
        Console.WriteLine("C8MAT_SFT_TEST");
        Console.WriteLine("  C8MAT_SFTF computes the forward slow Fourier transform.");
        Console.WriteLine("  C8MAT_SFTB computes the backward slow Fourier transform.");
        Console.WriteLine("");
        Console.WriteLine("  The data has dimensions N1 = " + n1 + ", N2 = " + n2 + "");

        seed = 123456789;

        x = UniformRNG.c8mat_uniform_01_new(n1, n2, ref seed);

        typeMethods.c8mat_print_some(n1, n2, x, 1, 1, 10, 10, "  The original data:");
        //
        //  Compute the slow Fourier transform of the data.
        //
        y = Slow.c8mat_sftf(n1, n2, x);

        typeMethods.c8mat_print_some(n1, n2, y, 1, 1, 10, 10, "  The Fourier coefficients:");
        //
        //  Now try to retrieve the data from the coefficients.
        //
        x2 = Slow.c8mat_sftb(n1, n2, y);

        typeMethods.c8mat_print_some(n1, n2, x2, 1, 1, 10, 10, "  The retrieved data:");
    }

    private static void c8vec_sft_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8VEC_SFT_TEST tests C8VEC_SFTB and C8VEC_SFTF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 June 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n = 36;
        int seed;
        Complex[] x;
        Complex[] x2;
        Complex[] y;

        Console.WriteLine("");
        Console.WriteLine("C8VEC_SFT_TEST");
        Console.WriteLine("  C8VEC_SFTF computes the forward slow Fourier transform.");
        Console.WriteLine("  C8VEC_SFTB computes the backward slow Fourier transform.");
        Console.WriteLine("");
        Console.WriteLine("  The number of data values, N = " + n + "");

        seed = 123456789;

        x = UniformRNG.c8vec_uniform_01_new(n, ref seed);

        typeMethods.c8vec_print_part(n, x, 10, "  The original data:");
        //
        //  Compute the slow Fourier transform of the data.
        //
        y = Slow.c8vec_sftf(n, x);

        typeMethods.c8vec_print_part(n, y, 10, "  The Fourier coefficients:");
        //
        //  Now try to retrieve the data from the coefficients.
        //
        x2 = Slow.c8vec_sftb(n, y);

        typeMethods.c8vec_print_part(n, x2, 10, "  The retrieved data:");
    }

    private static void r4vec_sft_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R4VEC_SFT_TEST tests R4VEC_SFTB and R4VEC_SFTF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 June 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        float[] a;
        float ahi = 5.0f;
        float alo = 0.0f;
        float azero = 0;
        float[] b;
        int i;
        int n = 36;
        int seed;
        float[] x;
        float[] z;

        Console.WriteLine("");
        Console.WriteLine("R4VEC_SFT_TEST");
        Console.WriteLine("  R4VEC_SFTF computes the forward slow Fourier transform.");
        Console.WriteLine("  R4VEC_SFTB computes the backward slow Fourier transform.");
        Console.WriteLine("");
        Console.WriteLine("  The number of data values, N = " + n + "");

        seed = 123456789;

        x = UniformRNG.r4vec_uniform_ab_new(n, alo, ahi, ref seed);

        typeMethods.r4vec_print_part(n, x, 10, "  The original data:");
        //
        //  Compute the slow Fourier transform of the data.
        //
        a = new float[n / 2];
        b = new float[n / 2];

        Slow.r4vec_sftf(n, x, ref azero, ref a, ref b);

        Console.WriteLine("");
        Console.WriteLine("  A (cosine) coefficients:");
        Console.WriteLine("");

        Console.WriteLine("  " + 0.ToString().PadLeft(4)
                               + "  " + azero.ToString().PadLeft(14) + "");

        for (i = 0; i < n / 2; i++)
        {
            Console.WriteLine("  " + i.ToString().PadLeft(4)
                                   + "  " + a[i].ToString().PadLeft(14) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  B (sine) coefficients:");
        Console.WriteLine("");

        for (i = 0; i < n / 2; i++)
        {
            Console.WriteLine("  " + i.ToString().PadLeft(4)
                                   + "  " + b[i].ToString().PadLeft(14) + "");
        }

        //
        //  Now try to retrieve the data from the coefficients.
        //
        z = Slow.r4vec_sftb(n, azero, a, b);

        typeMethods.r4vec_print_part(n, z, 10, "  The retrieved data:");
    }

    private static void r8vec_sct_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_SCT_TEST tests R8VEC_SCT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 February 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double ahi = 5.0;
        double alo = 0.0;
        double[] c;
        double[] d;
        double[] e;
        int i;
        int n = 256;
        int seed;

        Console.WriteLine("");
        Console.WriteLine("R8VEC_SCTTEST");
        Console.WriteLine("  R8VEC_SCT does a forward or backward slow cosine transform.");
        Console.WriteLine("");
        Console.WriteLine("  The number of data items is N = " + n + "");
        //
        //  Set the data values.
        //
        seed = 123456789;

        c = UniformRNG.r8vec_uniform_ab_new(n, alo, ahi, ref seed);

        typeMethods.r8vec_print_part(n, c, 10, "  The original data:");
        //
        //  Compute the coefficients.
        //
        d = Slow.r8vec_sct(n, c);

        typeMethods.r8vec_print_part(n, d, 10, "  The cosine coefficients:");
        //
        //  Now compute inverse transform of coefficients.  Should get back the
        //  original data.
        //
        e = Slow.r8vec_sct(n, d);

        for (i = 0; i < n; i++)
        {
            e[i] /= 2 * n;
        }

        typeMethods.r8vec_print_part(n, e, 10, "  The retrieved data:");
    }

    private static void r8vec_sft_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_SFT_TEST tests R8VEC_SFTB and R8VEC_SFTF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 February 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a;
        double ahi = 5.0;
        double alo = 0.0;
        double azero = 0;
        double[] b;
        int i;
        int n = 36;
        int seed;
        double[] x;
        double[] z;

        Console.WriteLine("");
        Console.WriteLine("R8VEC_SFT_TEST");
        Console.WriteLine("  R8VEC_SFTF computes the forward slow Fourier transform.");
        Console.WriteLine("  R8VEC_SFTB computes the backward slow Fourier transform.");
        Console.WriteLine("");
        Console.WriteLine("  The number of data values, N = " + n + "");

        seed = 123456789;

        x = UniformRNG.r8vec_uniform_ab_new(n, alo, ahi, ref seed);

        typeMethods.r8vec_print_part(n, x, 10, "  The original data:");
        //
        //  Compute the slow Fourier transform of the data.
        //
        a = new double[n / 2];
        b = new double[n / 2];

        Slow.r8vec_sftf(n, x, ref azero, ref a, ref b);

        Console.WriteLine("");
        Console.WriteLine("  A (cosine) coefficients:");
        Console.WriteLine("");

        Console.WriteLine("  " + 0.ToString().PadLeft(4)
                               + "  " + azero.ToString().PadLeft(14) + "");

        for (i = 0; i < n / 2; i++)
        {
            Console.WriteLine("  " + i.ToString().PadLeft(4)
                                   + "  " + a[i].ToString().PadLeft(14) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  B (sine) coefficients:");
        Console.WriteLine("");

        for (i = 0; i < n / 2; i++)
        {
            Console.WriteLine("  " + i.ToString().PadLeft(4)
                                   + "  " + b[i].ToString().PadLeft(14) + "");
        }

        //
        //  Now try to retrieve the data from the coefficients.
        //
        z = Slow.r8vec_sftb(n, azero, a, b);

        typeMethods.r8vec_print_part(n, z, 10, "  The retrieved data:");
    }

    private static void r8vec_sht_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_SHT_TEST tests R8VEC_SHT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 February 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double ahi = 5.0;
        double alo = 0.0;
        double[] c;
        double[] d;
        double[] e;
        int n = 17;
        int seed;

        Console.WriteLine("");
        Console.WriteLine("R8VEC_SHT_TEST");
        Console.WriteLine("  R8VEC_SHT does a forward or backward slow Hartley transform.");
        Console.WriteLine("");
        Console.WriteLine("  The number of data items is N = " + n + "");
        //
        //  Set the data values.
        //
        seed = 123456789;

        c = UniformRNG.r8vec_uniform_ab_new(n, alo, ahi, ref seed);

        typeMethods.r8vec_print_part(n, c, 10, "  The original data:");
        //
        //  Compute the coefficients.
        //
        d = Slow.r8vec_sht(n, c);

        typeMethods.r8vec_print_part(n, d, 10, "  The Hartley coefficients:");
        //
        //  Now compute inverse transform of coefficients.  Should get back the
        //  original data.
        //
        e = Slow.r8vec_sht(n, d);

        typeMethods.r8vec_print_part(n, e, 10, "  The retrieved data:");
    }

    private static void r8vec_sqct_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_SQCT_TEST tests R8VEC_SQCTB and R8VEC_SQCTF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 February 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double ahi = 5.0;
        double alo = 0.0;
        int n = 256;
        int seed;
        double[] x;
        double[] y;
        double[] z;

        Console.WriteLine("");
        Console.WriteLine("R8VEC_SQCT_TEST");
        Console.WriteLine("  R8VEC_SQCTF does a forward slow quarter wave cosine transform;");
        Console.WriteLine("  R8VEC_SQCTB does a backward slow quarter wave cosine transform.");
        Console.WriteLine("");
        Console.WriteLine("  The number of data items is N = " + n + "");
        //
        //  Set the data values.
        //
        seed = 123456789;

        x = UniformRNG.r8vec_uniform_ab_new(n, alo, ahi, ref seed);

        typeMethods.r8vec_print_part(n, x, 10, "  The original data:");
        //
        //  Compute the coefficients.
        //
        y = Slow.r8vec_sqctf(n, x);

        typeMethods.r8vec_print_part(n, y, 10, "  The cosine coefficients:");
        //
        //  Now compute inverse transform of coefficients.  Should get back the
        //  original data.
        //
        z = Slow.r8vec_sqctb(n, y);

        typeMethods.r8vec_print_part(n, z, 10, "  The retrieved data:");
    }

    private static void r8vec_sqst_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_SQST_TEST tests R8VEC_SQSTB and R8VEC_SQSTF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 February 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double ahi = 5.0;
        double alo = 0.0;
        int n = 256;
        int seed;
        double[] x;
        double[] y;
        double[] z;

        Console.WriteLine("");
        Console.WriteLine("R8VEC_SQST_TEST");
        Console.WriteLine("  R8VEC_SQSTF does a forward slow quarter wave sine transform;");
        Console.WriteLine("  R8VEC_SQSTB does a backward slow quarter wave sine transform.");
        Console.WriteLine("");
        Console.WriteLine("  The number of data items is N = " + n + "");
        //
        //  Set the data values.
        //
        seed = 123456789;

        x = UniformRNG.r8vec_uniform_ab_new(n, alo, ahi, ref seed);

        typeMethods.r8vec_print_part(n, x, 10, "  The original data:");
        //
        //  Compute the coefficients.
        //
        y = Slow.r8vec_sqstf(n, x);

        typeMethods.r8vec_print_part(n, y, 10, "  The sine coefficients:");
        //
        //  Now compute inverse transform of coefficients.  Should get back the
        //  original data.
        //
        z = Slow.r8vec_sqstb(n, y);

        typeMethods.r8vec_print_part(n, z, 10, "  The retrieved data:");
    }

    private static void r8vec_sst_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_SST_TEST tests R8VEC_SST.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 February 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double ahi = 5.0;
        double alo = 0.0;
        double[] c;
        double[] d;
        double[] e;
        int i;
        int n = 256;
        int seed;

        Console.WriteLine("");
        Console.WriteLine("R8VEC_SST_TEST");
        Console.WriteLine("  R8VEC_SST does a forward or backward slow sine transform.");
        Console.WriteLine("");
        Console.WriteLine("  The number of data items is N = " + n + "");
        //
        //  Set the data values.
        //
        seed = 123456789;

        c = UniformRNG.r8vec_uniform_ab_new(n, alo, ahi, ref seed);

        typeMethods.r8vec_print_part(n, c, 10, "  The original data:");
        //
        //  Compute the coefficients;
        //
        d = Slow.r8vec_sst(n, c);

        typeMethods.r8vec_print_part(n, d, 10, "  The sine coefficients:");
        //
        //  Now compute inverse transform of coefficients.  Should get back the
        //  original data.
        //
        e = Slow.r8vec_sst(n, d);

        for (i = 0; i < n; i++)
        {
            e[i] /= 2 * (n + 1);
        }

        typeMethods.r8vec_print_part(n, e, 10, "  The retrieved data:");
    }

    private static void r8vec_swt_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_SWT_TEST tests R8VEC_SWTB and R8VEC_SWTF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 June 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double ahi;
        double alo;
        double[] d;
        int i;
        int n = 36;
        int np1h;
        double[] s;
        int seed;
        double[] x;

        alo = 0.0;
        ahi = 5.0;

        Console.WriteLine("");
        Console.WriteLine("R8VEC_SWT_TEST");
        Console.WriteLine("  R8VEC_SWTF computes the forward slow wavelet transform.");
        Console.WriteLine("  R8VEC_SWTB computes the backward slow wavelet transform.");
        Console.WriteLine("");
        Console.WriteLine("  The number of data values, N = " + n + "");

        seed = 123456789;

        x = UniformRNG.r8vec_uniform_ab_new(n, alo, ahi, ref seed);

        typeMethods.r8vec_print_part(n, x, 10, "  The original data:");
        //
        //  Compute the slow wavelet transform of the data.
        //
        np1h = (n + 1) / 2;
        s = new double[np1h];
        d = new double[np1h];

        Slow.r8vec_swtf(n, x, ref s, ref d);

        Console.WriteLine("");
        Console.WriteLine("     I          S(I)            D(I)");
        Console.WriteLine("");
        for (i = 0; i < np1h; i++)
        {
            Console.WriteLine("  " + i.ToString().PadLeft(4)
                                   + "  " + s[i].ToString().PadLeft(14)
                                   + "  " + d[i].ToString().PadLeft(14) + "");
        }

        //
        //  Now try to retrieve the data from the coefficients.
        //

        x = Slow.r8vec_swtb(n, ref s, ref d);

        typeMethods.r8vec_print_part(n, x, 10, "  The retrieved data:");
    }
}