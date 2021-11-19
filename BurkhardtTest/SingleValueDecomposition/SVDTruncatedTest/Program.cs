using System;
using Burkardt;
using Burkardt.Types;
using Burkardt.Uniform;

namespace SVDTruncatedTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SVD_TRUNCATED.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 March 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int m;
        int n;

        Console.WriteLine("");
        Console.WriteLine("SVD_TRUNCATED");
        Console.WriteLine("  Demonstrate the use of the truncated or economy-size");
        Console.WriteLine("  Singular Value Decomposition (SVD) for cases where");
        Console.WriteLine("  the sizes of M and N are very different.");

        m = 4;
        n = 3;
        svd_truncated_u_test(m, n);

        m = 3;
        n = 4;
        svd_truncated_v_test(m, n);

        Console.WriteLine("");
        Console.WriteLine("SVD_TRUNCATED");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }


    private static void svd_truncated_u_test(int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SVD_TRUNCATED_U_TEST tests SVD_TRUNCATED_U.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 March 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a;
        double[] a_save;
        double err;
        int i;
        int j;
        int k;
        int seed;
        double[] sn;
        double[] un;
        double[] v;

        Console.WriteLine("");
        Console.WriteLine("SVD_TRUNCATED_U_TEST");
        Console.WriteLine("  M = " + m + "");
        Console.WriteLine("  N = " + n + "");

        a = new double[m * n];
        un = new double[m * n];
        sn = new double[n * n];
        v = new double[n * n];

        seed = 123456789;

        a_save = UniformRNG.r8mat_uniform_01_new(m, n, ref seed);

        typeMethods.r8mat_print(m, n, a_save, "  A:");

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m; i++)
            {
                a[i + j * m] = a_save[i + j * m];
            }
        }

        SingleValueDecomposition.svd_truncated_u(m, n, a, ref un, ref sn, ref v);

        typeMethods.r8mat_print(m, n, un, "  UN:");
        typeMethods.r8mat_print(n, n, sn, "  SN:");
        typeMethods.r8mat_print(n, n, v, "  V:");
        //
        //  Check the factorization by computing A = U * S * V'
        //
        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m; i++)
            {
                a[i + j * m] = 0.0;
                for (k = 0; k < n; k++)
                {
                    a[i + j * m] += un[i + k * m] * sn[k + k * n] * v[j + k * n];
                }
            }
        }

        err = 0.0;
        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m; i++)
            {
                err = Math.Max(err, Math.Abs(a[i + j * m] - a_save[i + j * m]));
            }
        }

        Console.WriteLine("");
        Console.WriteLine("  Maximum error |A - U*S*V'| = " + err + "");

        typeMethods.r8mat_print(m, n, a, "  Recomputed A = U * S * V':");
    }


    private static void svd_truncated_v_test(int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SVD_TRUNCATED_V_TEST tests SVD_TRUNCATED_V.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 March 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a;
        double[] a_save;
        double err;
        int i;
        int j;
        int k;
        int seed;
        double[] sm;
        double[] u;
        double[] vm;

        Console.WriteLine("");
        Console.WriteLine("SVD_TRUNCATED_V_TEST");
        Console.WriteLine("  M = " + m + "");
        Console.WriteLine("  N = " + n + "");

        a = new double[m * n];
        u = new double[m * m];
        sm = new double[m * m];
        vm = new double[n * m];

        seed = 123456789;

        a_save = UniformRNG.r8mat_uniform_01_new(m, n, ref seed);

        typeMethods.r8mat_print(m, n, a_save, "  A:");

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m; i++)
            {
                a[i + j * m] = a_save[i + j * m];
            }
        }

        SingleValueDecomposition.svd_truncated_v(m, n, a, ref u, ref sm, ref vm);

        typeMethods.r8mat_print(m, m, u, "  U:");
        typeMethods.r8mat_print(m, m, sm, "  SM:");
        typeMethods.r8mat_print(n, m, vm, "  VM:");
        //
        //  Check the factorization by computing A = U * S * V'
        //
        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m; i++)
            {
                a[i + j * m] = 0.0;
                for (k = 0; k < m; k++)
                {
                    a[i + j * m] += u[i + k * m] * sm[k + k * m] * vm[j + k * n];
                }
            }
        }

        err = 0.0;
        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m; i++)
            {
                err = Math.Max(err, Math.Abs(a[i + j * m] - a_save[i + j * m]));
            }
        }

        Console.WriteLine("");
        Console.WriteLine("  Maximum error |A - U*S*V'| = " + err + "");

        typeMethods.r8mat_print(m, n, a, "  Recomputed A = U * S * V':");

    }
}