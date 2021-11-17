using System;
using Burkardt.Types;

namespace Burkardt;

public static class ContinuedFraction
{
    public static void i4cf_evaluate(int n, int[] a, int[] b, ref int[] p, ref int[] q)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4CF_EVALUATE evaluates a continued fraction with I4 entries.
        //
        //  Discussion:
        //
        //    For convenience, we omit the parentheses or multiline display.
        //
        //    CF = A(0) + B(1) / A(1) + B(2) / A(2) + ... A(N-1) + B(N) / A(N).
        //
        //  Licensing:
        //
        //    I don't care what you do with this code.
        //
        //  Modified:
        //
        //    07 August 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, Charles Mesztenyi,
        //    John Rice, Henry Thatcher, Christoph Witzgall,
        //    Computer Approximations,
        //    Wiley, 1968.
        //
        //  Parameters:
        //
        //    Input, int N, the number of continued fraction
        //    coefficients.
        //
        //    Input, int A[N+1], B[N+1], the continued fraction 
        //    coefficients.
        //
        //    Output, int P[N+1], Q[N+1], the N+1 successive 
        //    approximations to the value of the continued fraction.
        //
    {
        int i;

        for (i = 0; i <= n; i++)
        {
            switch (i)
            {
                case 0:
                    p[i] = a[i] * 1 + 0;
                    q[i] = a[i] * 0 + 1;
                    break;
                case 1:
                    p[i] = a[i] * p[i - 1] + b[i] * 1;
                    q[i] = a[i] * q[i - 1] + b[i] * 0;
                    break;
                default:
                    p[i] = a[i] * p[i - 1] + b[i] * p[i - 2];
                    q[i] = a[i] * q[i - 1] + b[i] * q[i - 2];
                    break;
            }
        }

    }

    public static void i4cf_evaluate_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4CF_EVALUATE_TEST tests I4CF_EVALUATE.
        //
        //  Licensing:
        //
        //    I don't care what you do with this code.
        //
        //  Modified:
        //
        //    05 August 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] a =
        {
            3,
            6, 6, 6, 6, 6,
            6, 6, 6, 6, 6,
            6, 6, 6, 6, 6,
            6, 6, 6, 6
        };
        int[] b =
        {
            0,
            1, 9, 25, 49, 81,
            121, 169, 225, 289, 361,
            441, 529, 625, 729, 841,
            961, 1089, 1225, 1369
        };
        const int n = 19;
        int i;

        Console.WriteLine("");
        Console.WriteLine("I4CF_EVALUATE_TEST:");

        int[] p = new int[n + 1];
        int[] q = new int[n + 1];

        i4cf_evaluate(n, a, b, ref p, ref q);

        Console.WriteLine("");
        Console.WriteLine("  CF numerators, denominators, ratios:");
        Console.WriteLine("");

        for (i = 0; i <= n; i++)
        {
            double t = p[i] / (double) q[i];
            Console.WriteLine("  " + i.ToString().PadLeft(2)
                                   + "  " + p[i].ToString().PadLeft(20)
                                   + "  " + q[i].ToString().PadLeft(20)
                                   + "  " + t.ToString("0.################").PadLeft(24) + "");
        }

    }

    public static void i4scf_evaluate(int n, int[] a, ref int[] p, ref int[] q)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4SCF_EVALUATE evaluates a simple continued fraction with I4 entries.
        //
        //  Discussion:
        //
        //    The simple continued fraction with integer coefficients is:
        //
        //      SCF = A(0) + 1 / ( A(1) + 1 / ( A(2) ... + 1 / A(N) ) )
        //
        //    This routine returns the successive approximants P[i]/Q[i]
        //    to the value of the rational number represented by the continued
        //    fraction, with the value exactly equal to the final ratio P(N)/Q(N).
        //
        //  Licensing:
        //
        //    I don't care what you do with this code.
        //
        //  Modified:
        //
        //    04 August 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, Charles Mesztenyi,
        //    John Rice, Henry Thatcher, Christoph Witzgall,
        //    Computer Approximations,
        //    Wiley, 1968.
        //
        //  Parameters:
        //
        //    Input, int N, the number of continued fraction
        //    coefficients.
        //
        //    Input, int A[N+1], the continued fraction coefficients.
        //
        //    Output, int P[N+1], Q[N+1], the numerators and
        //    denominators of the successive approximations.
        //
    {
        int i;

        for (i = 0; i <= n; i++)
        {
            switch (i)
            {
                case 0:
                    p[i] = a[i] * 1 + 0;
                    q[i] = a[i] * 0 + 1;
                    break;
                case 1:
                    p[i] = a[i] * p[i - 1] + 1;
                    q[i] = a[i] * q[i - 1] + 0;
                    break;
                default:
                    p[i] = a[i] * p[i - 1] + p[i - 2];
                    q[i] = a[i] * q[i - 1] + q[i - 2];
                    break;
            }
        }

    }

    public static void i4scf_evaluate_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4SCF_EVALUATE_TEST tests I4SCF_EVALUATE.
        //
        //  Licensing:
        //
        //    I don't care what you do with this code.
        //
        //  Modified:
        //
        //    04 August 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] a =
        {
            3, 7, 15, 1, 292,
            1, 1, 1, 2, 1,
            3, 1, 14, 2, 1,
            1, 2, 2, 2, 2
        };
        const int n = 19;
        int i;

        Console.WriteLine("");
        Console.WriteLine("I4SCF_EVALUATE_TEST:");

        int[] p = new int[n + 1];
        int[] q = new int[n + 1];

        i4scf_evaluate(n, a, ref p, ref q);

        Console.WriteLine("");
        Console.WriteLine("  SCF numerators, denominators, ratios:");
        Console.WriteLine("");

        for (i = 0; i <= n; i++)
        {
            double t = p[i] / (double) q[i];
            Console.WriteLine("  " + i.ToString().PadLeft(2)
                                   + "  " + p[i].ToString().PadLeft(20)
                                   + "  " + q[i].ToString().PadLeft(20)
                                   + "  " + t.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void i4vec_print(int n, int[] a, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4VEC_PRINT prints an I4VEC.
        //
        //  Discussion:
        //
        //    An I4VEC is a vector of I4's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 November 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of components of the vector.
        //
        //    Input, int A[N], the vector to be printed.
        //
        //    Input, string TITLE, a title.
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine(title + "");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString().PadLeft(8)
                                   + ": " + a[i].ToString().PadLeft(8) + "");
        }
    }

    public static void i8cf_evaluate(int n, long[] a, long[] b, ref long[] p,
            ref long[] q)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I8CF_EVALUATE evaluates a continued fraction with I8 entries.
        //
        //  Discussion:
        //
        //    For convenience, we omit the parentheses or multiline display.
        //
        //    CF = A(0) + B(1) / A(1) + B(2) / A(2) + ... A(N-1) + B(N) / A(N).
        //
        //  Licensing:
        //
        //    I don't care what you do with this code.
        //
        //  Modified:
        //
        //    05 August 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, Charles Mesztenyi,
        //    John Rice, Henry Thatcher, Christoph Witzgall,
        //    Computer Approximations,
        //    Wiley, 1968.
        //
        //  Parameters:
        //
        //    Input, int N, the number of continued fraction
        //    coefficients.
        //
        //    Input, long long int A[N+1], B[N+1], the continued fraction 
        //    coefficients.
        //
        //    Output, long long int P[N+1], Q[N+1], the N successive 
        //    approximations to the value of the continued fraction.
        //
    {
        int i;

        for (i = 0; i <= n; i++)
        {
            switch (i)
            {
                case 0:
                    p[i] = a[i] * 1 + 0;
                    q[i] = a[i] * 0 + 1;
                    break;
                case 1:
                    p[i] = a[i] * p[i - 1] + b[i] * 1;
                    q[i] = a[i] * q[i - 1] + b[i] * 0;
                    break;
                default:
                    p[i] = a[i] * p[i - 1] + b[i] * p[i - 2];
                    q[i] = a[i] * q[i - 1] + b[i] * q[i - 2];
                    break;
            }
        }

    }

    public static void i8cf_evaluate_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I8CF_EVALUATE_TEST tests I8CF_EVALUATE.
        //
        //  Licensing:
        //
        //    I don't care what you do with this code.
        //
        //  Modified:
        //
        //    05 August 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        long[] a =
        {
            3,
            6, 6, 6, 6, 6,
            6, 6, 6, 6, 6,
            6, 6, 6, 6, 6,
            6, 6, 6, 6
        };
        long[] b =
        {
            0,
            1, 9, 25, 49, 81,
            121, 169, 225, 289, 361,
            441, 529, 625, 729, 841,
            961, 1089, 1225, 1369
        };
        int i;
        const int n = 19;


        Console.WriteLine("");
        Console.WriteLine("I8CF_EVALUATE_TEST:");

        long[] p = new long[n + 1];
        long[] q = new long[n + 1];

        i8cf_evaluate(n, a, b, ref p, ref q);

        Console.WriteLine("");
        Console.WriteLine("  CF numerators, denominators, ratios:");
        Console.WriteLine("");

        for (i = 0; i <= n; i++)
        {
            double t = p[i] / (double) q[i];
            Console.WriteLine("  " + i.ToString().PadLeft(2)
                                   + "  " + p[i].ToString().PadLeft(20)
                                   + "  " + q[i].ToString().PadLeft(20)
                                   + "  " + t.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void i8scf_evaluate(int n, long[] a, ref long[] p,
            ref long[] q)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I8SCF_EVALUATE evaluates a simple continued fraction with I8 entries.
        //
        //  Discussion:
        //
        //    The simple continued fraction with integer coefficients is:
        //
        //      SCF = A(0) + 1 / ( A(1) + 1 / ( A(2) ... + 1 / A(N) ) )
        //
        //    This routine returns the successive approximants P[i]/Q[i]
        //    to the value of the rational number represented by the continued
        //    fraction, with the value exactly equal to the final ratio P(N)/Q(N).
        //
        //  Licensing:
        //
        //    I don't care what you do with this code.
        //
        //  Modified:
        //
        //    04 August 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, Charles Mesztenyi,
        //    John Rice, Henry Thatcher, Christoph Witzgall,
        //    Computer Approximations,
        //    Wiley, 1968.
        //
        //  Parameters:
        //
        //    Input, int N, the number of continued fraction
        //    coefficients.
        //
        //    Input, long long int A[N+1], the continued fraction coefficients.
        //
        //    Output, long long int P[N+1], Q[N+1], the numerators and
        //    denominators of the successive approximations.
        //
    {
        int i;

        for (i = 0; i <= n; i++)
        {
            switch (i)
            {
                case 0:
                    p[i] = a[i] * 1 + 0;
                    q[i] = a[i] * 0 + 1;
                    break;
                case 1:
                    p[i] = a[i] * p[i - 1] + 1;
                    q[i] = a[i] * q[i - 1] + 0;
                    break;
                default:
                    p[i] = a[i] * p[i - 1] + p[i - 2];
                    q[i] = a[i] * q[i - 1] + q[i - 2];
                    break;
            }
        }

    }

    public static void i8scf_evaluate_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I8SCF_EVALUATE_TEST tests I8SCF_EVALUATE.
        //
        //  Licensing:
        //
        //    I don't care what you do with this code.
        //
        //  Modified:
        //
        //    07 August 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        long[] a =
        {
            3, 7, 15, 1, 292,
            1, 1, 1, 2, 1,
            3, 1, 14, 2, 1,
            1, 2, 2, 2, 2
        };
        int i;
        const int n = 19;

        Console.WriteLine("");
        Console.WriteLine("I8SCF_EVALUATE_TEST:");

        long[] p = new long[n + 1];
        long[] q = new long[n + 1];

        i8scf_evaluate(n, a, ref p, ref q);

        Console.WriteLine("");
        Console.WriteLine("  SCF numerators, denominators, ratios:");
        Console.WriteLine("");

        for (i = 0; i <= n; i++)
        {
            double t = p[i] / (double) q[i];
            Console.WriteLine("  " + i.ToString().PadLeft(2)
                                   + "  " + p[i].ToString().PadLeft(20)
                                   + "  " + q[i].ToString().PadLeft(20)
                                   + "  " + t.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void r8_to_i4scf(double r, int n, ref int[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_TO_I4SCF approximates an R8 with an I4 simple continued fraction.
        //
        //  Discussion:
        //
        //    The simple continued fraction with real coefficients is:
        //
        //      SCF = A(0) + 1 / ( A(1) + 1 / ( A(2) ... + 1 / A(N) ) )
        //
        //  Licensing:
        //
        //    I don't care what you do with this code.
        //
        //  Modified:
        //
        //    04 August 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Norman Richert,
        //    Strang's Strange Figures,
        //    American Mathematical Monthly,
        //    Volume 99, Number 2, February 1992, pages 101-107.
        //
        //  Parameters:
        //
        //    Input, double R, the real value.
        //
        //    Input, int N, the number of convergents to compute.
        //
        //    Output, int A[N+1], the partial quotients.
        //
    {
        int i;

        switch (r)
        {
            case 0.0:
            {
                for (i = 0; i <= n; i++)
                {
                    a[i] = 0;
                }

                return;
            }
        }

        double r2 = r;
        a[0] = (int) r2;

        for (i = 1; i <= n; i++)
        {
            r2 = 1.0 / (r2 - a[i - 1]);
            a[i] = (int) r2;
        }

    }

    public static void r8_to_i4scf_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_TO_I4SCF_TEST tests R8_TO_I4SCF.
        //
        //  Licensing:
        //
        //    I don't care what you do with this code.
        //
        //  Modified:
        //
        //    04 August 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        const int n = 19;

        Console.WriteLine("");
        Console.WriteLine("R8_TO_I4SCF_TEST:");

        double r = Math.PI;

        int[] a = new int[n + 1];

        r8_to_i4scf(r, n, ref a);

        i4vec_print(n + 1, a, "  SCF coefficients:");

        int[] p = new int[n + 1];
        int[] q = new int[n + 1];

        i4scf_evaluate(n, a, ref p, ref q);

        Console.WriteLine("");
        Console.WriteLine("  SCF numerators, denominators, ratios:");
        Console.WriteLine("");

        for (i = 0; i <= n; i++)
        {
            double t = p[i] / (double) q[i];
            Console.WriteLine("  " + i.ToString().PadLeft(2)
                                   + "  " + p[i].ToString().PadLeft(20)
                                   + "  " + q[i].ToString().PadLeft(20)
                                   + "  " + t.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void r8_to_i8scf(double r, int n, long[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_TO_I8SCF approximates an R8 with an I8 simple continued fraction.
        //
        //  Discussion:
        //
        //    The simple continued fraction with real coefficients is:
        //
        //      SCF = A(0) + 1 / ( A(1) + 1 / ( A(2) ... + 1 / A(N) ) )
        //
        //  Licensing:
        //
        //    I don't care what you do with this code.
        //
        //  Modified:
        //
        //    04 August 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Norman Richert,
        //    Strang's Strange Figures,
        //    American Mathematical Monthly,
        //    Volume 99, Number 2, February 1992, pages 101-107.
        //
        //  Parameters:
        //
        //    Input, double R, the real value.
        //
        //    Input, int N, the number of convergents to compute.
        //
        //    Output, long long int A[N+1], the partial quotients.
        //
    {
        int i;

        switch (r)
        {
            case 0.0:
            {
                for (i = 0; i <= n; i++)
                {
                    a[i] = 0;
                }

                return;
            }
        }

        double r2 = r;
        a[0] = (long) r2;

        for (i = 1; i <= n; i++)
        {
            r2 = 1.0 / (r2 - a[i - 1]);
            a[i] = (long) r2;
        }

    }

    public static void r8_to_i8scf_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_TO_I8SCF_EVALUATE_TEST tests R8_TO_I8SCF.
        //
        //  Licensing:
        //
        //    I don't care what you do with this code.
        //
        //  Modified:
        //
        //    04 August 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        const int n = 19;

        Console.WriteLine("");
        Console.WriteLine("R8_TO_I8SCF_TEST:");

        double r = Math.PI;

        long[] a = new long[n + 1];

        r8_to_i8scf(r, n, a);

        typeMethods.i8vec_print(n + 1, a, "  SCF coefficients:");

        long[] p = new long[n + 1];
        long[] q = new long[n + 1];

        i8scf_evaluate(n, a, ref p, ref q);

        Console.WriteLine("");
        Console.WriteLine("  SCF numerators, denominators, ratios:");
        Console.WriteLine("");

        for (i = 0; i <= n; i++)
        {
            double t = p[i] / (double) q[i];
            Console.WriteLine("  " + i.ToString().PadLeft(2)
                                   + "  " + p[i].ToString().PadLeft(20)
                                   + "  " + q[i].ToString().PadLeft(20)
                                   + "  " + t.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void r8cf_evaluate(int n, double[] a, double[] b, ref double[] p, ref double[] q)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CF_EVALUATE evaluates a continued fraction with R8 entries.
        //
        //  Discussion:
        //
        //    For convenience, we omit the parentheses or multiline display.
        //
        //    CF = A(0) + B(1) / A(1) + B(2) / A(2) + ... A(N-1) + B(N) / A(N).
        //
        //  Licensing:
        //
        //    I don't care what you do with this code.
        //
        //  Modified:
        //
        //    06 August 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of terms.
        //
        //    Input, double A[N+1], B[N+1], the continued fraction
        //    terms.
        //
        //    Output, double P[N+1], Q[N+1], the numerators
        //    and denominators of the successive partial sums of the continued
        //    fraction.
        //
    {
        int i;

        for (i = 0; i <= n; i++)
        {
            switch (i)
            {
                case 0:
                    p[i] = a[i] * 1.0 + 0.0;
                    q[i] = a[i] * 0.0 + 1.0;
                    break;
                case 1:
                    p[i] = a[i] * p[i - 1] + b[i] * 1.0;
                    q[i] = a[i] * q[i - 1] + b[i] * 0.0;
                    break;
                default:
                    p[i] = a[i] * p[i - 1] + b[i] * p[i - 2];
                    q[i] = a[i] * q[i - 1] + b[i] * q[i - 2];
                    break;
            }
        }

    }

    public static void r8cf_evaluate_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CF_EVALUATE_TEST tests R8CF_EVALUATE.
        //
        //  Licensing:
        //
        //    I don't care what you do with this code.
        //
        //  Modified:
        //
        //    03 August 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        const int n = 20;

        Console.WriteLine("");
        Console.WriteLine("R8CF_EVALUATE_TEST:");

        double[] a = new double[n + 1];
        double[] b = new double[n + 1];

        a[0] = 3.0;
        for (i = 1; i <= n; i++)
        {
            a[i] = 6.0;
        }

        b[0] = 0.0;
        for (i = 1; i <= n; i++)
        {
            double t = 2 * i - 1;
            b[i] = t * t;
        }

        double[] p = new double[n + 1];
        double[] q = new double[n + 1];

        r8cf_evaluate(n, a, b, ref p, ref q);

        Console.WriteLine("");
        Console.WriteLine("  CF numerators, denominators, ratios:");
        Console.WriteLine("");

        for (i = 0; i <= n; i++)
        {
            Console.WriteLine("  " + i.ToString().PadLeft(2)
                                   + "  " + p[i].ToString("0.################").PadLeft(24)
                                   + "  " + q[i].ToString("0.################").PadLeft(24)
                                   + "  " + (p[i] / q[i]).ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void r8scf_evaluate(int n, double[] a, ref double[] p, ref double[] q)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SCF_EVALUATE evaluates a simple continued fraction with R8 entries.
        //
        //  Discussion:
        //
        //    The simple continued fraction with real coefficients is:
        //
        //      SCF = A(0) + 1 / ( A(1) + 1 / ( A(2) ... + 1 / A(N) ) )
        //
        //    This routine returns the N successive approximants P[i]/Q[i]
        //    to the value of the rational number represented by the continued
        //    fraction, with the value exactly equal to the final ratio C(N)/D(N).
        //
        //  Licensing:
        //
        //    I don't care what you do with this code.
        //
        //  Modified:
        //
        //    03 August 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of continued fraction
        //    coefficients.
        //
        //    Input, double A[N+1], the continued fraction coefficients.
        //
        //    Output, double P[N+1], Q[N+1], the numerators and
        //    denominators of the successive approximations.
        //
    {
        int i;

        for (i = 0; i <= n; i++)
        {
            switch (i)
            {
                case 0:
                    p[i] = a[i] * 1.0 + 0.0;
                    q[i] = a[i] * 0.0 + 1.0;
                    break;
                case 1:
                    p[i] = a[i] * p[i - 1] + 1.0;
                    q[i] = a[i] * q[i - 1] + 0.0;
                    break;
                default:
                    p[i] = a[i] * p[i - 1] + p[i - 2];
                    q[i] = a[i] * q[i - 1] + q[i - 2];
                    break;
            }
        }

    }

    public static void r8scf_evaluate_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SCF_EVALUATE_TEST tests R8SCF_EVALUATE.
        //
        //  Licensing:
        //
        //    I don't care what you do with this code.
        //
        //  Modified:
        //
        //    03 August 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a =
        {
            3.0, 7.0, 15.0, 1.0, 292.0,
            1.0, 1.0, 1.0, 2.0, 1.0,
            3.0, 1.0, 14.0, 2.0, 1.0,
            1.0, 2.0, 2.0, 2.0, 2.0
        };
        int i;
        const int n = 19;

        Console.WriteLine("");
        Console.WriteLine("R8SCF_EVALUATE_TEST:");

        double[] p = new double[n + 1];
        double[] q = new double[n + 1];

        r8scf_evaluate(n, a, ref p, ref q);

        Console.WriteLine("");
        Console.WriteLine("  SCF numerators, denominators, ratios:");
        Console.WriteLine("");

        for (i = 0; i <= n; i++)
        {
            Console.WriteLine("  " + i.ToString().PadLeft(2)
                                   + "  " + p[i].ToString("0.################").PadLeft(24)
                                   + "  " + q[i].ToString("0.################").PadLeft(24)
                                   + "  " + (p[i] / q[i]).ToString("0.################").PadLeft(24) + "");
        }

    }
}