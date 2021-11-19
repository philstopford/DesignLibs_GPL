using System;
using Burkardt.Function;
using Burkardt.RationalNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace SubsetTestNS;

public static class RationalTest
{
    public static void rat_add_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RAT_ADD_TEST tests RAT_ADD.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 October 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int abot = 4;
        int atop = 3;
        int bbot = 7;
        int btop = 10;
        int cbot = 0;
        int ctop = 0;
        bool error = false;

        Console.WriteLine("");
        Console.WriteLine("RAT_ADD_TEST");
        Console.WriteLine("  RAT_ADD adds two rationals.");

        Rational.rat_add(atop, abot, btop, bbot, ref ctop, ref cbot, ref error);

        Console.WriteLine("");
        Console.WriteLine("  A = " + atop + "/" + abot + "");
        Console.WriteLine("  B = " + btop + "/" + bbot + "");
        Console.WriteLine("  C = A + B = " + ctop + "/" + cbot + "");
    }

    public static void rat_div_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RAT_DIV_TEST tests RAT_DIV.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int abot = 4;
        int atop = 3;
        int bbot = 7;
        int btop = 10;
        int cbot = 0;
        int ctop = 0;
        bool error = false;

        Console.WriteLine("");
        Console.WriteLine("RAT_DIV_TEST");
        Console.WriteLine("  RAT_DIV divides two rationals.");

        Rational.rat_div(atop, abot, btop, bbot, ref ctop, ref cbot, ref error);

        Console.WriteLine("");
        Console.WriteLine("  A = " + atop + "/" + abot + "");
        Console.WriteLine("  B = " + btop + "/" + bbot + "");
        Console.WriteLine("  C = A / B = " + ctop + "/" + cbot + "");
    }

    public static void rat_farey_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RAT_FAREY_TEST tests RAT_FAREY.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int FRAC_MAX = 20;

        int[] a = new int[FRAC_MAX];
        int[] b = new int[FRAC_MAX];
        int frac_num = 0;
        int i;
        int ihi;
        int ilo;
        int n;

        Console.WriteLine("");
        Console.WriteLine("RAT_FAREY_TEST");
        Console.WriteLine("  RAT_FAREY computes a row of the Farey fraction table.");

        for (n = 1; n <= 7; n++)
        {
            Rational.rat_farey(n, FRAC_MAX, ref frac_num, ref a, ref b);

            Console.WriteLine("");
            Console.WriteLine("  Row " + n + "");
            Console.WriteLine("  Number of fractions: " + frac_num + "");
            Console.WriteLine("");

            for (ilo = 0; ilo < frac_num; ilo += 20)
            {
                ihi = Math.Min(ilo + 20 - 1, frac_num - 1);
                string cout = "";
                for (i = ilo; i <= ihi; i++)
                {
                    cout += a[i].ToString().PadLeft(3);
                }

                Console.WriteLine(cout);
                cout = "";
                for (i = ilo; i <= ihi; i++)
                {
                    cout += b[i].ToString().PadLeft(3);
                }

                Console.WriteLine(cout);
            }

        }
    }

    public static void rat_farey2_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RAT_FAREY2_TEST tests RAT_FAREY2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N_MAX = 4;
        int TWO_POWER_MAX = 16;

        int[] a = new int[TWO_POWER_MAX + 1];
        int[] b = new int[TWO_POWER_MAX + 1];
        int i;
        int ihi;
        int ilo;
        int n;
        int two_power;

        Console.WriteLine("");
        Console.WriteLine("RAT_FAREY2_TEST");
        Console.WriteLine("  RAT_FAREY2 computes a row of the Farey fraction table.");

        for (n = 0; n <= N_MAX; n++)
        {
            Rational.rat_farey2(n, ref a, ref b);

            Console.WriteLine("");
            Console.WriteLine("  Row " + n + 1 + "");

            two_power = (int)Math.Pow(2, n);

            for (ilo = 0; ilo <= two_power; ilo += 20)

            {
                ihi = Math.Min(ilo + 20 - 1, two_power);
                Console.WriteLine("");
                string cout = "";
                for (i = ilo; i <= ihi; i++)
                {
                    cout += a[i].ToString().PadLeft(3);
                }

                Console.WriteLine(cout);
                cout = "";
                for (i = ilo; i <= ihi; i++)
                {
                    cout += b[i].ToString().PadLeft(3);
                }

                Console.WriteLine(cout);

            }

        }
    }

    public static void rat_mul_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RAT_MUL_TEST tests RAT_MUL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int abot = 4;
        int atop = 3;
        int bbot = 7;
        int btop = 10;
        int cbot = 0;
        int ctop = 0;
        bool error = false;

        Console.WriteLine("");
        Console.WriteLine("RAT_MUL_TEST");
        Console.WriteLine("  RAT_MUL multiplies two rationals.");

        Rational.rat_mul(atop, abot, btop, bbot, ref ctop, ref cbot, ref error);

        Console.WriteLine("");
        Console.WriteLine("  A = " + atop + "/" + abot + "");
        Console.WriteLine("  B = " + btop + "/" + bbot + "");
        Console.WriteLine("  C = A * B = " + ctop + "/" + cbot + "");
    }

    public static void rat_normalize_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RAT_NORMALIZE_TEST tests RAT_NORMALIZE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 June 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int a1;
        int a2;
        int b1;
        int b2;
        int i;
        int[] rat_top = { 3, 1, 20, 8, -10, 9, -11 };
        int[] rat_bot = { 4, 1000, 1, 4, 7, -15, -11 };

        Console.WriteLine("");
        Console.WriteLine("RAT_NORMALIZE_TEST");
        Console.WriteLine("  RAT_NORMALIZE normalizes a rational.");

        Console.WriteLine("");
        Console.WriteLine("       A       B         A       B");
        Console.WriteLine("                         Normalized");
        Console.WriteLine("");

        for (i = 0; i < 7; i++)
        {
            a1 = rat_top[i];
            b1 = rat_bot[i];
            a2 = a1;
            b2 = b1;
            Rational.rat_normalize(ref a2, ref b2);
            Console.WriteLine("  " + a1.ToString().PadLeft(6)
                                   + "  " + b1.ToString().PadLeft(6)
                                   + "  "
                                   + "  " + a2.ToString().PadLeft(6)
                                   + "  " + b2.ToString().PadLeft(6) + "");
        }
    }

    public static void rat_sum_formula_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RAT_SUM_FORMULA_TEST tests RAT_SUM_FORMULA.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 6;

        int[] a = new int[(N + 1) * (N + 1)];
        int[] b = new int[(N + 1) * (N + 1)];

        Console.WriteLine("");
        Console.WriteLine("RAT_SUM_FORMULA_TEST");
        Console.WriteLine("  RAT_SUM_FORMULA computes the coefficients for the");
        Console.WriteLine("  formulas for the sums of powers of integers.");

        Rational.rat_sum_formula(N, ref a, ref b);

        Rational.ratmat_print(N + 1, N + 1, a, b, "  Power Sum Coefficients:");

    }

    public static void rat_to_cfrac_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RAT_TO_CFRAC_TEST tests RAT_TO_CFRAC.
        //
        //  Discussion:
        //
        //    Compute the continued fraction form of 4096/15625.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 October 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int M = 10;

        int[] a = new int[M];
        int bot = 15625;
        bool error = false;
        int i;
        int n = 0;
        int[] p = new int[M];
        int[] q = new int[M];
        int top = 4096;

        Console.WriteLine("");
        Console.WriteLine("RAT_TO_CFRAC_TEST");
        Console.WriteLine("  RAT_TO_CFRAC fraction => continued fraction,");
        Console.WriteLine("");
        Console.WriteLine("  Regular fraction is " + top + "/" + bot + "");

        Rational.rat_to_cfrac(top, bot, M, ref n, ref a, ref error);

        typeMethods.i4vec1_print(n, a, "  Continued fraction coefficients:");

        Fraction.cfrac_to_rat(n, a, ref p, ref q);

        Console.WriteLine("");
        Console.WriteLine("  The continued fraction convergents.");
        Console.WriteLine("  The last row contains the value of the continued");
        Console.WriteLine("  fraction, written as a common fraction.");
        Console.WriteLine("");
        Console.WriteLine("  I, P(I), Q(I), P(I)/Q(I)");
        Console.WriteLine("");

        for (i = 0; i < n; i++)
        {
            Console.WriteLine(i.ToString().PadLeft(3) + "  "
                                                      + p[i].ToString().PadLeft(6) + "  "
                                                      + q[i].ToString().PadLeft(6) + "  "
                                                      + (p[i] / (double)q[i]).ToString().PadLeft(14) + "");
        }
    }

    public static void rat_to_dec_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RAT_TO_DEC_TEST tests RAT_TO_DEC.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 November 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int exponent = 0;
        int i;
        int mantissa = 0;
        double r1;
        double r2;
        double r3;
        int rat_bot = 0;
        int rat_bot2 = 0;
        int rat_top = 0;
        int rat_top2 = 0;
        int seed;

        Console.WriteLine("");
        Console.WriteLine("RAT_TO_DEC_TEST");
        Console.WriteLine("  RAT_TO_DEC fraction => decimal,");
        Console.WriteLine("");
        Console.WriteLine("  In this test, choose the top and bottom");
        Console.WriteLine("  of a rational at random, and compute the");
        Console.WriteLine("  equivalent real number.");
        Console.WriteLine("");
        Console.WriteLine("  Then convert to decimal, and the equivalent real.");
        Console.WriteLine("");
        Console.WriteLine("  Then convert back to rational and the equivalent real.");

        seed = 123456789;

        for (i = 1; i <= 10; i++)
        {
            rat_top = UniformRNG.i4_uniform_ab(-1000, 1000, ref seed);
            rat_bot = UniformRNG.i4_uniform_ab(1, 1000, ref seed);

            r1 = rat_top / (double)rat_bot;

            Rational.rat_to_dec(rat_top, rat_bot, ref mantissa, ref exponent);
            r2 = mantissa * Math.Pow(10.0, exponent);

            typeMethods.dec_to_rat(mantissa, exponent, ref rat_top2, ref rat_bot2);
            r3 = rat_top2 / (double)rat_bot2;

            Console.WriteLine("");
            Console.WriteLine("  " + r1 + " = " + rat_top + "/" + rat_bot + "");
            Console.WriteLine("  " + r2 + " = " + mantissa + "*10^(" + exponent + ")");
            Console.WriteLine("  " + r3 + " = " + rat_top2 + "/" + rat_bot2 + "");
        }
    }

    public static void rat_to_r8_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RAT_TO_R8_TEST tests RAT_TO_R8.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 October 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int a = 0;
        int b = 0;
        int i;
        int ndig = 4;
        double r;
        double r2;
        int seed;

        Console.WriteLine("");
        Console.WriteLine("RAT_TO_R8_TEST");
        Console.WriteLine("  RAT_TO_R8 converts a rational to a real number.");
        Console.WriteLine("");
        Console.WriteLine("  The maximum number of digits allowed is " + ndig + "");

        seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("     R   =>  A / B  =>  R2");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            r = UniformRNG.r8_uniform_01(ref seed);
            r = 10.0 * (r - 0.25);

            typeMethods.r8_to_rat(r, ndig, ref a, ref b);
            r2 = Rational.rat_to_r8(a, b);

            Console.WriteLine("  "
                              + r.ToString().PadLeft(10) + "  "
                              + a.ToString().PadLeft(6) + "  "
                              + b.ToString().PadLeft(6) + "  "
                              + r2.ToString().PadLeft(10) + "");
        }
    }

    public static void rat_to_s_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RAT_TO_S_TEST tests RAT_TO_S.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 June 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int a;
        int b;
        int i;
        int[] rat_top = { 3, 1, 20, 8, -10, 9, -11 };
        int[] rat_bot = { 4, 1000, 1, 4, 7, -15, -11 };
        string s;

        Console.WriteLine("");
        Console.WriteLine("RAT_TO_S_TEST");
        Console.WriteLine("  RAT_TO_S converts a rational to a string.");

        Console.WriteLine("");
        Console.WriteLine("       A       B      S");
        Console.WriteLine("");

        for (i = 0; i < 7; i++)
        {
            a = rat_top[i];
            b = rat_bot[i];
            s = Rational.rat_to_s(a, b);
            Console.WriteLine("  " + a.ToString().PadLeft(6)
                                   + "  " + b.ToString().PadLeft(6)
                                   + "      " + s + "");
        }
    }

    public static void rat_width_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RAT_WIDTH_TEST tests RAT_WIDTH.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N_TEST = 17;

        int a;
        int[] a_test =
        {
            1000, 1000, 1000, 1000, 1000, 1, -1, -10, -100, -1000,
            1, 10, 100, 1000, 10000, 17, 4000000
        };
        int b;
        int[] b_test =
        {
            3, 40, 500, 6000, 70000, 1, 200, 200, 200, 200,
            -200, -200, -200, -200, -200, 3000, 4000000
        };
        int i;
        int width;

        Console.WriteLine("");
        Console.WriteLine("RAT_WIDTH_TEST");
        Console.WriteLine("  RAT_WIDTH determines the \"width\" of a rational.");
        Console.WriteLine("");
        Console.WriteLine("  Top    Bottom  Width");
        Console.WriteLine("");

        for (i = 0; i < N_TEST; i++)
        {
            a = a_test[i];
            b = b_test[i];

            width = Rational.rat_width(a, b);

            Console.WriteLine("  "
                              + a.ToString().PadLeft(8) + "  "
                              + b.ToString().PadLeft(8) + "  "
                              + width.ToString().PadLeft(8) + "");
        }
    }

    public static void ratmat_det_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RATMAT_DET_TEST tests RATMAT_DET.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N3 = 3;

        int[] a3 = new int[N3 * N3];
        int[] b3 = new int[N3 * N3];
        bool error = false;
        int i;
        int idbot = 0;
        int idtop = 0;
        int j;
        int k;

        Console.WriteLine("");
        Console.WriteLine("RATMAT_DET_TEST");
        Console.WriteLine("  RATMAT_DET: determinant of a rational matrix.");
        Console.WriteLine("");

        k = 0;
        for (i = 0; i < N3; i++)
        {
            for (j = 0; j < N3; j++)
            {
                k += 1;
                a3[i + j * N3] = k;
            }
        }

        for (i = 0; i < N3; i++)
        {
            for (j = 0; j < N3; j++)
            {
                b3[i + j * N3] = 1;
            }
        }

        Rational.ratmat_print(N3, N3, a3, b3, "  The 123/456/789 matrix:");

        Rational.ratmat_det(N3, a3, b3, ref idtop, ref idbot, ref error);

        Console.WriteLine("");
        Console.WriteLine("  Determinant of the 123/456/789 matrix = "
                          + idtop + "/" + idbot + "");

        for (i = 0; i < N3; i++)
        {
            for (j = 0; j < N3; j++)
            {
                a3[i + j * N3] = 1;
            }
        }

        for (i = 0; i < N3; i++)
        {
            for (j = 0; j < N3; j++)
            {
                b3[i + j * N3] = i + j + 2;
            }
        }

        Rational.ratmat_print(N3, N3, a3, b3, "  The Hilbert matrix:");

        Rational.ratmat_det(N3, a3, b3, ref idtop, ref idbot, ref error);

        Console.WriteLine("");
        Console.WriteLine("  Determinant of the Hilbert matrix = "
                          + idtop + "/" + idbot + "");

        for (i = 0; i < N3; i++)
        {
            for (j = 0; j < N3; j++)
            {
                if (i == j)
                {
                    a3[i + j * N3] = 2;
                }
                else if (i == j + 1 || i == j - 1)
                {
                    a3[i + j * N3] = -1;
                }
                else
                {
                    a3[i + j * N3] = 0;
                }
            }
        }

        for (i = 0; i < N3; i++)
        {
            for (j = 0; j < N3; j++)
            {
                b3[i + j * N3] = 1;
            }
        }

        Rational.ratmat_print(N3, N3, a3, b3, "  The -1 2 -1 matrix:");

        Rational.ratmat_det(N3, a3, b3, ref idtop, ref idbot, ref error);

        Console.WriteLine("");
        Console.WriteLine("  Determinant of the -1,2,-1 matrix = "
                          + idtop + "/" + idbot + "");
    }

    public static void ratmat_print_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RATMAT_PRINT_TEST tests RATMAT_PRINT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] a = new int[4 * 3];
        int[] b = new int[4 * 3];
        int i;
        int j;
        int m = 4;
        int n = 3;

        Console.WriteLine("");
        Console.WriteLine("RATMAT_PRINT_TEST");
        Console.WriteLine("  RATMAT_PRINT prints a rational matrix.");
        Console.WriteLine("");

        for (i = 0; i < m; i++)
        {
            for (j = 0; j < n; j++)
            {
                a[i + j * m] = 1;
            }
        }

        for (i = 0; i < m; i++)
        {
            for (j = 0; j < n; j++)
            {
                b[i + j * m] = i + j + 2;
            }
        }

        Rational.ratmat_print(m, n, a, b, "  The Hilbert matrix:");

    }

}