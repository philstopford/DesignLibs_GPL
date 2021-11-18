using System;
using Burkardt.SphereNS;

namespace SphereLebedevRuleTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SPHERE_LEBEDEV_RULE_TEST.
        //
        //  Discussion:
        //
        //    SPHERE_LEBEDEV_RULE_TEST tests the SPHERE_LEBEDEV_RULE library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("SPHERE_LEBEDEV_RULE_TEST");

        Console.WriteLine("  Test the SPHERE_LEBEDEV_RULE library.");

        test01();
        test02();

        Console.WriteLine("");
        Console.WriteLine("SPHERE_LEBEDEV_RULE_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests AVAILABLE_TABLE, ORDER_TABLE, PRECISION_TABLE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //   12 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int available;
        int order;
        int precision;
        int rule;
        int rule_max = 65;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  List Lebedev rule properties.");
        Console.WriteLine("");
        Console.WriteLine("  Rule Avail Order  Prec");
        Console.WriteLine("");
        for (rule = 1; rule <= rule_max; rule++)
        {
            available = LebedevRule.available_table(rule);
            order = LebedevRule.order_table(rule);
            precision = LebedevRule.precision_table(rule);
            Console.WriteLine("  " + rule.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + available.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + order.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + precision.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "");
        }
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests the SPHERE_LEBEDEV_RULE functions.
        //
        //  Modified:
        //
        //    13 September 2010
        //
        //  Author:
        //
        //    Dmitri Laikov
        //
        //  Reference:
        //
        //    Vyacheslav Lebedev, Dmitri Laikov,
        //    A quadrature formula for the sphere of the 131st
        //    algebraic order of accuracy,
        //    Russian Academy of Sciences Doklady Mathematics,
        //    Volume 59, Number 3, 1999, pages 477-481.
        //
    {
        int nmax = 65;
        int mmax = (nmax * 2 + 3) * (nmax * 2 + 3) / 3;

        double alpha = 0;
        int available;
        double beta = 0;
        double err;
        double err_max;
        int i;
        double integral_exact;
        double integral_approx;
        int j;
        int k;
        int m;
        int n;
        int order;
        double[] w;
        double[] x;
        double[] y;
        double[] z;
        /*
        static double[] s = new double[nmax+2];
        static double[] xn = new double[mmax*(nmax+1)];
        static double[] yn = new double[mmax*(nmax+1)];
        static double[] zn = new double[mmax*(nmax+1)];
        */
        double[] s = new double[nmax + 2];
        double[] xn = new double[mmax * (nmax + 1)];
        double[] yn = new double[mmax * (nmax + 1)];
        double[] zn = new double[mmax * (nmax + 1)];

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  Generate each available rule and test for accuracy.");

        for (n = 1; n <= nmax; n++)
        {
            available = LebedevRule.available_table(n);

            switch (available)
            {
                case 1:
                {
                    order = LebedevRule.order_table(n);

                    w = new double[order];
                    x = new double[order];
                    y = new double[order];
                    z = new double[order];

                    LebedevRule.ld_by_order(order, ref x, ref y, ref z, ref w);

                    s[0] = 1.0;
                    for (k = 1; k <= n + 1; k++)
                    {
                        s[k] = (2 * k - 1) * s[k - 1];
                    }

                    //
                    //  For each abscissa X(M), compute the values 1, X(M)^2, X(M)^4, ..., X(M)^2*N.
                    //
                    for (m = 0; m < order; m++)
                    {
                        xn[m * (n + 1)] = 1.0;
                        yn[m * (n + 1)] = 1.0;
                        zn[m * (n + 1)] = 1.0;
                        for (k = 1; k <= n; k++)
                        {
                            xn[k + m * (n + 1)] = xn[k - 1 + m * (n + 1)] * x[m] * x[m];
                            yn[k + m * (n + 1)] = yn[k - 1 + m * (n + 1)] * y[m] * y[m];
                            zn[k + m * (n + 1)] = zn[k - 1 + m * (n + 1)] * z[m] * z[m];
                        }
                    }

                    err_max = 0.0;
                    for (i = 0; i <= n; i++)
                    {
                        for (j = 0; j <= n - i; j++)
                        {
                            k = n - i - j;
                            //
                            //  Apply Lebedev rule to x^2i y^2j z^2k.
                            //
                            integral_approx = 0.0;
                            for (m = 0; m < order; m++)
                            {
                                integral_approx += w[m] * xn[i + m * (n + 1)] * yn[j + m * (n + 1)] *
                                                   zn[k + m * (n + 1)];
                            }

                            //
                            //  Compute exact value of integral (aside from factor of 4 pi!).
                            //
                            integral_exact = s[i] * s[j] * s[k] / s[1 + i + j + k];
                            //
                            //  Record the maximum error for this rule.
                            //
                            err = Math.Abs((integral_approx - integral_exact) / integral_exact);
                            if (err_max < err)
                            {
                                err_max = err;
                            }
                        }
                    }

                    Console.WriteLine("");
                    ;
                    Console.WriteLine("  Order = " + order.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                                   + "  LMAXW = " + LebedevRule.precision_table(n)
                                                   + "  max error = " + err_max + "");
                    switch (order)
                    {
                        //
                        //  Convert (X,Y,Z) to (Theta,Phi) and print the data.
                        //
                        case <= 50:
                        {
                            for (m = 0; m < order; m++)
                            {
                                LebedevRule.xyz_to_tp(x[m], y[m], z[m], ref alpha, ref beta);
                                Console.WriteLine("  " + alpha.ToString("0.###############").PadLeft(20)
                                                       + "  " + beta.ToString("0.###############").PadLeft(20)
                                                       + "  " + w[m].ToString("0.###############").PadLeft(20) + "");
                            }

                            break;
                        }
                    }

                    break;
                }
            }
        }
    }
}