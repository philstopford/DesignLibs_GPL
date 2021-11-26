using System;
using System.Globalization;

namespace DiskRuleTest;

using QuadratureRule = Burkardt.Disk.QuadratureRule;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for DISK_RULE_TEST.
        //
        //  Discussion:
        //
        //    DISK01_RULE_TEST tests the DISK_RULE library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 March 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("DISK_RULE:");
        Console.WriteLine("  Test the DISK_RULE library.");

        test01();

        Console.WriteLine("");
        Console.WriteLine("DISK_RULE:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests DISK_RULE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 March 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int d;
        double[] exact =  {
                9.0,
                9.0, 18.0,
                117.0 / 4.0, 18.0, 225.0 / 4.0,
                279.0 / 4.0, 117.0 / 2.0, 225.0 / 4.0, 387.0 / 2.0,
                1773.0 / 8.0, 279.0 / 2.0, 1341.0 / 8.0, 387.0 / 2.0, 5769.0 / 8.0
            }
            ;
        int k;
        const int nr = 4;
        const int nt = 8;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  DISK_RULE can compute a rule Q(f) for the general disk");
        Console.WriteLine("  with center (XC,YC) and radius RC,");
        Console.WriteLine("  using NT equally spaced angles and NR radial distances.");
        Console.WriteLine("");
        Console.WriteLine("  NT = " + nt + "");
        Console.WriteLine("  NR = " + nr + "");
        Console.WriteLine("");
        Console.WriteLine("  Estimate integrals I(f) where f = x^ex * y^ey.");
        //
        //  Define the general disk.
        //
        const double xc = 1.0;
        const double yc = 2.0;
        const double rc = 3.0;
        //
        //  Put in the factor of PI in the exact values.
        //
        for (k = 0; k < 15; k++)
        {
            exact[k] *= Math.PI;
        }

        //
        //  Compute the quadrature rule.
        //
        double[] w = new double[nr * nt];
        double[] x = new double[nr * nt];
        double[] y = new double[nr * nt];

        QuadratureRule.disk_rule(nr, nt, xc, yc, rc, ref w, ref x, ref y);
        //
        //  Apply it to integrands.
        //
        Console.WriteLine("");
        Console.WriteLine("  EX    EY      I(f)            Q(f)");
        Console.WriteLine("");
        //
        //  Specify a monomial.
        //
        k = 0;

        for (d = 0; d <= 4; d++)
        {
            int ex;
            for (ex = d; 0 <= ex; ex--)
            {
                int ey = d - ex;

                double s = 0.0;
                int j;
                for (j = 0; j < nt; j++)
                {
                    int i;
                    for (i = 0; i < nr; i++)
                    {
                        s += w[i + j * nr] * Math.Pow(x[i + j * nr], ex) * Math.Pow(y[i + j * nr], ey);
                    }
                }

                double area = Math.PI * rc * rc;

                double q = area * s;

                Console.WriteLine("  " + ex.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                       + "  " + ey.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                       + "  " + exact[k].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + q.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

                k += 1;
            }
        }
    }
}