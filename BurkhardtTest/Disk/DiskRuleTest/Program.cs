﻿using System;

namespace DiskRuleTest
{
    using QuadratureRule = Burkardt.Disk.QuadratureRule;
    
    class Program
    {
        static void Main(string[] args)
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

        static void test01()

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
            double area;
            int d;
            int ex;
            int ey;
            double[] exact =  {
                9.0,
                9.0, 18.0,
                117.0 / 4.0, 18.0, 225.0 / 4.0,
                279.0 / 4.0, 117.0 / 2.0, 225.0 / 4.0, 387.0 / 2.0,
                1773.0 / 8.0, 279.0 / 2.0, 1341.0 / 8.0, 387.0 / 2.0, 5769.0 / 8.0
            }
            ;
            int i;
            int j;
            int k;
            int nr = 4;
            int nt = 8;
            double q;
            const double r8_pi = 3.141592653589793;
            double rc;
            double s;
            double[] w;
            double[] x;
            double xc;
            double[] y;
            double yc;

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
            xc = 1.0;
            yc = 2.0;
            rc = 3.0;
            //
            //  Put in the factor of PI in the exact values.
            //
            for (k = 0; k < 15; k++)
            {
                exact[k] = exact[k] * r8_pi;
            }

            //
            //  Compute the quadrature rule.
            //
            w = new double[nr * nt];
            x = new double[nr * nt];
            y = new double[nr * nt];

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
                for (ex = d; 0 <= ex; ex--)
                {
                    ey = d - ex;

                    s = 0.0;
                    for (j = 0; j < nt; j++)
                    {
                        for (i = 0; i < nr; i++)
                        {
                            s = s + w[i + j * nr] * Math.Pow(x[i + j * nr], ex) * Math.Pow(y[i + j * nr], ey);
                        }
                    }

                    area = r8_pi * rc * rc;

                    q = area * s;

                    Console.WriteLine("  " + ex.ToString().PadLeft(2)
                        + "  " + ey.ToString().PadLeft(2)
                        + "  " + exact[k].ToString().PadLeft(14)
                        + "  " + q.ToString().PadLeft(14) + "");

                    k = k + 1;
                }
            }
        }
    }
}