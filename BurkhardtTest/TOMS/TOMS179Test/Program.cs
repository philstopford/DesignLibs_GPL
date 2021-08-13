using System;
using Burkardt.AppliedStatistics;
using Burkardt.Probability;
using Burkardt.TOMSNS;

namespace TOMS179Test
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for TOMS179_TEST.
            //
            //  Discussion:
            //
            //    TOMS179_TEST tests the TOMS179 library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    30 January 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("TOMS179_TEST:");
            Console.WriteLine("  Test the TOMS179 library.");

            test01();
            test02();

            Console.WriteLine("");
            Console.WriteLine("TOMS179_TEST:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void test01()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 demonstrates the use of ALOGAM.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    30 January 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double fx = 0;
            double fx2;
            int ifault = 0;
            int n_data;
            double x = 0;

            Console.WriteLine("");
            Console.WriteLine("TEST01:");
            Console.WriteLine("  ALOGAM computes the logarithm of the ");
            Console.WriteLine("  Gamma function.  We compare the result");
            Console.WriteLine("  to tabulated values.");
            Console.WriteLine("");
            Console.WriteLine("          X                     "
                              + "FX                        FX2");
            Console.WriteLine("                                "
                              + "(Tabulated)               (ALOGAM)                DIFF");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                Burkardt.TestValues.Gamma.gamma_log_values(ref n_data, ref x, ref fx);

                if (n_data == 0)
                {
                    break;
                }

                fx2 = Algorithms.alogam(x, ref ifault);

                Console.WriteLine("  " + x.ToString("0.################").PadLeft(24)
                                       + "  " + fx.ToString("0.################").PadLeft(24)
                                       + "  " + fx2.ToString("0.################").PadLeft(24)
                                       + "  " + Math.Abs(fx - fx2).ToString("0.##########").PadLeft(10) + "");
            }

            return;
        }

        static void test02()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST02 demonstrates the use of MDBETA.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    30 January 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double fx = 0;
            double fx2;
            int ier = 0;
            int n_data;
            double p = 0;
            double q = 0;
            double x = 0;

            Console.WriteLine("");
            Console.WriteLine("TEST02:");
            Console.WriteLine("  MDBETA estimates the value of th modified Beta function.");
            Console.WriteLine("  Compare with tabulated values.");
            Console.WriteLine("");
            Console.WriteLine("         X         P         Q         "
                              + "Beta                       Beta                  DIFF");
            Console.WriteLine("                                       "
                              + "(Tabulated)                (MDBETA)");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                Beta.beta_cdf_values(ref n_data, ref p, ref q, ref x, ref fx);

                if (n_data == 0)
                {
                    break;
                }

                fx2 = TOMS.mdbeta(x, p, q, ref ier);

                Console.WriteLine("  " + x.ToString("0.####").PadLeft(8)
                                       + "  " + p.ToString("0.####").PadLeft(8)
                                       + "  " + q.ToString("0.####").PadLeft(8)
                                       + "  " + fx.ToString("0.################").PadLeft(24)
                                       + "  " + fx2.ToString("0.################").PadLeft(24)
                                       + "  " + Math.Abs(fx - fx2).ToString("0.####").PadLeft(10) + "");
            }

        }
    }
}