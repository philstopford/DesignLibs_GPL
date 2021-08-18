using System;
using Burkardt.Values;

namespace TestValuesTest
{
    public class StudentTest
    {
        public static void student_cdf_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    STUDENT_CDF_VALUES_TEST tests STUDENT_CDF_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    03 November 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double c = 0;
            double fx = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("STUDENT_CDF_VALUES_TEST:");
            Console.WriteLine("  STUDENT_CDF_VALUES returns values of");
            Console.WriteLine("  the Student T Cumulative Density Function.");
            Console.WriteLine("");
            Console.WriteLine("      C     X       CDF(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Student.student_cdf_values(ref n_data, ref c, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + c.ToString().PadLeft(16) + "  "
                                  + x.ToString().PadLeft(16) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void student_noncentral_cdf_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    STUDENT_NONCENTRAL_CDF_VALUES_TEST tests STUDENT_NONCENTRAL_CDF_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    08 February 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int df = 0;
            double fx = 0;
            double lambda = 0;
            ;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("STUDENT_NONCENTRAL_CDF_VALUES_TEST:");
            Console.WriteLine("  STUDENT_NONCENTRAL_CDF_VALUES returns values of");
            Console.WriteLine("  the noncentral Student T Cumulative Density Function.");
            Console.WriteLine("");
            Console.WriteLine("    DF     LAMBDA        X        CDF");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Student.student_noncentral_cdf_values(ref n_data, ref df, ref lambda, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + df.ToString().PadLeft(6) + "  "
                                  + lambda.ToString().PadLeft(8) + "  "
                                  + x.ToString().PadLeft(8) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }
    }
}