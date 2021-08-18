using System;
using Burkardt.Values;

namespace TestValuesTest
{
    public static class FactorialTest
    {

        public static void i4_factorial_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4_FACTORIAL_TEST tests I4_FACTORIAL_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    21 March 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int fn = 0;
            int n = 0;
            int n_data;
            Console.WriteLine("");
            Console.WriteLine("I4_FACTORIAL_TEST:");
            Console.WriteLine("   I4_FACTORIAL_VALUES returns values of");
            Console.WriteLine("   the factorial function.");
            Console.WriteLine("");
            Console.WriteLine("      N        Factorial(N)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Factorial.i4_factorial_values(ref n_data, ref n, ref fn);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + n.ToString().PadLeft(6) + "  "
                                  + fn.ToString().PadLeft(12) + "");
            }
        }

        public static void i4_factorial2_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4_FACTORIAL2_TEST tests I4_FACTORIAL2_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    12 June 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int fn = 0;
            int n = 0;
            int n_data;
            Console.WriteLine("");
            Console.WriteLine("I4_FACTORIAL2_TEST:");
            Console.WriteLine("   I4_FACTORIAL2_VALUES return;s values of");
            Console.WriteLine("   the double factorial function.");
            Console.WriteLine("");
            Console.WriteLine("      N         DoubleFactorial(N)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Factorial.i4_factorial2_values(ref n_data, ref n, ref fn);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + n.ToString().PadLeft(6) + "  "
                                  + fn.ToString().PadLeft(12) + "");
            }
        }

        public static void i4_fall_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4_FALL_VALUES_TEST tests I4_FALL_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 December 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int fmn = 0;
            int m = 0;
            int n = 0;
            int n_data;
            Console.WriteLine("");
            Console.WriteLine("I4_FALL_VALUES_TEST:");
            Console.WriteLine("  I4_FALL_VALUES returns some exact values");
            Console.WriteLine("  of the integer falling factorial function:");
            Console.WriteLine("");
            Console.WriteLine("     M     N      I4_FALL(M,N)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Factorial.i4_fall_values(ref n_data, ref m, ref n, ref fmn);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + m.ToString().PadLeft(6) + "  "
                                  + n.ToString().PadLeft(6) + "  "
                                  + fmn.ToString().PadLeft(12) + "");
            }
        }

        public static void i4_rise_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4_RISE_VALUES_TEST tests I4_RISE_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 December 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int fmn = 0;
            int m = 0;
            int n = 0;
            int n_data;
            Console.WriteLine("");
            Console.WriteLine("I4_RISE_VALUES_TEST:");
            Console.WriteLine("  I4_RISE_VALUES returns some exact values");
            Console.WriteLine("  of the integer rising factorial function:");
            Console.WriteLine("");
            Console.WriteLine("     M     N      I4_RISE(M,N)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Factorial.i4_rise_values(ref n_data, ref m, ref n, ref fmn);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + m.ToString().PadLeft(6) + "  "
                                  + n.ToString().PadLeft(6) + "  "
                                  + fmn.ToString().PadLeft(12) + "");
            }
        }

        public static void r8_factorial_values_test()
            //****************************************************************************80
            //
            //  Purpose: 
            //
            //    R8_FACTORIAL_VALUES_TEST tests R8_FACTORIAL_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    09 February 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double fn = 0;
            int n = 0;
            int n_data;
            Console.WriteLine("");
            Console.WriteLine("R8_FACTORIAL_VALUES_TEST:");
            Console.WriteLine("  R8_FACTORIAL_VALUES stores values of");
            Console.WriteLine("  the factorial function (using double arithmetic).");
            Console.WriteLine("");
            Console.WriteLine("      N       Factorial(N)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Factorial.r8_factorial_values(ref n_data, ref n, ref fn);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + n.ToString().PadLeft(6) + "  "
                                  + fn.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void r8_factorial_log_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_FACTORIAL_LOG_VALUES_TEST tests R8_FACTORIAL_LOG_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    09 February 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double fn = 0;
            int n = 0;
            int n_data;
            Console.WriteLine("");
            Console.WriteLine("R8_FACTORIAL_LOG_VALUES_TEST:");
            Console.WriteLine("  R8_FACTORIAL_LOG_VALUES stores values of");
            Console.WriteLine("  the logarithm of the factorial function");
            Console.WriteLine("  (using real arithmetic).");
            Console.WriteLine("");
            Console.WriteLine("      N       Log(Factorial(N))");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Factorial.r8_factorial_log_values(ref n_data, ref n, ref fn);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + n.ToString().PadLeft(6) + "  "
                                  + fn.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void r8_factorial2_values_test()
            //****************************************************************************80
            //
            //  Purpose: 
            //
            //    R8_FACTORIAL2_VALUES_TEST tests R8_FACTORIAL2_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 February 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double f = 0;
            int n = 0;
            int n_data;
            Console.WriteLine("");
            Console.WriteLine("R8_FACTORIAL2_VALUES_TEST:");
            Console.WriteLine("  R8_FACTORIAL2_VALUES stores values of");
            Console.WriteLine("  the double factorial function (using double arithmetic).");
            Console.WriteLine("");
            Console.WriteLine("      N       F");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Factorial.r8_factorial2_values(ref n_data, ref n, ref f);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + n.ToString().PadLeft(6) + "  "
                                  + f.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void r8_fall_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_FALL_VALUES_TEST tests R8_FALL_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 December 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double f = 0;
            int n = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("R8_FALL_VALUES_TEST:");
            Console.WriteLine("  R8_FALL_VALUES returns some exact values");
            Console.WriteLine("  of the falling factorial function:");
            Console.WriteLine("");
            Console.WriteLine("     X     N      R8_FALL(X,N)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Factorial.r8_fall_values(ref n_data, ref x, ref n, ref f);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString().PadLeft(8) + "  "
                                  + n.ToString().PadLeft(8) + "  "
                                  + f.ToString().PadLeft(12) + "");
            }
        }

        public static void r8_rise_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_RISE_VALUES_TEST tests R8_RISE_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 December 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double f = 0;
            int n = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("R8_RISE_VALUES_TEST:");
            Console.WriteLine("  R8_RISE_VALUES returns some exact values");
            Console.WriteLine("  of the  rising factorial function:");
            Console.WriteLine("");
            Console.WriteLine("     X     N      R8_RISE(X,N)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Factorial.r8_rise_values(ref n_data, ref x, ref n, ref f);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString().PadLeft(8) + "  "
                                  + n.ToString().PadLeft(8) + "  "
                                  + f.ToString().PadLeft(12) + "");
            }
        }

        public static void subfactorial_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SUBFACTORIAL_VALUES_TEST tests SUBFACTORIAL_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    21 March 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int fn = 0;
            int n = 0;
            int n_data;
            Console.WriteLine("");
            Console.WriteLine("SUBFACTORIAL_VALUES_TEST:");
            Console.WriteLine("  SUBFACTORIAL_VALUES returns values of");
            Console.WriteLine("  the subfactorial function.");
            Console.WriteLine("");
            Console.WriteLine("      N       Subfactorial[N]");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Factorial.subfactorial_values(ref n_data, ref n, ref fn);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + n.ToString().PadLeft(6) + "  "
                                  + fn.ToString().PadLeft(12) + "");
            }
        }

    }
}