using System;
using Burkardt.AppliedStatistics;
using Burkardt.RankingNS;
using Burkardt.Types;

namespace ComboTest
{
    partial class Program
    {
        static void r8_choose_test ( )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_CHOOSE_TEST tests R8_CHOOSE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    26 July 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double cnk;
            int k;
            int n;

            Console.WriteLine("");
            Console.WriteLine("R8_CHOOSE_TEST");
            Console.WriteLine("  R8_CHOOSE evaluates C(N,K).");
            Console.WriteLine("");
            Console.WriteLine("         N         K       CNK");
 
            for ( n = 0; n <= 5; n++ )
            {
                Console.WriteLine("");
                for ( k = 0; k <= n; k++ )
                {
                    cnk = typeMethods.r8_choose ( n, k );
                    Console.WriteLine(n.ToString().PadLeft(10) + "  "
                        + k.ToString().PadLeft(8) + "  "
                        + cnk.ToString().PadLeft(14) + "");
                }
            }
 
            return;
        }
        //****************************************************************************80

        static void r8_gamma_log_test ( )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_GAMMA_LOG_TEST tests R8_GAMMA_LOG.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    23 April 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double fx1 = 0;
            double fx2;
            int n_data;
            double x = 0;

            Console.WriteLine("");
            Console.WriteLine("R8_GAMMA_LOG_TEST:");
            Console.WriteLine("   R8_GAMMA_LOG evaluates the logarithm of the Gamma function.");
            Console.WriteLine("");
            Console.WriteLine("      X            GAMMA_LOG(X)     R8_GAMMA_LOG(X)");
            Console.WriteLine("");

            n_data = 0;

            for ( ; ; )
            {
                Algorithms.gamma_log_values ( ref n_data, ref x, ref fx1 );

                if ( n_data == 0 )
                {
                    break;
                }
                fx2 = typeMethods.r8_gamma_log ( x );
                Console.WriteLine("  " + x.ToString().PadLeft(12)
                    + "  " + fx1.ToString("0.################").PadLeft(24)
                    + "  " + fx2.ToString("0.################").PadLeft(24) + "");
            }
        }
    }
}