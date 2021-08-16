using System;
using Burkardt.Function;

namespace PolPakTest
{
    public static class gudermannianTest
    {
        public static void gud_test ( )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GUD_TEST tests GUD.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    23 May 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double fx = 0;
            double fx2 = 0;
            int n_data = 0;
            double x = 0;

            Console.WriteLine("");
            Console.WriteLine("GUD_TEST:");
            Console.WriteLine("  GUD evaluates the Gudermannian function.");
            Console.WriteLine("");
            Console.WriteLine("     X      Exact F       GUD(X)");
            Console.WriteLine("");

            n_data = 0;

            for ( ; ; )
            {
                Burkardt.TestValues.Gudermannian.gud_values ( ref n_data, ref x, ref fx );

                if ( n_data == 0 )
                {
                    break;
                }

                fx2 = Gudermannian.gud ( x );

                Console.WriteLine("  "
                                  + x.ToString().PadLeft(10)   + "  "
                                  + fx.ToString().PadLeft(10)  + "  "
                                  + fx2.ToString().PadLeft(10) + "");
            }

        }

        public static void agud_test ( )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    AGUD_TEST tests AGUD.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 June 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double g;
            int i;
            double x;
            double x2;

            Console.WriteLine("");
            Console.WriteLine("AGUD_TEST");
            Console.WriteLine("  AGUD computes the inverse Gudermannian;");
            Console.WriteLine("");
            Console.WriteLine("         X     GUD(X)     AGUD(GUD(X))");
            Console.WriteLine("");

            for ( i = 0; i <= 10; i++ )
            {
                x = 1.0 + ( ( double ) i ) / 5.0;
                g = Gudermannian.gud ( x );
                x2 = Gudermannian.agud ( g );

                Console.WriteLine("  " + x.ToString().PadLeft(10)
                                  + "  " + g.ToString().PadLeft(10)
                                  + "  " + x2.ToString().PadLeft(10)    + "");
            }

        }

    }
}