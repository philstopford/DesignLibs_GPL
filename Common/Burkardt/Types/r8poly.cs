using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static void r8poly_print ( int n, double[] a, string title )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8POLY_PRINT prints out a polynomial.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the dimension of A.
        //
        //    Input, double A[N+1], the polynomial coefficients.
        //    A(0) is the constant term and
        //    A(N) is the coefficient of X^N.
        //
        //    Input, string TITLE, a title.
        //
        {
            int i;
            double mag;
            char plus_minus;

            if ( 0 < title.Length )
            {
                Console.WriteLine("");
                Console.WriteLine(title + "");
            }
            Console.WriteLine("");

            if ( n < 0 )
            {
                Console.WriteLine("  p(x) = 0");
                return;
            }

            if ( a[n] < 0.0 )
            {
                plus_minus = '-';
            }
            else
            {
                plus_minus = ' ';
            }

            mag = Math.Abs ( a[n] );

            if ( 2 <= n )
            {
                Console.WriteLine("  p(x) = " + plus_minus
                    + mag.ToString().PadLeft(14) + " * x ^ " + n + "");
            }
            else if ( n == 1 )
            {
                Console.WriteLine("  p(x) = " + plus_minus
                    + mag.ToString().PadLeft(14) + " * x");
            }
            else if ( n == 0 )
            {
                Console.WriteLine("  p(x) = " + plus_minus
                    + mag.ToString().PadLeft(14) + "");
            }

            for ( i = n - 1; 0 <= i; i-- )
            {
                if ( a[i] < 0.0 )
                {
                    plus_minus = '-';
                }
                else
                {
                    plus_minus = '+';
                }

                mag = Math.Abs ( a[i] );

                if ( mag != 0.0 )
                {
                    if ( 2 <= i )
                    {
                        Console.WriteLine("         " + plus_minus
                            + mag.ToString().PadLeft(14) + " * x ^ " + i + "");
                    }
                    else if ( i == 1 )
                    {
                        Console.WriteLine("         " + plus_minus
                            + mag.ToString().PadLeft(14) + " * x");
                    }
                    else if ( i == 0 )
                    {
                        Console.WriteLine("         " + plus_minus
                            + mag.ToString().PadLeft(14) + "");
                    }
                }
            }
        }
    }
}