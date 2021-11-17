using System;

namespace Burkardt.PolynomialNS;

public static class Trinomial
{
    public static int trinomial ( int i, int j, int k )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRINOMIAL computes a trinomial coefficient.
        //
        //  Discussion:
        //
        //    The trinomial coefficient is a generalization of the binomial
        //    coefficient.  It may be interpreted as the number of combinations of
        //    N objects, where I objects are of type 1, J of type 2, and K of type 3.
        //    and N = I + J + K.
        //
        //    T(I,J,K) = N! / ( I! J! K! )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 April 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I, J, K, the factors.
        //    All should be nonnegative.
        //
        //    Output, int TRINOMIAL, the trinomial coefficient.
        //
    {
        int l;
        int t;
        int value;
        //
        //  Each factor must be nonnegative.
        //
        if ( i < 0 || j < 0 || k < 0 )
        {
            Console.WriteLine("");
            Console.WriteLine("TRINOMIAL - Fatal error!");
            Console.WriteLine("  Negative factor encountered.");
            return 1;
        }

        value = 1;

        t = 1;

        for ( l = 1; l <= i; l++ )
        {
            //
            //  value = value * t / l;
            //
            t += 1;
        }

        for ( l = 1; l <= j; l++ )
        {
            value = value * t / l;
            t += 1;
        }

        for ( l = 1; l <= k; l++ )
        {
            value = value * t / l;
            t += 1;
        }
  
        return value;
    }
}