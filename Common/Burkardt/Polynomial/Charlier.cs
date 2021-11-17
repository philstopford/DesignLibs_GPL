using System;

namespace Burkardt.PolynomialNS;

public static class Charlier
{
    public static void charlier ( int n, double a, double x, ref double[] value )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHARLIER evaluates Charlier polynomials at a point.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 March 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    J Simoes Pereira,
        //    Algorithm 234: Poisson-Charliers Polynomials,
        //    Communications of the ACM,
        //    Volume 7, Number 7, page 420, July 1964.
        //
        //    Walter Gautschi,
        //    Orthogonal Polynomials: Computation and Approximation,
        //    Oxford, 2004,
        //    ISBN: 0-19-850672-4,
        //    LC: QA404.5 G3555.
        //
        //    Gabor Szego,
        //    Orthogonal Polynomials,
        //    American Mathematical Society, 1975,
        //    ISBN: 0821810235,
        //    LC: QA3.A5.v23.
        //
        //    Eric Weisstein,
        //    CRC Concise Encyclopedia of Mathematics,
        //    CRC Press, 2002,
        //    Second edition,
        //    ISBN: 1584883472,
        //    LC: QA5.W45.
        //
        //  Parameters:
        //
        //    Input, int N, the maximum order of the polynomial.  
        //    N must be at least 0.
        //
        //    Input, double A, the parameter.  A must not be 0.
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double VALUE[0:N], the value of the polynomials at X.
        //
    {
        int i;

        switch (a)
        {
            case 0.0:
                Console.WriteLine("");
                Console.WriteLine("CHARLIER - Fatal error!");
                Console.WriteLine("  Parameter A cannot be zero.");
                return;
        }

        switch (n)
        {
            case < 0:
                Console.WriteLine("");
                Console.WriteLine("CHARLIER - Fatal error!");
                Console.WriteLine("  Parameter N must be nonnegative.");
                return;
        }

        value[0] = 1.0;

        switch (n)
        {
            case 0:
                return;
        }

        value[1] = - x / a;

        switch (n)
        {
            case 1:
                return;
        }

        for ( i = 1; i < n; i++ )
        {
            value[i+1] = ( ( i + a - x ) * value[i] - i * value[i-1] ) / a;
        }

    }

}