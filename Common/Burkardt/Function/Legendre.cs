using System;

namespace Burkardt.Function;

public static class Legendre
{
    public static void legendre_function_q(int n, double x, ref double[] cx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEGENDRE_FUNCTION_Q evaluates the Legendre Q functions.
        //
        //  Differential equation:
        //
        //    (1-X*X) Y'' - 2 X Y' + N (N+1) = 0
        //
        //  First terms:
        //
        //    Q(0,X) = 0.5 * log((1+X)/(1-X))
        //    Q(1,X) = Q(0,X)*X - 1
        //    Q(2,X) = Q(0,X)*(3*X*X-1)/4 - 1.5*X
        //    Q(3,X) = Q(0,X)*(5*X*X*X-3*X)/4 - 2.5*X^2 + 2/3
        //    Q(4,X) = Q(0,X)*(35*X^4-30*X^2+3)/16 - 35/8 * X^3 + 55/24 * X
        //    Q(5,X) = Q(0,X)*(63*X^5-70*X^3+15*X)/16 - 63/8*X^4 + 49/8*X^2 - 8/15
        //
        //  Recursion:
        //
        //    Q(0) = 0.5 * log ( (1+X) / (1-X) )
        //    Q(1) = 0.5 * X * log ( (1+X) / (1-X) ) - 1.0
        //
        //    Q(N) = ( (2*N-1) * X * Q(N-1) - (N-1) * Q(N-2) ) / N
        //
        //  Restrictions:
        //
        //    -1 < X < 1
        //
        //  Special values:
        //
        //    Note that the Legendre function Q(N,X) is equal to the
        //    associated Legendre function of the second kind,
        //    Q(N,M,X) with M = 0.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 February 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz, Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    National Bureau of Standards, 1964,
        //    ISBN: 0-486-61272-4,
        //    LC: QA47.A34.
        //
        //  Parameters:
        //
        //    Input, int N, the highest order function to evaluate.
        //
        //    Input, double X, the point at which the functions are to be
        //    evaluated.  X must satisfy -1 < X < 1.
        //
        //    Output, double CX[N+1], the values of the first N+1 Legendre
        //    functions at the point X.
        //
    {
        int i;
        switch (x)
        {
            //
            //  Check the value of X.
            //
            case <= -1.0:
            case >= 1.0:
                Console.WriteLine("");
                Console.WriteLine("LEGENDRE_FUNCTION_Q - Fatal error!");
                Console.WriteLine("  Illegal input value of X = " + x + "");
                Console.WriteLine("  But X must be between -1 and 1.");
                return;
        }

        switch (n)
        {
            case < 0:
                return;
        }

        cx[0] = 0.5 * Math.Log((1.0 + x) / (1.0 - x));

        switch (n)
        {
            case 0:
                return;
        }

        cx[1] = x * cx[0] - 1.0;

        for (i = 2; i <= n; i++)
        {
            cx[i] = ((2 * i - 1) * x * cx[i - 1]
                     + (-i + 1) * cx[i - 2])
                    / i;
        }

    }
}