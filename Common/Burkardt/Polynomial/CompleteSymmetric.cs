using System;

namespace Burkardt.PolynomialNS;

public static class CompleteSymmetric
{
    public static double complete_symmetric_poly(int n, int r, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COMPLETE_SYMMETRIC_POLY evaluates a complete symmetric polynomial.
        //
        //  Discussion:
        //
        //    N\R  0   1         2               3
        //      +--------------------------------------------------------
        //    0 |  1   0         0               0
        //    1 |  1   X1        X1^2            X1^3
        //    2 |  1   X1+X2     X1^2+X1X2+X2^2  X1^3+X1^2X2+X1X2^2+X2^3
        //    3 |  1   X1+X2+X3  ...
        //
        //    If X = ( 1, 2, 3, 4, 5, ... ) then
        //
        //    N\R  0     1     2     3     4 ...
        //      +--------------------------------------------------------
        //    0 |  1     0     0     0     0
        //    1 |  1     1     1     1     1
        //    2 |  1     3     7    15    31
        //    3 |  1     6    25    90   301
        //    4 |  1    10    65   350  1701
        //    5 |  1    15   140  1050  6951
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 November 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of variables.
        //    0 <= N.
        //
        //    Input, int R, the degree of the polynomial.
        //    0 <= R.
        //
        //    Input, double X[N], the value of the variables.
        //
        //    Output, double COMPLETE_SYMMETRIC_POLY, the value of TAU(N,R)(X).
        //
    {
        int i;
        int nn;

        switch (n)
        {
            case < 0:
                Console.WriteLine("");
                Console.WriteLine("COMPLETE_SYMMETRIC_POLY - Fatal error!");
                Console.WriteLine("  N < 0.");
                return 1;
        }

        switch (r)
        {
            case < 0:
                Console.WriteLine("");
                Console.WriteLine("COMPLETE_SYMMETRIC_POLY - Fatal error!");
                Console.WriteLine("  R < 0.");
                return 1;
        }

        double[] tau = new double[1 + Math.Max(n, r)];

        for (i = 0; i <= Math.Max(n, r); i++)
        {
            tau[i] = 0.0;
        }

        tau[0] = 1.0;
        for (nn = 1; nn <= n; nn++)
        {
            int rr;
            for (rr = 1; rr <= r; rr++)
            {
                tau[rr] += x[nn - 1] * tau[rr - 1];
            }
        }

        double value = tau[r];

        return value;
    }

}