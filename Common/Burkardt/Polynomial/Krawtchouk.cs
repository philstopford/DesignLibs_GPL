using System;

namespace Burkardt.PolynomialNS;

public static class Krawtchouk
{
    public static void krawtchouk(int n, double p, double x, int m, ref double[] v)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    KRAWTCHOUK evaluates the Krawtchouk polynomials at X.
        //
        //  Discussion:
        //
        //    The polynomial has a parameter P, which must be striclty between
        //    0 and 1, and a parameter M which must be a nonnegative integer.
        //
        //    The Krawtchouk polynomial of order N, with parameters P and M,
        //    evaluated at X, may be written K(N,P,X,M).
        //
        //    The first two terms are:
        //
        //      K(0,P,X,M) = 1
        //      K(1,P,X,M) = X - P * M
        //
        //    and the recursion, for fixed P and M is
        //
        //                             ( N + 1 ) * K(N+1,P,X,M) =
        //        ( X - ( N + P * ( M - 2 * N))) * K(N,  P,X,M)
        //       - ( M - N + 1 ) * P * ( 1 - P ) * K(N-1,P,X,M)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 March 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Walter Gautschi,
        //    Orthogonal Polynomials: Computation and Approximation,
        //    Oxford, 2004,
        //    ISBN: 0-19-850672-4,
        //    LC: QA404.5 G3555.
        //
        //  Parameters:
        //
        //    Input, int N, the highest order polynomial to evaluate.
        //    0 <= N.
        //
        //    Input, double P, the parameter.  0 < P < 1.
        //
        //    Input, double X, the evaluation parameter.
        //
        //    Input, int M, the parameter.  0 <= M.
        //
        //    Output, double V[N+1], the values of the Krawtchouk polynomials
        //    of orders 0 through N at X.
        //
    {
        int i;

        switch (n)
        {
            case < 0:
                Console.WriteLine("");
                Console.WriteLine("KRAWTCHOUK - Fatal error!");
                Console.WriteLine("  0 <= N is required.");
                return;
        }

        switch (p)
        {
            case <= 0.0:
            case >= 1.0:
                Console.WriteLine("");
                Console.WriteLine("KRAWTCHOUK - Fatal error!");
                Console.WriteLine("  0 < P < 1 is required.");
                return;
        }

        switch (m)
        {
            case < 0:
                Console.WriteLine("");
                Console.WriteLine("KRAWTCHOUK - Fatal error!");
                Console.WriteLine("  0 <= M is required.");
                return;
        }

        v[0] = 1.0;

        v[1] = n switch
        {
            >= 1 => x - p * m,
            _ => v[1]
        };

        for (i = 1; i < n; i++)
        {
            v[i + 1] = (
                (x - (i + p * (m - 2 * i))) * v[i]
                - (m - i + 1) * p * (1.0 - p) * v[i - 1]
            ) / (i + 1);
        }

    }
}