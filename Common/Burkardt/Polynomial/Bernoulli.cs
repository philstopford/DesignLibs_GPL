using Burkardt.SubsetNS;

namespace Burkardt.PolynomialNS
{
    public static class Bernoulli
    {
        public static double bernoulli_poly(int n, double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BERNOULLI_POLY evaluates the Bernoulli polynomial of order N at X.
            //
            //  Discussion:
            //
            //    Thanks to Bart Vandewoestyne for pointing out an error in the previous
            //    documentation, 31 January 2008.
            //
            //    Special values of the Bernoulli polynomial include:
            //
            //      B(N,0) = B(N,1) = B(N), the N-th Bernoulli number.
            //
            //      B'(N,X) = N * B(N-1,X)
            //
            //      B(N,X+1) - B(N,X) = N * X^(N-1)
            //      B(N,X) = (-1)^N * B(N,1-X)
            //
            //    A formula for the Bernoulli polynomial in terms of the Bernoulli
            //    numbers is:
            //
            //      B(N,X) = sum ( 0 <= K <= N ) B(K) * C(N,K) * X^(N-K)
            //
            //    The first few polynomials include:
            //
            //      B(0,X) = 1
            //      B(1,X) = X    - 1/2
            //      B(2,X) = X^2 -   X      +  1/6
            //      B(3,X) = X^3 - 3/2*X^2 +  1/2*X
            //      B(4,X) = X^4 - 2*X^3   +      X^2 - 1/30
            //      B(5,X) = X^5 - 5/2*X^4 +  5/3*X^3 - 1/6*X
            //      B(6,X) = X^6 - 3*X^5   +  5/2*X^4 - 1/2*X^2 + 1/42
            //      B(7,X) = X^7 - 7/2*X^6 +  7/2*X^5 - 7/6*X^3 + 1/6*X
            //      B(8,X) = X^8 - 4*X^7   + 14/3*X^6 - 7/3*X^4 + 2/3*X^2 - 1/30
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    31 January 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the Bernoulli polynomial to
            //    be evaluated.  N must be 0 or greater.
            //
            //    Input, double X, the value of X at which the polynomial is to
            //    be evaluated.
            //
            //    Output, double BERNOULLI_POLY, the value of B(N,X).
            //
        {
            int[] c;
            int i;
            double value;
            double[] work;

            work = new double[n + 1];
            Sequence.Bernoulli.bernoulli_number(n, ref work);
            //
            //  Get row N of Pascal's triangle.
            //
            c = new int[n + 1];
            for (i = 0; i <= n; i++)
            {
                Comb.comb_row_next(n, ref c);
            }

            value = 1.0;
            for (i = 1; i <= n; i++)
            {
                value = value * x + work[i] * (double)c[i];
            }

            return value;
        }

        public static double bernoulli_poly2(int n, double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BERNOULLI_POLY2 evaluates the N-th Bernoulli polynomial at X.
            //
            //  Discussion:
            //
            //    Thanks to Bart Vandewoestyne for pointing out an error in the previous
            //    documentation, 31 January 2008.
            //
            //    Special values of the Bernoulli polynomial include:
            //
            //      B(N,0) = B(N,1) = B(N), the N-th Bernoulli number.
            //
            //      B'(N,X) = N * B(N-1,X)
            //
            //      B(N,X+1) - B(N,X) = N * X^(N-1)
            //      B(N,X) = (-1)^N * B(N,1-X)
            //
            //    A formula for the Bernoulli polynomial in terms of the Bernoulli
            //    numbers is:
            //
            //      B(N,X) = sum ( 0 <= K <= N ) B(K) * C(N,K) * X^(N-K)
            //
            //    The first few polynomials include:
            //
            //      B(0,X) = 1
            //      B(1,X) = X    - 1/2
            //      B(2,X) = X^2 -   X      +  1/6
            //      B(3,X) = X^3 - 3/2*X^2 +  1/2*X
            //      B(4,X) = X^4 - 2*X^3   +      X^2 - 1/30
            //      B(5,X) = X^5 - 5/2*X^4 +  5/3*X^3 - 1/6*X
            //      B(6,X) = X^6 - 3*X^5   +  5/2*X^4 - 1/2*X^2 + 1/42
            //      B(7,X) = X^7 - 7/2*X^6 +  7/2*X^5 - 7/6*X^3 + 1/6*X
            //      B(8,X) = X^8 - 4*X^7   + 14/3*X^6 - 7/3*X^4 + 2/3*X^2 - 1/30
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    31 January 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the Bernoulli polynomial to
            //    be evaluated.  N must be 0 or greater.
            //
            //    Input, double X, the value at which the polynomial is to
            //    be evaluated.
            //
            //    Output, double BERNOULLI_POLY2, the value of B(N,X).
            //
        {
            double fact;
            int i;
            double value;

            fact = 1.0;

            value = Sequence.Bernoulli.bernoulli_number3(0);

            for (i = 1; i <= n; i++)
            {
                fact = fact * (double)(n + 1 - i) / (double)i;
                value = value * x + fact * Sequence.Bernoulli.bernoulli_number3(i);
            }

            return value;
        }
    }
}