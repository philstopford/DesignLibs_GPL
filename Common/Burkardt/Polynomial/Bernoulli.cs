using System;
using Burkardt.Sequence;
using Burkardt.SubsetNS;

namespace Burkardt.PolynomialNS;

public static class Bernoulli
{
    public static int poly_bernoulli(int n, int k)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLY_BERNOULLI evaluates the poly-Bernolli numbers with negative index.
        //
        //  Discussion:
        //
        //    The poly-Bernoulli numbers B_n^k were defined by M Kaneko
        //    formally as the coefficients of X^n/n! in a particular power
        //    series.  He also showed that, when the super-index is negative,
        //    we have
        //
        //      B_n^(-k) = Sum ( 0 <= j <= min ( n, k ) )
        //        (j!)^2 * S(n+1,j+1) * S(k+1,j+1)
        //
        //    where S(n,k) is the Stirling number of the second kind, the number of
        //    ways to partition a set of size n into k nonempty subset.
        //
        //    B_n^(-k) is also the number of "lonesum matrices", that is, 0-1
        //    matrices of n rows and k columns which are uniquely reconstructable
        //    from their row and column sums.
        //
        //    The poly-Bernoulli numbers get large very quickly.
        //
        //  Table:
        //
        //    \ K 0  1    2     3      4       5        6
        //    N
        //    0   1  1    1     1      1       1        1
        //    1   1  2    4     8     16      32       64
        //    2   1  4   14    46    146     454     1394
        //    3   1  8   46   230   1066    4718    20266
        //    4   1 16  146  1066   6902   41506   237686
        //    5   1 32  454  4718  41506  329462  2441314
        //    6   1 64 1394 20266 237686 2441314 22934774
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 March 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Chad Brewbaker,
        //    Lonesum (0,1) Matrices and Poly-Bernoulli Numbers of Negative Index,
        //    MS Thesis,
        //    Iowa State University, 2005.
        //
        //    M Kaneko,
        //    Poly-Bernoulli Numbers,
        //    Journal Theorie des Nombres Bordeaux,
        //    Volume 9, 1997, pages 221-228.
        //
        //  Parameters:
        //
        //    Input, int N, K, the indices.  N and K should be nonnegative.
        //
        //    Output, int POLY_BERNOULLI, the value of B_N^(-K).
        //
    {
        int b;
        int j;

        switch (n)
        {
            case < 0:
                b = 0;
                return b;
            case 0:
                b = 1;
                return b;
        }

        switch (k)
        {
            case <= 0:
                b = 0;
                return b;
        }

        switch (k)
        {
            case 0:
                b = 1;
                return b;
        }

        int jhi = Math.Min(n, k);
        int m = Math.Max(n, k) + 1;

        int[] s = Stirling.stirling2(m, m);

        int jfact = 1;
        b = 0;

        for (j = 0; j <= jhi; j++)
        {
            b += jfact * jfact * s[n + j * m] * s[k + j * m];

            jfact *= j + 1;
        }

        return b;
    }

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
        int i;

        double[] work = new double[n + 1];
        Sequence.Bernoulli.bernoulli_number(n, ref work);
        //
        //  Get row N of Pascal's triangle.
        //
        int[] c = new int[n + 1];
        for (i = 0; i <= n; i++)
        {
            Comb.comb_row_next(n, ref c);
        }

        double value = 1.0;
        for (i = 1; i <= n; i++)
        {
            value = value * x + work[i] * c[i];
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
        int i;

        double fact = 1.0;

        double value = Sequence.Bernoulli.bernoulli_number3(0);

        for (i = 1; i <= n; i++)
        {
            fact = fact * (n + 1 - i) / i;
            value = value * x + fact * Sequence.Bernoulli.bernoulli_number3(i);
        }

        return value;
    }
}