using Burkardt.SubsetNS;
using Burkardt.Types;

namespace Burkardt.MatrixNS;

public static class MoebiusMatrix
{
    public static void moebius_matrix ( int n, int[] a, ref int[] mu )

//****************************************************************************80
//
//  Purpose:
//
//    MOEBIUS_MATRIX finds the Moebius matrix from a covering relation.
//
//  Discussion:
//
//    This routine can be called with A and MU being the same matrix.
//    The routine will correctly compute the Moebius matrix, which
//    will, in this case, overwrite the input matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, number of elements in the partially ordered set.
//
//    Input, int A[N*N].  A(I,J) = 1 if I is covered by J,
//    0 otherwise.
//
//    Output, int MU[N*N], the Moebius matrix as computed by the routine.
//
    {
        int i;
        int j;
        int[] p1;
        int[] p2;

        p1 = new int[n];
//
//  Compute a reordering of the elements of the partially ordered matrix.
//
        Triang.triang ( n, a, ref p1 );
//
//  Copy the matrix.
//
        for ( i = 0; i < n; i++ )
        {
            for ( j = 0; j < n; j++ )
            {
                mu[i+j*n] = a[i+j*n];
            }
        }
//
//  Apply the reordering to MU.
//
        typeMethods.i4mat_2perm0 ( n, n, mu, p1, p1 );
//
//  Negate the (strict) upper triangular elements of MU.
//
        for ( i = 0; i < n-1; i++ )
        {
            for ( j = i+1; j < n; j++ )
            {
                mu[i+j*n] = - mu[i+j*n];
            }
        }
//
//  Compute the inverse of MU.
//
        typeMethods.i4mat_u1_inverse ( n, mu, ref mu );
//
//  All nonzero elements are reset to 1.
//
        for ( i = 0; i < n; i++ )
        {
            for ( j = i; j < n; j++ )
            {
                if ( mu[i+j*n] != 0 )
                {
                    mu[i+j*n] = 1;
                }
            }
        }
//
//  Invert the matrix again.
//
        typeMethods.i4mat_u1_inverse ( n, mu, ref mu );
//
//  Compute the inverse permutation.
//
        p2 = Permutation.perm0_inverse ( n, p1 );
//
//  Unpermute the rows and columns of MU.
//
        typeMethods.i4mat_2perm0 ( n, n, mu, p2, p2 );
    }

    public static void moebius_values ( ref int n_data, ref int n, ref int c )

//****************************************************************************80
//
//  Purpose:
//
//    MOEBIUS_VALUES returns some values of the Moebius function.
//
//  Discussion:
//
//    MU(N) is defined as follows:
//
//      MU(N) = 1 if N = 1;
//              0 if N is divisible by the square of a prime;
//              (-1)**K, if N is the product of K distinct primes.
//
//    In Mathematica, the function can be evaluated by:
//
//      MoebiusMu[n]
//
//  First values:
//
//     N  MU(N)
//
//     1    1
//     2   -1
//     3   -1
//     4    0
//     5   -1
//     6    1
//     7   -1
//     8    0
//     9    0
//    10    1
//    11   -1
//    12    0
//    13   -1
//    14    1
//    15    1
//    16    0
//    17   -1
//    18    0
//    19   -1
//    20    0
//
//  Note:
//
//    As special cases, MU(N) is -1 if N is a prime, and MU(N) is 0
//    if N is a square, cube, etc.
//
//  Formula:
//
//    The Moebius function is related to Euler's totient function:
//
//      PHI(N) = Sum ( D divides N ) MU(D) * ( N / D ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 February 2003
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
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the argument of the Moebius function.
//
//    Output, int &C, the value of the Moebius function.
//
    {
        const int N_MAX = 20;

        int[] c_vec = {
            1,  -1,  -1,   0,  -1,   1,  -1,   0,   0,   1,
            -1,   0,  -1,   1,   1,   0,  -1,   0,  -1,   0 };

        int[] n_vec = {
            1,   2,   3,   4,   5,   6,   7,   8,   9,  10,
            11,  12,  13,  14,  15,  16,  17,  18,  19,  20 };

        n_data = n_data switch
        {
            < 0 => 0,
            _ => n_data
        };

        n_data += 1;

        if ( N_MAX < n_data )
        {
            n_data = 0;
            n = 0;
            c = 0;
        }
        else
        {
            n = n_vec[n_data-1];
            c = c_vec[n_data-1];
        }
    }
}