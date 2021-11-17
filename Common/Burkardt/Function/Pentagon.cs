namespace Burkardt.Function;

public static class Pentagon
{
    public static int pentagon_num ( int n )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PENT_ENUM computes the N-th pentagonal number.
        //
        //  Discussion:
        //
        //    The pentagonal number P(N) counts the number of dots in a figure of
        //    N nested pentagons.  The pentagonal numbers are defined for both
        //    positive and negative N.
        //
        //  First values:
        //
        //    N   P
        //
        //   -5   40
        //   -4   26
        //   -3   15
        //   -2    7
        //   -1    2
        //    0    0
        //    1    1
        //    2    5
        //    3   12
        //    4   22
        //    5   35
        //    6   51
        //    7   70
        //    8   92
        //    9  117
        //   10  145
        //
        //    P(N) = ( N * ( 3 * N - 1 ) ) / 2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 May 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the index of the pentagonal number desired.
        //
        //    Output, int PENT_ENUM, the value of the N-th pentagonal number.
        //
    {
        int p;

        p = n * ( 3 * n - 1 ) / 2;

        return p;
    }
}