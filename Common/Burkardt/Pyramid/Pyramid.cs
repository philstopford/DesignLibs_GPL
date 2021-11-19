using Burkardt.Types;

namespace Burkardt.PyramidNS;

public static class Pyramid
{
    public static double pyra_unit_monomial(int[] expon)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PYRA_UNIT_MONOMIAL: monomial integral in a unit pyramid.
        //
        //  Discussion:
        //
        //    This function returns the value of the integral of X^ALPHA Y^BETA Z^GAMMA
        //    over the unit pyramid.
        //
        //    The integration region is:
        //
        //    - ( 1 - Z ) <= X <= 1 - Z
        //    - ( 1 - Z ) <= Y <= 1 - Z
        //              0 <= Z <= 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 March 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Arthur Stroud,
        //    Approximate Calculation of Multiple Integrals,
        //    Prentice Hall, 1971,
        //    ISBN: 0130438936,
        //    LC: QA311.S85.
        //
        //  Parameters:
        //
        //    Input, int EXPON[3], the exponents.
        //
        //    Output, double PYRA_UNIT_MONOMIAL, the integral of the monomial
        //    over the pyramid.
        //
    {
        double value = 0.0;

        switch (expon[0] % 2)
        {
            case 0 when expon[1] % 2 == 0:
            {
                int i_hi = 2 + expon[0] + expon[1];

                int i;
                for (i = 0; i <= i_hi; i++)
                {
                    value += typeMethods.r8_mop(i) * typeMethods.r8_choose(i_hi, i)
                             / (i + expon[2] + 1);
                }

                value = value
                    * 2.0 / (expon[0] + 1)
                    * 2.0 / (expon[1] + 1);
                break;
            }
        }

        return value;
    }

    public static double pyra_unit_volume()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PYRA_UNIT_VOLUME: volume of a unit pyramid with square base.
        //
        //  Discussion:
        //
        //    The volume of this unit pyramid is 4/3.
        //
        //    The integration region is:
        //
        //      - ( 1 - Z ) <= X <= 1 - Z
        //      - ( 1 - Z ) <= Y <= 1 - Z
        //                0 <= Z <= 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 March 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, double PYRA_UNIT_VOLUME, the volume of the pyramid.
        //
    {
        double volume = 4.0 / 3.0;

        return volume;
    }

    public static int pyramid_num(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PYRAMID_NUM returns the N-th pyramidal number.
        //
        //  Discussion:
        //
        //    The N-th pyramidal number P(N) is formed by the sum of the first
        //    N triangular numbers T(J):
        //
        //      T(J) = sum ( 1 <= J <= N ) J
        //
        //      P(N) = sum ( 1 <= I <= N ) T(I)
        //
        //    By convention, T(0) = 0.
        //
        //    The formula is:
        //
        //      P(N) = ( (N+1)^3 - (N+1) ) / 6
        //
        //    Note that this pyramid will have a triangular base.
        //
        //  First Values:
        //
        //      0
        //      1
        //      4
        //     10
        //     20
        //     35
        //     56
        //     84
        //    120
        //    165
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 May 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the index of the desired number, which must be
        //    at least 0.
        //
        //    Output, int PYRAMID_NUM, the N-th pyramidal number.
        //
    {
        int value = ((n + 1) * (n + 1) * (n + 1) - (n + 1)) / 6;

        return value;
    }

    public static int pyramid_square_num(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PYRAMID_SQUARE_NUM returns the N-th pyramidal square number.
        //
        //  Discussion:
        //
        //    The N-th pyramidal square number PS(N) is formed by the sum of the first
        //    N squares S:
        //
        //      S(I) = I^2
        //
        //      PS(N) = sum ( 1 <= I <= N ) S(I)
        //
        //    By convention, PS(0) = 0.
        //
        //    The formula is:
        //
        //      PS(N) = ( N * ( N + 1 ) * ( 2*N+1 ) ) / 6
        //
        //    Note that geometrically, this pyramid will have a square base.
        //
        //  Example:
        //
        //    0    0
        //    1    1
        //    2    5
        //    3   14
        //    4   30
        //    5   55
        //    6   91
        //    7  140
        //    8  204
        //    9  285
        //   10  385
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 August 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the index.
        //    0 <= N.
        //
        //    Output, int PYRAMID_SQUARE_NUM, the N-th 
        //    pyramid square number.
        //
    {
        int value = n * (n + 1) * (2 * n + 1) / 6;

        return value;
    }

}