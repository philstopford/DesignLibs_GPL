namespace Burkardt.TriangleNS
{
    public static class Monomial
    {
        public static double triangle_unit_monomial ( int ex, int ey )

            //****************************************************************************8
            //
            //  Purpose:
            //
            //    triangle_unit_monomial integrates a monomial over the unit triangle.
            //
            //  Discussion:
            //
            //    This routine integrates a monomial of the form
            //
            //      x^ex y^ey
            //
            //    where the exponents are nonnegative integers.  Note that
            //    if the combination 0^0 is encountered, it should be treated
            //    as 1.
            //
            //    Integral ( over unit triangle ) x^m y^n dx dy = m! * n! / ( m + n + 2 )!
            //
            //    The integration region is:
            //
            //      0 <= X
            //      0 <= Y
            //      X + Y <= 1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 April 2019
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int EX, EY, the exponents.
            //
            //    Output, double TRIANGLE_UNIT_MONOMIAL, the integral of the monomial.
            //
        {
            int i;
            int k;
            double value;

            value = 1.0;
            k = ex;

            for ( i = 1; i <= ey; i++ )
            {
                k = k + 1;
                value = value * ( double ) ( i ) / ( double ) ( k );
            }

            k = k + 1;
            value = value / ( double ) ( k );

            k = k + 1;
            value = value / ( double ) ( k );

            return value;
        }
    }
}