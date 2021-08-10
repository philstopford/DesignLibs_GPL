namespace Burkardt.TetrahedronNS
{
    public static class Integrals
    {
        public static double tet01_monomial_integral ( int dim_num, int[] expon )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TET01_MONOMIAL_INTEGRAL integrates a monomial over the unit tetrahedron.
            //
            //  Discussion:
            //
            //    This routine evaluates a monomial of the form
            //
            //      product ( 1 <= dim <= dim_num ) x(dim)^expon(dim)
            //
            //    where the exponents are nonnegative integers.  Note that
            //    if the combination 0^0 is encountered, it should be treated
            //    as 1.
            //
            //    Integral ( over unit tetrahedron ) x^l y^m z^n dx dy 
            //    = l! * m! * n! / ( l + m + n + 3 )!
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    04 July 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the spatial dimension.
            //
            //    Input, int EXPON[DIM_NUM], the exponents.
            //
            //    Output, double TET01_MONOMIAL_INTEGRAL, the value of the integral of the
            //    monomial.
            //
        {
            int i;
            int k;
            double value;
            //
            //  The first computation ends with VALUE = 1.0;
            //
            value = 1.0;

            k = 0;

            for ( i = 1; i <= expon[0]; i++ )
            {
                k = k + 1;
                //  value = value * ( double ) ( i ) / ( double ) ( k );
            }

            for ( i = 1; i <= expon[1]; i++ )
            {
                k = k + 1;
                value = value * ( double ) ( i ) / ( double ) ( k );
            }

            for ( i = 1; i <= expon[2]; i++ )
            {
                k = k + 1;
                value = value * ( double ) ( i ) / ( double ) ( k );
            }

            k = k + 1;
            value = value / ( double ) ( k );

            k = k + 1;
            value = value / ( double ) ( k );

            k = k + 1;
            value = value / ( double ) ( k );

            return value;
        }
    }
}