namespace Burkardt.LineNS;

public static class Monomial
{
    public static double[] line_monomial_moments ( double a, double b, int m )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_MONOMIAL_MOMENTS computes monomial moments in [A,B].
        //
        //  Discussion:
        //
        //    We use the uniform weight and the shifted and scaled monomial basis:
        //
        //      P(a,b,i;x) = xi(a,b;x)^(i-1)
        //       xi(a,b;x) = ( - ( b - x ) + ( x - a ) ) / ( b - a )
        //
        //    The i-th moment is
        //
        //      mom(i) = integral ( a <= x <= b ) P(a,b,i;x) dx
        //             = integral ( a <= x <= b ) xi(a,b;x)^(i-1) dx
        //             = 0.5 * ( b - a ) * integral ( -1 <= xi <= +1 ) xi^(i-1) dxi
        //             = 0.5 * ( b - a ) * xi^i / i | ( -1 <= xi <= +1 )
        //             = 0.5 * ( b - a ) * ( 1 - (-1)^i ) / i
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 April 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, the endpoints of the interval.
        //
        //    Input, int M, the number of basis polynomials.
        //
        //    Output, double LINE_MONOMIAL_MOMENTS[M], the moments.
        //
    {
        int i;

        double[] mom = new double[m];

        for ( i = 0; i < m; i++ )
        {
            mom[i] = ( b - a ) * (( i + 1 ) % 2) / (i + 1);
        }

        return mom;
    }
}