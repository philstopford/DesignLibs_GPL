using System;
using Burkardt.Graph;
using Burkardt.Types;

namespace Burkardt.TriangleNS;

public static class Integrals
{
    public static double triangle_area(double[] t)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_AREA returns the area of a triangle.
        //
        //  Discussion:
        //
        //    If the vertices are given in counter clockwise order, the area
        //    will be positive.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU GPL license.
        //
        //  Modified:
        //
        //    21 April 2015
        //
        //  Author:
        //
        //    John Burkardt.
        //
        //  Parameters:
        //
        //    Input, double T[2*3], the vertices of the triangle.
        //
        //    Output, double TRIANGLE_AREA, the area of the triangle.
        //
    {
        double value = 0.5 *
                       (
                           (t[0 + 1 * 2] - t[0 + 0 * 2]) * (t[1 + 2 * 2] - t[1 + 0 * 2])
                           - (t[0 + 2 * 2] - t[0 + 0 * 2]) * (t[1 + 1 * 2] - t[1 + 0 * 2])
                       );

        return value;
    }

    public static double triangle_monomial_integral(int i, int j, double[] t)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_MONOMIAL_INTEGRAL integrates a monomial over an arbitrary triangle.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 April 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I, J, the exponents of X and Y in the monomial.
        //    0 <= I, J.
        //
        //    Input, double T[2*3], the vertices of the triangle.
        //
        //    Output, double TRIANGLE_MONOMIAL_INTEGRAL, the integral 
        //    of X^I * Y^J over triangle T.
        //
    {
        double a = 0;
        double b = 0;
        double c = 0;
        double d = 0;
        double e = 0;
        double f = 0;
        //
        //  Get map coefficients from reference RS triangle to general XY triangle.
        //    R = a+b*X+c*Y
        //    S = d+e*X+f*Y
        //
        Map.rs_to_xy_map(t, ref a, ref b, ref c, ref d, ref e, ref f);
        //
        //  Set
        //    P1(R,S) = a+b*R+c*S
        //    P2(R,S) = d+e*R+f*S
        //
        const int d1 = 1;
        int m1 = (d1 + 1) * (d1 + 2) / 2;
        double[] p1 = new double[m1];
        p1[0] = a;
        p1[1] = b;
        p1[2] = c;

        const int d2 = 1;
        int m2 = (d2 + 1) * (d2 + 2) / 2;
        double[] p2 = new double[m2];
        p2[0] = d;
        p2[1] = e;
        p2[2] = f;
        //
        //  Exponentiate:
        //    P3(R,S) = P1(R,S)^i
        //    P4(R,S) = P2(R,S)^j
        //
        int d3 = i * d1;
        double[] p3 = typeMethods.poly_power_linear(d1, p1, i);

        int d4 = j * d2;
        double[] p4 = typeMethods.poly_power_linear(d2, p2, j);
        //
        //  Compute the product 
        //    P5(R,S) = P3(R,S) * P4(R,S)
        //
        int d5 = d3 + d4;
        double[] p5 = typeMethods.poly_product(d3, p3, d4, p4);
        //
        //  Compute the integral of P5(R,S) over the reference triangle.
        //
        double q = triangle01_poly_integral(d5, p5);
        //
        //  Multiply by the area of the physical triangle T(X,Y) divided by
        //  the area of the reference triangle.
        //
        q = q * triangle_area(t) / 0.5;

        return q;
    }

    public static double triangle_poly_integral(int d, double[] p, double[] t )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_POLY_INTEGRAL: polynomial integral over a triangle.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 April 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int D, the degree of the polynomial.
        //
        //    Input, double P[M], the polynomial coefficients.
        //    M = ((D+1)*(D+2))/2.
        //
        //    Input, double T[2*3], the vertices of the triangle.
        //
        //    Output, double TRIANGLE_POLY_INTEGRAL, the integral.
        //
    {
        int i = 0;
        int j = 0;
        int km1;

        int m = (d + 1) * (d + 2) / 2;

        double q = 0.0;
        for (km1 = 0; km1 < m; km1++)
        {
            int k = km1 + 1;
            typeMethods.i4_to_pascal(k, ref i, ref j);
            q += p[km1] * triangle_monomial_integral(i, j, t);
        }

        return q;
    }

    public static double triangle_xy_integral(double x1, double y1, double x2, double y2,
            double x3, double y3)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_XY_INTEGRAL computes the integral of XY over a triangle.
        //
        //  Discussion:
        //
        //    This function was written as a special test case for the general
        //    problem of integrating a monomial x^alpha * y^beta over a general 
        //    triangle.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 April 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X1, Y1, X2, Y2, X3, Y3, the coordinates of the
        //    triangle vertices.
        //
        //    Output, double TRIANGLE_XY_INTEGRAL, the integral of X*Y 
        //    over the triangle.
        //
    {
        //
        //  x = x1 * ( 1.0 - xi - eta )
        //    + x2 *         xi
        //    + x3 *              eta;
        //
        //  y = y1 * ( 1.0 - xi - eta )
        //    + y2 *         xi
        //    + y3 *              eta;
        //
        //  Rewrite as linear polynomials in (xi,eta):
        //
        //  x = x1 + ( x2 - x1 ) * xi + ( x3 - x1 ) * eta
        //  y = y1 + ( y2 - y1 ) * xi + ( y3 - y1 ) * eta
        //
        //  Jacobian:
        //
        //    J = [ ( x2 - x1 )  ( x3 - x1 ) ]
        //        [ ( y2 - y1 )  ( y3 - y1 ) ]
        //
        //    det J = ( x2 - x1 ) * ( y3 - y1 ) - ( y2 - y1 ) * ( x3 - x1 )
        //
        //  Integrand
        //
        //    x * y = ( x1 + ( x2 - x1 ) * xi + ( x3 - x1 ) * eta )
        //          * ( y1 + ( y2 - y1 ) * xi + ( y3 - y1 ) * eta )
        //
        //  Rewrite as linear combination of monomials:
        //
        //    x * y = 1      * x1 * y1
        //          + eta    * ( x1 * ( y3 - y1 ) + ( x3 - x1 ) * y1 )
        //          + xi     * ( x1 * ( y2 - y1 ) + ( x2 - x1 ) * y1 )
        //          + eta^2  * ( x3 - x1 ) * ( y3 - y1 )
        //          + xi*eta * ( ( x2 - x1 ) * ( y3 - y1 ) + ( x3 - x1 ) * ( y2 - y1 ) )
        //          + xi^2   * ( x2 - x1 ) * ( y2 - y1 )
        //
        double det = (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1);

        double p00 = x1 * y1;

        double p01 = x1 * (y3 - y1) + (x3 - x1) * y1;
        double p10 = x1 * (y2 - y1) + (x2 - x1) * y1;

        double p02 = (x3 - x1) * (y3 - y1);
        double p11 = (x2 - x1) * (y3 - y1) + (x3 - x1) * (y2 - y1);
        double p20 = (x2 - x1) * (y2 - y1);

        double q = 0.0;
        q += p00 * triangle01_monomial_integral(0, 0);
        q += p10 * triangle01_monomial_integral(1, 0);
        q += p01 * triangle01_monomial_integral(0, 1);
        q += p20 * triangle01_monomial_integral(2, 0);
        q += p11 * triangle01_monomial_integral(1, 1);
        q += p02 * triangle01_monomial_integral(0, 2);

        q *= det;

        return q;
    }

    public static double triangle01_monomial_integral(int i, int j)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE01_MONOMIAL_INTEGRAL: monomial integrals in the unit triangle.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 April 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I, J, the exponents.  
        //    Each exponent must be nonnegative.
        //
        //    Output, double TRIANGLE01_MONOMIAL_INTEGRAL, the integral.
        //
    {
        int l;

        int k = 0;
        double q = 1.0;

        for (l = 1; l <= i; l++)
        {
            k += 1;
            q = q * l / k;
        }

        for (l = 1; l <= j; l++)
        {
            k += 1;
            q = q * l / k;
        }

        for (l = 1; l <= 2; l++)
        {
            k += 1;
            q /= k;
        }

        return q;
    }
        
    public static double triangle01_monomial_integral ( int dim_num, int[]  expon)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE01_MONOMIAL_INTEGRAL integrates a monomial over the unit triangle.
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
        //    Integral ( over unit triangle ) x^m y^n dx dy = m! * n! / ( m + n + 2 )!
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
        //    Output, double TRIANGLE01_MONOMIAL_INTEGRAL, the value of the integral 
        //    of the monomial.
        //
    {
        int i;
        //
        //  The first computation ends with VALUE = 1.0;
        //
        double value = 1.0;

        int k = 0;

        for ( i = 1; i <= expon[0]; i++ )
        {
            k += 1;
            //  value = value * ( double ) ( i ) / ( double ) ( k );
        }

        for ( i = 1; i <= expon[1]; i++ )
        {
            k += 1;
            value = value * i / k;
        }

        k += 1;
        value /= k;

        k += 1;
        value /= k;

        return value;
    }

    public static double triangle01_monomial_quadrature ( int dim_num, int[] expon, 
            int point_num, double[] x, double[] weight )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRANGLE01_MONOMIAL_QUADRATURE applies quadrature to a monomial in a triangle.
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
        //    Input, int POINT_NUM, the number of points in the rule.
        //
        //    Input, double X[DIM_NUM*POINT_NUM], the quadrature points.
        //
        //    Input, double WEIGHT[POINT_NUM], the quadrature weights.
        //
        //    Output, double TRIANGLE01_MONOMIAL_QUADRATURE, the quadrature error.
        //
    {
        //
        //  Get the exact value of the integral of the unscaled monomial.
        //
        double scale = triangle01_monomial_integral ( dim_num, expon );
        //
        //  Evaluate the monomial at the quadrature points.
        //
        double[] value = MonomialNS.Monomial.monomial_value ( dim_num, point_num, expon, x );
        //
        //  Compute the weighted sum and divide by the exact value.
        //
        double area = 0.5;
        double quad = area * typeMethods.r8vec_dot_product ( point_num, weight, value ) / scale;
        //
        //  Error:
        //
        double exact = 1.0;
        double quad_error = Math.Abs ( quad - exact );
            
        return quad_error;
    }

    public static double triangle01_poly_integral(int d, double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE01_POLY_INTEGRAL: polynomial integral over the unit triangle.
        //
        //  Discussion:
        //
        //    The unit triangle is T = ( (0,0), (1,0), (0,1) ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 April 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, integer D, the degree of the polynomial.
        //
        //    Input, double P[M], the polynomial coefficients.
        //    M = ((D+1)*(D+2))/2.
        //
        //    Output, double TRIANGLE01_POLY_INTEGRAL, the integral.
        //
    {
        int i = 0;
        int j = 0;
        int km1;

        int m = (d + 1) * (d + 2) / 2;

        double q = 0.0;
        for (km1 = 0; km1 < m; km1++)
        {
            int k = km1 + 1;
            typeMethods.i4_to_pascal(k, ref i, ref j);
            q += p[km1] * triangle01_monomial_integral(i, j);
        }

        return q;
    }
}