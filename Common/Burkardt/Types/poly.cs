using System;
using Burkardt.PolynomialNS;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double[] poly_power(int d1, double[] p1, int n )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLY_POWER computes a power of a polynomial.
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
        //    Input, int D1, the degree of the polynomial.
        //
        //    Input, double P1(M1), the polynomial coefficients.
        //    M1 = ((D1+1)*(D1+2))/2.
        //
        //    Input, int N, the nonnegative integer power.
        //
        //    Output, double POLY_POWER[M2], the polynomial power.
        //    D2 = N * D1.
        //    M2 = ((D2+1)*(D2+2))/2.
        //
    {
        int i;
        //
        //  Set D2 to 0, to indicate that P2 currently contains only
        //  a constant term.
        //
        int d2 = 0;
        int m2 = (1) * (2) / 2;
        double[] p2 = new double[m2];
        p2[0] = 1.0;
        //
        //  Iterate N times:
        //    P3 <= P1 * P2
        //    P2 <= P3
        //
        for (i = 1; i <= n; i++)
        {
            int d3 = d1 + d2;
            double[] p3 = poly_product(d1, p1, d2, p2);
                
            d2 = d3;
            p2 = p3;
        }

        return p2;
    }

    public static double[] poly_power_linear(int d1, double[] p1, int n )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLY_POWER_LINEAR computes the polynomial ( a + b*x + c*y ) ^ n.
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
        //    Input, int D1, the degree of the linear polynomial,
        //    which should be 1 (or possibly 0).
        //
        //    Input, double P1(M1), the coefficients of the linear polynomial.
        //    M1 = ( (D1+1)*(D1+2) ) / 2, which should be 3.
        //
        //    Input, int N, the power to which the polynomial is to be 
        //    raised.  0 <= N.
        //
        //    Output, double P2(M2), the coefficients of the power polynomial.
        //    D2 = N * D1;
        //    M2 = ( (D2+1)*(D2+2) ) / 2.
        //
    {
        int i;

        switch (d1)
        {
            case < 0:
                Console.WriteLine("");
                Console.WriteLine("POLY_POWER_LINEAR - Fatal error!");
                Console.WriteLine("  D1 < 0.");
                return null;
        }

        switch (n)
        {
            case < 0:
                Console.WriteLine("");
                Console.WriteLine("POLY_POWER_LINEAR - Fatal error!");
                Console.WriteLine("  N < 0.");
                return null;
        }

        int d2 = n * d1;
        int m2 = (d2 + 1) * (d2 + 2) / 2;
        double[] p2 = new double[m2];

        switch (d1)
        {
            case 0:
                p2[0] = Math.Pow(p1[0], n);
                return p2;
        }

        switch (n)
        {
            case 0:
                p2[0] = 1.0;
                return p2;
        }

        //
        //  Use the Trinomial formula.
        //
        for (i = 0; i <= n; i++)
        {
            int j;
            for (j = 0; j <= n - i; j++)
            {
                int k;
                for (k = 0; k <= n - i - j; k++)
                {
                    //
                    //  We store X^J Y^K in location L.
                    //
                    int l = pascal_to_i4(j, k);
                    p2[l - 1] = Trinomial.trinomial(i, j, k)
                                * Math.Pow(p1[0], i) * Math.Pow(p1[1], j) * Math.Pow(p1[2], k);
                }
            }
        }

        return p2;
    }

    public static void poly_print(int d, ref double[] p, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLY_PRINT prints an XY polynomial.
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
        //    Input, int D, the degree of the polynomial.            //
        //    Output, double P[M], the coefficients of all monomials of 
        //    degree 0 through D.  P must contain ((D+1)*(D+2))/2 entries.
        //
        //    Input, string TITLE, a title string.
        //
    {
        int i = 0;
        int j = 0;
        int km1;

        int m = (d + 1) * (d + 2) / 2;

        for (km1 = 0; km1 < m; km1++)
        {
            if (p[km1] != 0.0)
            {
                break;
            }

            Console.WriteLine(title + " = 0");
            return;
        }

        Console.WriteLine(title + "");

        for (km1 = 0; km1 < m; km1++)
        {
            int k = km1 + 1;
            i4_to_pascal(k, ref i, ref j);

            string cout = "";

            if (p[km1] == 0.0)
            {
                continue;
            }

            cout += p[km1] switch
            {
                < 0.0 => "  -" + Math.Abs(p[km1]),
                _ => "  +" + p[km1]
            };

            if (i + j != 0)
            {
                cout += " ";
            }

            switch (i)
            {
                case 0:
                    break;
                case 1:
                    cout += "x";
                    break;
                default:
                    cout += "x^" + i;
                    break;
            }

            switch (j)
            {
                case 0:
                    break;
                case 1:
                    cout += "y";
                    break;
                default:
                    cout += "y^" + j;
                    break;
            }

            Console.WriteLine(cout);
        }
    }

    public static double[] poly_product(int d1, double[] p1, int d2, double[] p2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLY_PRODUCT computes P3(x,y) = P1(x,y) * P2(x,y) for polynomials.
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
        //    Input, int D1, the degree of factor 1.
        //
        //    Input, double P1[M1], the factor 1 coefficients.
        //    M1 = ((D1+1)*(D1+2))/2.
        //
        //    Input, int D2, the degree of factor 2.
        //
        //    Input, double P2[M2], the factor2 coefficients.
        //    M2 = ((D2+1)*(D2+2))/2.
        //
        //    Output, double POLY_PRODUCT[M3], the result coefficients.
        //    D3 = D1 + D2;
        //    M3 = ((D3+1)*(D3+2))/2.
        //
    {
        int i1 = 0;
        int i2 = 0;
        int j1 = 0;
        int j2 = 0;
        int k1m1;
        int k3m1;

        int m1 = (d1 + 1) * (d1 + 2) / 2;
        int m2 = (d2 + 1) * (d2 + 2) / 2;
        //
        //  Consider each entry in P1:
        //    P1(K1) * X^I1 * Y^J1
        //  and multiply it by each entry in P2:
        //    P2(K2) * X^I2 * Y^J2
        //  getting 
        //    P3(K3) = P3(K3) + P1(K1) * P2(X2) * X^(I1+I2) * Y(J1+J2)
        //
        int d3 = d1 + d2;
        int m3 = (d3 + 1) * (d3 + 2) / 2;
        double[] p3 = new double[m3];

        for (k3m1 = 0; k3m1 < m3; k3m1++)
        {
            p3[k3m1] = 0.0;
        }

        for (k1m1 = 0; k1m1 < m1; k1m1++)
        {
            int k1 = k1m1 + 1;
            i4_to_pascal(k1, ref i1, ref j1);
            int k2m1;
            for (k2m1 = 0; k2m1 < m2; k2m1++)
            {
                int k2 = k2m1 + 1;
                i4_to_pascal(k2, ref i2, ref j2);
                int i3 = i1 + i2;
                int j3 = j1 + j2;
                int k3 = pascal_to_i4(i3, j3);
                k3m1 = k3 - 1;
                p3[k3m1] += p1[k1m1] * p2[k2m1];
            }
        }

        return p3;
    }

    public static double[] polygon_centroid_2d(int n, double[] v)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_CENTROID_2D computes the centroid of a polygon in 2D.
        //
        //  Formula:
        //
        //    Denoting the centroid coordinates by CENTROID, then
        //
        //      CENTROID(1) = Integral ( Polygon interior ) x dx dy / Area ( Polygon )
        //      CENTROID(2) = Integral ( Polygon interior ) y dx dy / Area ( Polygon ).
        //
        //    Green's theorem states that
        //
        //      Integral ( Polygon boundary ) ( M dx + N dy ) =
        //      Integral ( Polygon interior ) ( dN/dx - dM/dy ) dx dy.
        //
        //    Using M = 0 and N = x * x / 2, we get:
        //
        //      CENTROID(1) = 0.5 * Integral ( Polygon boundary ) x * x dy,
        //
        //    which becomes
        //
        //      CENTROID(1) = 1/6 Sum ( 1 <= I <= N )
        //        ( X(I+1) + X(I) ) * ( X(I) * Y(I+1) - X(I+1) * Y(I))
        //
        //    where, when I = N, the index "I+1" is replaced by 1.
        //
        //    A similar calculation gives us a formula for CENTROID(2).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 August 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Gerard Bashein and Paul Detmer,
        //    Centroid of a Polygon,
        //    Graphics Gems IV, edited by Paul Heckbert,
        //    AP Professional, 1994.
        //
        //  Parameters:
        //
        //    Input, int N, the number of sides of the polygonal shape.
        //
        //    Input, double V[2*N], the coordinates of the vertices
        //    of the shape.
        //
        //    Output, double POLYGON_CENTROID_2D[2], the coordinates of the
        //    centroid of the shape.
        //
    {
        int i;
        //
        double area = 0.0;
        double[] centroid = new double[2];
        centroid[0] = 0.0;
        centroid[1] = 0.0;

        for (i = 0; i < n; i++)
        {
            int ip1;
            if (i < n - 1)
            {
                ip1 = i + 1;
            }
            else
            {
                ip1 = 0;
            }

            double temp = v[0 + i * 2] * v[1 + ip1 * 2] - v[0 + ip1 * 2] * v[1 + i * 2];

            area += temp;

            centroid[0] += (v[0 + ip1 * 2] + v[0 + i * 2]) * temp;
            centroid[1] += (v[1 + ip1 * 2] + v[1 + i * 2]) * temp;

        }

        area /= 2.0;

        centroid[0] /= 6.0 * area;
        centroid[1] /= 6.0 * area;

        return centroid;
    }
}