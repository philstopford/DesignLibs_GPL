using System;
using System.Globalization;
using Burkardt.IntegralNS;
using Burkardt.Types;

namespace Burkardt.ExactnessNS;

public static partial class Exactness
{
    public static void legendre_exactness(int n, double[] x, double[] w, int p_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEGENDRE_EXACTNESS investigates exactness of Legendre quadrature.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of points in the rule.
        //
        //    Input, double X[N], the quadrature points.
        //
        //    Input, double W[N], the quadrature weights.
        //
        //    Input, int P_MAX, the maximum exponent.
        //    0 <= P_MAX.
        //
    {
        double e;
        int i;
        int p;
        double q;
        double s;
        double[] v;

        Console.WriteLine("");
        Console.WriteLine("  Quadrature rule for the Legendre integral.");
        Console.WriteLine("  Rule of order N = " + n + "");
        Console.WriteLine("  Degree          Relative Error");
        Console.WriteLine("");

        v = new double[n];

        for (p = 0; p <= p_max; p++)
        {
            s = Integral.legendre_integral(p);

            for (i = 0; i < n; i++)
            {
                v[i] = Math.Pow(x[i], p);
            }

            q = typeMethods.r8vec_dot_product(n, w, v);

            e = s switch
            {
                0.0 => Math.Abs(q),
                _ => Math.Abs(q - s) / Math.Abs(s)
            };

            Console.WriteLine(p.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                                                      + e.ToString(CultureInfo.InvariantCulture).PadLeft(24) + "");
        }
    }

    public static void legendre_3d_exactness(double[] a, double[] b, int n, double[] x,
            double[] y, double[] z, double[] w, int t)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEGENDRE_3D_EXACTNESS: monomial exactness for the 3D Legendre integral.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 August 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A[3], the lower limits of integration.
        //
        //    Input, double B[3], the upper limits of integration.
        //
        //    Input, int N, the number of points in the rule.
        //
        //    Input, double X[N], Y[N], Z[N], the quadrature points.
        //
        //    Input, double W[N], the quadrature weights.
        //
        //    Input, int T, the maximum total degree.
        //    0 <= T.
        //
    {
        double e;
        int i;
        int j;
        int k;
        int l;
        int[] p = new int[3];
        double q;
        double s;
        int tt;
        double[] v;

        v = new double[n];

        Console.WriteLine("");
        Console.WriteLine("  Quadrature rule for the 3D Legendre integral.");
        Console.WriteLine("  Number of points in rule is " + n + "");
        Console.WriteLine("");
        Console.WriteLine("   D   I       J       K          Relative Error");

        for (tt = 0; tt <= t; tt++)
        {
            Console.WriteLine(tt.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "");

            for (k = 0; k <= tt; k++)
            {
                for (j = 0; j <= tt - k; j++)
                {
                    i = tt - j - k;

                    p[0] = i;
                    p[1] = j;
                    p[2] = k;

                    s = legendre_3d_monomial_integral(a, b, p);

                    for (l = 0; l < n; l++)
                    {
                        v[l] = Math.Pow(x[l], p[0]) * Math.Pow(y[l], p[1]) * Math.Pow(z[l], p[2]);
                    }

                    q = typeMethods.r8vec_dot_product(n, w, v);

                    e = s switch
                    {
                        0.0 => Math.Abs(q),
                        _ => Math.Abs(q - s) / Math.Abs(s)
                    };

                    Console.WriteLine(p[0].ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                                                                 + p[1].ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                                                                 + p[2].ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                                                                 + e.ToString(CultureInfo.InvariantCulture).PadLeft(24) + "");
                }
            }
        }
    }

    public static double legendre_3d_monomial_integral(double[] a, double[] b, int[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEGENDRE_3D_MONOMIAL_INTEGRAL the Legendre integral of a monomial.
        //
        //  Discussion:
        //
        //    The Legendre integral to be evaluated has the form
        //
        //      I(f) = integral ( z1 <= z <= z2 )
        //             integral ( y1 <= y <= y2 ) 
        //             integral ( x1 <= x <= x2 ) x^i y^j z^k dx dy dz
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 August 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A[3], the lower limits of integration.
        //
        //    Input, double B[3], the upper limits of integration.
        //
        //    Input, int P[3], the exponents of X and Y.
        //
        //    Output, double LEGENDRE_3D_MONOMIAL_INTEGRAL, the value of the exact integral.
        //
    {
        double value = 0;

        value = (Math.Pow(b[0], p[0] + 1) - Math.Pow(a[0], p[0] + 1)) / (p[0] + 1)
            * (Math.Pow(b[1], p[1] + 1) - Math.Pow(a[1], p[1] + 1)) / (p[1] + 1)
            * (Math.Pow(b[2], p[2] + 1) - Math.Pow(a[2], p[2] + 1)) / (p[2] + 1);

        return value;
    }

}