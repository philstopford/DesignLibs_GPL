using System;
using Burkardt.IntegralNS;
using Burkardt.Types;

namespace Burkardt.ExactnessNS;

public static partial class Exactness
{
    public static void gegenbauer_exactness(int n, double[] x, double[] w, int p_max,
            double lambda )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GEGENBAUER_EXACTNESS investigates exactness of Gegenbauer quadrature.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 January 2016
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
        //    Input, double LAMBDA, the parameter.
        //    -1/2 < LAMBDA.
        //
    {
        double e;
        int i;
        int p;
        double q;
        double s;
        double[] v;

        Console.WriteLine("");
        Console.WriteLine("  Quadrature rule for the Hermite integral.");
        Console.WriteLine("  Rule of order N = " + n + "");
        Console.WriteLine("  Lambda = " + lambda + "");
        Console.WriteLine("");
        Console.WriteLine("  Degree          Relative Error");
        Console.WriteLine("");

        v = new double[n];

        for (p = 0; p <= p_max; p++)
        {
            s = Integral.gegenbauer_integral(p, lambda);

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

            Console.WriteLine(p.ToString().PadLeft(6) + "  "
                                                      + e.ToString().PadLeft(24) + "");
        }
    }
}