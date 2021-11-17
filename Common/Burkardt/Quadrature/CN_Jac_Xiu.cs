using System;
using Burkardt.IntegralNS;
using Burkardt.Types;

namespace Burkardt.Quadrature;

public static class CN_Jac_Xiu
{
    public static void cn_jac_02_xiu(int n, double alpha, double beta, int o, ref double[] x,
            ref double[] w)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CN_JAC_02_XIU implements the Xiu precision 2 rule for region CN_JAC.
        //
        //  Discussion:
        //
        //    The rule has order 
        //
        //      O = N + 1.
        //
        //    The rule has precision P = 2.
        //
        //    CN is the cube [-1,+1]^N with the Jacobi (beta) weight function
        //
        //      w(alpha,beta;x) = product ( 1 <= i <= n ) (1-x(i))^beta (1+x(i))^alpha.
        //
        //    with -1 < alpha, -1 < beta.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 March 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Dongbin Xiu,
        //    Numerical integration formulas of degree two,
        //    Applied Numerical Mathematics,
        //    Volume 58, 2008, pages 1515-1520.
        //
        //  Parameters:
        //
        //    Input, int N, the spatial dimension.
        //
        //    Input, double ALPHA, BETA, the parameters.
        //    -1.0 < ALPHA, -1.0 < BETA.
        //
        //    Input, int O, the order.
        //
        //    Output, double X[N*O], the abscissas.
        //
        //    Output, double W[O], the weights.
        //
    {
        double arg;
        double c1;
        double delta0;
        int expon;
        double gamma0;
        int i;
        int j;
            
        int r;
        double volume;
        double volume_1d;

        switch (alpha)
        {
            case <= -1.0:
                Console.WriteLine("");
                Console.WriteLine("CN_JAC_02_XIU - Fatal error!");
                Console.WriteLine("  ALPHA <= -1.0");
                return;
        }

        switch (beta)
        {
            case <= -1.0:
                Console.WriteLine("");
                Console.WriteLine("CN_JAC_02_XIU - Fatal error!");
                Console.WriteLine("  BETA <= -1.0");
                return;
        }

        for (j = 0; j < o; j++)
        {
            i = 0;
            for (r = 1; r <= n / 2; r++)
            {
                arg = 2 * r * j * Math.PI / (n + 1);

                x[i + j * n] = Math.Sqrt(2.0) * Math.Cos(arg);
                i += 1;
                x[i + j * n] = Math.Sqrt(2.0) * Math.Sin(arg);
                i += 1;
            }

            if (i < n)
            {
                x[i + j * n] = typeMethods.r8_mop(j);
                i += 1;
            }
        }

        gamma0 = (alpha + beta + 2.0) / 2.0;
        delta0 = (alpha - beta) / 2.0;
        c1 = 2.0 * (alpha + 1.0) * (beta + 1.0) / (alpha + beta + 3.0)
                                                / (alpha + beta + 2.0);

        for (j = 0; j < o; j++)
        {
            for (i = 0; i < n; i++)
            {
                x[i + j * n] = (Math.Sqrt(gamma0 * c1) * x[i + j * n] - delta0) / gamma0;
            }
        }

        expon = 0;
        volume_1d = C1.c1_jac_monomial_integral(alpha, beta, expon);
        volume = Math.Pow(volume_1d, n);

        for (j = 0; j < o; j++)
        {
            w[j] = volume / o;
        }
    }

    public static int cn_jac_02_xiu_size(int n, double alpha, double beta)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CN_JAC_02_XIU_SIZE sizes the Xiu rule for region CN_JAC.
        //
        //  Discussion:
        //
        //    The rule has order 
        //
        //      O = N + 1.
        //
        //    The rule has precision P = 2.
        //
        //    CN is the cube [-1,+1]^N with the Jacobi (beta) weight function
        //
        //      w(alpha,beta;x) = product ( 1 <= i <= n ) (1-x(i))^beta (1+x(i))^alpha.
        //
        //    with -1 < alpha, -1 < beta.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 January 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Dongbin Xiu,
        //    Numerical integration formulas of degree two,
        //    Applied Numerical Mathematics,
        //    Volume 58, 2008, pages 1515-1520.
        //
        //  Parameters:
        //
        //    Input, int N, the spatial dimension.
        //
        //    Input, double ALPHA, BETA, the parameters.
        //    -1.0 < ALPHA, -1.0 < BETA.
        //
        //    Output, int CN_JAC_02_XIU_SIZE, the order.
        //
    {
        int o;

        switch (alpha)
        {
            case <= -1.0:
                Console.WriteLine("");
                Console.WriteLine("CN_JAC_02_XIU_SIZE - Fatal error!");
                Console.WriteLine("  ALPHA <= -1.0");
                return 1;
        }

        switch (beta)
        {
            case <= -1.0:
                Console.WriteLine("");
                Console.WriteLine("CN_JAC_02_XIU_SIZE - Fatal error!");
                Console.WriteLine("  BETA <= -1.0");
                return 1;
            default:
                o = n + 1;

                return o;
        }
    }

}