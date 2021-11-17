using System;
using Burkardt.IntegralNS;
using Burkardt.Types;

namespace Burkardt.Quadrature;

public static class CN_Jac
{
    public static void cn_jac_00_1(int n, double alpha, double beta, int o, ref double[] x,
            ref double[] w)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CN_JAC_00_1 implements the midpoint rule for region CN_JAC.
        //
        //  Discussion:
        //
        //    The rule has order O = 1.
        //
        //    The rule has precision P = 0.
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
        //    26 January 2010
        //
        //  Author:
        //
        //    John Burkardt
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
        int expon;
        int k;
        double volume;

        switch (alpha)
        {
            case <= -1.0:
                Console.WriteLine("");
                Console.WriteLine("CN_JAC_00_1 - Fatal error!");
                Console.WriteLine("  ALPHA <= -1.0");
                return;
        }

        switch (beta)
        {
            case <= -1.0:
                Console.WriteLine("");
                Console.WriteLine("CN_JAC_00_1 - Fatal error!");
                Console.WriteLine("  BETA <= -1.0");
                return;
        }

        expon = 0;
        volume = C1.c1_jac_monomial_integral(alpha, beta, expon);
        volume = Math.Pow(volume, n);

        typeMethods.r8vec_zero(n * o, ref x);

        k = -1;
        //
        //  1 point.
        //
        k += 1;
        w[k] = volume;

    }

    public static int cn_jac_00_1_size(int n, double alpha, double beta)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CN_JAC_00_1_SIZE sizes the midpoint rule for region CN_JAC.
        //
        //  Discussion:
        //
        //    The rule has order O = 1.
        //
        //    The rule has precision P = 0.
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
        //  Parameters:
        //
        //    Input, int N, the spatial dimension.
        //
        //    Input, double ALPHA, BETA, the parameters.
        //    -1.0 < ALPHA, -1.0 < BETA.
        //
        //    Output, int CN_JAC_00_1_SIZE the order.
        //
    {
        int o;

        switch (alpha)
        {
            case <= -1.0:
                Console.WriteLine("");
                Console.WriteLine("CN_JAC_00_1_SIZE - Fatal error!");
                Console.WriteLine("  ALPHA <= -1.0");
                return 1;
        }

        switch (beta)
        {
            case <= -1.0:
                Console.WriteLine("");
                Console.WriteLine("CN_JAC_00_1_SIZE - Fatal error!");
                Console.WriteLine("  BETA <= -1.0");
                return 1;
            default:
                o = 1;

                return o;
        }
    }

    public static void cn_jac_01_1(int n, double alpha, double beta, int o, ref double[] x,
            ref double[] w)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CN_JAC_01_1 implements a precision 1 rule for region CN_JAC.
        //
        //  Discussion:
        //
        //    The rule has order O = 1.
        //
        //    The rule has precision P = 1.
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
        //    26 January 2010
        //
        //  Author:
        //
        //    John Burkardt
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
        int expon;
        int i;
        int k;
        double value1;
        double value2;
        double volume;

        switch (alpha)
        {
            case <= -1.0:
                Console.WriteLine("");
                Console.WriteLine("CN_JAC_01_1 - Fatal error!");
                Console.WriteLine("  ALPHA <= -1.0");
                return;
        }

        switch (beta)
        {
            case <= -1.0:
                Console.WriteLine("");
                Console.WriteLine("CN_JAC_01_1 - Fatal error!");
                Console.WriteLine("  BETA <= -1.0");
                return;
        }

        expon = 0;
        value1 = C1.c1_jac_monomial_integral(alpha, beta, expon);
        volume = Math.Pow(value1, n);

        expon = 1;
        value2 = C1.c1_jac_monomial_integral(alpha, beta, expon);

        typeMethods.r8vec_zero(n * o, ref x);

        k = -1;
        //
        //  1 point.
        //
        k += 1;
        for (i = 0; i < n; i++)
        {
            x[i + k * n] = value2 / value1;
        }

        w[k] = volume;
    }

    public static int cn_jac_01_1_size(int n, double alpha, double beta)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CN_JAC_01_1_SIZE sizes a precision 1 rule for region CN_JAC.
        //
        //  Discussion:
        //
        //    The rule has order O = 1.
        //
        //    The rule has precision P = 1.
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
        //  Parameters:
        //
        //    Input, int N, the spatial dimension.
        //
        //    Input, double ALPHA, BETA, the parameters.
        //    -1.0 < ALPHA, -1.0 < BETA.
        //
        //    Output, int CN_JAC_01_1_SIZE, the order.
        //
    {
        int o;

        switch (alpha)
        {
            case <= -1.0:
                Console.WriteLine("");
                Console.WriteLine("CN_JAC_01_1_SIZE - Fatal error!");
                Console.WriteLine("  ALPHA <= -1.0");
                return 1;
        }

        switch (beta)
        {
            case <= -1.0:
                Console.WriteLine("");
                Console.WriteLine("CN_JAC_01_1_SIZE - Fatal error!");
                Console.WriteLine("  BETA <= -1.0");
                return 1;
            default:
                o = 1;

                return o;
        }
    }

}