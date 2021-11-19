using System;
using Burkardt.Types;

namespace Burkardt.Quadrature;

public static class Epn_Lag
{
    public static void epn_glg_00_1(int n, double alpha, int o, ref double[] x, ref double[] w)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EPN_GLG_00_1 implements the "midpoint rule" for region EPN_GLG.
        //
        //  Discussion:
        //
        //    The rule has order O = 1.
        //
        //    The rule has precision P = 0.
        //
        //    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
        //    Laguerre weight function:
        //
        //      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 January 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the spatial dimension.
        //
        //    Input, double ALPHA, the exponent of X in the weight function.
        //    -1.0 < ALPHA.
        //
        //    Input, int O, the order.
        //
        //    Output, double X[N*O], the abscissas.
        //
        //    Output, double W[O], the weights.
        //
    {
        int i;

        switch (alpha)
        {
            case <= -1.0:
                Console.WriteLine("");
                Console.WriteLine("EPN_GLG_00_1 - Fatal error!");
                Console.WriteLine("  ALPHA <= -1.0");
                return;
        }

        const int expon = 0;
        double volume = IntegralNS.Epn_Lag.ep1_glg_monomial_integral(expon, alpha);
        volume = Math.Pow(volume, n);

        typeMethods.r8vec_zero(n * o, ref x);

        int k = -1;
        //
        //  1 point.
        //
        k += 1;
        for (i = 0; i < n; i++)
        {
            x[i + k * n] = 1.0;
        }

        w[k] = volume;
    }

    public static int epn_glg_00_1_size(int n, double alpha)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EPN_GLG_00_1_SIZE sizes the midpoint rule for region EPN_GLG.
        //
        //  Discussion:
        //
        //    The rule has order O = 1.
        //
        //    The rule has precision P = 0.
        //
        //    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
        //    Laguerre weight function:
        //
        //      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 January 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the spatial dimension.
        //
        //    Input, double ALPHA, the exponent of X in the weight function.
        //    -1.0 < ALPHA.
        //
        //    Output, int EPN_GLG_00_1_SIZE, the order.
        //
    {
        switch (alpha)
        {
            case <= -1.0:
                Console.WriteLine("");
                Console.WriteLine("EPN_GLG_00_1_SIZE - Fatal error!");
                Console.WriteLine("  ALPHA <= -1.0");
                return 1;
            default:
                const int o = 1;

                return o;
        }
    }

    public static void epn_glg_01_1(int n, double alpha, int o, ref double[] x, ref double[] w)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EPN_GLG_01_1 implements a precision 1 rule for region EPN_GLG.
        //
        //  Discussion:
        //
        //    The rule has order O = 1.
        //
        //    The rule has precision P = 1.
        //
        //    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
        //    Laguerre weight function:
        //
        //      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 January 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the spatial dimension.
        //
        //    Input, double ALPHA, the exponent of X in the weight function.
        //    -1.0 < ALPHA.
        //
        //    Input, int O, the order.
        //
        //    Output, double X[N*O], the abscissas.
        //
        //    Output, double W[O], the weights.
        //
    {
        int i;

        switch (alpha)
        {
            case <= -1.0:
                Console.WriteLine("");
                Console.WriteLine("EPN_GLG_01_1 - Fatal error!");
                Console.WriteLine("  ALPHA <= -1.0");
                return;
        }

        int expon = 0;
        double value1 = IntegralNS.Epn_Lag.ep1_glg_monomial_integral(expon, alpha);
        double volume = Math.Pow(value1, n);

        expon = 1;
        double value2 = IntegralNS.Epn_Lag.ep1_glg_monomial_integral(expon, alpha);

        typeMethods.r8vec_zero(n * o, ref x);

        int k = -1;
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

    public static int epn_glg_01_1_size(int n, double alpha)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EPN_GLG_01_1_SIZE sizes a precision 1 rule for region EPN_GLG.
        //
        //  Discussion:
        //
        //    The rule has order O = 1.
        //
        //    The rule has precision P = 1.
        //
        //    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
        //    Laguerre weight function:
        //
        //      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 January 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the spatial dimension.
        //
        //    Input, double ALPHA, the exponent of X in the weight function.
        //    -1.0 < ALPHA.
        //
        //    Output, int EPN_GLG_01_1_SIZE, the order.
        //
    {
        switch (alpha)
        {
            case <= -1.0:
                Console.WriteLine("");
                Console.WriteLine("EPN_GLG_01_1_SIZE - Fatal error!");
                Console.WriteLine("  ALPHA <= -1.0");
                return 1;
            default:
                const int o = 1;

                return o;
        }
    }

    public static void epn_glg_02_xiu(int n, double alpha, int o, ref double[] x, ref double[] w)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EPN_GLG_02_XIU implements the Xiu precision 2 rule for region EPN_GLG.
        //
        //  Discussion:
        //
        //    The rule has order 
        //
        //      O = N + 1.
        //
        //    The rule has precision P = 2.
        //
        //    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
        //    Laguerre weight function:
        //
        //      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
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
        //    Input, double ALPHA, the exponent of X in the weight function.
        //    -1.0 < ALPHA.
        //
        //    Input, int O, the order.
        //
        //    Output, double X[N*O], the abscissas.
        //
        //    Output, double W[O], the weights.
        //
    {
        int i;
        int j;

        switch (alpha)
        {
            case <= -1.0:
                Console.WriteLine("");
                Console.WriteLine("EPN_GLG_02_XIU - Fatal error!");
                Console.WriteLine("  ALPHA <= -1.0");
                return;
        }

        for (j = 0; j < o; j++)
        {
            i = 0;
            int r;
            for (r = 1; r <= n / 2; r++)
            {
                double arg = 2 * r * j * Math.PI / (n + 1);

                x[i + j * n] = Math.Sqrt(2.0) * Math.Cos(arg);
                i += 1;
                x[i + j * n] = Math.Sqrt(2.0) * Math.Sin(arg);
                i += 1;
            }

            if (i >= n)
            {
                continue;
            }

            x[i + j * n] = typeMethods.r8_mop(j);
        }

        double gamma0 = -1.0;
        double delta0 = alpha + 1.0;
        double c1 = -alpha - 1.0;

        for (j = 0; j < o; j++)
        {
            for (i = 0; i < n; i++)
            {
                x[i + j * n] = (Math.Sqrt(gamma0 * c1) * x[i + j * n] - delta0) / gamma0;
            }
        }

        int expon = 0;
        double volume_1d = IntegralNS.Epn_Lag.ep1_glg_monomial_integral(expon, alpha);
        double volume = Math.Pow(volume_1d, n);

        for (j = 0; j < o; j++)
        {
            w[j] = volume / o;
        }
    }

    public static int epn_glg_02_xiu_size(int n, double alpha)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EPN_GLG_02_XIU_SIZE sizes the Xiu rule for region EPN_GLG.
        //
        //  Discussion:
        //
        //    The rule has order 
        //
        //      O = N + 1.
        //
        //    The rule has precision P = 2.
        //
        //    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
        //    Laguerre weight function:
        //
        //      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 January 2010
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
        //    Input, double ALPHA, the exponent of X in the weight function.
        //    -1.0 < ALPHA.
        //
        //    Output, int EPN_GLG_02_XIU_SIZE, the order.
        //
    {
        switch (alpha)
        {
            case <= -1.0:
                Console.WriteLine("");
                Console.WriteLine("EPN_GLG_02_XIUI_SIZE - Fatal error!");
                Console.WriteLine("  ALPHA <= -1.0");
                return 1;
            default:
                int o = n + 1;

                return o;
        }
    }

    public static void epn_lag_00_1(int n, int o, ref double[] x, ref double[] w)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EPN_LAG_00_1 implements the "midpoint rule" for region EPN_LAG.
        //
        //  Discussion:
        //
        //    The rule has order O = 1.
        //
        //    The rule has precision P = 0.
        //
        //    EPN is the N-dimensional positive space [0,+oo)^N with exponential 
        //    or Laguerre weight function:
        //
        //      w(x(1:n)) = exp ( - sum ( x(1:n) ) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 January 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the spatial dimension.
        //
        //    Input, int O, the order.
        //
        //    Output, double X[N*O], the abscissas.
        //
        //    Output, double W[O], the weights.
        //
    {
        int i;

        const int expon = 0;
        double volume = IntegralNS.Epn_Lag.ep1_lag_monomial_integral(expon);
        volume = Math.Pow(volume, n);

        typeMethods.r8vec_zero(n * o, ref x);

        int k = -1;
        //
        //  1 point.
        //
        k += 1;
        for (i = 0; i < n; i++)
        {
            x[i + k * n] = 1.0;
        }

        w[k] = volume;
    }

    public static int epn_lag_00_1_size(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EPN_LAG_00_1_SIZE sizes the midpoint rule for region EPN_LAG.
        //
        //  Discussion:
        //
        //    The rule has order O = 1.
        //
        //    The rule has precision P = 0.
        //
        //    EPN is the N-dimensional positive space [0,+oo)^N with exponential 
        //    or Laguerre weight function:
        //
        //      w(x(1:n)) = exp ( - sum ( x(1:n) ) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 January 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the spatial dimension.
        //
        //    Output, int EPN_LAG_00_1_SIZE, the order.
        //
    {
        int o = 1;

        return o;
    }

    public static void epn_lag_01_1(int n, int o, ref double[] x, ref double[] w)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EPN_LAG_01_1 implements a precision 1 rule for region EPN_LAG.
        //
        //  Discussion:
        //
        //    The rule has order O = 1.
        //
        //    The rule has precision P = 1.
        //
        //    EPN is the N-dimensional positive space [0,+oo)^N with exponential 
        //    or Laguerre weight function:
        //
        //      w(x(1:n)) = exp ( - sum ( x(1:n) ) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 January 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the spatial dimension.
        //
        //    Input, int O, the order.
        //
        //    Output, double X[N*O], the abscissas.
        //
        //    Output, double W[O], the weights.
        //
    {
        int i;

        int expon = 0;
        double value1 = IntegralNS.Epn_Lag.ep1_lag_monomial_integral(expon);
        double volume = Math.Pow(value1, n);

        expon = 1;
        double value2 = IntegralNS.Epn_Lag.ep1_lag_monomial_integral(expon);

        typeMethods.r8vec_zero(n * o, ref x);

        int k = -1;
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

    public static int epn_lag_01_1_size(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EPN_LAG_01_1_SIZE sizes a precision 1 rule for region EPN_LAG.
        //
        //  Discussion:
        //
        //    The rule has order O = 1.
        //
        //    The rule has precision P = 1.
        //
        //    EPN is the N-dimensional positive space [0,+oo)^N with exponential 
        //    or Laguerre weight function:
        //
        //      w(x(1:n)) = exp ( - sum ( x(1:n) ) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 January 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the spatial dimension.
        //
        //    Output, int EPN_LOG_01_1_SIZE, the order.
        //
    {
        const int o = 1;

        return o;
    }

    public static void epn_lag_02_xiu(int n, int o, ref double[] x, ref double[] w)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EPN_LAG_02_XIU implements the Xiu precision 2 rule for region EPN_LAG.
        //
        //  Discussion:
        //
        //    The rule has order 
        //
        //      O = N + 1.
        //
        //    The rule has precision P = 2.
        //
        //    EPN is the N-dimensional positive space [0,+oo)^N with exponential 
        //    or Laguerre weight function:
        //
        //      w(x(1:n)) = exp ( - sum ( x(1:n) ) )
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
        //    Input, int O, the order.
        //
        //    Output, double X[N*O], the abscissas.
        //
        //    Output, double W[O], the weights.
        //
    {
        int i;
        int j;

        for (j = 0; j < o; j++)
        {
            i = 0;
            int r;
            for (r = 1; r <= n / 2; r++)
            {
                double arg = 2 * r * j * Math.PI / (n + 1);

                x[i + j * n] = Math.Sqrt(2.0) * Math.Cos(arg);
                i += 1;
                x[i + j * n] = Math.Sqrt(2.0) * Math.Sin(arg);
                i += 1;
            }

            if (i < n)
            {
                x[i + j * n] = typeMethods.r8_mop(j);
            }
        }

        const double gamma0 = -1.0;
        const double delta0 = 1.0;
        const double c1 = -1.0;

        for (j = 0; j < o; j++)
        {
            for (i = 0; i < n; i++)
            {
                x[i + j * n] = (Math.Sqrt(gamma0 * c1) * x[i + j * n] - delta0) / gamma0;
            }
        }

        int expon = 0;
        double volume_1d = IntegralNS.Epn_Lag.ep1_lag_monomial_integral(expon);
        double volume = Math.Pow(volume_1d, n);

        for (j = 0; j < o; j++)
        {
            w[j] = volume / o;
        }
    }

    public static int epn_lag_02_xiu_size(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EPN_LAG_02_XIU_SIZE sizes the Xiu rule for region EPN_LAG.
        //
        //  Discussion:
        //
        //    The rule has order 
        //
        //      O = N + 1.
        //
        //    The rule has precision P = 2.
        //
        //    EPN is the N-dimensional positive space [0,+oo)^N with exponential 
        //    or Laguerre weight function:
        //
        //      w(x(1:n)) = exp ( - sum ( x(1:n) ) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 January 2010
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
        //    Output, int EPN_LAG_02_XIU_SIZE, the order.
        //
    {
        int o = n + 1;

        return o;
    }

}