using System;

namespace Burkardt.Chebyshev;

public static class Chebyshev2
{
    public static void chebyshev2_compute_np ( int n, int np, double[] p, ref double[] x,
            ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHEBYSHEV2_COMPUTE_NP computes a Chebyshev type 2 quadrature rule.
        //
        //  Discussion:
        //
        //    The integral:
        //
        //      Integral ( -1 <= X <= 1 ) F(X)  sqrt ( 1 - x^2 )  dX
        //
        //    The quadrature rule:
        //
        //      Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Philip Davis, Philip Rabinowitz,
        //    Methods of Numerical Integration,
        //    Second Edition,
        //    Dover, 2007,
        //    ISBN: 0486453391,
        //    LC: QA299.3.D28.
        //
        //  Parameters:
        //
        //    Input, int N, the order.
        //    1 <= N.
        //
        //    Input, int NP, the number of parameters.
        //
        //    Input, double P[NP], parameters which are not needed by this function.
        //
        //    Output, double X[N], the abscissas.
        //
        //    Output, double W[N], the weights.
        //
    {
        chebyshev2_compute ( n, ref x, ref w );
    }
        
    public static void chebyshev2_compute ( int n, ref double[] x, ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHEBYSHEV2_COMPUTE computes a Chebyshev type 2 quadrature rule.
        //
        //  Discussion:
        //
        //    The integral:
        //
        //      integral ( -1 <= x <= 1 ) f(x)  sqrt ( 1 - x^2 )  dx
        //
        //    The quadrature rule:
        //
        //      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Philip Davis, Philip Rabinowitz,
        //    Methods of Numerical Integration,
        //    Second Edition,
        //    Dover, 2007,
        //    ISBN: 0486453391,
        //    LC: QA299.3.D28.
        //
        //  Parameters:
        //
        //    Input, int N, the order.
        //    1 <= N.
        //
        //    Output, double X[N], the abscissas.
        //
        //    Output, double W[N], the weights.
        //
    {
        double angle;
        int i;
            

        switch (n)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("CHEBYSHEV2_COMPUTE - Fatal error!");
                Console.WriteLine("  Illegal value of N = " + n + "");
                return;
        }

        for ( i = 0; i < n; i++ )
        {
            angle = Math.PI * (n - i) / (n + 1);
            w[i] = Math.PI / (n + 1) * Math.Pow ( Math.Sin ( angle ), 2 );
            x[i] = Math.Cos ( angle );
        }

        x[(n - 1) / 2] = (n % 2) switch
        {
            1 => 0.0,
            _ => x[(n - 1) / 2]
        };
    }
}