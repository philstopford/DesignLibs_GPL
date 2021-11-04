using System;
using Burkardt.Types;

namespace Burkardt.Quadrature
{
    public static class Fejer2
    {
        public static double f2_abscissa(int order, int i)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    F2_ABSCISSA returns the I-th abscissa for the Fejer type 2 rule.
            //
            //  Discussion:
            //
            //    Our convention is that the abscissas are numbered from left to
            //    right.
            //
            //    This rule is defined on [-1,1].
            //
            //  Modified:
            //
            //    31 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int ORDER, the order of the Fejer type 2 rule.
            //    1 <= ORDER.
            //
            //    Input, int I, the index of the desired abscissa.  1 <= I <= ORDER.
            //
            //    Output, double F2_ABSCISSA, the value of the I-th 
            //    abscissa in the Fejer type 2 rule of order ORDER.
            //
        {
            double pi = 3.141592653589793;
            double value;

            if (order < 1)
            {
                value = -typeMethods.r8_huge();
                return value;
            }

            if (i < 1 || order < i)
            {
                Console.WriteLine("");
                Console.WriteLine("F2_ABSCISSA - Fatal error!");
                Console.WriteLine("  1 <= I <= ORDER is required.");
                Console.WriteLine("  I = " + i + "");
                Console.WriteLine("  ORDER = " + order + "");
                return(1);
            }

            if (order == 1)
            {
                value = 0.0;
            }
            else if (2 * (order + 1 - i) == order + 1)
            {
                value = 0.0;
            }
            else
            {
                value = Math.Cos((double) (order + 1 - i) * pi
                            / (double) (order + 1));
            }

            return value;
        }

        public static double[] f2_weights(int order)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    F2_WEIGHTS computes weights for a Fejer type 2 rule.
            //
            //  Modified:
            //
            //    28 May 2007
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
            //    Walter Gautschi,
            //    Numerical Quadrature in the Presence of a Singularity,
            //    SIAM Journal on Numerical Analysis,
            //    Volume 4, Number 3, 1967, pages 357-362.
            //
            //    Joerg Waldvogel,
            //    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
            //    BIT Numerical Mathematics,
            //    Volume 43, Number 1, 2003, pages 1-18.
            //
            //  Parameters:
            //
            //    Input, int ORDER, the order.
            //
            //    Output, double F2_WEIGHTS[ORDER], the weights.
            //
        {
            int i;
            int j;
            double p;
            double[] theta;
            double[] w;

            if (order < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("F2_WEIGHTS - Fatal error!");
                Console.WriteLine("  ORDER < 1.");
                return(null);
            }

            w = new double[order];

            if (order == 1)
            {
                w[0] = 2.0;
                return w;
            }
            else if (order == 2)
            {
                w[0] = 1.0;
                w[1] = 1.0;
                return w;
            }

            theta = new double[order];

            for (i = 1; i <= order; i++)
            {
                theta[i - 1] = (double) (order + 1 - i) * Math.PI
                               / (double) (order + 1);
            }

            for (i = 1; i <= order; i++)
            {
                w[i - 1] = 1.0;

                for (j = 1; j <= ((order - 1) / 2); j++)
                {
                    w[i - 1] = w[i - 1] - 2.0 * Math.Cos(2.0 * (double) (j) * theta[i - 1])
                        / (double) (4 * j * j - 1);
                }

                if (2 < order)
                {
                    p = 2.0 * (double) (((order + 1) / 2)) - 1.0;
                    w[i - 1] = w[i - 1] - Math.Cos((p + 1.0) * theta[i - 1]) / p;
                }

            }

            for (i = 0; i < order; i++)
            {
                w[i] = 2.0 * w[i] / (double) (order + 1);
            }

            return w;
        }

        public static void fejer2_compute_np ( int n, int np, double[] p, ref double[] x, ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FEJER2_COMPUTE_NP computes a Fejer type 2 rule.
        //
        //  Discussion:
        //
        //    Our convention is that the abscissas are numbered from left to right.
        //
        //    The rule is defined on [-1,1].
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
            fejer2_compute ( n, ref x, ref w );
        }
        
        public static void fejer2_compute ( int n, ref double[] x, ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FEJER2_COMPUTE computes a Fejer type 2 rule.
        //
        //  Discussion:
        //
        //    Our convention is that the abscissas are numbered from left to right.
        //
        //    The rule is defined on [-1,1].
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
            int i;
            int j;
            double p;
            double pi = 3.141592653589793;
            double theta;

            if ( n < 1 )
            {
                Console.WriteLine("");
                Console.WriteLine("FEJER2_COMPUTE - Fatal error!");
                Console.WriteLine("  Illegal value of N = " + n + "");
                return;
            }
            else if ( n == 1 )
            {
                x[0] = 0.0;
                w[0] = 2.0;
                return;
            }

            for ( i = 0; i < n; i++ )
            {
                x[i] =  Math.Cos ( ( double ) ( n - i ) * pi
                                   / ( double ) ( n + 1 ) );
            }
            if ( ( n % 2 ) == 1 )
            {
                x[(n-1)/2] = 0.0;
            }

            if ( n == 2 )
            {
                w[0] = 1.0;
                w[1] = 1.0;
            }
            else
            {
                for ( i = 0; i < n; i++ )
                {
                    theta = ( double ) ( n - i ) * pi
                            / ( double ) ( n + 1 );

                    w[i] = 1.0;

                    for ( j = 1; j <= ( ( n - 1 ) / 2 ); j++ )
                    {
                        w[i] = w[i] - 2.0 *  Math.Cos ( 2.0 * ( double ) ( j ) * theta )
                            / ( double ) ( 4 * j * j - 1 );
                    }
                    p = 2.0 * ( double ) ( ( ( n + 1 ) / 2 ) ) - 1.0;
                    w[i] = w[i] -  Math.Cos ( ( p + 1.0 ) * theta ) / p;
                }
                for ( i = 0; i < n; i++ )
                {
                    w[i] = 2.0 * w[i] / ( double ) ( n + 1 );
                }
            }
        }
    }
}