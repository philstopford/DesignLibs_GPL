using System;
using Burkardt.Table;
using Burkardt.Types;

namespace Burkardt
{
    public static class ClenshawCurtis
    {
        public static void clenshaw_curtis_compute(int order, ref double[] x, ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CLENSHAW_CURTIS_COMPUTE computes a Clenshaw Curtis quadrature rule.
        //
        //  Discussion:
        //
        //    The integration interval is [ -1, 1 ].
        //
        //    The weight function is w(x) = 1.0.
        //
        //    The integral to approximate:
        //
        //      Integral ( -1 <= X <= 1 ) F(X) dX
        //
        //    The quadrature rule:
        //
        //      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 March 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ORDER, the order of the rule.
        //    1 <= ORDER.
        //
        //    Output, double X[ORDER], the abscissas.
        //
        //    Output, double W[ORDER], the weights.
        //
        {
            double b;
            int i;
            int j;
            double theta;

            if (order < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("CLENSHAW_CURTIS_COMPUTE - Fatal error!");
                Console.WriteLine("  Illegal value of ORDER = " + order + "");
                return;
            }
            else if (order == 1)
            {
                x[0] = 0.0;
                w[0] = 2.0;
            }
            else
            {
                for (i = 0; i < order; i++)
                {
                    x[i] = Math.Cos((double) (order - 1 - i) * Math.PI
                                    / (double) (order - 1));
                }

                x[0] = -1.0;
                if ((order % 2) == 1)
                {
                    x[(order - 1) / 2] = 0.0;
                }

                x[order - 1] = +1.0;

                for (i = 0; i < order; i++)
                {
                    theta = (double) (i) * Math.PI / (double) (order - 1);

                    w[i] = 1.0;

                    for (j = 1; j <= (order - 1) / 2; j++)
                    {
                        if (2 * j == (order - 1))
                        {
                            b = 1.0;
                        }
                        else
                        {
                            b = 2.0;
                        }

                        w[i] = w[i] - b * Math.Cos(2.0 * (double) (j) * theta)
                            / (double) (4 * j * j - 1);
                    }
                }

                w[0] = w[0] / (double) (order - 1);
                for (i = 1; i < order - 1; i++)
                {
                    w[i] = 2.0 * w[i] / (double) (order - 1);
                }

                w[order - 1] = w[order - 1] / (double) (order - 1);
            }

            return;
        }

        public static double[] cc_compute_points ( int n )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CC_COMPUTE_POINTS computes Clenshaw Curtis quadrature points.
            //
            //  Discussion:
            //
            //    Our convention is that the abscissas are numbered from left to right.
            //
            //    This rule is defined on [-1,1].
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    30 September 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the rule.
            //
            //    Output, double CC_COMPUTE_POINTS[N], the abscissas.
            //
        {
            int i;
            double[] x;

            if ( n < 1 )
            {
                Console.WriteLine("");
                Console.WriteLine("CC_COMPUTE_POINTS - Fatal error!");
                Console.WriteLine("  N < 1.");
                return null;
            }

            x = new double[n];

            if ( n == 1 )
            {
                x[0] = 0.0;
            }
            else
            {
                for ( i = 1; i <= n; i++ )
                {
                    x[i-1] =  Math.Cos ( ( double ) ( n - i ) * Math.PI 
                                    / ( double ) ( n - 1     ) );
                }
                x[0] = -1.0;
                if ( ( n % 2 ) == 1 )
                {
                    x[(n-1)/2] = 0.0;
                }
                x[n-1] = +1.0;
            }
            return x;
        }
        
        public static double[] ccn_compute_points_new(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CCN_COMPUTE_POINTS: compute Clenshaw Curtis Nested points.
            //
            //  Discussion:
            //
            //    We want to compute the following sequence:
            //
            //    1/2,
            //    0, 1
            //    1/4, 3/4
            //    1/8, 3/8, 5/8, 7/8,
            //    1/16, 3/16, 5/16, 7/16, 9/16, 11/16, 13/16, 15/16, and so on.
            //
            //    But we would prefer that the numbers in each row be regrouped in pairs
            //    that are symmetric about 1/2, with the number above 1/2 coming first.
            //    Thus, the last row might become:
            //    (9/16, 7/16), (11/16, 5/16), ..., (15/16, 1/16).
            //
            //    Once we have our sequence, we apply the Chebyshev transformation
            //    which maps [0,1] to [-1,+1].
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 March 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of elements to compute.
            //
            //    Output, double CCN_COMPUTE_POINTS_NEW[N], the elements of the sequence.
            //
        {
            int d;
            int i;
            int k;
            int m;
            int td;
            int tu;

            double[] x = new double[n];
            //
            //  Handle first three entries specially.
            //
            if (1 <= n)
            {
                x[0] = 0.5;
            }

            if (2 <= n)
            {
                x[1] = 1.0;
            }

            if (3 <= n)
            {
                x[2] = 0.0;
            }

            m = 3;
            d = 2;

            while (m < n)
            {
                tu = d + 1;
                td = d - 1;

                k = Math.Min(d, n - m);

                for (i = 1; i <= k; i++)
                {
                    if ((i % 2) == 1)
                    {
                        x[m + i - 1] = tu / 2.0 / (double) (k);
                        tu = tu + 2;
                    }
                    else
                    {
                        x[m + i - 1] = td / 2.0 / (double) (k);
                        td = td - 2;
                    }
                }

                m = m + k;
                d = d * 2;
            }

            //
            //  Apply the Chebyshev transformation.
            //
            for (i = 0; i < n; i++)
            {
                x[i] = Math.Cos(x[i] * Math.PI);
            }

            x[0] = 0.0;

            if (2 <= n)
            {
                x[1] = -1.0;
            }

            if (3 <= n)
            {
                x[2] = +1.0;
            }

            return x;
        }

        public static double[] nc_compute_new(int n, double x_min, double x_max, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NC_COMPUTE_NEW computes a Newton-Cotes quadrature rule.
            //
            //  Discussion:
            //
            //    For the interval [X_MIN,X_MAX], the Newton-Cotes quadrature rule
            //    estimates
            //
            //      Integral ( X_MIN <= X <= X_MAX ) F(X) dX
            //
            //    using N abscissas X and weights W:
            //
            //      Sum ( 1 <= I <= N ) W(I) * F ( X(I) ).
            //
            //    For the CLOSED rule, the abscissas include the end points.
            //    For the OPEN rule, the abscissas do not include the end points.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 November 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order.
            //
            //    Input, double X_MIN, X_MAX, the endpoints of the interval.
            //
            //    Input, double X[N], the abscissas.
            //
            //    Output, double NC_COMPUTE_NEW[N], the weights.
            //
        {
            int i;
            int j;
            int k;
            double yvala;
            double yvalb;

            double[] d = new double[n];
            double[] w = new double[n];

            for (i = 0; i < n; i++)
            {
                //
                //  Compute the Lagrange basis polynomial which is 1 at XTAB(I),
                //  and zero at the other nodes.
                //
                for (j = 0; j < n; j++)
                {
                    d[j] = 0.0;
                }

                d[i] = 1.0;

                for (j = 2; j <= n; j++)
                {
                    for (k = j; k <= n; k++)
                    {
                        d[n + j - k - 1] = (d[n + j - k - 1 - 1] - d[n + j - k - 1]) /
                                           (x[n + 1 - k - 1] - x[n + j - k - 1]);
                    }
                }

                for (j = 1; j <= n - 1; j++)
                {
                    for (k = 1; k <= n - j; k++)
                    {
                        d[n - k - 1] = d[n - k - 1] - x[n - k - j] * d[n - k];
                    }
                }

                //
                //  Evaluate the antiderivative of the polynomial at the left and
                //  right endpoints.
                //
                yvala = d[n - 1] / (double) (n);
                for (j = n - 2; 0 <= j; j--)
                {
                    yvala = yvala * x_min + d[j] / (double) (j + 1);
                }

                yvala = yvala * x_min;

                yvalb = d[n - 1] / (double) (n);
                for (j = n - 2; 0 <= j; j--)
                {
                    yvalb = yvalb * x_max + d[j] / (double) (j + 1);
                }

                yvalb = yvalb * x_max;

                w[i] = yvalb - yvala;
            }

            return w;
        }

        public static void rescale(double a, double b, int n, ref double[] x, ref double[] w)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    RESCALE rescales a Legendre quadrature rule from [-1,+1] to [A,B].
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 October 2009
            //
            //  Author:
            //
            //    John Burkardt.
            //
            //  Reference:
            //
            //    Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin,
            //    A fast algorithm for the calculation of the roots of special functions,
            //    SIAM Journal on Scientific Computing,
            //    Volume 29, Number 4, pages 1420-1438, 2007.
            //
            //  Parameters:
            //
            //    Input, double A, B, the endpoints of the new interval.
            //
            //    Input, int N, the order.
            //
            //    Input/output, double X[N], on input, the abscissas for [-1,+1].
            //    On output, the abscissas for [A,B].
            //
            //    Input/output, double W[N], on input, the weights for [-1,+1].
            //    On output, the weights for [A,B].
            //
        {
            int i;

            for (i = 0; i < n; i++)
            {
                x[i] = ((a + b) + (b - a) * x[i]) / 2.0;
            }

            for (i = 0; i < n; i++)
            {
                w[i] = (b - a) * w[i] / 2.0;
            }

            return;
        }

        public static void rule_write(int order, string filename, double[] x, double[] w,
                double[] r)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    RULE_WRITE writes a quadrature rule to three files.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 February 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int ORDER, the order of the rule.
            //
            //    Input, double A, the left endpoint.
            //
            //    Input, double B, the right endpoint.
            //
            //    Input, string FILENAME, specifies the output filenames.
            //    "filename_w.txt", "filename_x.txt", "filename_r.txt"
            //    defining weights, abscissas, and region.
            //
        {
            string filename_r;
            string filename_w;
            string filename_x;

            filename_w = filename + "_w.txt";
            filename_x = filename + "_x.txt";
            filename_r = filename + "_r.txt";

            Console.WriteLine("");
            Console.WriteLine("  Creating quadrature files.");
            Console.WriteLine("");
            Console.WriteLine("  Root file name is     \"" + filename + "\".");
            Console.WriteLine("");
            Console.WriteLine("  Weight file will be   \"" + filename_w + "\".");
            Console.WriteLine("  Abscissa file will be \"" + filename_x + "\".");
            Console.WriteLine("  Region file will be   \"" + filename_r + "\".");

            typeMethods.r8mat_write(filename_w, 1, order, w);
            typeMethods.r8mat_write(filename_x, 1, order, x);
            typeMethods.r8mat_write(filename_r, 1, 2, r);

        }
    }
}