using System;
using Burkardt.Quadrature;
using Burkardt.Table;
using Burkardt.Types;

namespace Burkardt.ClenshawCurtisNS
{
    public static class ClenshawCurtis
    {
        public static int cce_order(int l)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CCE_ORDER: order of a Clenshaw-Curtis Exponential rule from the level.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    22 February 2014
            //
            //  Author:
            //
            //    John Burkardt.
            //
            //  Parameters:
            //
            //    Input, int L, the level of the rule.  
            //    1 <= L.
            //
            //    Output, int CCE_ORDER, the order of the rule.
            //
        {
            int n;

            if (l < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("CCE_ORDER - Fatal error!");
                Console.WriteLine("  1 <= L required.");
                Console.WriteLine("  Input L = " + l + "");
                return (1);
            }
            else if (l == 1)
            {
                n = 1;
            }
            else
            {
                n = (int)Math.Pow(2, l - 1) + 1;
            }

            return n;
        }

        public static int ccl_order(int l)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CCL_ORDER computes the order of a CCL rule from the level.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    19 February 2014
            //
            //  Author:
            //
            //    John Burkardt.
            //
            //  Parameters:
            //
            //    Input, int L, the level of the rule.  
            //    1 <= L.
            //
            //    Output, int CCL_ORDER, the order of the rule.
            //
        {
            int n;

            if (l < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("CCL_ORDER - Fatal error!");
                Console.WriteLine("  1 <= L required.");
                Console.WriteLine("  Input L = " + l + "");
                return (1);
            }

            n = 2 * l - 1;

            return n;
        }

        public static int ccs_order(int l)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CCS_ORDER: order of a "slow growth" Clenshaw Curtis quadrature rule.
            //
            //  Discussion:
            //
            //    Our convention is that the abscissas are numbered from left to right.
            //
            //    The rule is defined on [0,1].
            //
            //    The integral to approximate:
            //
            //      Integral ( 0 <= X <= 1 ) F(X) dX
            //
            //    The quadrature rule:
            //
            //      Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
            //
            //    The input value L requests a rule of precision at least 2*L-1.
            //
            //    In order to preserve nestedness, this function returns the order
            //    of a rule which is the smallest value of the form 1+2^E which
            //    is greater than or equal to 2*L-1.
            //
            //     L  2*L-1   N
            //    --  -----  --
            //     1      1   1
            //     2      3   3
            //     3      5   5
            //     4      7   9
            //     5      9   9
            //     6     11  17
            //     7     13  17
            //     8     15  17
            //     9     17  17
            //    10     19  33
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    08 December 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int L, the level of the rule.
            //    1 <= L.
            //
            //    Output, int CCS_ORDER, the appropriate order.
            //
        {
            int n;

            if (l < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("CCS_ORDER - Fatal error!");
                Console.WriteLine("  Illegal value of L = " + l + "");
                return (1);
            }

            //
            //  Find the order N that satisfies the precision requirement.
            //
            if (l == 1)
            {
                n = 1;
            }
            else
            {
                n = 3;
                while (n < 2 * l - 1)
                {
                    n = 2 * n - 1;
                }
            }

            return n;
        }

        public class ccResult
        {
            public double[] x;
            public double[] w;
        }

        public static ccResult cc(int n, double[] x_, double[] w_)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CC computes a Clenshaw Curtis quadrature rule based on order.
            //
            //  Discussion:
            //
            //    Our convention is that the abscissas are numbered from left to right.
            //
            //    The rule is defined on [0,1].
            //
            //    The integral to approximate:
            //
            //      Integral ( 0 <= X <= 1 ) F(X) dX
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
            //    08 December 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the rule.
            //    1 <= N.
            //
            //    Output, double X[N], the abscissas.
            //
            //    Output, double W[N], the weights.
            //
        {
            double b;
            int i;
            int j;
            double pi = 3.141592653589793;
            double theta;
            ccResult result = new ccResult();
            result.w = w_;
            result.x = x_;

            if (n < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("CC - Fatal error!");
                Console.WriteLine("  Illegal value of N = " + n + "");
                return result;
            }

            if (n == 1)
            {
                result.x[0] = 0.0;
                result.w[0] = 2.0;
            }
            else
            {
                for (i = 0; i < n; i++)
                {
                    result.x[i] = Math.Cos((double)(n - 1 - i) * pi
                                           / (double)(n - 1));
                }

                result.x[0] = -1.0;
                if ((n % 2) == 1)
                {
                    result.x[(n + 1) / 2 - 1] = 0.0;
                }

                result.x[n - 1] = +1.0;

                for (i = 0; i < n; i++)
                {
                    theta = (double)(i) * pi
                            / (double)(n - 1);

                    result.w[i] = 1.0;

                    for (j = 1; j <= (n - 1) / 2; j++)
                    {
                        if (2 * j == (n - 1))
                        {
                            b = 1.0;
                        }
                        else
                        {
                            b = 2.0;
                        }

                        result.w[i] = result.w[i] - b * Math.Cos(2.0 * (double)(j) * theta)
                            / (double)(4 * j * j - 1);
                    }
                }

                result.w[0] = result.w[0] / (double)(n - 1);
                for (j = 1; j < n - 1; j++)
                {
                    result.w[j] = 2.0 * result.w[j] / (double)(n - 1);
                }

                result.w[n - 1] = result.w[n - 1] / (double)(n - 1);
            }

            //
            //  Transform from [-1,+1] to [0,1].
            //
            QuadratureRule.rule_adjust(-1.0, +1.0, 0.0, 1.0, n, ref result.x, ref result.w);

            return result;
        }

        public static double cc_abscissa(int order, int i)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CC_ABSCISSA returns the I-th abscissa of the Clenshaw Curtis rule.
            //
            //  Discussion:
            //
            //    Our convention is that the abscissas are numbered from left to
            //    right.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    16 October 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int ORDER, the order of the rule.
            //
            //    Input, int I, the index of the desired abscissa.  1 <= I <= ORDER.
            //
            //    Output, double CC_ABSCISSA, the value of the I-th 
            //    abscissa in the rule of order ORDER.
            //
        {
            double pi = 3.141592653589793;
            double value;

            if (order < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("CC_ABSCISSA - Fatal error!");
                Console.WriteLine("  Input value of ORDER < 1.");
                Console.WriteLine("  Input value of ORDER = " + order + "");
                return (1);
            }

            if (i < 1 || order < i)
            {
                Console.WriteLine("");
                Console.WriteLine("CC_ABSCISSA - Fatal error!");
                Console.WriteLine("  1 <= I <= ORDER is required.");
                Console.WriteLine("  I = " + i + "");
                Console.WriteLine("  ORDER = " + order + "");
                return (1);
            }

            if (order == 1)
            {
                value = 0.0;
                return value;
            }

            value = Math.Cos((double)(order - i) * pi
                             / (double)(order - 1));

            if (2 * i - 1 == order)
            {
                value = 0.0;
            }

            return value;
        }

        public static double[] product_weights_cc(int dim_num, int[] order_1d, int order_nd)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PRODUCT_WEIGHTS_CC computes weights for a Clenshaw Curtis product rule.
            //
            //  Discussion:
            //
            //    This routine computes the weights for a quadrature rule which is
            //    a product of closed rules of varying order.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    09 November 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the spatial dimension.
            //
            //    Input, int ORDER_1D[DIM_NUM], the order of the 1D rules.
            //
            //    Input, int ORDER_ND, the order of the product rule.
            //
            //    Output, double PRODUCT_WEIGHTS_CC[DIM_NUM*ORDER_ND], 
            //    the product rule weights.
            //
        {
            int dim;
            int order;
            double[] w_1d;
            double[] w_nd;

            typeMethods.r8vecDPData data = new typeMethods.r8vecDPData();

            w_nd = new double[order_nd];

            for (order = 0; order < order_nd; order++)
            {
                w_nd[order] = 1.0;
            }

            for (dim = 0; dim < dim_num; dim++)
            {
                w_1d = cc_weights(order_1d[dim]);

                typeMethods.r8vec_direct_product2(ref data, dim, order_1d[dim], w_1d, dim_num,
                    order_nd, ref w_nd);
            }

            return w_nd;
        }

        public static double[] cc_weights(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CC_WEIGHTS computes Clenshaw Curtis weights.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    09 November 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Charles Clenshaw, Alan Curtis,
            //    A Method for Numerical Integration on an Automatic Computer,
            //    Numerische Mathematik,
            //    Volume 2, Number 1, December 1960, pages 197-205.
            //
            //  Parameters:
            //
            //    Input, int N, the order of the rule.
            //
            //    Output, double CC_WEIGHTS[N], the weights of the rule.
            //
        {
            double b;
            int i;
            int j;
            double pi = 3.141592653589793;
            double theta;
            double[] w;

            w = new double[n];

            if (n == 1)
            {
                w[0] = 2.0;
                return w;
            }

            for (i = 1; i <= n; i++)
            {
                theta = (double)(i - 1) * pi / (double)(n - 1);

                w[i - 1] = 1.0;

                for (j = 1; j <= (n - 1) / 2; j++)
                {
                    if (2 * j == (n - 1))
                    {
                        b = 1.0;
                    }
                    else
                    {
                        b = 2.0;
                    }

                    w[i - 1] = w[i - 1] - b * Math.Cos(2.0 * (double)(j) * theta)
                        / (double)(4 * j * j - 1);
                }
            }

            w[0] = w[0] / (double)(n - 1);
            for (i = 1; i < n - 1; i++)
            {
                w[i] = 2.0 * w[i] / (double)(n - 1);
            }

            w[n - 1] = w[n - 1] / (double)(n - 1);

            return w;
        }

        public static int index_to_level_closed(int dim_num, int[] t, int order, int level_max, int tIndex = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    INDEX_TO_LEVEL_CLOSED determines the level of a point given its index.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    09 November 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Fabio Nobile, Raul Tempone, Clayton Webster,
            //    A Sparse Grid Stochastic Collocation Method for Partial Differential
            //    Equations with Random Input Data,
            //    SIAM Journal on Numerical Analysis,
            //    Volume 46, Number 5, 2008, pages 2309-2345.
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the spatial dimension.
            //
            //    Input, int T[DIM_NUM], the grid indices of a point in a 1D closed rule.
            //    0 <= T[I] <= ORDER.
            //
            //    Input, int ORDER, the order of the rule.
            //
            //    Input, int LEVEL_MAX, the level with respect to which the
            //    index applies.
            //
            //    Output, int INDEX_TO_LEVEL_CLOSED, the first level on which
            //    the point associated with the given index will appear.
            //
        {
            int dim;
            int level;
            int s;
            int value;

            value = 0;

            for (dim = 0; dim < dim_num; dim++)
            {
                s = t[tIndex + dim];

                s = typeMethods.i4_modp(s, order);

                if (s == 0)
                {
                    level = 0;
                }
                else
                {
                    level = level_max;

                    while ((s % 2) == 0)
                    {
                        s = s / 2;
                        level = level - 1;
                    }
                }

                if (level == 0)
                {
                    level = 1;
                }
                else if (level == 1)
                {
                    level = 0;
                }

                value = value + level;
            }

            return value;
        }

        public static void level_to_order_closed(int dim_num, int[] level, ref int[] order)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LEVEL_TO_ORDER_CLOSED converts a level to an order for closed rules.
            //
            //  Discussion:
            //
            //    Sparse grids can naturally be nested.  A natural scheme is to use
            //    a series of one-dimensional rules arranged in a series of "levels"
            //    whose order roughly doubles with each step.
            //
            //    The arrangement described here works naturally for the Clenshaw Curtis
            //    and Newton Cotes closed rules.  
            //
            //    The idea is that we start with LEVEL = 0, ORDER = 1 indicating the single 
            //    point at the center, and for all values afterwards, we use the 
            //    relationship
            //
            //      ORDER = 2^LEVEL + 1
            //
            //    The following table shows how the growth will occur:
            //
            //    Level    Order
            //
            //    0          1
            //    1          3 =  2 + 1
            //    2          5 =  4 + 1
            //    3          9 =  8 + 1
            //    4         17 = 16 + 1
            //    5         33 = 32 + 1
            //
            //    For the Clenshaw Curtis and Newton Cotes Closed rules, the point growth
            //    is nested.  If we have ORDER points on a particular LEVEL, the next
            //    level includes all these old points, plus ORDER-1 new points, formed
            //    in the gaps between successive pairs of old points.
            //
            //    Level    Order = New + Old
            //
            //    0          1   =  1  +  0
            //    1          3   =  2  +  1
            //    2          5   =  2  +  3
            //    3          9   =  4  +  5
            //    4         17   =  8  +  9
            //    5         33   = 16  + 17
            //
            //    In this routine, we assume that a vector of levels is given,
            //    and the corresponding orders are desired.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    09 November 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Fabio Nobile, Raul Tempone, Clayton Webster,
            //    A Sparse Grid Stochastic Collocation Method for Partial Differential
            //    Equations with Random Input Data,
            //    SIAM Journal on Numerical Analysis,
            //    Volume 46, Number 5, 2008, pages 2309-2345.
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the spatial dimension.
            //
            //    Input, int LEVEL[DIM_NUM], the nesting level.
            //
            //    Output, int ORDER[DIM_NUM], the order (number of points) 
            //    of the rule.
            //
        {
            int dim;

            for (dim = 0; dim < dim_num; dim++)
            {
                if (level[dim] < 0)
                {
                    order[dim] = -1;
                }
                else if (level[dim] == 0)
                {
                    order[dim] = 1;
                }
                else
                {
                    order[dim] = (int)Math.Pow(2, level[dim]) + 1;
                }
            }
        }

        public static void clenshaw_curtis_compute(int order, ref double[] x, ref double[] w)

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
                    x[i] = Math.Cos((double)(order - 1 - i) * Math.PI
                                    / (double)(order - 1));
                }

                x[0] = -1.0;
                if ((order % 2) == 1)
                {
                    x[(order - 1) / 2] = 0.0;
                }

                x[order - 1] = +1.0;

                for (i = 0; i < order; i++)
                {
                    theta = (double)(i) * Math.PI / (double)(order - 1);

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

                        w[i] = w[i] - b * Math.Cos(2.0 * (double)(j) * theta)
                            / (double)(4 * j * j - 1);
                    }
                }

                w[0] = w[0] / (double)(order - 1);
                for (i = 1; i < order - 1; i++)
                {
                    w[i] = 2.0 * w[i] / (double)(order - 1);
                }

                w[order - 1] = w[order - 1] / (double)(order - 1);
            }

            return;
        }

        public static double[] cc_compute_points(int n)

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

            if (n < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("CC_COMPUTE_POINTS - Fatal error!");
                Console.WriteLine("  N < 1.");
                return null;
            }

            x = new double[n];

            if (n == 1)
            {
                x[0] = 0.0;
            }
            else
            {
                for (i = 1; i <= n; i++)
                {
                    x[i - 1] = Math.Cos((double)(n - i) * Math.PI
                                        / (double)(n - 1));
                }

                x[0] = -1.0;
                if ((n % 2) == 1)
                {
                    x[(n - 1) / 2] = 0.0;
                }

                x[n - 1] = +1.0;
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
                        x[m + i - 1] = tu / 2.0 / (double)(k);
                        tu = tu + 2;
                    }
                    else
                    {
                        x[m + i - 1] = td / 2.0 / (double)(k);
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


        public static void level_to_order_open(int dim_num, int[] level, ref int[] order)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LEVEL_TO_ORDER_OPEN converts a level to an order for open rules.
            //
            //  Discussion:
            //
            //    Sparse grids can naturally be nested.  A natural scheme is to use
            //    a series of one-dimensional rules arranged in a series of "levels"
            //    whose order roughly doubles with each step.
            //
            //    The arrangement described here works naturally for the Fejer Type 1,
            //    Fejer Type 2, Newton Cotes Open, Newton Cotes Half Open,
            //    and Gauss-Patterson rules.  It also can be used, partially, to describe
            //    the growth of Gauss-Legendre rules.
            //
            //    The idea is that we start with LEVEL = 0, ORDER = 1 indicating the single 
            //    point at the center, and for all values afterwards, we use the relationship
            //
            //      ORDER = 2**(LEVEL+1) - 1.
            //
            //    The following table shows how the growth will occur:
            //
            //    Level    Order
            //
            //    0          1
            //    1          3 =  4 - 1
            //    2          7 =  8 - 1
            //    3         15 = 16 - 1
            //    4         31 = 32 - 1
            //    5         63 = 64 - 1
            //
            //    For the Fejer Type 1, Fejer Type 2, Newton Cotes Open, 
            //    Newton Cotes Open Half, and Gauss-Patterson rules, the point growth is
            //    nested.  If we have ORDER points on a particular LEVEL, the next level 
            //    includes all these old points, plus ORDER+1 new points, formed in the 
            //    gaps between successive pairs of old points plus an extra point at each 
            //    end.
            //
            //    Level    Order = New + Old
            //
            //    0          1   =  1  +  0
            //    1          3   =  2  +  1
            //    2          7   =  4  +  3
            //    3         15   =  8  +  7
            //    4         31   = 16  + 15
            //    5         63   = 32  + 31
            //
            //    If we use a series of Gauss-Legendre rules, then there is almost no 
            //    nesting, except that the central point is shared.  If we insist on 
            //    producing a comparable series of such points, then the "nesting" behavior
            //    is as follows:
            //
            //    Level    Order = New + Old
            //
            //    0          1   =  1  +  0
            //    1          3   =  2  +  1
            //    2          7   =  6  +  1
            //    3         15   = 14  +  1
            //    4         31   = 30  +  1
            //    5         63   = 62  +  1
            //
            //    Moreover, if we consider ALL the points used in such a set of "nested" 
            //    Gauss-Legendre rules, then we must sum the "NEW" column, and we see that
            //    we get roughly twice as many points as for the truly nested rules.
            //
            //    In this routine, we assume that a vector of levels is given,
            //    and the corresponding orders are desired.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 April 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Fabio Nobile, Raul Tempone, Clayton Webster,
            //    A Sparse Grid Stochastic Collocation Method for Partial Differential
            //    Equations with Random Input Data,
            //    SIAM Journal on Numerical Analysis,
            //    Volume 46, Number 5, 2008, pages 2309-2345.
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the spatial dimension.
            //
            //    Input, int LEVEL[DIM_NUM], the nesting level.
            //
            //    Output, int ORDER[DIM_NUM], the order (number of points) 
            //    of the rule.
            //
        {
            int dim;

            for (dim = 0; dim < dim_num; dim++)
            {
                if (level[dim] < 0)
                {
                    order[dim] = -1;
                }
                else if (level[dim] == 0)
                {
                    order[dim] = 1;
                }
                else
                {
                    order[dim] = (int)Math.Pow(2, level[dim] + 1) - 1;
                }
            }
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
                yvala = d[n - 1] / (double)(n);
                for (j = n - 2; 0 <= j; j--)
                {
                    yvala = yvala * x_min + d[j] / (double)(j + 1);
                }

                yvala = yvala * x_min;

                yvalb = d[n - 1] / (double)(n);
                for (j = n - 2; 0 <= j; j--)
                {
                    yvalb = yvalb * x_max + d[j] / (double)(j + 1);
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