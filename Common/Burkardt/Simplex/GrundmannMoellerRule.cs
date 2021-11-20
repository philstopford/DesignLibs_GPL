using System;
using Burkardt.Composition;
using Burkardt.Types;

namespace Burkardt.SimplexNS;

public static class GrundmannMoellerRule
{
    public static void gm_general_rule_set(int rule, int m, int n, double[] t, ref double[] w,
            ref double[] x)

        //*****************************************************************************/
        //
        //  Purpose:
        //
        //    GM_GENERAL_RULE_SET sets a Grundmann-Moeller rule for a general simplex.
        //
        //  Discussion:
        //
        //    The vertices of the simplex are given by the array T.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 March 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Axel Grundmann, Michael Moeller,
        //    Invariant Integration Formulas for the N-Simplex 
        //    by Combinatorial Methods,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 15, Number 2, April 1978, pages 282-290.
        //
        //  Parameters:
        //
        //    Input, int RULE, the index of the rule.
        //    0 <= RULE.
        //
        //    Input, int M, the spatial dimension.
        //    1 <= M.
        //
        //    Input, int N, the number of points in the rule.
        //
        //    Input, double T[M*(M+1)], the vertices of the simplex.
        //
        //    Output, double W[N], the weights.
        //
        //    Output, double X[M*N], the abscissas.
        //
    {
        int j;
        //
        //  Get the unit rule.
        //
        double[] w1 = new double[n];
        double[] x1 = new double[m * n];

        gm_unit_rule_set(rule, m, n, ref w1, ref x1);
        //
        //  Compute the volume of the unit simplex.
        //
        double volume1 = Simplex.simplex_unit_volume(m);
        //
        //  Compute the volume of the general simplex.
        //
        double volume = Simplex.simplex_general_volume(m, t);
        //
        //  Convert the points.
        //
        Simplex.simplex_unit_to_general(m, n, t, x1, ref x);
        //
        //  Convert the weights.
        //
        for (j = 0; j < n; j++)
        {
            w[j] = w1[j] * volume / volume1;
        }
    }

    public static int gm_rule_size(int rule, int m)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GM_RULE_SIZE determines the size of a Grundmann-Moeller rule.
        //
        //  Discussion:
        //
        //    This rule returns the value of N, the number of points associated
        //    with a GM rule of given index.
        //
        //    After calling this rule, the user can use the value of N to
        //    allocate space for the weight vector as W(N) and the abscissa
        //    vector as X(M,N), and then call GM_RULE_SET.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Axel Grundmann, Michael Moeller,
        //    Invariant Integration Formulas for the N-Simplex
        //    by Combinatorial Methods,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 15, Number 2, April 1978, pages 282-290.
        //
        //  Parameters:
        //
        //    Input, int RULE, the index of the rule.
        //    0 <= RULE.
        //
        //    Input, int M, the spatial dimension.
        //    1 <= M.
        //
        //    Output, int GM_RULE_SIZE, the number of points in the rule.
        //
    {
        int arg1 = m + rule + 1;

        int n = typeMethods.i4_choose(arg1, rule);

        return n;
    }

    public static void gm_unit_rule_set(int rule, int m, int n, ref double[] w, ref double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GM_UNIT_RULE_SET sets a Grundmann-Moeller rule.
        //
        //  Discussion:
        //
        //    This is a revised version of the calculation which seeks to compute
        //    the value of the weight in a cautious way that apublic static voids intermediate
        //    overflow.  Thanks to John Peterson for pointing out the problem on
        //    26 June 2008.
        //
        //    This rule returns weights and abscissas of a Grundmann-Moeller
        //    quadrature rule for the M-dimensional unit simplex.
        //
        //    The dimension N can be determined by calling GM_RULE_SIZE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 June 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Axel Grundmann, Michael Moeller,
        //    Invariant Integration Formulas for the N-Simplex
        //    by Combinatorial Methods,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 15, Number 2, April 1978, pages 282-290.
        //
        //  Parameters:
        //
        //    Input, int RULE, the index of the rule.
        //    0 <= RULE.
        //
        //    Input, int M, the spatial dimension.
        //    1 <= M.
        //
        //    Input, int N, the number of points in the rule.
        //
        //    Output, double W[N], the weights.
        //
        //    Output, double X[M*N], the abscissas.
        //
    {
        int i;

        int d = 2 * rule + 1;
        int k = 0;
        double one_pm = 1.0;

        int[] beta = new int[m + 1];

        for (i = 0; i <= rule; i++)
        {
            double weight = one_pm;

            int j_hi = Math.Max(m, Math.Max(d, d + m - i));

            int j;
            for (j = 1; j <= j_hi; j++)
            {
                if (j <= m)
                {
                    weight *= j;
                }

                if (j <= d)
                {
                    weight *= d + m - 2 * i;
                }

                if (j <= 2 * rule)
                {
                    weight /= 2.0;
                }

                if (j <= i)
                {
                    weight /= j;
                }

                if (j <= d + m - i)
                {
                    weight /= j;
                }
            }

            one_pm = -one_pm;

            int beta_sum = rule - i;
            bool more = false;
            int h = 0;
            int t = 0;

            for (;;)
            {
                Comp.comp_next(beta_sum, m + 1, ref beta, ref more, ref h, ref t);

                w[k] = weight;
                int dim;
                for (dim = 0; dim < m; dim++)
                {
                    x[dim + k * m] = (2 * beta[dim + 1] + 1)
                                     / (double)(d + m - 2 * i);
                }

                k += 1;

                if (!more)
                {
                    break;
                }
            }
        }

        //
        //  Normalize.
        //
        double volume1 = Simplex.simplex_unit_volume(m);
        for (i = 0; i < n; i++)
        {
            w[i] *= volume1;
        }
    }
}