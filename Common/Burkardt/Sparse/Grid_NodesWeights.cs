using System;
using Burkardt.ClenshawCurtisNS;
using Burkardt.Composition;
using Burkardt.Quadrature;
using Burkardt.Types;

namespace Burkardt.Sparse;

public static class Grid_NodesWeights
{
    public static void nwspgr(Func<int, double[], double[], ClenshawCurtis.ccResult> rule,
            Func<int, int> rule_order, int dim, int k, int r_size, ref int s_size,
            ref double[] nodes, ref double[] weights)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NWSPGR generates nodes and weights for sparse grid integration.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 April 2013
        //
        //  Author:
        //
        //    Original MATLAB version by Florian Heiss, Viktor Winschel.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Florian Heiss, Viktor Winschel,
        //    Likelihood approximation by numerical integration on sparse grids,
        //    Journal of Econometrics,
        //    Volume 144, 2008, pages 62-80.
        //
        //  Parameters:
        //
        //    Input, void RULE ( int n, double x[], double w[] ), the name of a function
        //    which is given the order N and returns the points X and weights W of the
        //    corresponding 1D quadrature rule.
        //
        //    Input, int RULE_ORDER ( int l ), the name of a function which
        //    is given the level L and returns the order N of the corresponding 1D rule.
        //
        //    Input, int DIM, the spatial dimension.
        //
        //    Input, int K, the level of the sparse rule.
        //
        //    Input, int R_SIZE, the "size" of the sparse rule.
        //
        //    Output, int &S_SIZE, the size of the sparse rule, after
        //    duplicate points have been merged.
        //
        //    Output, double NODES[DIM*R_SIZE], the nodes of the sparse rule.
        //
        //    Output, double WEIGHTS[R_SIZE], the weights of the sparse rule.
        //
    {
        int i;
        int j;
        int level;
        int n;
        int q;

        for (j = 0; j < r_size; j++)
        {
            for (i = 0; i < dim; i++)
            {
                nodes[i + j * dim] = 0.0;
            }
        }

        for (j = 0; j < r_size; j++)
        {
            weights[j] = 0.0;
        }

        //
        //  Create cell arrays that will contain the points and weights 
        //  for levels 1 through K.
        //
        int[] n1d = new int[k];
        int[] x1d_off = new int[k + 1];
        int[] w1d_off = new int[k + 1];

        x1d_off[0] = 0;
        w1d_off[0] = 0;

        for (level = 1; level <= k; level++)
        {
            n = rule_order(level);
            n1d[level - 1] = n;
            x1d_off[level] = x1d_off[level - 1] + n;
            w1d_off[level] = w1d_off[level - 1] + n;
        }

        int n1d_total = x1d_off[k];
        //
        //  Calculate all the 1D rules needed.
        //
        double[] x1d = new double[n1d_total];
        double[] w1d = new double[n1d_total];

        for (level = 1; level <= k; level++)
        {
            n = n1d[level - 1];

            double[] x = new double[n];
            double[] w = new double[n];

            ClenshawCurtis.ccResult result = rule(n, x, w);
            x = result.x;
            w = result.w;
                
            typeMethods.r8cvv_rset(n1d_total, x1d, k, x1d_off, level - 1, x);
            typeMethods.r8cvv_rset(n1d_total, w1d, k, w1d_off, level - 1, w);

        }

        //
        //  Construct the sparse grid.
        //
        int minq = Math.Max(0, k - dim);
        int maxq = k - 1;
        //
        //  Q is the level total.
        //
        int[] lr = new int[dim];
        int[] nr = new int[dim];

        int r = 0;

        for (q = minq; q <= maxq; q++)
        {
            //
            //  BQ is the combinatorial coefficient applied to the component
            //  product rules which have level Q.
            //
            int bq = typeMethods.i4_mop(maxq - q) * typeMethods.i4_choose(dim - 1, dim + q - k);
            //
            //  Compute the D-dimensional row vectors that sum to DIM+Q.
            //
            int seq_num = Comp.num_seq(q, dim);

            int[] is_ = Comp.get_seq(dim, q + dim, seq_num);
            //
            //  Allocate new rows for nodes and weights.
            //
            int[] rq = new int[seq_num];

            for (j = 0; j < seq_num; j++)
            {
                rq[j] = 1;
                for (i = 0; i < dim; i++)
                {
                    level = is_[j + i * seq_num] - 1;
                    rq[j] *= n1d[level];
                }
            }

            //
            //  Generate every product rule whose level total is Q.
            //
            int j2;
            for (j2 = 0; j2 < seq_num; j2++)
            {
                for (i = 0; i < dim; i++)
                {
                    lr[i] = is_[j2 + i * seq_num];
                }

                for (i = 0; i < dim; i++)
                {
                    nr[i] = rule_order(lr[i]);
                }

                int[] roff = typeMethods.r8cvv_offset(dim, nr);

                int nc = typeMethods.i4vec_sum(dim, nr);
                double[] wc = new double[nc];
                double[] xc = new double[nc];
                //
                //  Corrected first argument in calls to R8CVV to N1D_TOTAL,
                //  19 April 2013.
                //
                for (i = 0; i < dim; i++)
                {
                    double[] xr = typeMethods.r8cvv_rget_new(n1d_total, x1d, k, x1d_off, lr[i] - 1);
                    double[] wr = typeMethods.r8cvv_rget_new(n1d_total, w1d, k, w1d_off, lr[i] - 1);
                    typeMethods.r8cvv_rset(nc, xc, dim, roff, i, xr);
                    typeMethods.r8cvv_rset(nc, wc, dim, roff, i, wr);
                }

                int np = rq[j2];
                double[] wp = new double[np];
                double[] xp = new double[dim * np];

                typeMethods.tensor_product_cell(nc, xc, wc, dim, nr, roff, np, ref xp, ref wp);
                //
                //  Append the new nodes and weights to the arrays.
                //
                for (j = 0; j < np; j++)
                {
                    for (i = 0; i < dim; i++)
                    {
                        nodes[i + (r + j) * dim] = xp[i + j * dim];
                    }
                }

                for (j = 0; j < np; j++)
                {
                    weights[r + j] = bq * wp[j];
                }

                //
                //  Update the count.
                //
                r += rq[j2];

            }
        }

        //
        //  Reorder the rule so the points are in ascending lexicographic order.
        //
        QuadratureRule.rule_sort(dim, r_size, ref nodes, ref weights);
        //
        //  Suppress duplicate points and merge weights.
        //
        r = 0;
        for (j = 1; j < r_size; j++)
        {
            bool equal = true;
            for (i = 0; i < dim; i++)
            {
                if (!(Math.Abs(nodes[i + r * dim] - nodes[i + j * dim]) > typeMethods.r8_epsilon()))
                {
                    continue;
                }

                equal = false;
                break;
            }

            switch (equal)
            {
                case true:
                    weights[r] += weights[j];
                    break;
                default:
                {
                    r += 1;
                    weights[r] = weights[j];
                    for (i = 0; i < dim; i++)
                    {
                        nodes[i + r * dim] = nodes[i + j * dim];
                    }

                    break;
                }
            }
        }

        r += 1;
        s_size = r;
        //
        //  Zero out unneeded entries.
        //
        for (j = r; j < r_size; j++)
        {
            for (i = 0; i < dim; i++)
            {
                nodes[i + j * dim] = 0.0;
            }
        }

        for (j = r; j < r_size; j++)
        {
            weights[j] = 0.0;
        }

        //
        //  Normalize the weights to sum to 1.
        //
        double t = typeMethods.r8vec_sum(r, weights);

        for (j = 0; j < r; j++)
        {
            weights[j] /= t;
        }
    }

    public static int nwspgr_size(Func<int, int> rule_order, int dim, int k)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NWSPGR_SIZE determines the size of a sparse grid rule.
        //
        //  Discussion:
        //
        //    This routine does a "raw" count, that is, it does not notice that many
        //    points may be repeated, in which case, the size of the rule could be
        //    reduced by merging repeated points and combining the corresponding weights.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 December 2012
        //
        //  Author:
        //
        //    John Burkardt.
        //
        //  Reference:
        //
        //    Florian Heiss, Viktor Winschel,
        //    Likelihood approximation by numerical integration on sparse grids,
        //    Journal of Econometrics,
        //    Volume 144, 2008, pages 62-80.
        //
        //  Parameters:
        //
        //    Input, int RULE_ORDER ( int l ), the name of a function which
        //    is given the level L and returns the order N of the corresponding rule.
        //
        //    Input, int DIM, the dimension of the integration problem.
        //
        //    Input, int K, the level.  When using the built in 1D 
        //    rules, the resulting sparse grid will be exact for polynomials up to total
        //    order 2*K-1.  When using the built-in quadrature rules, the maximum value 
        //    of K that is available is 25.
        //
        //    Output, int NWSPGR_SIZE, the "size" of the rule, that is,
        //    the number of weights and multidimensional quadrature points that will
        //    be needed.  The size of the rule will be reduced when duplicate points
        //    are merged.
        //
    {
        int level;
        int q;
        //
        //  Determine the size of each 1D rule.
        //
        int[] n1d = new int[k];

        for (level = 1; level <= k; level++)
        {
            int n = rule_order(level);
            n1d[level - 1] = n;
        }

        //
        //  Go through the motions of generating the rules.
        //
        int minq = Math.Max(0, k - dim);
        int maxq = k - 1;
        int r_size = 0;

        for (q = minq; q <= maxq; q++)
        {
            //
            //  Compute the D-dimensional vectors that sum to Q+DIM.
            //
            int seq_num = Comp.num_seq(q, dim);

            int[] is_ = Comp.get_seq(dim, q + dim, seq_num);
            //
            //  Determine the size of each rule.
            //
            int[] rq = new int[seq_num];

            int j;
            for (j = 0; j < seq_num; j++)
            {
                rq[j] = 1;
                int i;
                for (i = 0; i < dim; i++)
                {
                    rq[j] *= n1d[is_[j + i * seq_num] - 1];
                }
            }

            //
            //  Add the sizes to the total.
            //
            r_size += typeMethods.i4vec_sum(seq_num, rq);

        }

        return r_size;
    }
}