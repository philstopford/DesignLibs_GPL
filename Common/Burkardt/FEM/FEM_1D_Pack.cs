using System;
using Burkardt.Types;

namespace Burkardt.FEM
{
    public static class FEM_1D_Pack
    {

        public static void bandwidth_mesh(int element_order, int element_num, int[] element_node,
        ref int ml, ref int mu, ref int m )
//****************************************************************************80
//
//  Purpose:
//
//    BANDWIDTH_MESH determines the bandwidth of the coefficient matrix.
//
//  Discussion:
//
//    The quantity computed here is the "geometric" bandwidth determined
//    by the finite element mesh alone.
//
//    If a single finite element variable is associated with each node
//    of the mesh, and if the nodes and variables are numbered in the
//    same way, then the geometric bandwidth is the same as the bandwidth
//    of a typical finite element matrix.
//
//    The bandwidth M is defined in terms of the lower and upper bandwidths:
//
//      M = ML + 1 + MU
//
//    where 
//
//      ML = maximum distance from any diagonal entry to a nonzero
//      entry in the same row, but earlier column,
//
//      MU = maximum distance from any diagonal entry to a nonzero
//      entry in the same row, but later column.
//
//    Because the finite element node adjacency relationship is symmetric,
//    we are guaranteed that ML = MU.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 January 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ELEMENT_ORDER, the order of the elements.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM];
//    ELEMENT_NODE(I,J) is the global index of local node I in element J.
//
//    Output, int *ML, *MU, the lower and upper bandwidths of the matrix.
//
//    Output, int *M, the bandwidth of the matrix.
//
        {
            int element;
            int global_i;
            int global_j;
            int local_i;
            int local_j;

            ml = 0;
            mu = 0;

            for (element = 0; element < element_num; element++)
            {
                for (local_i = 0; local_i < element_order; local_i++)
                {
                    global_i = element_node[local_i + element * element_order];

                    for (local_j = 0; local_j < element_order; local_j++)
                    {
                        global_j = element_node[local_j + element * element_order];

                        mu = Math.Max(mu, global_j - global_i);
                        ml = Math.Max(ml, global_i - global_j);
                    }
                }
            }

            m = ml + 1 + mu;
        }


        public static void legendre_com(int order, ref double[] xtab, ref double[] weight )
//****************************************************************************80
//
//  Purpose: 
//
//    LEGENDRE_COM computes abscissas and weights for Gauss-Legendre quadrature.
//
//  Integration interval:
//
//    [ -1, 1 ]
//
//  Weight function:
//
//    1.
//
//  Integral to approximate:
//
//    Integral ( -1 <= X <= 1 ) F(X) dX.
//
//  Approximate integral:
//
//    sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NORDER, the order of the rule.
//    NORDER must be greater than 0.
//
//    Output, double XTAB[NORDER], the abscissas of the rule.
//
//    Output, double WEIGHT[NORDER], the weights of the rule.
//    The weights are positive, symmetric, and should sum to 2.
//
        {
            double d1;
            double d2pn;
            double d3pn;
            double d4pn;
            double dp;
            double dpn;
            double e1;
            double fx;
            double h;
            int i;
            int iback;
            int k;
            int m;
            int mp1mi;
            int ncopy;
            int nmove;
            double p;
            double pk;
            double pkm1;
            double pkp1;
            double t;
            double u;
            double v;
            double x0;
            double xtemp;

            if (order < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("LEGENDRE_COM - Fatal error!");
                Console.WriteLine("  Illegal value of NORDER = " + order + "");
                return;
            }

            e1 = (double) (order * (order + 1));

            m = (order + 1) / 2;

            for (i = 1; i <= (order + 1) / 2; i++)
            {
                mp1mi = m + 1 - i;
                t = Math.PI * (double) (4 * i - 1) / (double) (4 * order + 2);
                x0 = Math.Cos(t) * (1.0 - (1.0 - 1.0 /
                    (double) (order)) / (double) (8 * order * order));

                pkm1 = 1.0;
                pk = x0;

                for (k = 2; k <= order; k++)
                {
                    pkp1 = 2.0 * x0 * pk - pkm1 - (x0 * pk - pkm1) / (double) (k);
                    pkm1 = pk;
                    pk = pkp1;
                }

                d1 = (double) (order) * (pkm1 - x0 * pk);

                dpn = d1 / (1.0 - x0 * x0);

                d2pn = (2.0 * x0 * dpn - e1 * pk) / (1.0 - x0 * x0);

                d3pn = (4.0 * x0 * d2pn + (2.0 - e1) * dpn) / (1.0 - x0 * x0);

                d4pn = (6.0 * x0 * d3pn + (6.0 - e1) * d2pn) / (1.0 - x0 * x0);

                u = pk / dpn;
                v = d2pn / dpn;
//
//  Initial approximation H:
//
                h = -u * (1.0 + 0.5 * u * (v + u * (v * v - d3pn
                    / (3.0 * dpn))));
//
//  Refine H using one step of Newton's method:
//
                p = pk + h * (dpn + 0.5 * h * (d2pn + h / 3.0
                    * (d3pn + 0.25 * h * d4pn)));

                dp = dpn + h * (d2pn + 0.5 * h * (d3pn + h * d4pn / 3.0));

                h = h - p / dp;

                xtemp = x0 + h;

                xtab[mp1mi - 1] = xtemp;

                fx = d1 - h * e1 * (pk + 0.5 * h * (dpn + h / 3.0
                    * (d2pn + 0.25 * h * (d3pn + 0.2 * h * d4pn))));

                weight[mp1mi - 1] = 2.0 * (1.0 - xtemp * xtemp) / (fx * fx);
            }

            if ((order % 2) == 1)
            {
                xtab[0] = 0.0;
            }

//
//  Shift the data up.
//
            nmove = (order + 1) / 2;
            ncopy = order - nmove;

            for (i = 1; i <= nmove; i++)
            {
                iback = order + 1 - i;
                xtab[iback - 1] = xtab[iback - ncopy - 1];
                weight[iback - 1] = weight[iback - ncopy - 1];
            }

//
//  Reflect values for the negative abscissas.
//
            for (i = 0; i < order - nmove; i++)
            {
                xtab[i] = -xtab[order - 1 - i];
                weight[i] = weight[order - 1 - i];
            }

        }

        public static double[] local_basis_1d(int order, double[] node_x, double x )
//****************************************************************************80
//
//  Purpose:
//
//    LOCAL_BASIS_1D evaluates the basis functions in an element.
//
//  Discussion:
//
//    PHI(I)(X) = product ( J ~= I ) ( X         - NODE_X(I) ) 
//                                 / ( NODE_X(J) - NODE_X(I) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 March 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ORDER, the order of the element.
//    0 <= ORDER.  ORDER = 1 means piecewise linear.
//
//    Input, double NODE_X[ORDER], the element nodes.  
//    These must be distinct.  Basis function I is 1 when X = NODE_X(I) 
//    and 0 when X is equal to any other node.
//
//    Input, double X, the point at which the basis functions are to 
//    be evaluated.
//
//    Output, double LOCAL_BASIS_1D[ORDER], the basis functions.
//
        {
            double[] phi = new double[order];

            for (int j = 0; j < order; j++)
            {
                phi[j] = 1.0;
                for (int i = 0; i < order; i++)
                {
                    if (j != i)
                    {
                        phi[j] = (phi[j] * (x - node_x[i])) / (node_x[j] - node_x[i]);
                    }
                }
            }

            return phi;
        }

        public static double[] local_basis_prime_1d(int order, double[] node_x, double x )
//****************************************************************************80
//
//  Purpose:
//
//    LOCAL_BASIS_PRIME_1D evaluates the basis function derivatives in an element.
//
//  Discussion:
//
//    PHI(I)(X) = product ( J ~= I ) ( X - NODE_X(I) ) 
//                                 / ( NODE_X(J) - NODE_X(I) )
//
//    dPHIdx(I)(X) = sum ( J ~= I ) ( 1 / ( NODE_X(J) - NODE_X(I) ) *
//      product ( K ~= ( J, I ) ) ( X - NODE_X(I) ) / ( NODE_X(J) - NODE_X(I) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ORDER, the order of the element.
//    0 <= ORDER.  ORDER = 1 means piecewise linear.
//
//    Input, double NODE_X[ORDER], the element nodes.  
//    These must be distinct.  Basis function I is 1 when X = NODE_X(I) 
//    and 0 when X is equal to any other node.
//
//    Input, double X, the point at which the basis functions are to 
//    be evaluated.
//
//    Output, double LOCAL_BASIS_PRIME_1D[ORDER], the basis functions.
//
        {
            double[] dphidx = new double[order];

            for (int i = 0; i < order; i++)
            {
                dphidx[i] = 0.0;
                for (int j = 0; j < order; j++)
                {
                    if (j != i)
                    {
                        double term = 1.0 / (node_x[j] - node_x[i]);
                        for (int k = 0; k < order; k++)
                        {
                            if (k != i && k != j)
                            {
                                term = term * (x - node_x[i]) / (node_x[k] - node_x[i]);
                            }
                        }

                        dphidx[i] = dphidx[i] + term;
                    }
                }
            }

            return dphidx;
        }

        public static double[] local_fem_1d(int order, double[] node_x, double[] node_v,
        int sample_num, double[] sample_x )
//****************************************************************************80
//
//  Purpose:
//
//    LOCAL_FEM_1D evaluates a local finite element function.
//
//  Discussion:
//
//    A local finite element function is a finite element function
//    defined over a single element.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 March 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ORDER, the order of the element.
//    0 <= ORDER.  ORDER = 1 means piecewise linear.
//
//    Input, double NODE_X[ORDER], the element nodes.  
//    These must be distinct.  Basis function I is 1 when X = NODE_X(I) and 0 
//    when X is equal to any other node.
//
//    Input, double NODE_V[ORDER], the value of the finite element 
//    function at each node.
//
//    Input, int SAMPLE_NUM, the number of sample points.
//
//    Input, double SAMPLE_X[SAMPLE_NUM], the sample points at which 
//    the local finite element function is to be evaluated.
//
//    Output, double LOCAL_FEM_1D[SAMPLE_NUM], the values of the local 
//    finite element basis functions.
//
        {
            double x;

            double[] sample_v = new double[sample_num];

            for (int sample = 0; sample < sample_num; sample++)
            {
                x = sample_x[sample];
                double[] phi = local_basis_1d(order, node_x, x);
                sample_v[sample] = typeMethods.r8vec_dot_product(order, node_v, phi);
            }

            return sample_v;
        }

    }
}