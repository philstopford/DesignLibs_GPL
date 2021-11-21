using System;
using Burkardt.Types;

namespace Burkardt.FEM;

public static class FEM_1D_Pack
{



    public static void legendre_com(int order, ref double[] xtab, ref double[] weight)
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
        int i;

        switch (order)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("LEGENDRE_COM - Fatal error!");
                Console.WriteLine("  Illegal value of NORDER = " + order + "");
                return;
        }

        double e1 = order * (order + 1);

        int m = (order + 1) / 2;

        for (i = 1; i <= (order + 1) / 2; i++)
        {
            int mp1mi = m + 1 - i;
            double t = Math.PI * (4 * i - 1) / (4 * order + 2);
            double x0 = Math.Cos(t) * (1.0 - (1.0 - 1.0 /
                order) / (8 * order * order));

            double pkm1 = 1.0;
            double pk = x0;

            int k;
            for (k = 2; k <= order; k++)
            {
                double pkp1 = 2.0 * x0 * pk - pkm1 - (x0 * pk - pkm1) / k;
                pkm1 = pk;
                pk = pkp1;
            }

            double d1 = order * (pkm1 - x0 * pk);

            double dpn = d1 / (1.0 - x0 * x0);

            double d2pn = (2.0 * x0 * dpn - e1 * pk) / (1.0 - x0 * x0);

            double d3pn = (4.0 * x0 * d2pn + (2.0 - e1) * dpn) / (1.0 - x0 * x0);

            double d4pn = (6.0 * x0 * d3pn + (6.0 - e1) * d2pn) / (1.0 - x0 * x0);

            double u = pk / dpn;
            double v = d2pn / dpn;
            //
            //  Initial approximation H:
            //
            double h = -u * (1.0 + 0.5 * u * (v + u * (v * v - d3pn
                / (3.0 * dpn))));
            //
            //  Refine H using one step of Newton's method:
            //
            double p = pk + h * (dpn + 0.5 * h * (d2pn + h / 3.0
                * (d3pn + 0.25 * h * d4pn)));

            double dp = dpn + h * (d2pn + 0.5 * h * (d3pn + h * d4pn / 3.0));

            h -= p / dp;

            double xtemp = x0 + h;

            xtab[mp1mi - 1] = xtemp;

            double fx = d1 - h * e1 * (pk + 0.5 * h * (dpn + h / 3.0
                * (d2pn + 0.25 * h * (d3pn + 0.2 * h * d4pn))));

            weight[mp1mi - 1] = 2.0 * (1.0 - xtemp * xtemp) / (fx * fx);
        }

        xtab[0] = (order % 2) switch
        {
            1 => 0.0,
            _ => xtab[0]
        };

        //
        //  Shift the data up.
        //
        int nmove = (order + 1) / 2;
        int ncopy = order - nmove;

        for (i = 1; i <= nmove; i++)
        {
            int iback = order + 1 - i;
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


    public static double[] local_fem_1d(int order, double[] node_x, double[] node_v,
            int sample_num, double[] sample_x)
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
        double[] sample_v = new double[sample_num];

        for (int sample = 0; sample < sample_num; sample++)
        {
            double x = sample_x[sample];
            double[] phi = LocalBasis.local_basis_1d(order, node_x, x);
            sample_v[sample] = typeMethods.r8vec_dot_product(order, node_v, phi);
        }

        return sample_v;
    }

}