namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static void tensor_product(int d, int[] order1d, int n1d, double[] x1d,
        double[] w1d, int n, ref double[] xnd, ref double[] wnd )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TENSOR_PRODUCT generates a tensor product quadrature rule.
        //
        //  Discussion:
        //
        //    The Kronecker product of an K by L matrix A and an M by N matrix B
        //    is the K*M by L*N matrix formed by
        //
        //      a(1,1) * B,  a(1,2) * B,  ..., a(1,l) * B
        //      a(2,1) * B,  a(2,2) * B,  ..., a(2,l) * B
        //      ..........   ..........   .... ..........
        //      a(k,1) * B,  a(k,2) * B,  ..., a(k,l) * B
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 December 2012
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
        //    Input, int D, the spatial dimension.
        //
        //    Input, int ORDER1D[D], the order of each 1D rule.
        //
        //    Input, int N1D, the number of 1D items.
        //
        //    Input, double X1D[N1D], the 1D nodes.
        //
        //    Input, double W1D[N1D], the 1D weights.
        //
        //    Input, int N, the number of N-dimensional items.
        //
        //    Output, double XND[D*N], the nodes.
        //
        //    Output, double WND[N], the weights.
        //
        {
            int i;
            int i1;
            int i2;
            r8vecDPData data1 = new r8vecDPData();
            r8vecDPData data2 = new r8vecDPData();
            
            //
            //  Compute the weights.
            //
            i2 = -1;
            for (i = 0; i < d; i++)
            {
                i1 = i2 + 1;
                i2 = i2 + order1d[i];
                r8vec_direct_product2(ref data2, i, order1d[i], w1d, d, n, ref wnd, factorValueIndex: + i1);
            }

            //
            //  Compute the points.
            //
            i2 = -1;
            for (i = 0; i < d; i++)
            {
                i1 = i2 + 1;
                i2 = i2 + order1d[i];
                r8vec_direct_product(ref data1, i, order1d[i], x1d, d, n, ref xnd, factorValueIndex: + i1);
            }
        }

        public static void tensor_product_cell(int nc, double[] xc, double[] wc, int dim, int[] nr,
        int[] roff, int np, ref double[] xp, ref double[] wp )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TENSOR_PRODUCT_CELL generates a tensor product quadrature rule.
        //
        //  Discussion:
        //
        //    The Kronecker product of an K by L matrix A and an M by N matrix B
        //    is the K*M by L*N matrix formed by
        //
        //      a(1,1) * B,  a(1,2) * B,  ..., a(1,l) * B
        //      a(2,1) * B,  a(2,2) * B,  ..., a(2,l) * B
        //      ..........   ..........   .... ..........
        //      a(k,1) * B,  a(k,2) * B,  ..., a(k,l) * B
        //
        //    The 1D factors are stored in a kind of cell array structure.
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
        //  Parameters:
        //
        //    Input, int NC, the number of items in the cell arrays.
        //
        //    Input, double XC[NC], a cell array containing points for 
        //    1D rules.
        //
        //    Input, double WC[NC], a cell array containing weights for
        //    1D rules.
        //
        //    Input, int DIM, the spatial dimension.
        //
        //    Input, int NR[DIM], the length of each row of the 
        //    cell array.
        //
        //    Input, int ROFF[DIM+1], offsets for the cell arrays.
        //
        //    Input, int NP, the number of points in the product rule.
        //
        //    Output, double XP[DIM*NP], the nodes.
        //
        //    Output, double WP[NP], the weights.
        //
        {
            int i;
            int n1d;
            double[] w1d;
            double[] x1d;
            r8vecDPData data1 = new r8vecDPData();
            r8vecDPData data2 = new r8vecDPData();
            //
            //  Compute the weights.
            //
            for (i = 0; i < dim; i++)
            {
                n1d = nr[i];
                w1d = r8cvv_rget_new(nc, wc, dim, roff, i);
                r8vec_direct_product2(ref data2, i, n1d, w1d, dim, np, ref wp);
            }

            //
            //  Compute the points.
            //
            for (i = 0; i < dim; i++)
            {
                n1d = nr[i];
                x1d = r8cvv_rget_new(nc, xc, dim, roff, i);
                r8vec_direct_product(ref data1, i, n1d, x1d, dim, np, ref xp);
            }
        }
    }
}