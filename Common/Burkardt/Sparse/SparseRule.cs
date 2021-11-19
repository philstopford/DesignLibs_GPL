using System;

namespace Burkardt.Sparse;

public static class SparseRule
{
    public static int symmetric_sparse_size ( int nr, int dim, double[] nodes, double x0 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SYMMETRIC_SPARSE_SIZE sizes a symmetric sparse rule.
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
        //    Input, int DIM, the dimension.
        //
        //    Input, int NR, the dimension of the rule in the 
        //    positive orthant.
        //
        //    Input, double NODES[NR*DIM], the nodes for the positive orthant.
        //
        //    Input, double X0, the point of symmetry for the 1D rule, 
        //    typically 0.
        //
        //    Output, int SYMMETRIC_SPARSE_SIZE, the dimension of the rule 
        //    when "unfolded" to the full space.
        //
    {
        int j;
        int r;
        //
        //  Count the size of the full rule.
        //
        int nr2 = 0;

        for ( r = 0; r < nr; r++ )
        {
            int count = 1;
            for ( j = 0; j < dim; j++ )
            {
                if ( Math.Abs(nodes[r+j*nr] - x0) > double.Epsilon )
                {
                    count = 2 * count;
                }
            }
            nr2 += count;
        }

        return nr2;
    }
}