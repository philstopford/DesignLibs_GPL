using System;
using Burkardt.Types;

namespace Burkardt;

public static class Diaphony
{
    public static double diaphony_compute ( int dim_num, int point_num, double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIAPHONY_COMPUTE evaluates the diaphony of a N-dimensional point set.
        //
        //  Discussion:
        //
        //    The diaphony is analogous to, and related to, the discrepancy,
        //    and is a measure of how well spread a set of point is.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 January 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Peter Heelekalek, Harald Niederreiter,
        //    The Weighted Spectral Test: Diaphony,
        //    ACM Transactions on Modeling and Computer Simulation,
        //    Volume 8, Number 1, January 1998, pages 43-60.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int POINT_NUM, the number of points.
        //
        //    Input, double X[DIM_NUM*POINT_NUM], the point set, which is
        //    presumed to lie in the DIM_NUM dimensional unit hypercube.
        //
        //    Output, double DIAPHONY_COMPUTE, the value of the diaphony.
        //
    {
        int i;

        double d = 0.0;

        for ( i = 0; i < point_num; i++ )
        {
            int j;
            for ( j = 0; j < point_num; j++ )
            {
                double prod = 1.0;
                int k;
                for ( k = 0; k < dim_num; k++ )
                {
                    double z = typeMethods.r8_modp ( x[k+i*dim_num] - x[k+j*dim_num], 1.0 );
                    prod *= 1.0 + 2.0 * Math.PI * Math.PI * ( z * z - z + 1.0 / 6.0 );
                }
                d = d + prod - 1.0;
            }
        }

        double bot = Math.Pow ( point_num, 2 )
                     * ( Math.Pow ( 1.0 + Math.PI * Math.PI / 3.0, dim_num ) - 1.0 );

        d /= bot;

        d = Math.Sqrt ( d );

        return d;
    }
}