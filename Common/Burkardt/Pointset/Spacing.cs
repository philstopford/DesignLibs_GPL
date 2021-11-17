using System;
using Burkardt.Types;

namespace Burkardt.Pointset;

public static class Spacing
{
    public static double[] pointset_spacing ( int dim_num, int n, double[] z )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POINTSET_SPACING determines the minimum spacing between points in the set.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Max Gunzburger and John Burkardt,
        //    Uniformity Measures for Point Samples in Hypercubes.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int N, the number of points.
        //
        //    Input, double Z[DIM_NUM*N], the point distribution.
        //
        //    Output, double POINTSET_SPACING(N), the minimum distance between each
        //    point and a distinct point in the set.
        //
    {
        double dist;
        int i;
        int j1;
        int j2;
        double[] gamma;

        gamma = new double[n];

        for ( j1 = 0; j1 < n; j1++ )
        {
            gamma[j1] = typeMethods.r8_huge ( );

            for ( j2 = 0; j2 < n; j2++ )
            {
                if ( j2 != j1 )
                {
                    dist = 0.0;
                    for ( i = 0; i < dim_num; i++ )
                    {
                        dist += Math.Pow ( z[i+j1*dim_num] - z[i+j2*dim_num], 2 );
                    }
                    gamma[j1] = Math.Min ( gamma[j1], dist );
                }
            }
        }

        for ( j1 = 0; j1 < n; j1++ )
        {
            gamma[j1] = Math.Sqrt ( gamma[j1] );
        }

        return gamma;
    }
    //*****
}