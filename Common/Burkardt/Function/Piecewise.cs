namespace Burkardt.Function;

public static class Piecewise
{
    public static double[] piecewise_linear ( int nd, double[] xd, double[] yd, int nv, double[] xv )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PIECEWISE_LINEAR evaluates a piecewise linear spline.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 January 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ND, the number of data points.
        //
        //    Input, double XD[ND], YD[ND], the data values.
        //
        //    Input, int NV, the number of evaluation points.
        //
        //    Input, double XV[NV], the evaluation arguments.
        //
        //    Output, double PIECEWISE_LINEAR[NV], the values.
        //
    {
        int id;
        int iv;
        double[] yv;

        yv = new double[nv];

        for ( iv = 0; iv < nv; iv++ )
        {
            if ( xv[iv] < xd[0] )
            {
                yv[iv] = yd[0];
            }
            else if ( xd[nd-1] < xv[iv] )
            {
                yv[iv] = yd[nd-1];
            }
            else
            {
                for ( id = 1; id < nd; id++ )
                {
                    if ( xv[iv] < xd[id] )
                    {
                        yv[iv] = ( ( xd[id] - xv[iv]            ) * yd[id-1] 
                                   + (          xv[iv] - xd[id-1] ) * yd[id] ) 
                                 / ( xd[id]          - xd[id-1] );
                        break;
                    }
                }
            }
        }
        return yv;
    }
}