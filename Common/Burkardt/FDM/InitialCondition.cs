using System;

namespace Burkardt.FDM
{
    public static class InitialCondition
    {
        public         static double[] initial_condition ( int nx, double[] x )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    INITIAL_CONDITION sets the initial condition.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 December 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NX, the number of nodes.
            //
            //    Input, double X[NX], the coordinates of the nodes.
            //
            //    Output, double INITIAL_CONDITION[NX], the value of the initial condition.
            //
        {
            int i;
            double[] u;

            u = new double[nx];

            for ( i = 0; i < nx; i++ )
            {
                if  ( 0.4 <= x[i] && x[i] <= 0.6 )
                {
                    u[i] = Math.Pow ( 10.0 * x[i] - 4.0, 2 )
                           * Math.Pow ( 6.0 - 10.0 * x[i], 2 );
                }
                else
                {
                    u[i] = 0.0;
                }
            }
            return u;
        }
    }
}