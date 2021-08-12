namespace Burkardt.ODE
{
    public static class Lorenz
    {
        public static double[] lorenz_rhs ( double t, int m, double[] x, int index )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LORENZ_RHS evaluates the right hand side of the Lorenz ODE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    10 October 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double T, the value of the independent variable.
            //
            //    Input, int M, the spatial dimension.
            //
            //    Input, double X[M], the values of the dependent variables
            //    at time T.
            //
            //    Output, double DXDT[M], the values of the derivatives
            //    of the dependent variables at time T.
            //
        {
            double beta = 8.0 / 3.0;
            double[] dxdt;
            double rho = 28.0;
            double sigma = 10.0;

            dxdt = new double[m];

            dxdt[0] = sigma * ( x[index + 1] - x[index + 0] );
            dxdt[1] = x[index + 0] * ( rho - x[index + 2] ) - x[index + 1];
            dxdt[2] = x[index + 0] * x[index + 1] - beta * x[index + 2];

            return dxdt;
        }
    }
}