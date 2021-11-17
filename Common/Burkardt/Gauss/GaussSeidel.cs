using System;

namespace Burkardt;

public static class GaussSeidel
{
    public static void gauss_seidel ( int n, double[] r, ref double[] u, ref double dif_l1, int rIndex = 0, int uIndex = 0 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GAUSS_SEIDEL carries out one step of a Gauss-Seidel iteration.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 December 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    William Hager,
        //    Applied Numerical Linear Algebra,
        //    Prentice-Hall, 1988,
        //    ISBN13: 978-0130412942,
        //    LC: QA184.H33.
        //
        //  Parameters:
        //
        //    Input, int N, the number of unknowns.
        //
        //    Input, double R[N], the right hand side.
        //
        //    Input/output, double U[N], the estimated solution.
        //
        //    Output, double &DIF_L1, the L1 norm of the difference between the
        //    input and output solution estimates.
        //
    {
        int i;
        double u_old;

        dif_l1 = 0.0;

        for ( i = 1; i < n - 1; i++ )
        {
            u_old = u[uIndex + i];
            u[uIndex + i] = 0.5 * ( u[uIndex + i-1] + u[uIndex + i+1] + r[rIndex + i] );
            dif_l1 += Math.Abs ( u[uIndex + i] - u_old );
        }
    }
}