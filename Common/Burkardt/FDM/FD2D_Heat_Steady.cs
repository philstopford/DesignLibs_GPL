using System;
using Burkardt.Types;

namespace Burkardt.FDM;

public static class FD2D_Heat_Steady
{
    public static double[] fd2d_heat_steady(int nx, int ny, double[] x, double[] y,
            Func<double, double, double> d, Func<double, double, double> f)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FD2D_HEAT_STEADY solves the steady 2D heat equation.
        //
        //  Discussion:
        //
        //    Nodes are assigned a single index K, which increases as:
        //
        //    (NY-1)*NX+1  (NY-1)*NX+2  ...  NY * NX
        //           ....         ....  ...    .....
        //           NX+1         NX+2  ...   2 * NX
        //              1            2  ...       NX
        //
        //    Therefore, the neighbors of an interior node numbered C are
        //
        //             C+NY
        //              |
        //      C-1 --- C --- C+1
        //              |
        //             C-NY
        //
        //    Nodes on the lower boundary satisfy:
        //      1 <= K <= NX
        //    Nodes on the upper boundary satisfy:
        //      (NY-1)*NX+1 <= K <= NY * NX
        //    Nodes on the left boundary satisfy:
        //      mod ( K, NX ) = 1
        //    Nodes on the right boundary satisfy:
        //      mod ( K, NX ) = 0
        //
        //    If we number rows from bottom I = 1 to top I = NY
        //    and columns from left J = 1 to right J = NX, we have
        //      K = ( I - 1 ) * NX + J
        //    and
        //      J = 1 + mod ( K - 1, NX )
        //      I = 1 + ( K - J ) / NX
        //      
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NX, NY, the number of grid points in X and Y.
        //
        //    Input, double X[NX], Y[NY], the coordinates of grid lines.
        //
        //    Input, double D ( X, Y ), evaluates the thermal
        //    conductivity.
        //
        //    Input, double F ( X, Y ), evaluates the heat 
        //    source term.
        //
        //    Output, double FD2D_HEAT_STEADY[NX*NY], the approximation to the solution
        //    at the grid points.
        //
    {
        //
        //  Set the total number of unknowns.
        //
        int n = nx * ny;
        //
        //  Set up the matrix and right hand side.
        //
        double[] a = new double[n * n];
        double[] u = new double[n];
        //
        //  Define the matrix at interior points.
        //
        interior(nx, ny, x, y, d, f, n, ref a, ref u);
        //
        //  Handle boundary conditions.
        //
        boundary(nx, ny, x, y, n, a, u);
        //
        //  Solve the linear system.
        //
        typeMethods.r8mat_fs(n, a, ref u);

        return u;
    }

    public static void interior(int nx, int ny, double[] x, double[] y,
            Func<double, double, double> d, Func<double, double, double> f, int n,
            ref double[] a, ref double[] rhs)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INTERIOR sets up the matrix and right hand side at interior nodes.
        //
        //  Discussion:
        //
        //    Nodes are assigned a single index K, which increases as:
        //
        //    (NY-1)*NX+1  (NY-1)*NX+2  ...  NY * NX
        //           ....         ....  ...    .....
        //           NX+1         NX+2  ...   2 * NX
        //              1            2  ...       NX
        //
        //    Therefore, the neighbors of an interior node numbered C are
        //
        //             C+NY
        //              |
        //      C-1 --- C --- C+1
        //              |
        //             C-NY
        //
        //    If we number rows from bottom I = 1 to top I = NY
        //    and columns from left J = 1 to right J = NX, then the relationship
        //    between the single index K and the row and column indices I and J is:
        //      K = ( I - 1 ) * NX + J
        //    and
        //      J = 1 + mod ( K - 1, NX )
        //      I = 1 + ( K - J ) / NX
        //      
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NX, NY, the number of grid points in X and Y.
        //
        //    Input, double X[NX], Y[NY], the coordinates of grid lines.
        //
        //    Input, double D ( double X, double Y ), evaluates the thermal
        //    conductivity.
        //
        //    Input, double function F ( double X, double Y ), evaluates the heat 
        //    source term.
        //
        //    Input, int N, the number of nodes.
        //
        //    Output, double A[N*N], the system matrix, with the entries for 
        //    the interior nodes filled in.
        //
        //    Output, double RHS[N], the system right hand side, with the 
        //    entries for the interior nodes filled in.
        //
    {
        int ic;
        //
        //  For now, assume X and Y are equally spaced.
        //
        double dx = x[1] - x[0];
        double dy = y[1] - y[0];

        for (ic = 1; ic < ny - 1; ic++)
        {
            int jc;
            for (jc = 1; jc < nx - 1; jc++)
            {
                int in_ = ic + 1;
                int is_ = ic - 1;
                int je = jc + 1;
                int jw = jc - 1;

                int kc = ic * nx + jc;
                int ke = kc + 1;
                int kw = kc - 1;
                int kn = kc + nx;
                int ks = kc - nx;

                double dce = d(0.5 * (x[jc] + x[je]), y[ic]);
                double dcw = d(0.5 * (x[jc] + x[jw]), y[ic]);
                double dcn = d(x[jc], 0.5 * (y[ic] + y[in_]));
                double dcs = d(x[jc], 0.5 * (y[ic] + y[is_]));

                a[kc + kc * n] = (dce + dcw) / dx / dx + (dcn + dcs) / dy / dy;
                a[kc + ke * n] = -dce / dx / dx;
                a[kc + kw * n] = -dcw / dx / dx;
                a[kc + kn * n] = -dcn / dy / dy;
                a[kc + ks * n] = -dcs / dy / dy;

                rhs[kc] = f(x[jc], y[ic]);
            }
        }
    }

    public static void boundary(int nx, int ny, double[] x, double[] y, int n, double[] a,
            double[] rhs )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BOUNDARY sets up the matrix and right hand side at boundary nodes.
        //
        //  Discussion:
        //
        //    For this simple problem, the boundary conditions specify that the solution
        //    is 100 on the left side, and insulated on the right, top and bottom.
        //
        //    Nodes are assigned a single index K, which increases as:
        //
        //    (NY-1)*NX+1  (NY-1)*NX+2  ...  NY * NX
        //           ....         ....  ...    .....
        //           NX+1         NX+2  ...   2 * NX
        //              1            2  ...       NX
        //
        //    The index K of a node on the lower boundary satisfies:
        //      1 <= K <= NX
        //    The index K of a node on the upper boundary satisfies:
        //      (NY-1)*NX+1 <= K <= NY * NX
        //    The index K of a node on the left boundary satisfies:
        //      mod ( K, NX ) = 1
        //    The index K of a node on the right boundary satisfies:
        //      mod ( K, NX ) = 0
        //
        //    If we number rows from bottom I = 1 to top I = NY
        //    and columns from left J = 1 to right J = NX, then the relationship
        //    between the single index K and the row and column indices I and J is:
        //      K = ( I - 1 ) * NX + J
        //    and
        //      J = 1 + mod ( K - 1, NX )
        //      I = 1 + ( K - J ) / NX
        //      
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NX, NY, the number of grid points in X and Y.
        //
        //    Input, double X[NX], Y[NY], the coordinates of grid lines.
        //
        //    Input, int N, the number of nodes.
        //
        //    Input/output, double A[N*N].  On input, the system matrix, with the 
        //    entries for the interior nodes filled in.  On output, the entries for
        //    the boundary nodes have been set as well.
        //
        //    Input, double RHS[N], on input, the system right hand side, 
        //    with the entries for the interior nodes filled in.  On output, the entries for
        //    the boundary nodes have been set as well.
        //
    {
        int i;
        int kc;
        //
        //  Left boundary.
        //
        int j = 0;
        for (i = 1; i < ny - 1; i++)
        {
            kc = i * nx + j;
            a[kc + kc * n] += 1.0;
            rhs[kc] = 10.0;
        }

        //
        //  Right boundary.
        //
        j = nx - 1;
        for (i = 1; i < ny - 1; i++)
        {
            kc = i * nx + j;
            a[kc + kc * n] += 1.0;
            rhs[kc] = 100.0;
        }

        //
        //  Lower boundary.
        //
        i = 0;
        for (j = 0; j < nx; j++)
        {
            kc = i * nx + j;
            a[kc + kc * n] += 1.0;
            rhs[kc] = 0.0;
        }

        //
        //  Upper boundary.
        //
        i = ny - 1;
        for (j = 0; j < nx; j++)
        {
            kc = i * nx + j;
            a[kc + kc * n] += 1.0;
            rhs[kc] = 0.0;
        }
    }
}