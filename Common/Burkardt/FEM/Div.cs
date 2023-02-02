using System;
using Burkardt.Types;

namespace Burkardt.FEM;

public class Div
{
    public static void div_q4(int m, int n, double[] u, double[] v, double xlo, double xhi,
            double ylo, double yhi, ref double[] div, ref double[] vort )

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    DIV_Q4 estimates the divergence and vorticity of a discrete field.
        //
        //  Discussion:
        //
        //    The routine is given the values of a vector field ( U(X,Y), V(X,Y) ) at
        //    an array of points ( X(1:M), Y(1:N) ).
        //
        //    The routine models the vector field over the interior of this region using
        //    a bilinear interpolant.  It then uses the interpolant to estimate the
        //    value of the divergence:
        //
        //      DIV(X,Y) = dU/dX + dV/dY
        //
        //    and the vorticity:
        //
        //      VORT(X,Y) = dV/dX - dU/dY
        //
        //    at the center point of each of the bilinear elements.
        //
        //        |       |       |
        //      (3,1)---(3,2)---(3,3)---
        //        |       |       |
        //        | [2,1] | [2,2] |
        //        |       |       |
        //      (2,1)---(2,2)---(2,3)---
        //        |       |       |
        //        | [1,1] | [1,2] |
        //        |       |       |
        //      (1,1)---(1,2)---(1,3)---
        //
        //    Here, the nodes labeled with parentheses represent the points at
        //    which the original (U,V) data is given, while the nodes labeled
        //    with square brackets represent the centers of the bilinear
        //    elements, where the approximations to the divergence and vorticity
        //    are made.
        //
        //    The reason for evaluating the divergence and vorticity in this way
        //    is that the bilinear interpolant to the (U,V) data is not
        //    differentiable at the boundaries of the elements, nor especially at
        //    the nodes, but is an (infinitely differentiable) bilinear function
        //    in the interior of each element.  If a value at the original nodes
        //    is strongly desired, then the average at the four surrounding
        //    central nodes may be taken.
        //
        //  Element Q4:
        //
        //    |
        //    1  4-----3
        //    |  |     |
        //    |  |     |
        //    S  |     |
        //    |  |     |
        //    |  |     |
        //    0  1-----2
        //    |
        //    +--0--R--1-->
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 February 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of data rows.  M must be at least 2.
        //
        //    Input, int N, the number of data columns.  N must be at least 2.
        //
        //    Input, double U[M*N], V[M*N], the value of the components 
        //    of a vector quantity whose divergence and vorticity are desired. 
        //    A common example would be that U and V are the horizontal and 
        //    vertical velocity component of a flow field.
        //
        //    Input, double XLO, XHI, the minimum and maximum X coordinates.
        //
        //    Input, double YLO, YHI, the minimum and maximum Y coordinates.
        //
        //    Output, double DIV[(M-1)*(N-1)], an estimate for
        //    the divergence in the bilinear element that lies between
        //    data rows I and I+1, and data columns J and J+1.
        //
        //    Output, double VORT[(M-1)*(N-1)], an estimate for
        //    the vorticity in the bilinear element that lies between
        //    data rows I and I+1, and data columns J and J+1.
        //
    {
        double[] dphidx = new double[4];
        double[] dphidy = new double[4];
        int i;
        double[] p = new double[2];
        double[] phi = new double[4];
        double[] q = new double[2 * 4];

        switch (m)
        {
            case <= 1:
                Console.WriteLine("");
                Console.WriteLine("DIV_Q4 - Fatal error!");
                Console.WriteLine("  M must be at least 2,");
                Console.WriteLine("  but the input value of M is " + m + "");
                return;
        }

        switch (n)
        {
            case <= 1:
                Console.WriteLine("");
                Console.WriteLine("DIV_Q4 - Fatal error!");
                Console.WriteLine("  N must be at least 2,");
                Console.WriteLine("  but the input value of N is " + n + "");
                return;
        }

        if (Math.Abs(xhi - xlo) <= typeMethods.r8_epsilon())
        {
            Console.WriteLine("");
            Console.WriteLine("DIV_Q4 - Fatal error!");
            Console.WriteLine("  XHI and XLO must be distinct,");
            Console.WriteLine("  but the input value of XLO is " + xlo + "");
            Console.WriteLine("  and the input value of XHI is " + xhi + "");
            return;
        }

        if (Math.Abs(yhi - ylo) <= typeMethods.r8_epsilon())
        {
            Console.WriteLine("");
            Console.WriteLine("DIV_Q4 - Fatal error!");
            Console.WriteLine("  YHI and YLO must be distinct,");
            Console.WriteLine("  but the input value of YLO is " + ylo + "");
            Console.WriteLine("  and the input value of YHI is " + yhi + "");
            return;
        }

        for (i = 1; i <= m - 1; i++)
        {
            double yb = ((2 * m - 2 * i) * ylo
                         + (2 * i - 2) * yhi)
                        / (2 * m - 2);
            p[1] = ((2 * m - 2 * i - 1) * ylo
                    + (2 * i - 1) * yhi)
                   / (2 * m - 2);
            double yt = ((2 * m - 2 * i - 2) * ylo
                         + 2 * i * yhi)
                        / (2 * m - 2);

            q[1 + 0 * 2] = yb;
            q[1 + 1 * 2] = yb;
            q[1 + 2 * 2] = yt;
            q[1 + 3 * 2] = yt;

            int j;
            for (j = 1; j <= n - 1; j++)
            {
                double xl = ((2 * n - 2 * j) * xlo
                             + (2 * j - 2) * xhi)
                            / (2 * n - 2);
                p[0] = ((2 * n - 2 * j - 1) * xlo
                        + (2 * j - 1) * xhi)
                       / (2 * n - 2);
                double xr = ((2 * n - 2 * j - 2) * xlo
                             + 2 * j * xhi)
                            / (2 * n - 2);

                q[0 + 0 * 2] = xl;
                q[0 + 1 * 2] = xr;
                q[0 + 2 * 2] = xr;
                q[0 + 3 * 2] = xl;
                //
                //  Evaluate the basis function and derivatives at the center of the element.
                //
                Basis_mn.basis_mn_q4(q, 1, p, ref phi, ref dphidx, ref dphidy);
                //
                //  Note the following formula for the value of U and V at the same
                //  point that the divergence and vorticity are being evaluated.
                //
                //         umid =  u(i  ,j  ) * phi[0] &
                //               + u(i  ,j+1) * phi[1] &
                //               + u(i+1,j+1) * phi[2] &
                //               + u(i+1,j  ) * phi[3] 
                //
                //         vmid =  v(i  ,j  ) * phi[0] &
                //               + v(i  ,j+1) * phi[1] &
                //               + v(i+1,j+1) * phi[2] &
                //               + v(i+1,j  ) * phi[3] 
                //
                div[i - 1 + (j - 1) * (m - 1)] =
                    u[i - 1 + (j - 1) * m] * dphidx[0] + v[i - 1 + (j - 1) * m] * dphidy[0]
                                                       + u[i - 1 + j * m] * dphidx[1] +
                                                       v[i - 1 + j * m] * dphidy[1]
                                                       + u[i + j * m] * dphidx[2] + v[i + j * m] * dphidy[2]
                                                       + u[i + (j - 1) * m] * dphidx[3] +
                                                       v[i + (j - 1) * m] * dphidy[3];

                vort[i - 1 + (j - 1) * (m - 1)] =
                    v[i - 1 + (j - 1) * m] * dphidx[0] - u[i - 1 + (j - 1) * m] * dphidy[0]
                    + v[i - 1 + j * m] * dphidx[1] - u[i - 1 + j * m] * dphidy[1]
                    + v[i + j * m] * dphidx[2] - u[i + j * m] * dphidy[2]
                    + v[i + (j - 1) * m] * dphidx[3] - u[i + (j - 1) * m] * dphidy[3];
            }
        }

    }
}