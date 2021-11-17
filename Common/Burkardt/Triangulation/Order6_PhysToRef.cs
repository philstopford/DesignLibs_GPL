using System;

namespace Burkardt.TriangulationNS;

public static partial class Triangulation
{
    public static void triangle_order6_physical_to_reference(double[] t, int n,
            double[] phy, ref double[] ref_ )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_ORDER6_PHYSICAL_TO_REFERENCE maps a physical point to a reference point.
        //
        //  Discussion:
        //
        //    Given the vertices of an order 6 physical triangle and a point
        //    (X,Y) in the physical triangle, the routine computes the value
        //    of the corresponding image point (R,S) in reference space.
        //
        //    The mapping from (R,S) to (X,Y) has the form:
        //
        //      X(R,S) = A1 * R * R + B1 * R * S + C1 * S * S
        //             + D1 * R     + E1 * S     + F1
        //
        //      Y(R,S) = A2 * R * R + B2 * R * S + C2 * S * S
        //             + D2 * R     + E2 * S     + F2
        //
        //  Reference Element T3:
        //
        //    |
        //    1  3
        //    |  |.
        //    |  | .
        //    S  6  5
        //    |  |   .
        //    |  |    .
        //    0  1--4--2
        //    |
        //    +--0--R--1-->
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double T(2,6), the coordinates of the vertices.
        //    The vertices are assumed to be the images of (0,0), (1,0), (0,1),
        //    (1/2,0), (1/2,1/2) and (0,1/2), in that order.
        //
        //    Input, int N, the number of points to transform.
        //
        //    Input, double PHY(2,N), the coordinates of points in the
        //    physical space.
        //
        //    Output, double REF(2,N), the coordinates of the corresponding
        //    points in the reference space.
        //
    {
        double[] a = new double[2];
        double[] b = new double[2];
        double[] c = new double[2];
        double[] d = new double[2];
        double det;
        double[] dx = new double[2];
        double[] e = new double[2];
        double[] f = new double[2];
        double[] fun = new double[2];
        double fun_norm;
        int i;
        int it;
        int j;
        double[] jac = new double[2 * 2];
        int it_max = 10;
        double it_tol = 0.000001;
        //
        //  Set iteration parameters.
        //
        for (i = 0; i < 2; i++)
        {
            a[i] = 2.0 * t[i + 0 * 2] + 2.0 * t[i + 1 * 2] - 4.0 * t[i + 3 * 2];
            b[i] = 4.0 * t[i + 0 * 2] - 4.0 * t[i + 3 * 2] + 4.0 * t[i + 4 * 2] - 4.0 * t[i + 5 * 2];
            c[i] = 2.0 * t[i + 0 * 2] + 2.0 * t[i + 2 * 2] - 4.0 * t[i + 5 * 2];

            d[i] = -3.0 * t[i + 0 * 2] - t[i + 1 * 2] + 4.0 * t[i + 3 * 2];
            e[i] = -3.0 * t[i + 0 * 2] - t[i + 2 * 2] + 4.0 * t[i + 5 * 2];

            f[i] = t[i + 0 * 2];
        }

        //
        //  Initialize the points by inverting the linear map.
        //
        triangle_order3_physical_to_reference(t, n, phy, ref ref_);
        //
        //  Carry out the Newton iteration.
        //
        for (j = 0; j < n; j++)
        {
            for (it = 0; it < it_max; it++)
            {
                for (i = 0; i < 2; i++)
                {
                    fun[i] = a[i] * ref_[
                                 0 + j * 2] * ref_[
                                 0 + j * 2]
                             +b[i] * ref_[
                                 0 + j * 2] * ref_[
                                 1 + j * 2]
                             +c[i] * ref_[
                                 1 + j * 2] * ref_[
                                 1 + j * 2]
                             +d[i] * ref_[
                                 0 + j * 2]
                             +e[i] * ref_[
                                 1 + j * 2]
                             +f[i]
                             - phy[i + j * 2];
                }

                fun_norm = Math.Sqrt(Math.Pow(fun[0], 2) + Math.Pow(fun[1], 2));

                if (fun_norm <= it_tol)
                {
                    break;
                }

                jac[0 + 0 * 2] = 2.0 * a[0] * ref_[
                    0 + j * 2] +b[0] * ref_[
                    1 + j * 2] +d[0];
                jac[1 + 0 * 2] = 2.0 * a[1] * ref_[
                    0 + j * 2] +b[1] * ref_[
                    1 + j * 2] +d[1];
                jac[0 + 1 * 2] = b[0] * ref_[
                    0 + j * 2] +2.0 * c[0] * ref_[
                    1 + j * 2] +e[0];
                jac[1 + 1 * 2] = b[1] * ref_[
                    0 + j * 2] +2.0 * c[1] * ref_[
                    1 + j * 2] +e[1];

                det = jac[0 + 0 * 2] * jac[1 + 1 * 2] - jac[0 + 1 * 2] * jac[1 + 0 * 2];

                switch (det)
                {
                    case 0.0:
                        Console.WriteLine("");
                        Console.WriteLine("TRIANGLE_ORDER6_PHYSICAL_TO_REFERENCE - Fatal error!");
                        Console.WriteLine("  The jacobian of the mapping is singular.");
                        break;
                }

                dx[0] = (jac[1 + 1 * 2] * fun[0] - jac[0 + 1 * 2] * fun[1]) / det;
                dx[1] = (-jac[1 + 0 * 2] * fun[0] + jac[0 + 0 * 2] * fun[1]) / det;

                ref_[
                    0 + j * 2] -= dx[0];
                ref_[
                    1 + j * 2] -= dx[1];
            }
        }
    }
}