using System;
using Burkardt.Types;

namespace Burkardt
{
    public class Tetrahedron
    {
        public static double[] tetrahedron_barycentric(double[] tetra, double[] p )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TETRAHEDRON_BARYCENTRIC returns the barycentric coordinates of a point.
        //
        //  Discussion:
        //
        //    The barycentric coordinates of a point P with respect to
        //    a tetrahedron are a set of four values C(1:4), each associated
        //    with a vertex of the tetrahedron.  The values must sum to 1.
        //    If all the values are between 0 and 1, the point is contained
        //    within the tetrahedron.
        //
        //    The barycentric coordinate of point X related to vertex A can be
        //    interpreted as the ratio of the volume of the tetrahedron with 
        //    vertex A replaced by vertex X to the volume of the original 
        //    tetrahedron.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double TETRA[3*4], the vertices of the tetrahedron.
        //
        //    Input, double P[3], the point to be checked.
        //
        //    Output, double C[4], the barycentric coordinates of the point with
        //    respect to the tetrahedron.
        //
        {
            int N = 3;
            int RHS_NUM = 1;

            double[] a = new double[N * (N + RHS_NUM)];
            double[] c;
            int info;
            //
            //  Set up the linear system
            //
            //    ( X2-X1  X3-X1  X4-X1 ) C1    X - X1
            //    ( Y2-Y1  Y3-Y1  Y4-Y1 ) C2  = Y - Y1
            //    ( Z2-Z1  Z3-Z1  Z4-Z1 ) C3    Z - Z1
            //
            //  which is satisfied by the barycentric coordinates.
            //

            a[0 + 0 * N] = tetra[0 + 1 * 3] - tetra[0 + 0 * 3];
            a[1 + 0 * N] = tetra[1 + 1 * 3] - tetra[1 + 0 * 3];
            a[2 + 0 * N] = tetra[2 + 1 * 3] - tetra[2 + 0 * 3];

            a[0 + 1 * N] = tetra[0 + 2 * 3] - tetra[0 + 0 * 3];
            a[1 + 1 * N] = tetra[1 + 2 * 3] - tetra[1 + 0 * 3];
            a[2 + 1 * N] = tetra[2 + 2 * 3] - tetra[2 + 0 * 3];

            a[0 + 2 * N] = tetra[0 + 3 * 3] - tetra[0 + 0 * 3];
            a[1 + 2 * N] = tetra[1 + 3 * 3] - tetra[1 + 0 * 3];
            a[2 + 2 * N] = tetra[2 + 3 * 3] - tetra[2 + 0 * 3];

            a[0 + 3 * N] = p[0] - tetra[0 + 0 * 3];
            a[1 + 3 * N] = p[1] - tetra[1 + 0 * 3];
            a[2 + 3 * N] = p[2] - tetra[2 + 0 * 3];
            //
            //  Solve the linear system.
            //
            info = typeMethods.r8mat_solve(N, RHS_NUM, ref a);

            if (info != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("TETRAHEDRON_BARYCENTRIC - Fatal error!");
                Console.WriteLine("  The linear system is singular.");
                Console.WriteLine("  The input data does not form a proper tetrahedron.");
                return null;
            }

            c = new double[4];

            c[1] = a[0 + 3 * N];
            c[2] = a[1 + 3 * N];
            c[3] = a[2 + 3 * N];

            c[0] = 1.0 - c[1] - c[2] - c[3];

            return c;
        }

        public static double tetrahedron_volume(double[] tetra )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TETRAHEDRON_VOLUME computes the volume of a tetrahedron in 3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double TETRA[3*4], the coordinates of the vertices.
        //
        //    Output, double TETRAHEDRON_VOLUME, the volume of the tetrahedron.
        //
        {
            double[] a = new double[4 * 4];
            int i;
            int j;
            double volume;

            for (i = 0; i < 3; i++)
            {
                for (j = 0; j < 4; j++)
                {
                    a[i + j * 4] = tetra[i + j * 3];
                }
            }

            i = 3;
            for (j = 0; j < 4; j++)
            {
                a[i + j * 4] = 1.0;
            }

            volume = Math.Abs(typeMethods.r8mat_det_4d(a)) / 6.0;

            return volume;
        }
    }
}