using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.TetrahedronNS
{
    public class Tetrahedron
    {
        public static void reference_to_physical_t4 ( double[] t, int n, double[] ref_, ref double[] phy )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    REFERENCE_TO_PHYSICAL_T4 maps T4 reference points to physical points.
        //
        //  Discussion:
        //
        //    Given the vertices of an order 4 physical tetrahedron and a point
        //    (R,S,T) in the reference tetrahedron, the routine computes the value
        //    of the corresponding image point (X,Y,Z) in physical space.
        //
        //    This routine will also be correct for an order 10 tetrahedron,
        //    if the mapping between reference and physical space
        //    is linear.  This implies, in particular, that the sides of the
        //    image tetrahedron are straight, the faces are flat, and 
        //    the "midside" nodes in the physical tetrahedron are
        //    halfway along the edges of the physical tetrahedron.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double T[3*4], the coordinates of the vertices.
        //    The vertices are assumed to be the images of (0,0,0), (1,0,0),
        //    (0,1,0) and (0,0,1) respectively.
        //
        //    Input, int N, the number of objects to transform.
        //
        //    Input, double REF[3*N], points in the reference tetrahedron.
        //
        //    Output, double PHY[3*N], corresponding points in the
        //    physical tetrahedron.
        //
        {
            int i;
            int j;

            for ( i = 0; i < 3; i++ )
            {
                for ( j = 0; j < n; j++ )
                {
                    phy[i+j*3] = t[i+0*3] * ( 1.0 - ref_[0+j*3] - ref_[1+j*3] - ref_[2+j*3] ) 
                    + t[i+1*3] *       + ref_[0+j*3] 
                    + t[i+2*3] *                    + ref_[1+j*3]                  
                    + t[i+3*3] *                                 + ref_[2+j*3];
                }
            }
        }
        public static double[] reference_to_physical_tet4(double[] t, int n, double[] ref_)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    REFERENCE_TO_PHYSICAL_TET4 maps TET4 reference points to physical points.
            //
            //  Discussion:
            //
            //    Given the vertices of an order 4 physical tetrahedron and a point 
            //    (R,S,T) in the reference tetrahedron, the routine computes the value 
            //    of the corresponding image point (X,Y,Z) in physical space.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    10 August 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double T[3*4], the coordinates of the vertices.  
            //    The vertices are assumed to be the images of (1,0,0), (0,1,0),
            //    (0,0,1) and (0,0,0) respectively.
            //
            //    Input, int N, the number of points to transform.
            //
            //    Input, double REF[3*N], points in the reference element.
            //
            //    Output, double REFERENCE_TO_PHYSICAL_TET4[3*N], corresponding points in the
            //    physical element.
            //
        {
            int i;
            int j;
            double[] phy;

            phy = new double[3 * n];

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < 3; i++)
                {
                    phy[i + j * 3] =
                        t[i + 0 * 3] * ref_[0 + j * 3]
                        + t[i + 1 * 3] * ref_[1 + j * 3]
                        + t[i + 2 * 3] * ref_[2 + j * 3]
                        + t[i + 3 * 3] * (1.0 - ref_[0 + j * 3] - ref_[1 + j * 3] - ref_[2 + j * 3]);
                }
            }

            return phy;
        }

        public static double[] tetrahedron_barycentric(double[] tetra, double[] p)

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

        public static void tetrahedron_order4_physical_to_reference(double[] tetra, int n,
        double[] phy, ref double[] ref_ )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TETRAHEDRON_ORDER4_PHYSICAL_TO_REFERENCE maps physical points to reference points.
        //
        //  Discussion:
        //
        //    Given the vertices of an order 4 physical tetrahedron and a point
        //    (X,Y,Z) in the physical tetrahedron, the routine computes the value
        //    of the corresponding image point (R,S,T) in reference space.
        //
        //    This routine may be appropriate for an order 10 tetrahedron,
        //    if the mapping between reference and physical space is linear.  
        //    This implies, in particular, that the edges of the image tetrahedron 
        //    are straight, the faces are flat, and the "midside" nodes in the
        //    physical tetrahedron are halfway along the sides of the physical 
        //    tetrahedron.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double TETRA[3*4], the coordinates of the vertices.  
        //    The vertices are assumed to be the images of
        //    (0,0,0), (1,0,0), (0,1,0) and (0,0,1) respectively.
        //
        //    Input, int N, the number of points to transform.
        //
        //    Input, double PHY[3*N], the coordinates of physical points
        //    to be transformed.
        //
        //    Output, double REF[3*N], the coordinates of the corresponding
        //    points in the reference space.
        //
        {
            double[] a = new double[3 * 3];
            double det;
            int i;
            int j;
            //
            //  Set up the matrix.
            //
            for (i = 0; i < 3; i++)
            {
                a[i + 0 * 3] = tetra[i + 1 * 3] - tetra[i + 0 * 3];
                a[i + 1 * 3] = tetra[i + 2 * 3] - tetra[i + 0 * 3];
                a[i + 2 * 3] = tetra[i + 3 * 3] - tetra[i + 0 * 3];
            }

            //
            //  Compute the determinant.
            //
            det = a[0 + 0 * 3] * (a[1 + 1 * 3] * a[2 + 2 * 3] - a[1 + 2 * 3] * a[2 + 1 * 3])
                  + a[0 + 1 * 3] * (a[1 + 2 * 3] * a[2 + 0 * 3] - a[1 + 0 * 3] * a[2 + 2 * 3])
                  + a[0 + 2 * 3] * (a[1 + 0 * 3] * a[2 + 1 * 3] - a[1 + 1 * 3] * a[2 + 0 * 3]);
            //
            //  If the determinant is zero, bail out.
            //
            if (det == 0.0)
            {
                for (j = 0; j < n; j++)
                {
                    for (i = 0; i < 3; i++)
                    {
                        ref_[i
                        +j * 3] = 0.0;
                    }
                }

                return;
            }

            //
            //  Compute the solution.
            //
            for (j = 0; j < n; j++)
            {
                ref_[
                0 + j * 3] = ((a[1 + 1 * 3] * a[2 + 2 * 3] - a[1 + 2 * 3] * a[2 + 1 * 3])
                              * (phy[0 + j * 3] - tetra[0 + 0 * 3])
                              - (a[0 + 1 * 3] * a[2 + 2 * 3] - a[0 + 2 * 3] * a[2 + 1 * 3])
                              * (phy[1 + j * 3] - tetra[1 + 0 * 3])
                              + (a[0 + 1 * 3] * a[1 + 2 * 3] - a[0 + 2 * 3] * a[1 + 1 * 3])
                              * (phy[2 + j * 3] - tetra[2 + 0 * 3])
                    ) / det;

                ref_[
                1 + j * 3] = (-(a[1 + 0 * 3] * a[2 + 2 * 3] - a[1 + 2 * 3] * a[2 + 0 * 3])
                              * (phy[0 + j * 3] - tetra[0 + 0 * 3])
                              + (a[0 + 0 * 3] * a[2 + 2 * 3] - a[0 + 2 * 3] * a[2 + 0 * 3])
                              * (phy[1 + j * 3] - tetra[1 + 0 * 3])
                              - (a[0 + 0 * 3] * a[1 + 2 * 3] - a[0 + 2 * 3] * a[1 + 0 * 3])
                              * (phy[2 + j * 3] - tetra[2 + 0 * 3])
                    ) / det;

                ref_[
                2 + j * 3] = ((a[1 + 0 * 3] * a[2 + 1 * 3] - a[1 + 1 * 3] * a[2 + 0 * 3])
                              * (phy[0 + j * 3] - tetra[0 + 0 * 3])
                              - (a[0 + 0 * 3] * a[2 + 1 * 3] - a[0 + 1 * 3] * a[2 + 0 * 3])
                              * (phy[1 + j * 3] - tetra[1 + 0 * 3])
                              + (a[0 + 0 * 3] * a[1 + 1 * 3] - a[0 + 1 * 3] * a[1 + 0 * 3])
                              * (phy[2 + j * 3] - tetra[2 + 0 * 3])
                    ) / det;
            }
        }

                public static void tetrahedron_reference_sample(int n, ref int seed, ref double[] p)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TETRAHEDRON_REFERENCE_SAMPLE samples points in the reference tetrahedron.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    08 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the  number of points to sample.
            //
            //    Input/output, int *SEED, a seed for the random number generator.
            //
            //    Output, double P[3*N], random points in the tetrahedron.
            //
        {
            double alpha;
            double beta;
            double gamma;
            int j;
            double r;

            for (j = 0; j < n; j++)
            {
                r = UniformRNG.r8_uniform_01(ref seed);
                //
                //  Interpret R as a percentage of the tetrahedron's volume.
                //
                //  Imagine a plane, parallel to face 1, so that the volume between
                //  vertex 1 and the plane is R percent of the full tetrahedron volume.
                //
                //  The plane will intersect sides 12, 13, and 14 at a fraction
                //  ALPHA = R^1/3 of the distance from vertex 1 to vertices 2, 3, and 4.
                //
                alpha = Math.Pow(r, 1.0 / 3.0);
                //
                //  Determine the coordinates of the points on sides 12, 13 and 14 intersected
                //  by the plane, which form a triangle TR.
                //
                //  Now choose, uniformly at random, a point in this triangle.
                //
                r = UniformRNG.r8_uniform_01(ref seed);
                //
                //  Interpret R as a percentage of the triangle's area.
                //
                //  Imagine a line L, parallel to side 1, so that the area between
                //  vertex 1 and line L is R percent of the full triangle's area.
                //
                //  The line L will intersect sides 2 and 3 at a fraction
                //  ALPHA = SQRT ( R ) of the distance from vertex 1 to vertices 2 and 3.
                //
                beta = Math.Sqrt(r);
                //
                //  Determine the coordinates of the points on sides 2 and 3 intersected
                //  by line L.
                //
                //  Now choose, uniformly at random, a point on the line L.
                //
                gamma = UniformRNG.r8_uniform_01(ref seed);

                p[0 + j * 3] = alpha * (1.0 - beta) * gamma;
                p[1 + j * 3] = alpha * beta * (1.0 - gamma);
                p[2 + j * 3] = alpha * beta * gamma;
            }
        }

        public static void tetrahedron_sample(double[] tetra, int n, ref int seed, ref double[] p)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TETRAHEDRON_SAMPLE returns random points in a tetrahedron.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    06 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double TETRA[3*4], the coordinates of the vertices.
            //
            //    Input/output, int *SEED, a seed for the random number generator.
            //
            //    Output, double P[3*N], random points in the tetrahedron.
            //
        {
            int DIM_NUM = 3;

            double alpha;
            double beta;
            double gamma;
            int i;
            int j;
            int k;
            double[] p12;
            double[] p13;
            double r;
            double[] t;

            p12 = new double[DIM_NUM];
            p13 = new double[DIM_NUM];
            t = new double[DIM_NUM * 3];

            for (k = 0; k < n; k++)
            {
                r = UniformRNG.r8_uniform_01(ref seed);
                //
                //  Interpret R as a percentage of the tetrahedron's volume.
                //
                //  Imagine a plane, parallel to face 1, so that the volume between
                //  vertex 1 and the plane is R percent of the full tetrahedron volume.
                //
                //  The plane will intersect sides 12, 13, and 14 at a fraction
                //  ALPHA = R^1/3 of the distance from vertex 1 to vertices 2, 3, and 4.
                //
                alpha = Math.Pow(r, 1.0 / 3.0);
                //
                //  Determine the coordinates of the points on sides 12, 13 and 14 intersected
                //  by the plane, which form a triangle TR.
                //
                for (i = 0; i < DIM_NUM; i++)
                {
                    for (j = 0; j < 3; j++)
                    {
                        t[i + j * 3] = (1.0 - alpha) * tetra[i + 0 * 3]
                                       + alpha * tetra[i + (j + 1) * 3];
                    }
                }

                //
                //  Now choose, uniformly at random, a point in this triangle.
                //
                r = UniformRNG.r8_uniform_01(ref seed);
                //
                //  Interpret R as a percentage of the triangle's area.
                //
                //  Imagine a line L, parallel to side 1, so that the area between
                //  vertex 1 and line L is R percent of the full triangle's area.
                //
                //  The line L will intersect sides 2 and 3 at a fraction
                //  ALPHA = SQRT ( R ) of the distance from vertex 1 to vertices 2 and 3.
                //
                beta = Math.Sqrt(r);
                //
                //  Determine the coordinates of the points on sides 2 and 3 intersected
                //  by line L.
                //
                for (i = 0; i < DIM_NUM; i++)
                {
                    p12[i] = (1.0 - beta) * t[i + 0 * 3]
                             + beta * t[i + 1 * 3];

                    p13[i] = (1.0 - beta) * t[i + 0 * 3]
                             + beta * t[i + 2 * 3];
                }

                //
                //  Now choose, uniformly at random, a point on the line L.
                //
                gamma = UniformRNG.r8_uniform_01(ref seed);

                for (i = 0; i < DIM_NUM; i++)
                {
                    p[i + k * 3] = gamma * p12[i] + (1.0 - gamma) * p13[i];
                }
            }
        }

        public static void tetrahedron_order4_reference_to_physical(double[] tetra, int n,
                double[] ref_, ref double[] phy)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TETRAHEDRON_ORDER4_REFERENCE_TO_PHYSICAL maps reference points to physical points.
            //
            //  Discussion:
            //
            //    Given the vertices of an order 4 physical tetrahedron and a point
            //    (R,S,T) in the reference triangle, the routine computes the value
            //    of the corresponding image point (X,Y,Z) in physical space.
            //
            //    This routine will also be correct for an order 10 tetrahedron,
            //    if the mapping between reference and physical space
            //    is linear.  This implies, in particular, that the sides of the
            //    image tetrahedron are straight, the faces are flat, and 
            //    the "midside" nodes in the physical tetrahedron are
            //    halfway along the edges of the physical tetrahedron.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    06 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double TETRA[3*4], the coordinates of the vertices.
            //    The vertices are assumed to be the images of (0,0,0), (1,0,0),
            //    (0,1,0) and (0,0,1) respectively.
            //
            //    Input, int N, the number of points to transform.
            //
            //    Input, double REF[3*N], points in the reference tetrahedron
            //
            //    Output, double PHY[3*N], corresponding points in the
            //    physical tetrahedron.
            //
        {
            int i;
            int j;

            for (i = 0; i < 3; i++)
            {
                for (j = 0; j < n; j++)
                {
                    phy[i + j * 3] = tetra[i + 0 * 3] * (1.0 - ref_[0 + j * 3] - ref_[1 + j * 3] - ref_[2 + j * 3])
                                     + tetra[i + 1 * 3] * +ref_[0 + j * 3]
                                     + tetra[i + 2 * 3] * +ref_[1 + j * 3]
                                     + tetra[i + 3 * 3] * +ref_[2 + j * 3];
                }
            }
        }

        public static double tetrahedron_volume(double[] tetra)

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