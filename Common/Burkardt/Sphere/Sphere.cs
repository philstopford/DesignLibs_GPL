using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.SphereNS
{
    public static class Sphere
    {
        public static double[] sphere_spiralpoints ( double r, double[] pc, int n )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_SPIRALPOINTS produces spiral points on an implicit sphere.
        //
        //  Discussion:
        //
        //    The points should be arranged on the sphere in a pleasing design.
        //
        //    The implicit form of a sphere in 3D is:
        //
        //        pow ( P[0] - PC[0], 2 ) 
        //      + pow ( P[1] - PC[1], 2 ) 
        //      + pow ( P[2] - PC[2], 2 ) = pow ( R, 2 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Edward Saff, Arno Kuijlaars,
        //    Distributing Many Points on a Sphere,
        //    The Mathematical Intelligencer,
        //    Volume 19, Number 1, 1997, pages 5-11.
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the sphere.
        //
        //    Input, double PC[3], the coordinates of the center of the sphere.
        //
        //    Input, int N, the number of points to create.
        //
        //    Output, double SPHERE_SPIRALPOINTS[3*N], the coordinates of the grid points.
        //
        {
            double cosphi = 0.0;
            int i;
            double[] p;
            
            double sinphi = 0.0;
            double theta = 0.0;

            p = new double[3*n];

            for ( i = 0; i < n; i++ )
            {
                cosphi = ( ( double ) ( n - i - 1 ) * ( -1.0 ) 
                           + ( double ) (     i     ) * (  1.0 ) ) 
                         / ( double ) ( n     - 1 );

                sinphi = Math.Sqrt ( 1.0 - cosphi * cosphi );

                if ( i == 0 || i == n - 1 )
                {
                    theta = 0.0;
                }
                else
                {
                    theta = theta + 3.6 / ( sinphi * Math.Sqrt ( ( double ) n ) );
                    theta = typeMethods.r8_modp ( theta, 2.0 * Math.PI );
                }
                p[0+i*3] = pc[0] + r * sinphi * Math.Cos ( theta );
                p[1+i*3] = pc[1] + r * sinphi * Math.Sin ( theta );
                p[2+i*3] = pc[2] + r * cosphi;
            }

            return p;
        }
        
        public static double[] sphere_unit_sample ( int n, ref int seed )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_UNIT_SAMPLE picks a random point on the unit sphere.
            //
            //  Discussion:
            //
            //    The unit sphere in 3D satisfies the equation:
            //
            //      X * X + Y * Y + Z * Z = 1
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    29 August 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of points to generate.
            //
            //    Input/output, int *SEED, a seed for the random number generator.
            //
            //    Output, double SPHERE_UNIT_SAMPLE[3*N], the sample point.
            //
        {
            int j;
            double phi;
            
            double theta;
            double vdot;
            double[] x;

            x = new double[3*n];

            for ( j = 0; j < n; j++ )
            {
                //
                //  Pick a uniformly random VDOT, which must be between -1 and 1.
                //  This represents the dot product of the random vector with the Z unit vector.
                //
                //   this works because the surface area of the sphere between
                //  Z and Z + dZ is independent of Z.  So choosing Z uniformly chooses
                //  a patch of area uniformly.
                //
                vdot = 2.0 * UniformRNG.r8_uniform_01 ( ref seed ) - 1.0;

                phi = Helpers.arc_cosine ( vdot );
                //
                //  Pick a uniformly random rotation between 0 and 2 Pi around the
                //  axis of the Z vector.
                //
                theta = 2.0 * Math.PI * UniformRNG.r8_uniform_01 ( ref seed );

                x[0+j*3] = Math.Cos ( theta ) * Math.Sin ( phi );
                x[1+j*3] = Math.Sin ( theta ) * Math.Sin ( phi );
                x[2+j*3] = Math.Cos ( phi );
            }

            return x;
        }
        
        public static double sphere01_monomial_quadrature ( int[] expon, int point_num, double[] xyz, 
                double[] w )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE01_MONOMIAL_QUADRATURE applies quadrature to a monomial in a sphere.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    05 July 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the spatial dimension.
            //
            //    Input, int EXPON[DIM_NUM], the exponents.
            //
            //    Input, int POINT_NUM, the number of points in the rule.
            //
            //    Input, double XYZ[DIM_NUM*POINT_NUM], the quadrature points.
            //
            //    Input, double W[POINT_NUM], the quadrature weights.
            //
            //    Output, double SPHERE01_MONOMIAL_QUADRATURE, the quadrature error.
            //
        {
            double exact;
            double quad;
            double quad_error;
            double[] value;
            //
            //  Get the exact value of the integral.
            //
            exact = Integrals.sphere01_monomial_integral ( expon );
            //
            //  Evaluate the monomial at the quadrature points.
            //
            value = Burkardt.MonomialNS.Monomial.monomial_value ( 3, point_num, expon, xyz );
            //
            //  Compute the weighted sum.
            //
            quad = typeMethods.r8vec_dot ( point_num, w, value );
            //
            //  Error:
            //
            quad_error = Math.Abs ( quad - exact );

            return quad_error;
        }

    }
}