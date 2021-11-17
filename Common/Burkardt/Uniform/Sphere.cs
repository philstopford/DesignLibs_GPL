using System;
using Burkardt.RandomNS;
using Burkardt.Types;

namespace Burkardt.Uniform;

public static class Sphere
{
    public static double[] uniform_in_sphere01_map(int m, int n, ref typeMethods.r8vecNormalData data, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UNIFORM_IN_SPHERE01_MAP maps uniform points into the unit sphere.
        //
        //  Discussion:
        //
        //    The sphere has center 0 and radius 1.
        //
        //    We first generate a point ON the sphere, and then distribute it
        //    IN the sphere.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 August 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Russell Cheng,
        //    Random Variate Generation,
        //    in Handbook of Simulation,
        //    edited by Jerry Banks,
        //    Wiley, 1998, pages 168.
        //
        //    Reuven Rubinstein,
        //    Monte Carlo Optimization, Simulation, and Sensitivity
        //    of Queueing Networks,
        //    Krieger, 1992,
        //    ISBN: 0894647644,
        //    LC: QA298.R79.
        //
        //  Parameters:
        //
        //    Input, int M, the dimension of the space.
        //
        //    Input, int N, the number of points.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double X[M*N], the points.
        //
    {
        int j;

        double exponent = 1.0 / m;

        double[] x = new double[m * n];

        for (j = 0; j < n; j++)
        {
            //
            //  Fill a vector with normally distributed values.
            //
            double[] v = typeMethods.r8vec_normal_01_new(m, ref data, ref seed);
            //
            //  Compute the length of the vector.
            //
            double norm = typeMethods.r8vec_norm(m, v);
            //
            //  Normalize the vector.
            //
            int i;
            for (i = 0; i < m; i++)
            {
                v[i] /= norm;
            }

            //
            //  Now compute a value to map the point ON the sphere INTO the sphere.
            //
            double r = UniformRNG.r8_uniform_01(ref seed);
            r = Math.Pow(r, exponent);

            for (i = 0; i < m; i++)
            {
                x[i + j * m] = r * v[i];
            }

        }

        return x;
    }

    public static double[] uniform_on_hemisphere01_phong(int n, int m, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UNIFORM_ON_HEMISPHERE01_PHONG maps uniform points onto the unit hemisphere.
        //
        //  Discussion:
        //
        //    The sphere has center 0 and radius 1.
        //
        //    The Phong density is used, with exponent M:
        //
        //    rho ( theta, phi; m ) = ( m + 1 ) * cos ( phi )**M / ( 2 * Math.PI )
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
        //  Reference:
        //
        //    Peter Shirley,
        //    Nonuniform Random Point Sets Via Warping,
        //    Graphics Gems, Volume III,
        //    edited by David Kirk,
        //    AP Professional, 1992,
        //    ISBN: 0122861663,
        //    LC: T385.G6973.
        //
        //  Parameters:
        //
        //    Input, int N, the number of points.
        //
        //    Input, int M, the Phong exponent.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double UNIFORM_ON_HEMISPHERE01_PHONG[3*N], the points.
        //
    {
        const int DIM_NUM = 3;

        int j;

        double[] x = new double[DIM_NUM * n];

        double power = 1.0 / (m + 1);

        for (j = 0; j < n; j++)
        {
            double phi = UniformRNG.r8_uniform_01(ref seed);

            phi = Math.Acos(Math.Pow(1.0 - phi, power));

            double theta = UniformRNG.r8_uniform_01(ref seed);

            theta = 2.0 * Math.PI * theta;

            x[0 + j * DIM_NUM] = Math.Cos(theta) * Math.Sin(phi);
            x[1 + j * DIM_NUM] = Math.Sin(theta) * Math.Sin(phi);
            x[2 + j * DIM_NUM] = Math.Cos(phi);
        }

        return x;
    }

    public static double[] uniform_on_sphere01_map(int dim_num, int n, ref typeMethods.r8vecNormalData data, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UNIFORM_ON_SPHERE01_MAP maps uniform points onto the unit sphere.
        //
        //  Discussion:
        //
        //    The sphere has center 0 and radius 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 November 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Russell Cheng,
        //    Random Variate Generation,
        //    in Handbook of Simulation,
        //    edited by Jerry Banks,
        //    Wiley, 1998, pages 168.
        //
        //    Reuven Rubinstein,
        //    Monte Carlo Optimization, Simulation, and Sensitivity
        //    of Queueing Networks,
        //    Krieger, 1992,
        //    ISBN: 0894647644,
        //    LC: QA298.R79.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the dimension of the space.
        //
        //    Input, int N, the number of points.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double UNIFORM_ON_SPHERE01_MAP[DIM_NUM*N], the points.
        //
    {
        int j;

        double[] x = typeMethods.r8mat_normal_01_new(dim_num, n, ref data, ref seed);

        for (j = 0; j < n; j++)
        {
            //
            //  Compute the length of the vector.
            //
            double norm = 0.0;
            int i;
            for (i = 0; i < dim_num; i++)
            {
                norm += x[i + j * dim_num] * x[i + j * dim_num];
            }

            norm = Math.Sqrt(norm);
            //
            //  Normalize the vector.
            //
            for (i = 0; i < dim_num; i++)
            {
                x[i + j * dim_num] /= norm;
            }

        }

        return x;
    }

    public static double[] uniform_on_sphere01_patch_tp(int n, double phi1, double phi2,
            double theta1, double theta2, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UNIFORM_ON_SPHERE01_PATCH_TP maps uniform points to a spherical TP patch.
        //
        //  Discussion:
        //
        //    The sphere has center 0 and radius 1.
        //
        //    A sphere TP patch on the surface of the unit sphere contains those
        //    points with radius R = 1 and angles (THETA,PHI) such that
        //
        //      0.0 <= THETA1 <= THETA <= THETA2 <= 2 * PI
        //      0.0 <= PHI1   <= PHI   <= PHI2   <=     PI
        //
        //    mapped to
        //
        //      X = cos ( THETA ) * sin ( PHI )
        //      Y = sin ( THETA ) * sin ( PHI )
        //      Z =                 cos ( PHI )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 August 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Peter Shirley,
        //    Nonuniform Random Point Sets Via Warping,
        //    Graphics Gems, Volume III,
        //    edited by David Kirk,
        //    AP Professional, 1992,
        //    ISBN: 0122861663,
        //    LC: T385.G6973.
        //
        //  Parameters:
        //
        //    Input, int N, the number of points.
        //
        //    Input, double PHI1, PHI2, the latitudinal angle range.
        //
        //    Input, double THETA1, THETA2, the longitudinal angle range.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double UNIFORM_ON_SPHERE01_PATCH_TP[2*N], the TP points.
        //
    {
        int j;

        double[] tp = UniformRNG.r8mat_uniform_01_new(2, n, ref seed);

        for (j = 0; j < n; j++)
        {
            tp[0 + j * 2] = (1.0 - tp[0 + j * 2]) * theta1
                            + tp[0 + j * 2] * theta2;

            tp[1 + j * 2] = Math.Acos((1.0 - tp[1 + j * 2]) * Math.Cos(phi1)
                                      + tp[1 + j * 2] * Math.Cos(phi2));
        }

        return tp;
    }

    public static double[] uniform_on_sphere01_patch_xyz(int n, double phi1, double phi2,
            double theta1, double theta2, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UNIFORM_ON_SPHERE01_PATCH_XYZ maps uniform points to a spherical XYZ patch.
        //
        //  Discussion:
        //
        //    The sphere has center 0 and radius 1.
        //
        //    A sphere XYZ patch on the surface of the unit sphere contains those
        //    points with radius R = 1 and angles (THETA,PHI) such that
        //
        //      0.0 <= THETA1 <= THETA <= THETA2 <= 2 * PI
        //      0.0 <= PHI1   <= PHI   <= PHI2   <=     PI
        //
        //    mapped to
        //
        //      X = cos ( THETA ) * sin ( PHI )
        //      Y = sin ( THETA ) * sin ( PHI )
        //      Z =                 cos ( PHI )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Peter Shirley,
        //    Nonuniform Random Point Sets Via Warping,
        //    Graphics Gems, Volume III,
        //    edited by David Kirk,
        //    AP Professional, 1992,
        //    ISBN: 0122861663,
        //    LC: T385.G6973.
        //
        //  Parameters:
        //
        //    Input, int N, the number of points.
        //
        //    Input, double PHI1, PHI2, the latitudinal angle range.
        //
        //    Input, double THETA1, THETA2, the longitudinal angle range.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double UNIFORM_ON_SPHERE01_PATCH_XYZ[3*N], the points.
        //
    {
        const int DIM_NUM = 3;

        int j;

        double[] x = new double[DIM_NUM * n];

        for (j = 0; j < n; j++)
        {
            double phi = UniformRNG.r8_uniform_01(ref seed);

            phi = Math.Acos((1.0 - phi) * Math.Cos(phi1)
                            + phi * Math.Cos(phi2));

            double theta = UniformRNG.r8_uniform_01(ref seed);

            theta = (1.0 - theta) * theta1
                    + theta * theta2;

            x[0 + j * DIM_NUM] = Math.Cos(theta) * Math.Sin(phi);
            x[1 + j * DIM_NUM] = Math.Sin(theta) * Math.Sin(phi);
            x[2 + j * DIM_NUM] = Math.Cos(phi);
        }

        return x;
    }

    public static double[] uniform_on_sphere01_triangle_xyz(int n, double[] v1, double[] v2,
            double[] v3, ref int seed )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UNIFORM_ON_SPHERE01_TRIANGLE_XYZ: sample spherical triangle, XYZ coordinates.
        //
        //  Discussion:
        //
        //    The sphere has center 0 and radius 1.
        //
        //    A spherical triangle on the surface of the unit sphere contains those
        //    points with radius R = 1, bounded by the vertices V1, V2, V3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 August 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    James Arvo,
        //    Stratified sampling of spherical triangles,
        //    Computer Graphics Proceedings, Annual Conference Series,
        //    ACM SIGGRAPH '95, pages 437-438, 1995.
        //
        //  Parameters:
        //
        //    Input, int N, the number of points.
        //
        //    Input, double V1[3], V2[3], V3[3], the XYZ coordinates of
        //    the vertices of the spherical triangle.
        //
        //    Input/output, int &SEED, a seed for the random
        //    number generator.
        //
        //    Output, double UNIFORM_ON_SPHERE01_TRIANGLE_XYZ[3*N], the XYZ
        //    coordinates of the sample points.
        //
    {
        double a = 0;
        double alpha = 0;
        double b = 0;
        double beta = 0;
        double c = 0;
        double gamma = 0;
        int j;
        //
        //  Compute the sides, angles, and area of the spherical triangle;
        //  for now, we assume R = 1.
        //
        double r = 1.0;

        Stri.stri_vertices_to_sides(r, v1, v2, v3, ref a, ref b, ref c);

        Stri.stri_sides_to_angles(r, a, b, c, ref alpha, ref beta, ref gamma);

        double area = Stri.stri_angles_to_area(r, alpha, beta, gamma);

        double[] x = new double[3 * n];
        double[] v31 = new double[3];
        double[] v4 = new double[3];
        double[] v42 = new double[3];

        for (j = 0; j < n; j++)
        {
            //
            //  Select the new area.
            //
            double xsi1 = UniformRNG.r8_uniform_01(ref seed);
            double area_hat = xsi1 * area;
            //
            //  Compute the sine and cosine of the angle phi.
            //
            double s = Math.Sin(area_hat - alpha);
            double t = Math.Cos(area_hat - alpha);
            //
            //  Compute the pair that determines beta_hat.
            //
            double u = t - Math.Cos(alpha);
            double v = s + Math.Sin(alpha) * Math.Cos(c);
            //
            //  Q is the cosine of the new edge length b_hat.
            //
            double q = ((v * t - u * s) * Math.Cos(alpha) - v)
                       / ((v * s + u * t) * Math.Sin(alpha));
            //
            //  V31 = normalized ( V3 - ( V3 dot V1 ) * V1 )
            //
            double w = typeMethods.r8vec_dot_product(3, v3, v1);

            int i;
            for (i = 0; i < 3; i++)
            {
                v31[i] = v3[i] - w * v1[i];
            }

            double temp = typeMethods.r8vec_norm(3, v31);
            for (i = 0; i < 3; i++)
            {
                v31[i] /= temp;
            }

            //
            //  V4 is the third vertex of the subtriangle V1, V2, V4.
            //
            for (i = 0; i < 3; i++)
            {
                v4[i] = q * v1[i] + Math.Sqrt(1.0 - q * q) * v31[i];
            }

            //
            //  Select cos theta, which will sample along the edge from V2 to V4.
            //
            double xsi2 = UniformRNG.r8_uniform_01(ref seed);
            double z = 1.0 - xsi2 * (1.0 - typeMethods.r8vec_dot_product(3, v4, v2));
            //
            //  V42 = normalized ( V4 - ( V4 dot V2 ) * V2 )
            //
            w = typeMethods.r8vec_dot_product(3, v4, v2);
            for (i = 0; i < 3; i++)
            {
                v42[i] = v4[i] - w * v2[i];
            }

            temp = typeMethods.r8vec_norm(3, v42);
            for (i = 0; i < 3; i++)
            {
                v42[i] /= temp;
            }

            //
            //  Construct the point.
            //
            for (i = 0; i < 3; i++)
            {
                x[i + j * 3] = z * v2[i] + Math.Sqrt(1.0 - z * z) * v42[i];
            }
        }

        return x;
    }

}