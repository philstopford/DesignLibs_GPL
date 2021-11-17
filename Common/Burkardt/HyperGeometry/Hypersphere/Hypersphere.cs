using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.HyperGeometry.HypersphereNS;

public static class Hypersphere
{
    public static void cartesian_to_hypersphere(int m, int n, double[] c, double[] x,
            ref double[] r, ref double[] theta )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CARTESIAN_TO_HYPERSPHERE: Cartesian to hypersphere coordinate transform.
        //
        //  Discussion:
        //
        //    We allow the trivial case M = 1; in that case alone, the value R
        //    must be assumed to have a sign.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 December 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the spatial dimension.
        //    1 <= M.
        //
        //    Input, int N, the number of points to transform.
        //
        //    Input, double C[M], the center of the hypersphere.
        //
        //    Input, double X[M*N], the Cartesian coordinates of the points.
        //
        //    Output, double R[N], the radius of the points on the 
        //    hypersphere.  Except for the trivial case M = 1, R is assumed nonnegative.
        //
        //    Output, double THETA[(M-1)*N], the coordinate angles of the 
        //    points, measured in radians.
        //
    {
        int i;
        int i1;
        int j;
        double t;
        double[] top;
        double[] x2;
        switch (m)
        {
            //
            //  Handle special case of M = 1.
            //
            case 1:
            {
                for (j = 0; j < n; j++)
                {
                    r[j] = x[0 + j * m] - c[0];
                }

                return;
            }
        }

        //
        //  Subtract the center.
        //
        x2 = new double[m * n];
        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m; i++)
            {
                x2[i + j * m] = x[i + j * m] - c[i];
            }
        }

        //
        //  Compute R.
        //
        for (j = 0; j < n; j++)
        {
            t = 0.0;
            for (i = 0; i < m; i++)
            {
                t += Math.Pow(x2[i + m * j], 2);
            }

            r[j] = Math.Sqrt(t);
        }

        //
        //  Compute M-2 components of THETA.
        //
        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m - 1; i++)
            {
                theta[i + j * (m - 1)] = 0.0;
            }
        }

        for (i = 1; i < m - 1; i++)
        {
            for (i1 = 0; i1 <= i - 1; i1++)
            {
                for (j = 0; j < n; j++)
                {
                    theta[i1 + j * (m - 1)] += Math.Pow(x2[i + j * m], 2);
                }
            }
        }

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m - 2; i++)
            {
                theta[i + j * (m - 1)] += Math.Pow(x2[m - 1 + j * m], 2);
            }
        }

        for (i = 0; i < m - 2; i++)
        {
            for (j = 0; j < n; j++)
            {
                theta[i + j * (m - 1)] = Math.Atan2(Math.Sqrt(theta[i + j * (m - 1)]), x2[i + j * m]);
            }
        }

        //
        //  Compute last component of THETA.
        //
        top = new double[n];

        for (j = 0; j < n; j++)
        {
            top[j] = Math.Sqrt(Math.Pow(x2[m - 1 + j * m], 2) + Math.Pow(x2[m - 2 + j * m], 2)) + x2[m - 2 + j * m];
        }

        for (j = 0; j < n; j++)
        {
            theta[m - 2 + j * (m - 1)] = 2.0 * Math.Atan2(x2[m - 1 + j * m], top[j]);
        }
    }

    public static double hypersphere_01_area(int m)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HYPERSPHERE_01_AREA computes the surface area of a unit hypersphere.
        //
        //  Discussion:
        //
        //    The unit hypersphere satisfies the equation:
        //
        //      Sum ( 1 <= I <= M ) X(I) * X(I) = 1
        //
        //    M   Area
        //
        //     2    2        * PI
        //     3    4        * PI
        //     4  ( 2 /   1) * PI^2
        //     5  ( 8 /   3) * PI^2
        //     6  ( 1 /   1) * PI^3
        //     7  (16 /  15) * PI^3
        //     8  ( 1 /   3) * PI^4
        //     9  (32 / 105) * PI^4
        //    10  ( 1 /  12) * PI^5
        //
        //    For the unit hypersphere, Area(M) = M * Volume(M)
        //
        //    Sphere_Unit_Area ( M ) = 2 * PI^(M/2) / Gamma ( M / 2 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 December 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the dimension of the space.
        //
        //    Output, double HYPERSPHERE_01_AREA, the area.
        //
    {
        double area;
        int i;
        int m2;
            

        switch (m % 2)
        {
            case 0:
            {
                m2 = m / 2;
                area = 2.0 * Math.Pow(Math.PI, m2);
                for (i = 1; i <= m2 - 1; i++)
                {
                    area /= i;
                }

                break;
            }
            default:
            {
                m2 = (m - 1) / 2;
                area = Math.Pow(2.0, m) * Math.Pow(Math.PI, m2);
                for (i = m2 + 1; i <= 2 * m2; i++)
                {
                    area /= i;
                }

                break;
            }
        }

        return area;
    }

    public static void hypersphere_01_area_values(ref int n_data, ref int m, ref double area )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HYPERSPHERE_01_AREA_VALUES returns some areas of the unit hypersphere.
        //
        //  Discussion:
        //
        //    The unit hypersphere satisfies the equation:
        //
        //      Sum ( 1 <= I <= M ) X(I) * X(I) = 1
        //
        //     M         Area
        //
        //     2    2        * PI
        //     3  ( 4 /    ) * PI
        //     4  ( 2 /   1) * PI^2
        //     5  ( 8 /   3) * PI^2
        //     6  ( 1 /   1) * PI^3
        //     7  (16 /  15) * PI^3
        //     8  ( 1 /   3) * PI^4
        //     9  (32 / 105) * PI^4
        //    10  ( 1 /  12) * PI^5
        //
        //    For the unit hypersphere, Area(M) = M * Volume(M)
        //
        //    Sphere_Unit_Area ( M ) = 2 * PI^(M/2) / Gamma ( M / 2 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 December 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input/output, int &N_DATA.
        //    On input, if N_DATA is 0, the first test data is returned, and
        //    N_DATA is set to the index of the test data.  On each subsequent
        //    call, N_DATA is incremented and that test data is returned.  When
        //    there is no more test data, N_DATA is set to 0.
        //
        //    Output, int &M, the spatial dimension.
        //
        //    Output, double &AREA, the area.
        //
    {
        int NMAX = 20;

        double[] area_vec =  {
                0.2000000000000000E+01,
                0.6283185307179586E+01,
                0.1256637061435917E+02,
                0.1973920880217872E+02,
                0.2631894506957162E+02,
                0.3100627668029982E+02,
                0.3307336179231981E+02,
                0.3246969701133415E+02,
                0.2968658012464836E+02,
                0.2550164039877345E+02,
                0.2072514267328890E+02,
                0.1602315322625507E+02,
                0.1183817381218268E+02,
                0.8389703410491089E+01,
                0.5721649212349567E+01,
                0.3765290085742291E+01,
                0.2396678817591364E+01,
                0.1478625959000308E+01,
                0.8858104195716824E+00,
                0.5161378278002812E+00
            }
            ;
        int[] m_vec =  {
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8,
                9,
                10,
                11,
                12,
                13,
                14,
                15,
                16,
                17,
                18,
                19,
                20
            }
            ;

        n_data = n_data switch
        {
            < 0 => 0,
            _ => n_data
        };

        if (NMAX <= n_data)
        {
            n_data = 0;
            m = 0;
            area = 0.0;
        }
        else
        {
            m = m_vec[n_data];
            area = area_vec[n_data];
            n_data += 1;
        }
    }

    public static double[] hypersphere_01_interior_uniform(int m, int n, ref typeMethods.r8vecNormalData data, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HYPERSPHERE_01_INTERIOR_UNIFORM: uniform points inside unit hypersphere.
        //
        //  Discussion:
        //
        //    The hypersphere has center 0 and radius 1.
        //
        //    This routine is valid for any spatial dimension.
        //
        //    We first generate a point ON the hypersphere, and then distribute it
        //    IN the hypersphere.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 December 2013
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
        //    Input/output, int &SEED, a seed for the random
        //    number generator.
        //
        //    Output, double HYPERSPHERE_01_INTERIOR_UNIFORM[M*N], the points.
        //
    {
        double exponent;
        int i;
        int j;
        double norm;
        double r;
        double[] v;
        double[] x;

        x = new double[m * n];

        exponent = 1.0 / m;

        for (j = 0; j < n; j++)
        {
            //
            //  Fill a vector with normally distributed values.
            //
            v = typeMethods.r8vec_normal_01_new(m, ref data, ref seed);
            //
            //  Compute the length of the vector.
            //
            norm = 0.0;
            for (i = 0; i < m; i++)
            {
                norm += Math.Pow(v[i], 2);
            }

            norm = Math.Sqrt(norm);
            //
            //  Normalize the vector.
            //
            for (i = 0; i < m; i++)
            {
                x[i + j * m] = v[i] / norm;
            }

            //
            //  Now compute a value to map the point ON the hypersphere INTO the hypersphere.
            //
            r = UniformRNG.r8_uniform_01(ref seed);

            for (i = 0; i < m; i++)
            {
                x[i + j * m] = Math.Pow(r, exponent) * x[i + j * m];
            }
        }

        return x;
    }

    public static double[] hypersphere_01_surface_uniform(int m, int n, ref typeMethods.r8vecNormalData data, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HYPERSPHERE_01_SURFACE_UNIFORM: uniform points on unit hypersphere surface.
        //
        //  Discussion:
        //
        //    The hypersphere has center 0 and radius 1.
        //
        //    This procedure is valid for any spatial dimension DIM_NUM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 December 2013
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
        //    George Marsaglia,
        //    Choosing a point from the surface of a sphere,
        //    Annals of Mathematical Statistics,
        //    Volume 43, Number 2, April 1972, pages 645-646.
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
        //    Input/output, int &SEED, a seed for the random
        //    number generator.
        //
        //    Output, double HYPERSPHERE_01_UNIFORM_SURFACE[M*N], the points.
        //
    {
        int i;
        int j;
        double norm;
        double[] x;
        //
        //  Fill a matrix with normally distributed values.
        //
        x = typeMethods.r8mat_normal_01_new(m, n, ref data, ref seed);
        //
        //  Normalize each column.
        //
        for (j = 0; j < n; j++)
        {
            //
            //  Compute the length of the vector.
            //
            norm = 0.0;
            for (i = 0; i < m; i++)
            {
                norm += Math.Pow(x[i + j * m], 2);
            }

            norm = Math.Sqrt(norm);
            //
            //  Normalize the vector.
            //
            for (i = 0; i < m; i++)
            {
                x[i + j * m] /= norm;
            }
        }

        return x;
    }

    public static double hypersphere_01_volume(int m)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HYPERSPHERE_01_VOLUME computes the volume of a unit hypersphere.
        //
        //  Discussion:
        //
        //    The unit hypersphere satisfies the equation:
        //
        //      Sum ( 1 <= I <= M ) X(I) * X(I) = 1
        //
        //     M  Volume
        //
        //     1    2
        //     2    1        * PI
        //     3  ( 4 /   3) * PI
        //     4  ( 1 /   2) * PI^2
        //     5  ( 8 /  15) * PI^2
        //     6  ( 1 /   6) * PI^3
        //     7  (16 / 105) * PI^3
        //     8  ( 1 /  24) * PI^4
        //     9  (32 / 945) * PI^4
        //    10  ( 1 / 120) * PI^5
        //
        //    For the unit hypersphere, Volume(M) = 2 * PI * Volume(M-2)/ M
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 December 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the dimension of the space.
        //
        //    Output, double HYPERSPHERE_01_VOLUME, the volume.
        //
    {
        int i;
        int m2;
            
        double volume;

        switch (m % 2)
        {
            case 0:
            {
                m2 = m / 2;
                volume = 1.0;
                for (i = 1; i <= m2; i++)
                {
                    volume = volume * Math.PI / i;
                }

                break;
            }
            default:
            {
                m2 = (m - 1) / 2;
                volume = Math.Pow(Math.PI, m2) * Math.Pow(2.0, m);
                for (i = m2 + 1; i <= 2 * m2 + 1; i++)
                {
                    volume /= i;
                }

                break;
            }
        }

        return volume;
    }

    public static void hypersphere_01_volume_values(ref int n_data, ref int m, ref double volume )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HYPERSPHERE_01_VOLUME_VALUES returns some volumes of the unit hypersphere.
        //
        //  Discussion:
        //
        //    The unit hypersphere satisfies the equation:
        //
        //      Sum ( 1 <= I <= M ) X(I) * X(I) = 1
        //
        //     M  Volume
        //
        //     1    1
        //     2    1        * PI
        //     3  ( 4 /   3) * PI
        //     4  ( 1 /   2) * PI^2
        //     5  ( 8 /  15) * PI^2
        //     6  ( 1 /   6) * PI^3
        //     7  (16 / 105) * PI^3
        //     8  ( 1 /  24) * PI^4
        //     9  (32 / 945) * PI^4
        //    10  ( 1 / 120) * PI^5
        //
        //    For the unit hypersphere, Volume(M) = 2 * PI * Volume(M-2) / M
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 December 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input/output, int &N_DATA.
        //    On input, if N_DATA is 0, the first test data is returned, and
        //    N_DATA is set to the index of the test data.  On each subsequent
        //    call, N_DATA is incremented and that test data is returned.  When
        //    there is no more test data, N_DATA is set to 0.
        //
        //    Output, int &M, the spatial dimension.
        //
        //    Output, double &VOLUME, the volume.
        //
    {
        const int N_MAX = 20;

        int[] m_vec =  {
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8,
                9,
                10,
                11,
                12,
                13,
                14,
                15,
                16,
                17,
                18,
                19,
                20
            }
            ;
        double[] volume_vec =  {
                0.2000000000000000E+01,
                0.3141592653589793E+01,
                0.4188790204786391E+01,
                0.4934802200544679E+01,
                0.5263789013914325E+01,
                0.5167712780049970E+01,
                0.4724765970331401E+01,
                0.4058712126416768E+01,
                0.3298508902738707E+01,
                0.2550164039877345E+01,
                0.1884103879389900E+01,
                0.1335262768854589E+01,
                0.9106287547832831E+00,
                0.5992645293207921E+00,
                0.3814432808233045E+00,
                0.2353306303588932E+00,
                0.1409811069171390E+00,
                0.8214588661112823E-01,
                0.4662160103008855E-01,
                0.2580689139001406E-01
            }
            ;

        n_data = n_data switch
        {
            < 0 => 0,
            _ => n_data
        };

        if (N_MAX <= n_data)
        {
            n_data = 0;
            m = 0;
            volume = 0.0;
        }
        else
        {
            m = m_vec[n_data];
            volume = volume_vec[n_data];
            n_data += 1;
        }
    }

    public static double hypersphere_area(int m, double r)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HYPERSPHERE_AREA computes the surface area of a hypersphere.
        //
        //  Discussion:
        //
        //    A hypersphere satisfies the equation:
        //
        //      sum ( ( P(1:M) - C(1:M) )^2 ) = R^2
        //
        //    M   Area
        //
        //    2      2       * PI   * R
        //    3      4       * PI   * R^2
        //    4      2       * PI^2 * R^3
        //    5      (8/3)   * PI^2 * R^4
        //    6                PI^3 * R^5
        //    7      (16/15) * PI^3 * R^6
        //
        //    Sphere_Area ( M, R ) = 2 * PI^(M/2) * R^(M-1) / Gamma ( M / 2 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 December 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the dimension of the space.
        //
        //    Input, double R, the radius.
        //
        //    Output, double HYPERSPHERE_AREA, the area.
        //
    {
        double value = 0;

        value = Math.Pow(r, m - 1) * hypersphere_01_area(m);

        return value;
    }

    public static double[] hypersphere_stereograph(int m, int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HYPERSPHERE_STEREOGRAPH: stereographic mapping of points on a hypersphere.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //   16 December 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the spatial dimension.
        //    M must be at least 2.
        //
        //    Input, int N, the number of points.
        //
        //    Input, double X[M*N], the points to be mapped.
        //
        //    Output, double HYPERSPHERE_STEREOGRAPH[M-1)*N], the stereographically 
        //    mapped points.
        //
    {
        int i;
        int j;
        double[] x2;

        x2 = new double[(m - 1) * n];

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m - 1; i++)
            {
                x2[i + j * (m - 1)] = x[i + j * m];
            }
        }

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m - 1; i++)
            {
                x2[i + j * (m - 1)] /= (1.0 - x[m - 1 + j * m]);
            }
        }

        return x2;
    }

    public static double[] hypersphere_stereograph_inverse(int m, int n, double[] x2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HYPERSPHERE_STEREOGRAPH_INVERSE inverts a stereographic map.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 December 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the spatial dimension.
        //    M must be at least 2.
        //
        //    Input, int N, the number of points.
        //
        //    Input, double X2[(M-1)*N], points in the plane.
        //
        //    Input, double HYPERSPHERE_STEREOGRAPH_INVERSE[M*N], points mapped 
        //    back to the hypersphere.
        //
    {
        double[] d;
        int i;
        int j;
        double[] x;

        x = new double[m * n];

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m - 1; i++)
            {
                x[i + j * m] = 2.0 * x2[i + j * (m - 1)];
            }
        }

        d = new double[n];

        for (j = 0; j < n; j++)
        {
            d[j] = 0.0;
            for (i = 0; i < m - 1; i++)
            {
                d[j] += Math.Pow(x2[i + j * (m - 1)], 2);
            }
        }

        for (j = 0; j < n; j++)
        {
            x[m - 1 + j * m] = d[j] - 1.0;
        }

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m; i++)
            {
                x[i + j * m] /= (d[j] + 1.0);
            }
        }
            
        return x;
    }

    public static double[] hypersphere_surface_uniform(int m, int n, double r, double[] c,
            ref typeMethods.r8vecNormalData data, ref int seed )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HYPERSPHERE_SURFACE_UNIFORM: uniform hypersphere surface samples
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 December 2013
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
        //    George Marsaglia,
        //    Choosing a point from the surface of a sphere,
        //    Annals of Mathematical Statistics,
        //    Volume 43, Number 2, April 1972, pages 645-646.
        //
        //    Reuven Rubinstein,
        //    Monte Carlo Optimization, Simulation, and Sensitivity
        //    of Queueing Networks,
        //    Wiley, 1986, page 234.
        //
        //  Parameters:
        //
        //    Input, int M, the dimension of the space.
        //
        //    Input, int N, the number of points.
        //
        //    Input, real R, the radius.
        //
        //    Input, real C[M], the center.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, real HYPERSPHERE_SURFACE_UNIFORM[M*N], the points.
        //
    {
        int i;
        int j;
        double[] x;

        x = hypersphere_01_surface_uniform(m, n, ref data, ref seed);
        //
        //  Scale by the radius.
        //
        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m; i++)
            {
                x[i + j * m] = r * x[i + j * m];
            }
        }

        //
        //  Shift to the center.
        //
        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m; i++)
            {
                x[i + j * m] += c[i];
            }
        }

        return x;
    }

    public static double[] hypersphere_to_cartesian(int m, int n, double[] c, double[] r,
            double[] theta )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HYPERSPHERE_TO_CARTESIAN: hypersphere to Cartesian coordinate transform.
        //
        //  Discussion:
        //
        //    We allow the trivial case M = 1; in that case alone, the value R
        //    must be assumed to have a sign.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 December 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, integer M, the spatial dimension.
        //    1 <= M.
        //
        //    Input, integer N, the number of points to transform.
        //
        //    Input, real C[M], the center of the hypersphere.
        //
        //    Input, real R[N], the radius of the points on the hypersphere.
        //    Except for the trivial case M = 1, R is assumed nonnegative.
        //
        //    Input, real THETA[(M-1)*N], the coordinate angles of the points,
        //    measured in radians.
        //
        //    Output, real HYPERSPHERE_TO_CARTESIAN[M*N], the Cartesian 
        //    coordinates of the points.
        //
    {
        int i;
        int i1;
        int i2;
        int j;
        double[] x;

        x = new double[m * n];

        switch (m)
        {
            case 1:
            {
                for (j = 0; j < n; j++)
                {
                    x[0 + j * m] = r[j];
                }

                break;
            }
            default:
            {
                for (j = 0; j < n; j++)
                {
                    for (i = 0; i < m; i++)
                    {
                        x[i + j * m] = r[j];
                    }
                }

                for (j = 0; j < n; j++)
                {
                    for (i1 = 0; i1 < m - 1; i1++)
                    {
                        x[i1 + j * m] *= Math.Cos(theta[i1 + j * (m - 1)]);
                        for (i2 = i1 + 1; i2 < m; i2++)
                        {
                            x[i2 + j * m] *= Math.Sin(theta[i1 + j * (m - 1)]);
                        }
                    }
                }

                break;
            }
        }

        //
        //  Add the center.
        //
        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m; i++)
            {
                x[i + j * m] += c[i];
            }
        }

        return x;
    }

    public static double hypersphere_volume(int m, double r)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HYPERSPHERE_VOLUME computes the volume of a hypersphere.
        //
        //  Discussion:
        //
        //    A hypersphere satisfies the equation:
        //
        //      sum ( ( X(1:M) - PC(1:M) )^2 ) = R^2
        //
        //    where R is the radius and PC is the center.
        //
        //    Results include:
        //
        //    M     Volume
        //    -     -----------------------
        //    2                PI   * R^2
        //    3     (4/3)    * PI   * R^3
        //    4     (1/2)    * PI^2 * R^4
        //    5     (8/15)   * PI^2 * R^5
        //    6     (1/6)    * PI^3 * R^6
        //    7     (16/105) * PI^3 * R^7
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 December 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the dimension of the space.
        //
        //    Input, double R, the radius.
        //
        //    Output, double HYPERSPHERE_VOLUME, the volume.
        //
    {
        double value = Math.Pow(r, m) * hypersphere_01_volume(m);

        return value;
    }

    public static double hypersphere_unit_volume(int m)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HYPERSPHERE_UNIT_VOLUME: volume of a unit hypersphere in M dimensions.
        //
        //  Discussion:
        //
        //    The unit sphere in M dimensions satisfies the equation:
        //
        //      Sum ( 1 <= I <= M ) X(I) * X(I) = 1
        //
        //     M  Volume
        //
        //     1    2
        //     2    1        * PI
        //     3  ( 4 /   3) * PI
        //     4  ( 1 /   2) * PI^2
        //     5  ( 8 /  15) * PI^2
        //     6  ( 1 /   6) * PI^3
        //     7  (16 / 105) * PI^3
        //     8  ( 1 /  24) * PI^4
        //     9  (32 / 945) * PI^4
        //    10  ( 1 / 120) * PI^5
        //
        //    For the unit sphere, Volume(M) = 2 * PI * Volume(M-2)/ M
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
        //  Parameters:
        //
        //    Input, int M, the dimension of the space.
        //
        //    Output, double HYPERSPHERE_UNIT_VOLUME, the volume of the sphere.
        //
    {
        int i;
        int m2;
            
        double volume;

        switch (m % 2)
        {
            case 0:
            {
                m2 = m / 2;
                volume = 1.0;
                for (i = 1; i <= m2; i++)
                {
                    volume = volume * Math.PI / i;
                }

                break;
            }
            default:
            {
                m2 = (m - 1) / 2;
                volume = Math.Pow(Math.PI, m2) * Math.Pow(2.0, m);
                for (i = m2 + 1; i <= 2 * m2 + 1; i++)
                {
                    volume /= i;
                }

                break;
            }
        }

        return volume;
    }
}