﻿using System;
using Burkardt.IntegralNS;
using Burkardt.Quadrature;
using Burkardt.SubsetNS;
using Burkardt.Types;

namespace Burkardt.Stroud
{
    public static class Sphere
    {
        public static double sphere_05_nd(Func<int, double[], double> func, int n, double[] center,
                double r)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_05_ND approximates an integral on the surface of a sphere in ND.
            //
            //  Integration region:
            //
            //    R1*R1 <= sum ( X(1:N) - CENTER(1:N) )^2 <= R2*R2
            //
            //  Discussion:
            //
            //    A 2*N+2^N points 5-th degree formula is used, Stroud number UN:5-2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Arthur Stroud,
            //    Approximate Calculation of Multiple Integrals,
            //    Prentice Hall, 1971,
            //    ISBN: 0130438936,
            //    LC: QA311.S85.
            //
            //  Parameters:
            //
            //    Input, Func < int, double[], double> func, the name of the 
            //    user supplied function to be integrated.
            //
            //    Input, int N, the dimension of the space.
            //
            //    Input, double CENTER[N], the center of the sphere.
            //
            //    Input, double R, the radius of the sphere.
            //
            //    Output, double SPHERE_05_ND, the approximate integral of the function.
            //
        {
            int i;
            int iadd = 0;
            int ihi;
            int[] ix;
            bool more;
            int ncard = 0;
            double quad;
            double result;
            double volume;
            double w1;
            double w2;
            double[] x;
            double x1;
            double x2;

            x1 = 1.0;
            x2 = 1.0 / Math.Sqrt((double)(n));

            w1 = 1.0 / (double)(n * (n + 2));
            w2 = (double)(n) / (double)((n + 2) * (int)Math.Pow(2, n));

            x = new double[n];

            for (i = 0; i < n; i++)
            {
                x[i] = center[i];
            }

            quad = 0.0;

            for (i = 0; i < n; i++)
            {
                x[i] = center[i] + r * x1;
                quad = quad + w1 * func(n, x);
                x[i] = center[i] - r * x1;
                quad = quad + w1 * func(n, x);
                x[i] = center[i];
            }

            more = false;
            ihi = (int)Math.Pow(2, n);

            for (i = 0; i < n; i++)
            {
                x[i] = center[i] - r * x2;
            }

            ix = new int[n];

            for (i = 0; i < ihi; i++)
            {
                Subset.subset_gray_next(n, ref ix, ref more, ref ncard, ref iadd);

                if (iadd != 0)
                {
                    x[iadd - 1] = center[iadd - 1] - (x[iadd - 1] - center[iadd - 1]);
                }

                quad = quad + w2 * func(n, x);
            }

            volume = sphere_area_nd(n, r);
            result = quad * volume;

            return result;
        }

        public static double sphere_07_1_nd(Func<int, double[], double> func, int n,
                double[] center, double r)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_07_1_ND approximates an integral on the surface of a sphere in ND.
            //
            //  Integration region:
            //
            //    sum ( X(1:N) - CENTER(1:N) )^2 = R * R.
            //
            //  Discussion:
            //
            //    A 2^N + 2*N*N point 7th degree formula is used, Stroud number UN:7-1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    17 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Arthur Stroud,
            //    Approximate Calculation of Multiple Integrals,
            //    Prentice Hall, 1971,
            //    ISBN: 0130438936,
            //    LC: QA311.S85.
            //
            //  Parameters:
            //
            //    Input, Func < int, double[], double> func, the name of the 
            //    user supplied function to be integrated.
            //
            //    Input, int N, the dimension of the space.
            //
            //    Input, double CENTER[N], the center of the sphere.
            //
            //    Input, double R, the radius of the sphere.
            //
            //    Output, double SPHERE_07_1_ND, the approximate integral of the function.
            //
        {
            int i;
            int iadd = 0;
            int[] ix;
            int j;
            int jhi;
            bool more;
            int ncard = 0;
            double quad;
            double result;
            double volume;
            double w1;
            double w2;
            double w3;
            double[] x;
            double x1;
            double x2;
            double x3;

            ix = new int[n];
            x = new double[n];

            for (i = 0; i < n; i++)
            {
                x[i] = center[i];
            }

            w1 = (double)(8 - n)
                 / (double)(n * (n + 2) * (n + 4));

            w2 = (double)(n * n * n)
                 / (double)((int)Math.Pow(2, n) * n * (n + 2) * (n + 4));

            w3 = 4.0 / (double)(n * (n + 2) * (n + 4));

            x1 = 1.0;
            x2 = 1.0 / Math.Sqrt((double)(n));
            x3 = 1.0 / Math.Sqrt(2.0);

            quad = 0.0;
            //
            //  First term.
            //
            for (i = 0; i < n; i++)
            {
                x[i] = center[i] + r * x1;
                quad = quad + w1 * func(n, x);
                x[i] = center[i] - r * x1;
                quad = quad + w1 * func(n, x);
                x[i] = center[i];
            }

            //
            //  Second term.
            //
            for (i = 0; i < n; i++)
            {
                x[i] = center[i] - r * x2;
            }

            more = false;
            jhi = (int)Math.Pow(2, n);

            for (j = 0; j < jhi; j++)
            {
                Subset.subset_gray_next(n, ref ix, ref more, ref ncard, ref iadd);

                if (iadd != 0)
                {
                    x[iadd - 1] = center[iadd - 1] - (x[iadd - 1] - center[iadd - 1]);
                }

                quad = quad + w2 * func(n, x);
            }

            //
            //  Third term.
            //
            for (i = 0; i < n; i++)
            {
                x[i] = center[i];
            }

            for (i = 0; i < n - 1; i++)
            {
                for (j = i + 1; j < n; j++)
                {
                    x[i] = center[i] + r * x3;
                    x[j] = center[j] + r * x3;
                    quad = quad + w3 * func(n, x);
                    x[i] = center[i] - r * x3;
                    x[j] = center[j] + r * x3;
                    quad = quad + w3 * func(n, x);
                    x[i] = center[i] + r * x3;
                    x[j] = center[j] - r * x3;
                    quad = quad + w3 * func(n, x);
                    x[i] = center[i] - r * x3;
                    x[j] = center[j] - r * x3;
                    quad = quad + w3 * func(n, x);
                    x[i] = center[i];
                    x[j] = center[j];
                }
            }

            volume = sphere_area_nd(n, r);
            result = quad * volume;

            return result;
        }

        public static double sphere_area_3d(double r)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_AREA_3D computes the area of a sphere in 3D.
            //
            //  Integration region:
            //
            //    X*X + Y*Y + Z*Z = R * R
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    11 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double R, the radius of the sphere.
            //
            //    Output, double SPHERE_AREA_3D, the area of the sphere.
            //
        {
            double pi = 3.141592653589793;
            double value;

            value = 4.0 * pi * r * r;

            return value;
        }

        public static double sphere_area_nd(int n, double r)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_AREA_ND computes the area of a sphere in ND.
            //
            //  Integration region:
            //
            //    sum ( X(1:N)^2 ) = R * R
            //
            //  Discussion:
            //
            //    N   Area
            //
            //    2   2       * PI   * R
            //    3   4       * PI   * R^2
            //    4   2       * PI^2 * R^3
            //    5   (8/3)   * PI^2 * R^4
            //    6             PI^3 * R^5
            //    7   (16/15) * PI^3 * R^6
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    11 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the dimension of the space.
            //
            //    Input, double R, the radius of the sphere.
            //
            //    Output, double SPHERE_AREA_ND, the area of the sphere.
            //
        {
            double value;

            value = sphere_unit_area_nd(n) * Math.Pow(r, n - 1);

            return value;
        }

        public static double sphere_cap_area_2d(double r, double h)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_CAP_AREA_2D computes the surface area of a spherical cap in 2D.
            //
            //  Discussion:
            //
            //    Draw any radius of the sphere and note the point P where the radius
            //    intersects the sphere.  Consider the point on the radius line which is
            //    H units from P.  Draw the circle that lies in the plane perpendicular to
            //    the radius, and which intersects the sphere.  The circle divides the sphere
            //    into two pieces, and the corresponding disk divides the solid sphere into
            //    two pieces.  The spherical cap is the part of the solid sphere that
            //    includes the point P.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    12 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double R, the radius of the sphere.
            //
            //    Input, double H, the "height" of the spherical cap. 
            //
            //    Output, double SPHERE_CAP_AREA_2D, the area of the spherical cap.
            //
        {
            double area;
            double pi = 3.141592653589793;
            double theta;

            if (h <= 0.0)
            {
                area = 0.0;
            }
            else if (2.0 * r <= h)
            {
                area = 2.0 * pi * r;
            }
            else
            {
                theta = 2.0 * Math.Asin(Math.Sqrt(r * r - (r - h) * (r - h)) / r);

                area = r * theta;

                if (r <= h)
                {
                    area = 2.0 * pi * r - area;
                }
            }

            return area;
        }

        public static double sphere_cap_area_3d(double r, double h)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_CAP_AREA_3D computes the surface area of a spherical cap in 3D.
            //
            //  Discussion:
            //
            //    Draw any radius of the sphere and note the point P where the radius
            //    intersects the sphere.  Consider the point on the radius line which is
            //    H units from P.  Draw the circle that lies in the plane perpendicular to
            //    the radius, and which intersects the sphere.  The circle divides the sphere
            //    into two pieces, and the corresponding disk divides the solid sphere into
            //    two pieces.  The spherical cap is the part of the solid sphere that
            //    includes the point P.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    12 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double R, the radius of the sphere.
            //
            //    Input, double H, the "height" of the spherical cap. 
            //
            //    Output, double SPHERE_CAP_AREA_3D, the area of the spherical cap.
            //
        {
            double area;
            double pi = 3.141592653589793;

            if (h <= 0.0)
            {
                area = 0.0;
            }
            else if (2.0 * r <= h)
            {
                area = 4.0 * pi * r * r;
            }
            else
            {
                area = 2.0 * pi * r * h;
            }

            return area;
        }

        public static double sphere_cap_area_nd(int dim_num, double r, double h)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_CAP_AREA_ND computes the area of a spherical cap in ND.
            //
            //  Discussion:
            //
            //    The spherical cap is a portion of the surface of the sphere:
            //
            //      sum ( X(1:N)^2 ) = R * R
            //
            //    which is no more than H units from the uppermost point on the sphere.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    12 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Thomas Ericson, Victor Zinoviev,
            //    Codes on Euclidean Spheres,
            //    Elsevier, 2001
            //    QA166.7 E75
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the dimension of the space.
            //
            //    Input, double R, the radius of the sphere.
            //
            //    Input, double H, the "thickness" of the spherical cap,
            //    which is normally between 0 and 2 * R.
            //
            //    Output, double SPHERE_CAP_AREA_ND, the area of the spherical cap.
            //
        {
            double area;
            double haver_sine;
            int i;
            double theta;
            double ti;
            double tj;
            double tk;

            if (h <= 0.0)
            {
                area = 0.0;
                return area;
            }

            if (2.0 * r <= h)
            {
                area = sphere_area_nd(dim_num, r);
                return area;
            }

            //
            //  For cases where R < H < 2 * R, work with the complementary region.
            //
            haver_sine = Math.Sqrt((2.0 * r - h) * h);

            theta = Math.Asin(haver_sine / r);

            if (dim_num < 1)
            {
                area = -1.0;
            }
            else if (dim_num == 1)
            {
                area = 0.0;
            }
            else if (dim_num == 2)
            {
                area = 2.0 * theta * r;
            }
            else
            {
                ti = theta;

                tj = ti;
                ti = 1.0 - Math.Cos(theta);

                for (i = 2; i <= dim_num - 2; i++)
                {
                    tk = tj;
                    tj = ti;
                    ti = ((double)(i - 1) * tk
                          - Math.Cos(theta) * Math.Pow(Math.Sin(theta), i - 1))
                         / (double)(i);
                }

                area = sphere_k(dim_num - 1) * ti * Math.Pow(r, dim_num - 1);
            }

            //
            //  Adjust for cases where R < H < 2R.
            //
            if (r < h)
            {
                area = sphere_area_nd(dim_num, r) - area;
            }

            return area;
        }

        public static double sphere_cap_volume_2d(double r, double h)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_CAP_VOLUME_2D computes the volume of a spherical cap in 2D.
            //
            //  Discussion:
            //
            //    Draw any radius R of the circle and denote as P the point where the
            //    radius intersects the circle.  Now consider the point Q which lies
            //    on the radius and which is H units from P.  The line which is
            //    perpendicular to the radius R and passes through Q divides the
            //    circle into two pieces.  The piece including the point P is the
            //    spherical (circular) cap of height (or thickness) H.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    13 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double R, the radius of the sphere.
            //
            //    Input, double H, the "height" of the spherical cap.
            //
            //    Output, double VOLUME, the volume (area) of the spherical cap.
            //
        {
            double pi = 3.141592653589793;
            double theta;
            double volume;

            if (h <= 0.0)
            {
                volume = 0.0;
            }
            else if (2.0 * r <= h)
            {
                volume = pi * r * r;
            }
            else
            {
                theta = 2.0 * Math.Asin(Math.Sqrt(r * r - (r - h) * (r - h)) / r);
                volume = r * r * (theta - Math.Sin(theta)) / 2.0;
                if (r < h)
                {
                    volume = pi * r * r - volume;
                }
            }

            return volume;
        }

        public static double sphere_cap_volume_3d(double r, double h)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_CAP_VOLUME_3D computes the volume of a spherical cap in 3D.
            //
            //  Discussion:
            //
            //    Draw any radius of the sphere and note the point P where the radius
            //    intersects the sphere.  Consider the point on the radius line which is
            //    H units from P.  Draw the circle that lies in the plane perpendicular to
            //    the radius, and which intersects the sphere.  The circle divides the sphere
            //    into two pieces, and the corresponding disk divides the solid sphere into
            //    two pieces.  The part of the solid sphere that includes the point P
            //    is the spherical cap of height (or thickness) H.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    16 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double R, the radius of the sphere.
            //
            //    Input, double H, the "height" of the spherical cap.
            //
            //    Output, double SPHERE_CAP_VOLUME_3D, the volume of the spherical cap.
            //
        {
            double pi = 3.141592653589793;
            double volume;

            if (h <= 0.0)
            {
                volume = 0.0;
            }
            else if (2.0 * r <= h)
            {
                volume = (4.0 / 3.0) * pi * r * r * r;
            }
            else
            {
                volume = (1.0 / 3.0) * pi * h * h * (3.0 * r - h);
            }

            return volume;
        }

        public static double sphere_cap_volume_nd(int dim_num, double r, double h)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_CAP_VOLUME_ND computes the volume of a spherical cap in ND.
            //
            //  Discussion:
            //
            //    The spherical cap is a portion of the surface and interior of the sphere:
            //
            //      sum ( X(1:N)^2 ) <= R * R
            //
            //    which is no more than H units from some point P on the sphere.
            //
            //
            //    The algorithm proceeds from the observation that the N-dimensional
            //    sphere can be parameterized by a quantity RC that runs along the
            //    radius from the center to the point P.  The value of RC at the
            //    base of the spherical cap is (R-H) and at P it is R.  We intend to
            //    use RC as our integration parameeter.
            //
            //    The volume of the spherical cap is then the integral, as RC goes
            //    from (R-H) to R, of the N-1 dimensional volume of the sphere
            //    of radius RS, where RC * RC + RS * RS = R * R.
            //
            //    The volume of the N-1 dimensional sphere of radius RS is simply 
            //    some constants times RS**(N-1).
            // 
            //    After factoring out the constant terms, and writing RC = R * Math.Cos ( T ),
            //    and RS = R * Math.Sin ( T ), and letting 
            //      T_MAX = Math.Asin ( Math.Sqrt ( ( 2.0D+00 * r - h ) * h / r ) ),
            //    the "interesting part" of our integral becomes
            //
            //      constants * R**N * Integral ( T = 0 to T_MAX ) sin**N ( T ) dT
            //
            //    The integral of sin**N ( T ) dT can be handled by recursion.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    16 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the space.
            //
            //    Input, double R, the radius of the sphere.
            //
            //    Input, double H, the "thickness" of the spherical cap,
            //    which is normally between 0 and 2 * R.
            //
            //    Output, double SPHERE_CAP_VOLUME_ND, the volume of the spherical cap.
            //
        {
            double angle;
            double factor1;
            double factor2;
            double volume;
            double volume2;

            if (h <= 0.0)
            {
                volume = 0.0;
                return volume;
            }

            if (2.0 * r <= h)
            {
                volume = sphere_volume_nd(dim_num, r);
                return volume;
            }

            if (dim_num < 1)
            {
                volume = -1.0;
            }
            else if (dim_num == 1)
            {
                volume = h;
            }
            else
            {
                factor1 = sphere_unit_volume_nd(dim_num - 1);

                angle = Math.Asin(Math.Sqrt((2.0 * r - h) * h / r));

                factor2 = SinPower.sin_power_int(0.0, angle, dim_num);

                volume = factor1 * factor2 * Math.Pow(r, dim_num);

                if (r < h)
                {
                    volume2 = sphere_volume_nd(dim_num, r);
                    volume = volume2 - volume;
                }
            }

            return volume;
        }

        public static double sphere_k(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_K computes a factor useful for spherical computations.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    12 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Thomas Ericson, Victor Zinoviev,
            //    Codes on Euclidean Spheres,
            //    Elsevier, 2001
            //    QA166.7 E75
            //
            //  Parameters:
            //
            //    Input, int N, the dimension of the space.
            //
            //    Output, double SPHERE_K, the factor.
            //
        {
            double pi = 3.141592653589793;
            double value;

            if ((n % 2) == 0)
            {
                value = Math.Pow(2.0 * pi, n / 2);
            }
            else
            {
                value = 2.0 * Math.Pow(2.0 * pi, (n - 1) / 2);
            }

            value = value / (double)(typeMethods.i4_factorial2(n - 2));

            return value;
        }

        public static double sphere_monomial_int_nd(int n, double r, int[] e)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_MONOMIAL_INT_ND integrates a monomial on the surface of a sphere in ND.
            //
            //  Integration region:
            //
            //    sum ( X(1:N)^2 ) = R * R.
            //
            //  Discussion:
            //
            //    The sphere may have nonunit radius, but it must be centered at 0.
            //
            //    The monomial is F(X) = X(1)^E(1) * X(2)^E(2) * ... * X(N)^E(N).
            //
            //    This routine is useful for testing the accuracy of quadrature
            //    rules on the sphere.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    18 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Philip Davis, Philip Rabinowitz,
            //    Methods of Numerical Integration,
            //    Second Edition,
            //    Dover, 2007,
            //    ISBN: 0486453391,
            //    LC: QA299.3.D28.
            //
            //  Parameters:
            //
            //    Input, int N, the dimension of the space.
            //
            //    Input, double R, the radius of the sphere.
            //
            //    Input, int E[N], the exponents of X, Y and Z in the monomial.
            //    Each exponent must be nonnegative.
            //
            //    Output, double SPHERE_MONOMIAL_INT_ND, the integral.
            //
        {
            bool all_zero;
            bool any_odd;
            int e_sum;
            int i;
            double integral;
            double pi = 3.141592653589793;

            integral = 0.0;

            for (i = 0; i < n; i++)
            {
                if (e[i] < 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("SPHERE_MONOMIAL_INT_ND - Fatal error!");
                    Console.WriteLine("  All exponents must be nonnegative.");
                    return (1);
                }
            }

            all_zero = true;
            for (i = 0; i < n; i++)
            {
                if (e[i] != 0)
                {
                    all_zero = false;
                    break;
                }
            }

            any_odd = false;
            for (i = 0; i < n; i++)
            {
                if ((e[i] % 2) == 1)
                {
                    any_odd = true;
                    break;
                }
            }

            e_sum = 0;
            for (i = 0; i < n; i++)
            {
                e_sum = e_sum + e[i];
            }

            if (all_zero)
            {
                integral = 2.0 * Math.Sqrt(Math.Pow(pi, n))
                           / typeMethods.r8_gamma(0.5 * (double)(n));
            }
            else if (any_odd)
            {
                integral = 0.0;
            }
            else
            {
                integral = 2.0;

                for (i = 0; i < n; i++)
                {
                    integral = integral * typeMethods.r8_gamma(0.5 * (double)(e[i] + 1));
                }

                integral = integral / typeMethods.r8_gamma(0.5 * ((double)(e_sum + n)));
            }

            integral = integral * Math.Pow(r, e_sum + 2);

            return integral;
        }

        public static double sphere_shell_03_nd(Func<int, double[], double> func, int n,
                double[] center, double r1, double r2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_SHELL_03_ND approximates an integral inside a spherical shell in ND.
            //
            //  Integration region:
            //
            //    R1*R1 <= sum ( X(1:N) - CENTER(1:N) )^2 <= R2*R2.
            //
            //  Discussion:
            //
            //    An 2*N point 3-rd degree formula is used, Stroud number SN-Shell:3-1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Arthur Stroud,
            //    Approximate Calculation of Multiple Integrals,
            //    Prentice Hall, 1971,
            //    ISBN: 0130438936,
            //    LC: QA311.S85.
            //
            //  Parameters:
            //
            //    Input, Func < int, double[], double> func, the name of the 
            //    user supplied function to be integrated.
            //
            //    Input, int N, the dimension of the space.
            //
            //    Input, double CENTER[N], the center of the spheres.
            //
            //    Input, double R1, R2, the inner and outer radiuses that
            //    define the spherical shell.
            //
            //    Output, double SPHERE_SHELL_03_ND, the approximate integral of the function.
            //
        {
            int i;
            double quad;
            double r;
            double result;
            double rho;
            double volume;
            double w;
            double[] x;

            if (r1 == r2)
            {
                result = 0.0;
                return result;
            }

            rho = r1 / r2;

            r = (double)(n) * (1.0 - Math.Pow(rho, n + 2))
                / ((double)(n + 2) * (1.0 - Math.Pow(rho, n)));
            r = Math.Sqrt(r);
            w = 1.0 / (double)(2 * n);

            x = new double[n];

            for (i = 0; i < n; i++)
            {
                x[i] = center[i];
            }

            quad = 0.0;
            for (i = 0; i < n; i++)
            {
                x[i] = center[i] + r * r2;
                quad = quad + w * func(n, x);
                x[i] = center[i] - r * r2;
                quad = quad + w * func(n, x);
                x[i] = center[i];
            }

            volume = sphere_shell_volume_nd(n, r1, r2);
            result = quad * volume;

            return result;
        }

        public static double sphere_shell_volume_nd(int n, double r1, double r2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_SHELL_VOLUME_ND computes the volume of a spherical shell in ND.
            //
            //  Integration region:
            //
            //    R1*R1 <= sum ( X(1:N) - CENTER(1:N) )^2 <= R2*R2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the dimension of the space.
            //
            //    Input, double R1, R2, the radiuses of the inner and 
            //    outer spheres.
            //
            //    Output, double SPHERE_SHELL_VOLUME_ND, the volume of the
            //    spherical shell.
            //
        {
            double volume;

            volume = Ball.ball_volume_nd(n, r2) - Ball.ball_volume_nd(n, r1);

            return volume;
        }

        public static double sphere_unit_03_nd(Func<int, double[], double> func, int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_UNIT_03_ND approximates an integral on the surface of the unit sphere in ND.
            //
            //  Integration region:
            //
            //    sum ( X(1:N)^2 ) = 1.
            //
            //  Discussion:
            //
            //    A 2*N point 3rd degree formula is used, Stroud number UN:3-1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    23 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Arthur Stroud,
            //    Approximate Calculation of Multiple Integrals,
            //    Prentice Hall, 1971,
            //    ISBN: 0130438936,
            //    LC: QA311.S85.
            //
            //  Parameters:
            //
            //    Input, Func < int, double[], double> func, the name of the
            //    user supplied function to be integrated.
            //
            //    Input, int N, the dimension of the space.
            //
            //    Output, double SPHERE_UNIT_03_ND, the approximate integral of the function.
            //
        {
            int i;
            double quad;
            double result;
            double volume;
            double w;
            double[] x;

            x = new double[n];

            for (i = 0; i < n; i++)
            {
                x[i] = 0.0;
            }

            w = 1.0 / (double)(2 * n);

            quad = 0.0;
            for (i = 0; i < n; i++)
            {
                x[i] = 1.0;
                quad = quad + w * func(n, x);
                x[i] = -1.0;
                quad = quad + w * func(n, x);
                x[i] = 0.0;
            }

            volume = sphere_unit_area_nd(n);
            result = quad * volume;

            return result;
        }

        public static double sphere_unit_04_nd(Func<int, double[], double> func, int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_UNIT_04_ND approximates an integral on the surface of the unit sphere in ND.
            //
            //  Integration region:
            //
            //    sum ( X(1:N)^2 ) = 1.
            //
            //  Discussion:
            //
            //    A 2*N*N point 5th degree formula is used, Stroud number UN:5-1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    23 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Arthur Stroud,
            //    Approximate Calculation of Multiple Integrals,
            //    Prentice Hall, 1971,
            //    ISBN: 0130438936,
            //    LC: QA311.S85.
            //
            //  Parameters:
            //
            //    Input, Func < int, double[], double> func, the name of the
            //    user supplied function to be integrated.
            //
            //    Input, int N, the dimension of the space.
            //
            //    Output, double SPHERE_UNIT_04_ND, the approximate integral of the function.
            //
        {
            int i;
            int j;
            double quad;
            double result;
            double s;
            double volume;
            double w1;
            double w2;
            double[] x;

            x = new double[n];

            for (i = 0; i < n; i++)
            {
                x[i] = 0.0;
            }

            w1 = (double)(4 - n) / (double)(2 * n * (n + 2));

            quad = 0.0;

            for (i = 0; i < n; i++)
            {
                x[i] = 1.0;
                quad = quad + w1 * func(n, x);
                x[i] = -1.0;
                quad = quad + w1 * func(n, x);
                x[i] = 0.0;
            }

            s = 1.0 / Math.Sqrt(2.0);
            w2 = 1.0 / (double)(n * (n + 2));

            for (i = 0; i < n; i++)
            {
                x[i] = s;

                for (j = i + 1; j < n; j++)
                {
                    x[j] = s;
                    quad = quad + w2 * func(n, x);
                    x[j] = -s;
                    quad = quad + w2 * func(n, x);
                    x[j] = 0.0;
                }

                x[i] = -s;

                for (j = i + 1; j < n; j++)
                {
                    x[j] = s;
                    quad = quad + w2 * func(n, x);
                    x[j] = -s;
                    quad = quad + w2 * func(n, x);
                    x[j] = 0.0;
                }

                x[i] = 0.0;
            }

            volume = sphere_unit_area_nd(n);
            result = quad * volume;

            return result;
        }

        public static double sphere_unit_05_nd(Func<int, double[], double> func, int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_UNIT_05_ND approximates an integral on the surface of the unit sphere in ND.
            //
            //  Integration region:
            //
            //    sum ( X(1:N)^2 ) = 1.
            //
            //  Discussion:
            //
            //    A 2*N+2**N points 5-th degree formula is used, Stroud number UN:5-2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    23 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Arthur Stroud,
            //    Approximate Calculation of Multiple Integrals,
            //    Prentice Hall, 1971,
            //    ISBN: 0130438936,
            //    LC: QA311.S85.
            //
            //  Parameters:
            //
            //    Input, Func < int, double[], double> func, the name of the
            //    user supplied function to be integrated.
            //
            //    Input, int N, the dimension of the space.
            //
            //    Output, double SPHERE_UNIT_05_ND, the approximate integral of the function.
            //
        {
            int i;
            int iadd = 0;
            int ihi;
            int[] ix;
            bool more;
            int ncard = 0;
            double quad;
            double result;
            double volume;
            double w1;
            double w2;
            double[] x;
            double x1;
            double x2;

            x1 = 1.0;
            x2 = 1.0 / Math.Sqrt((double)(n));

            w1 = 1.0 / (double)(n * (n + 2));
            w2 = (double)(n) / (double)((n + 2) * (int)Math.Pow(2, n));

            ix = new int[n];
            x = new double[n];

            for (i = 0; i < n; i++)
            {
                x[i] = 0.0;
            }

            quad = 0.0;

            for (i = 0; i < n; i++)
            {
                x[i] = x1;
                quad = quad + w1 * func(n, x);
                x[i] = -x1;
                quad = quad + w1 * func(n, x);
                x[i] = 0.0;
            }

            more = false;
            ihi = (int)Math.Pow(2, n);

            for (i = 0; i < n; i++)
            {
                x[i] = -x2;
            }

            for (i = 0; i < ihi; i++)
            {
                Subset.subset_gray_next(n, ref ix, ref more, ref ncard, ref iadd);

                if (iadd != 0)
                {
                    x[iadd - 1] = -x[iadd - 1];
                }

                quad = quad + w2 * func(n, x);
            }

            volume = sphere_unit_area_nd(n);
            result = quad * volume;

            return result;
        }

        public static double sphere_unit_07_3d(Func<double, double, double, double> func)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_UNIT_07_3D approximates an integral on the surface of the unit sphere in 3D.
            //
            //  Integration region:
            //
            //    X*X + Y*Y + Z*Z = 1.
            //
            //  Discussion:
            //
            //    A 32 point 7-th degree formula is used.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Arthur Stroud,
            //    Approximate Calculation of Multiple Integrals,
            //    Prentice Hall, 1971,
            //    ISBN: 0130438936,
            //    LC: QA311.S85.
            //
            //  Parameters:
            //
            //    Input, Func< double, double, double, double > func, the name of the 
            //    user supplied function to be integrated.
            //
            //    Output, double SPHERE_UNIT_07_3D, the approximate integral of the function.
            //
        {
            double angle;
            int i;
            int j;
            int k;
            int order1 = 2;
            int order2 = 4;
            int order3 = 4;
            double pi = 3.141592653589793;
            double quad;
            double result;
            double volume;
            double[] weight1 = new double[2];
            double[] weight2 = new double[4];
            double[] weight3 = new double[4];
            double x;
            double[] xtab1 = new double[2];
            double[] xtab2 = new double[4];
            double[] xtab3 = new double[4];
            double y;
            double z;
            //
            //  Set XTAB1 and WATE1.
            //
            xtab1[0] = -1.0;
            xtab1[1] = 1.0;
            weight1[0] = 1.0;
            weight1[1] = 1.0;
            //
            //  Set XTAB2 and WATE2.
            //
            for (j = 0; j < order2; j++)
            {
                angle = pi * (double)(2 * j + 1) / (double)(2 * order2);
                xtab2[j] = Math.Cos(angle);
            }

            for (j = 0; j < order2; j++)
            {
                weight2[j] = 1.0 / (double)(4 * order2);
            }

            //
            //  Set XTAB3 and WATE3.
            //
            LegendreQuadrature.legendre_set(order3, ref xtab3, ref weight3);

            quad = 0.0;
            for (i = 0; i < order1; i++)
            {
                for (j = 0; j < order2; j++)
                {
                    for (k = 0; k < order3; k++)
                    {
                        x = xtab1[i] * Math.Sqrt(1.0 - xtab2[j] * xtab2[j])
                                     * Math.Sqrt(1.0 - xtab3[k] * xtab3[k]);
                        y = xtab1[i] * xtab2[j] * Math.Sqrt(1.0 - xtab3[k] * xtab3[k]);
                        z = xtab1[i] * xtab3[k];

                        quad = quad + weight1[i] * weight2[j] * weight3[k] * func(x, y, z);
                    }
                }
            }

            volume = sphere_unit_area_3d();
            result = quad * volume;

            return result;
        }

        public static double sphere_unit_07_1_nd(Func<int, double[], double> func, int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_UNIT_07_1_ND approximates an integral on the surface of the unit sphere in ND.
            //
            //  Integration region:
            //
            //    sum ( X(1:N)^2 ) = 1.
            //
            //  Discussion:
            //
            //    A 2**N + 2*N*N point 7th degree formula is used, Stroud number UN:7-1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Arthur Stroud,
            //    Approximate Calculation of Multiple Integrals,
            //    Prentice Hall, 1971,
            //    ISBN: 0130438936,
            //    LC: QA311.S85.
            //
            //  Parameters:
            //
            //    Input, Func < int, double[], double> func, the name of the 
            //    user supplied function to be integrated.
            //
            //    Input, int N, the dimension of the space.
            //
            //    Output, double SPHERE_UNIT_07_1_ND, the approximate integral of the function.
            //
        {
            int i;
            int iadd = 0;
            int[] ix;
            int j;
            int jhi;
            bool more;
            int ncard = 0;
            double quad;
            double result;
            double volume;
            double w1;
            double w2;
            double w3;
            double[] x;
            double x1;
            double x2;
            double x3;

            ix = new int[n];
            x = new double[n];

            w1 = (double)(8 - n) / (double)(n * (n + 2) * (n + 4));
            w2 = (double)(n * n * n)
                 / (double)((int)Math.Pow(2, n) * n * (n + 2) * (n + 4));
            w3 = 4.0 / (double)(n * (n + 2) * (n + 4));

            x1 = 1.0;
            x2 = 1.0 / Math.Sqrt((double)(n));
            x3 = 1.0 / Math.Sqrt(2.0);

            for (i = 0; i < n; i++)
            {
                x[i] = 0.0;
            }

            quad = 0.0;
            //
            //  First term.
            //
            for (i = 0; i < n; i++)
            {
                x[i] = x1;
                quad = quad + w1 * func(n, x);
                x[i] = -x1;
                quad = quad + w1 * func(n, x);
                x[i] = 0.0;
            }

            //
            //  Second term.
            //
            for (i = 0; i < n; i++)
            {
                x[i] = -x2;
            }

            more = false;
            jhi = (int)Math.Pow(2, n);

            for (j = 0; j < jhi; j++)
            {
                Subset.subset_gray_next(n, ref ix, ref more, ref ncard, ref iadd);


                if (iadd != 0)
                {
                    x[iadd - 1] = -x[iadd - 1];
                }

                quad = quad + w2 * func(n, x);
            }

            //
            //  Third term.
            //
            for (i = 0; i < n; i++)
            {
                x[i] = 0.0;
            }

            for (i = 0; i < n - 1; i++)
            {
                for (j = i + 1; j < n; j++)
                {
                    x[i] = x3;
                    x[j] = x3;
                    quad = quad + w3 * func(n, x);
                    x[i] = -x3;
                    x[j] = x3;
                    quad = quad + w3 * func(n, x);
                    x[i] = x3;
                    x[j] = -x3;
                    quad = quad + w3 * func(n, x);
                    x[i] = -x3;
                    x[j] = -x3;
                    quad = quad + w3 * func(n, x);
                    x[i] = 0.0;
                    x[j] = 0.0;
                }
            }

            volume = sphere_unit_area_nd(n);
            result = quad * volume;

            return result;
        }

        public static double sphere_unit_07_2_nd(Func<int, double[], double> func, int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_UNIT_07_2_ND approximates an integral on the surface of the unit sphere in ND.
            //
            //  Integration region:
            //
            //    sum ( X(1:N)^2 ) = 1.
            //
            //  Discussion:
            //
            //    A 2^N * ( N + 1 ) point 7th degree formula is used, Stroud number UN:7-2.
            //
            //    Some of the weights in this quadrature formula are negative.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Arthur Stroud,
            //    Approximate Calculation of Multiple Integrals,
            //    Prentice Hall, 1971,
            //    ISBN: 0130438936,
            //    LC: QA311.S85.
            //
            //  Parameters:
            //
            //    Input, Func < int, double[], double> func, the name of the 
            //    user supplied function to be integrated.
            //
            //    Input, int N, the dimension of the space.
            //
            //    Output, double SPHERE_UNIT_07_2_ND, the approximate integral of the function.
            //
        {
            int iadd = 0;
            int i;
            int[] ix;
            int j;
            int jhi;
            bool more;
            int ncard = 0;
            double quad;
            double result;
            double volume;
            double w1;
            double w2;
            double[] x;
            double x1;
            double x2;
            double x3;

            ix = new int[n];
            x = new double[n];

            for (i = 0; i < n; i++)
            {
                x[i] = 0.0;
            }

            w1 = -(double)(n * n)
                 / (double)((int)Math.Pow(2, n + 3) * (n + 2));
            w2 = (double)((n + 4) * (n + 4))
                 / (double)((int)Math.Pow(2, n + 3) * n * (n + 2));
            x1 = 1.0 / Math.Sqrt((double)(n));
            x2 = Math.Sqrt(5.0 / (double)(n + 4));
            x3 = 1.0 / Math.Sqrt((double)(n + 4));

            quad = 0.0;

            for (j = 0; j < n; j++)
            {
                x[j] = -x1;
            }

            more = false;
            jhi = (int)Math.Pow(2, n);

            for (j = 0; j < jhi; j++)
            {
                Subset.subset_gray_next(n, ref ix, ref more, ref ncard, ref iadd);

                if (iadd != 0)
                {
                    x[iadd - 1] = -x[iadd - 1];
                }

                quad = quad + w1 * func(n, x);
            }

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    x[j] = -x3;
                }

                x[i] = -x2;
                more = false;

                for (j = 0; j < jhi; j++)
                {
                    Subset.subset_gray_next(n, ref ix, ref more, ref ncard, ref iadd);

                    if (iadd != 0)
                    {
                        x[iadd - 1] = -x[iadd - 1];
                    }

                    quad = quad + w2 * func(n, x);
                }
            }

            volume = sphere_unit_area_nd(n);
            result = quad * volume;

            return result;
        }

        public static double sphere_unit_11_3d(Func<double, double, double, double> func)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_UNIT_11_3D approximates an integral on the surface of the unit sphere in 3D.
            //
            //  Integration region:
            //
            //    X*X + Y*Y + Z*Z = 1.
            //
            //  Discussion:
            //
            //    A 50 point 11-th degree formula is used, Stroud number U3:11-1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    AD McLaren,
            //    Mathematics of Computation,
            //    Volume 17, pages 361-383, 1963.
            //
            //    Arthur Stroud,
            //    Approximate Calculation of Multiple Integrals,
            //    Prentice Hall, 1971,
            //    ISBN: 0130438936,
            //    LC: QA311.S85.
            //
            //  Parameters:
            //
            //    Input, Func< double, double, double, double > func, the name of the 
            //    user supplied function to be integrated.
            //
            //    Output, double SPHERE_UNIT_11_3D, the approximate integral of the function.
            //
        {
            int i;
            int j;
            int k;
            int l;
            double quad;
            double result;
            double volume;
            double w1;
            double w2;
            double w3;
            double w4;
            double x;
            double y;
            double z;

            quad = 0.0;

            w1 = 9216.0 / 725760.0;
            x = 1.0;
            y = 0.0;
            z = 0.0;

            for (i = 0; i < 2; i++)
            {
                x = -x;
                for (j = 0; j < 3; j++)
                {
                    typeMethods.r8_swap3(ref x, ref y, ref z);
                    quad = quad + w1 * func(x, y, z);
                }
            }

            w2 = 16384.0 / 725760.0;
            x = Math.Sqrt(0.5);
            y = Math.Sqrt(0.5);
            z = 0.0;

            for (i = 0; i < 2; i++)
            {
                x = -x;
                for (j = 0; j < 2; j++)
                {
                    y = -y;
                    for (k = 0; k < 3; k++)
                    {
                        typeMethods.r8_swap3(ref x, ref y, ref z);
                        quad = quad + w2 * func(x, y, z);
                    }
                }
            }

            w3 = 15309.0 / 725760.0;
            x = Math.Sqrt(1.0 / 3.0);
            y = Math.Sqrt(1.0 / 3.0);
            z = Math.Sqrt(1.0 / 3.0);

            for (i = 0; i < 2; i++)
            {
                x = -x;
                for (j = 0; j < 2; j++)
                {
                    y = -y;
                    for (k = 0; k < 2; k++)
                    {
                        z = -z;
                        quad = quad + w3 * func(x, y, z);
                    }
                }
            }

            w4 = 14641.0 / 725760.0;
            x = Math.Sqrt(1.0 / 11.0);
            y = Math.Sqrt(1.0 / 11.0);
            z = 3.0 * Math.Sqrt(1.0 / 11.0);

            for (i = 0; i < 2; i++)
            {
                x = -x;
                for (j = 0; j < 2; j++)
                {
                    y = -y;
                    for (k = 0; k < 2; k++)
                    {
                        z = -z;
                        for (l = 0; l < 3; l++)
                        {
                            typeMethods.r8_swap3(ref x, ref y, ref z);
                            quad = quad + w4 * func(x, y, z);
                        }
                    }
                }
            }

            volume = sphere_unit_area_3d();
            result = quad * volume;

            return result;
        }

        public static double sphere_unit_11_nd(Func<int, double[], double> func, int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_UNIT_11_ND: integral on the surface of the unit sphere in ND.
            //
            //  Discussion:
            //
            //    The integration region:
            //
            //      sum ( X(1:N)^2 ) = 1
            //
            //    An 2^N * ( N^2 + N + 1 ) point formula of degree 5 is used.
            //
            //    (For N = 3, the number of points is actually only 56, and
            //     for N = 4, the number of points is actually only 240.)
            //
            //    One element of COEF31 was changed from
            //      0.0236339091329 to
            //      0.0236639091329
            //    by Stroud, when going from his paper to his later textbook.
            //    This correction was pointed out by David Wright, 16 February 2010.
            //
            //    One element of COEF21 was incorrectly transcribed.  The correct
            //    value of COEF21(7) is 0.0337329118818D+00, as pointed out by
            //    John Nolan, 23 April 2013.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    23 April 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Arthur Stroud,
            //    A Fifth Degree Integration Formula for the N-Simplex,
            //    SIAM Journal on Numerical Analysis,
            //    Volume 6, Number 1, March 1969.
            //
            //    Arthur Stroud,
            //    Approximate Calculation of Multiple Integrals,
            //    Prentice Hall, 1971,
            //    ISBN: 0130438936,
            //    LC: QA311.S85.
            //
            //  Parameters:
            //
            //    Input, Func < int, double[], double> func, the name of the 
            //    user supplied function to be integrated.
            //
            //    Input, int N, the dimension of the space.  For this routine,
            //    it must be the case that 3 <= N <= 16.
            //
            //    Output, double SPHERE_UNIT_11_ND, the approximate integral.
            //
        {
            double area;
            double[] coef1 =
            {
                0.0,
                0.0,
                0.128571428571,
                0.0518518518518,
                0.0211979378646,
                0.281250000000,
                1.11934731935,
                2.82751322751,
                5.68266145619,
                9.93785824515,
                15.8196616478,
                23.5285714285,
                33.2409299392,
                45.1113811729,
                59.2754264177,
                75.8518518518
            };
            double[] coef21 =
            {
                0.0,
                0.0,
                0.163795782462,
                0.0967270533860,
                0.0638253880175,
                0.0452340041459,
                0.0337329118818,
                0.0261275095270,
                0.0208331595340,
                0.0169937111647,
                0.0141147212492,
                0.0118949128383,
                0.0101424250926,
                0.00873046796644,
                0.00757257014768,
                0.00660813369775
            };
            double[] coef22 =
            {
                0.0,
                0.0,
                0.126680408014,
                0.0514210947621,
                0.0213579471658,
                -0.108726067638,
                -0.371589499738,
                -0.786048144448,
                -1.36034060198,
                -2.09547695631,
                -2.98784764467,
                -4.03107480702,
                -5.21726499521,
                -6.53783099707,
                -7.98401677102,
                -9.54722261180
            };
            double[] coef31 =
            {
                0.0,
                0.0,
                0.0,
                0.0592592592592,
                0.0236639091329,
                0.0525940190875,
                0.0925052768546,
                0.141316953438,
                0.196818580052,
                0.257027634179,
                0.320299222258,
                0.385326226441,
                0.451098131789,
                0.516849445559,
                0.582010515746,
                0.646165210110
            };
            double[] coef32 =
            {
                0.0,
                0.0,
                0.0,
                0.0,
                0.0316246294890,
                0.0207194729760,
                0.0144303800811,
                0.0105348984135,
                0.00798435122193,
                0.00623845929545,
                0.00499896882962,
                0.00409176297655,
                0.00341037426698,
                0.00288710646943,
                0.00247745182907,
                0.00215128820597
            };
            int i;
            int iadd = 0;
            int[] ix;
            int j;
            int k;
            bool more;
            int ncard = 0;
            double quad;
            double r1;
            double r2;
            double result;
            double s1;
            double s2;
            double u1;
            double u2;
            double v1;
            double v2;
            double[] x;

            result = 0.0;

            if (n < 3 || 16 < n)
            {
                Console.WriteLine("");
                Console.WriteLine("SPHERE_UNIT_11_ND - Fatal error!");
                Console.WriteLine("  Input spatial dimension N out of range.");
                Console.WriteLine("  N = " + n + "");
                return (1);
            }

            ix = new int[n];
            x = new double[n];

            quad = 0.0;
            //
            //  S1
            //
            for (i = 0; i < n; i++)
            {
                x[i] = 1.0 / Math.Sqrt((double)(n));
            }

            more = false;

            for (;;)
            {
                Subset.subset_gray_next(n, ref ix, ref more, ref ncard, ref iadd);

                if (iadd != 0)
                {
                    x[iadd - 1] = -x[iadd - 1];
                }

                quad = quad + coef1[n - 1] * func(n, x);

                if (!more)
                {
                    break;
                }
            }

            //
            //  S21
            //
            r1 = ((double)(n + 6) - 4.0 * Math.Sqrt(3.0))
                 / (double)(n * n + 12 * n - 12);
            r1 = Math.Sqrt(r1);

            s1 = ((double)(7 * n - 6)
                  + (double)(4 * (n - 1)) * Math.Sqrt(3.0))
                 / (double)(n * n + 12 * n - 12);
            s1 = Math.Sqrt(s1);

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    x[j] = r1;
                }

                x[i] = s1;

                more = false;

                for (;;)
                {
                    Subset.subset_gray_next(n, ref ix, ref more, ref ncard, ref iadd);

                    if (iadd != 0)
                    {
                        x[iadd - 1] = -x[iadd - 1];
                    }

                    quad = quad + coef21[n - 1] * func(n, x);

                    if (!more)
                    {
                        break;
                    }
                }
            }

            //
            //  S22
            //
            r2 = ((double)(n + 6) + 4.0 * Math.Sqrt(3.0))
                 / (double)(n * n + 12 * n - 12);
            r2 = Math.Sqrt(r2);

            s2 = ((double)(7 * n - 6)
                  - (double)(4 * (n - 1)) * Math.Sqrt(3.0))
                 / (double)(n * n + 12 * n - 12);
            s2 = Math.Sqrt(s2);

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    x[j] = r2;
                }

                x[i] = s2;

                more = false;

                for (;;)
                {
                    Subset.subset_gray_next(n, ref ix, ref more, ref ncard, ref iadd);

                    if (iadd != 0)
                    {
                        x[iadd - 1] = -x[iadd - 1];
                    }

                    quad = quad + coef22[n - 1] * func(n, x);

                    if (!more)
                    {
                        break;
                    }
                }
            }

            //
            //  S31
            //
            u1 = ((double)(n + 12) + 8.0 * Math.Sqrt(3.0))
                 / (double)(n * n + 24 * n - 48);
            u1 = Math.Sqrt(u1);

            v1 = ((double)(7 * n - 12)
                  - (double)(4 * n - 8) * Math.Sqrt(3.0))
                 / (double)(n * n + 24 * n - 48);
            v1 = Math.Sqrt(v1);

            for (i = 0; i < n - 1; i++)
            {
                for (j = i + 1; j < n; j++)
                {
                    for (k = 0; k < n; k++)
                    {
                        x[k] = u1;
                    }

                    x[i] = v1;
                    x[j] = v1;

                    more = false;

                    for (;;)
                    {
                        Subset.subset_gray_next(n, ref ix, ref more, ref ncard, ref iadd);

                        if (iadd != 0)
                        {
                            x[iadd - 1] = -x[iadd - 1];
                        }

                        quad = quad + coef31[n - 1] * func(n, x);

                        if (!more)
                        {
                            break;
                        }
                    }
                }
            }

            //
            //  S32
            //
            u2 = ((double)(n + 12) - 8.0 * Math.Sqrt(3.0))
                 / (double)(n * n + 24 * n - 48);
            u2 = Math.Sqrt(u2);

            v2 = ((double)(7 * n - 12)
                  + (double)(4 * n - 8) * Math.Sqrt(3.0))
                 / (double)(n * n + 24 * n - 48);
            v2 = Math.Sqrt(v2);

            for (i = 0; i < n - 1; i++)
            {
                for (j = i + 1; j < n; j++)
                {
                    for (k = 0; k < n; k++)
                    {
                        x[k] = u2;
                    }

                    x[i] = v2;
                    x[j] = v2;

                    more = false;

                    for (;;)
                    {
                        Subset.subset_gray_next(n, ref ix, ref more, ref ncard, ref iadd);

                        if (iadd != 0)
                        {
                            x[iadd - 1] = -x[iadd - 1];
                        }

                        quad = quad + coef32[n - 1] * func(n, x);

                        if (!more)
                        {
                            break;
                        }
                    }
                }
            }

            area = sphere_unit_area_nd(n);
            result = quad * area / Math.Pow(2.0, n);

            return result;
        }

        public static double sphere_unit_14_3d(Func<double, double, double, double> func)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_UNIT_14_3D: integral on the surface of the unit sphere in 3D.
            //
            //  Integration region:
            //
            //    X*X + Y*Y + Z*Z = 1.
            //
            //  Discussion:
            //
            //    A 72 point 14-th degree formula is used, Stroud number U3:14-1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    25 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    AD McLaren,
            //    Mathematics of Computation,
            //    Volume 17, pages 361-383, 1963.
            //
            //    Arthur Stroud,
            //    Approximate Calculation of Multiple Integrals,
            //    Prentice Hall, 1971,
            //    ISBN: 0130438936,
            //    LC: QA311.S85.
            //
            //  Parameters:
            //
            //    Input, Func< double, double, double, double > func, the name of the 
            //    user supplied function to be integrated.
            //
            //    Output, double SPHERE_UNIT_14_3D, the approximate integral of the function.
            //
        {
            int i;
            int j;
            int k;
            double quad;
            double result;
            double temp;
            double volume;
            double w1;
            double w2;
            double x;
            double[] xtab =
            {
                -0.151108275, 0.315838353, 0.346307112, -0.101808787, -0.409228403
            };
            double y;
            double[] ytab =
            {
                0.155240600, 0.257049387, 0.666277790, 0.817386065, 0.501547712
            };
            double z;
            double[] ztab =
            {
                0.976251323, 0.913330032, 0.660412970, 0.567022920, 0.762221757
            };

            quad = 0.0;

            w1 = 125.0 / 10080.0;
            x = 0.525731112;
            y = 0.850650808;
            z = 0.0;

            for (i = 0; i < 2; i++)
            {
                x = -x;
                for (j = 0; j < 2; j++)
                {
                    y = -y;
                    for (k = 0; k < 3; k++)
                    {
                        typeMethods.r8_swap3(ref x, ref y, ref z);
                        quad = quad + w1 * func(x, y, z);
                    }
                }
            }

            w2 = 143.0 / 10080.0;

            for (i = 0; i < 5; i++)
            {
                x = xtab[i];
                y = ytab[i];
                z = ztab[i];

                for (j = 0; j < 3; j++)
                {
                    temp = x;
                    x = z;
                    z = -y;
                    y = -temp;

                    for (k = 0; k < 3; k++)
                    {
                        typeMethods.r8_swap3(ref x, ref y, ref z);
                        quad = quad + w2 * func(x, y, z);
                    }

                    y = -y;
                    z = -z;
                    quad = quad + w2 * func(x, y, z);
                }
            }

            volume = sphere_unit_area_3d();
            result = quad * volume;

            return result;
        }

        public static double sphere_unit_15_3d(Func<double, double, double, double> func)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_UNIT_15_3D approximates an integral on the surface of the unit sphere in 3D.
            //
            //  Integration region:
            //
            //    X*X + Y*Y + Z*Z = 1.
            //
            //  Discussion:
            //
            //    A 128 point 15-th degree spherical product Gauss formula is used.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    25 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Arthur Stroud,
            //    Approximate Calculation of Multiple Integrals,
            //    Prentice Hall, 1971,
            //    ISBN: 0130438936,
            //    LC: QA311.S85.
            //
            //  Parameters:
            //
            //    Input, Func< double, double, double, double > func, the name of the 
            //    user supplied function to be integrated.
            //
            //    Output, double SPHERE_UNIT_15_3D, the approximate integral of the function.
            //
        {
            double angle;
            int i;
            int j;
            int k;
            int order = 8;
            double pi = 3.141592653589793;
            double quad;
            double result;
            double volume;
            double[] weight;
            double x;
            double[] xtab;
            double y;
            double z;

            weight = new double[order];
            xtab = new double[order];

            LegendreQuadrature.legendre_set(order, ref xtab, ref weight);

            for (i = 0; i < order; i++)
            {
                weight[i] = weight[i] / 32.0;
            }

            quad = 0.0;

            for (j = 0; j < order; j++)
            {
                for (k = 0; k < 16; k++)
                {
                    angle = (double)(k) * pi / 8.0;
                    x = Math.Sqrt(1.0 - xtab[j] * xtab[j]) * Math.Cos(angle);
                    y = Math.Sqrt(1.0 - xtab[j] * xtab[j]) * Math.Sin(angle);
                    z = xtab[j];

                    quad = quad + weight[j] * func(x, y, z);
                }
            }

            volume = sphere_unit_area_3d();
            result = quad * volume;

            return result;
        }

        public static double sphere_unit_area_3d()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_UNIT_AREA_3D computes the surface area of the unit sphere in 3D.
            //
            //  Integration region:
            //
            //    X*X + Y*Y + Z*Z = 1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    15 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, double SPHERE_UNIT_AREA_3D, the area of the sphere.
            //
        {
            double area;
            double pi = 3.141592653589793;

            area = 4.0 * pi;

            return area;
        }

        public static double sphere_unit_area_nd(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_UNIT_AREA_ND computes the surface area of the unit sphere in ND.
            //
            //  Integration region:
            //
            //    sum ( ( X(1:N) - CENTER(1:N) )^2 ) = R * R.
            //
            //  Discussion:
            //
            //    N   Area
            //
            //    2   2       * PI
            //    3   4       * PI
            //    4   2       * PI^2
            //    5   (8/3)   * PI^2
            //    6             PI^3
            //    7   (16/15) * PI^3
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    12 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the dimension of the space.
            //
            //    Output, double SPHERE_UNIT_AREA_ND, the area of the sphere.
            //
        {
            double area;
            int i;
            int m;
            double pi = 3.141592653589793;

            if ((n % 2) == 0)
            {
                m = n / 2;
                area = 2.0 * Math.Pow(pi, m);
                for (i = 1; i <= m - 1; i++)
                {
                    area = area / (double)(i);
                }
            }
            else
            {
                m = (n - 1) / 2;
                area = Math.Pow(2.0, n) * Math.Pow(pi, m);
                for (i = m + 1; i <= 2 * m; i++)
                {
                    area = area / (double)(i);
                }
            }

            return area;
        }

        public static void sphere_unit_area_values(ref int n_data, ref int n, ref double area)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_UNIT_AREA_VALUES returns some areas of the unit sphere in ND.
            //
            //  Discussion:
            //
            //    The formula for the surface area of the unit sphere in N dimensions is:
            //
            //      Sphere_Unit_Area ( N ) = 2 * PI**(N/2) / Gamma ( N / 2 )
            //
            //    Some values of the function include:
            //
            //       N   Area
            //
            //       2    2        * PI
            //       3  ( 4 /    ) * PI
            //       4  ( 2 /   1) * PI^2
            //       5  ( 8 /   3) * PI^2
            //       6  ( 1 /   1) * PI^3
            //       7  (16 /  15) * PI^3
            //       8  ( 1 /   3) * PI^4
            //       9  (32 / 105) * PI^4
            //      10  ( 1 /  12) * PI^5
            //
            //    For the unit sphere, Area(N) = N * Volume(N)
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      2 * Pi^(n/2) / Gamma[n/2]
            //
            //  Modified:
            //
            //    20 August 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Stephen Wolfram,
            //    The Mathematica Book,
            //    Fourth Edition,
            //    Cambridge University Press, 1999,
            //    ISBN: 0-521-64314-7,
            //    LC: QA76.95.W65.
            //
            //  Parameters:
            //
            //    Input/output, int *N_DATA.
            //    On input, if N_DATA is 0, the first test data is returned, and
            //    N_DATA is set to the index of the test data.  On each subsequent
            //    call, N_DATA is incremented and that test data is returned.  When
            //    there is no more test data, N_DATA is set to 0.
            //
            //    Output, int *N, the spatial dimension.
            //
            //    Output, double *AREA, the area of the unit sphere 
            //    in that dimension.
            //
        {
            int N_MAX = 20;

            double[] area_vec =
            {
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
            };

            int[] n_vec =
            {
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
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                n = 0;
                area = 0.0;
            }
            else
            {
                n = n_vec[n_data - 1];
                area = area_vec[n_data - 1];
            }

        }

        public static double sphere_unit_monomial_nd(int n, int[] p)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_UNIT_MONOMIAL_ND integrates a monomial on the surface of the unit sphere in ND.
            //
            //  Integration region:
            //
            //    sum ( X(1:N)^2 ) == 1
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Gerald Folland,
            //    How to Integrate a Polynomial Over a Sphere,
            //    American Mathematical Monthly,
            //    Volume 108, May 2001, pages 446-448.
            //
            //  Parameters:
            //
            //    Input, int N, the dimension of the space.
            //
            //    Input, int P[N], the exponents of X(1) through X(N) in the monomial.
            //    The exponents P(N) must be nonnegative.
            //
            //    Output, double SPHERE_UNIT_MONOMIAL_ND, the integral of
            //    X1**P(1)*X2**P(2)*...*XN**P(N) over the unit sphere.
            //
        {
            double arg1;
            double arg2;
            int i;
            double temp;
            double value;

            for (i = 0; i < n; i++)
            {
                if (p[i] % 2 == 1)
                {
                    value = 0.0;
                    return value;
                }
            }

            temp = 0.0;
            arg2 = 0.0;

            for (i = 0; i < n; i++)
            {
                arg1 = (double)(p[i] + 1) / 2.0;
                temp = temp + typeMethods.r8_gamma_log(arg1);
                arg2 = arg2 + arg1;
            }

            temp = temp - typeMethods.r8_gamma_log(arg2);

            value = 2.0 * Math.Exp(temp);

            return value;
        }

        public static double sphere_unit_volume_nd(int dim_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_UNIT_VOLUME_ND computes the volume of a unit sphere in ND.
            //
            //  Discussion:
            //
            //    The unit sphere in ND satisfies:
            //
            //      sum ( 1 <= I <= DIM_NUM ) X(I) * X(I) = 1
            //
            //    Results for the first few values of DIM_NUM are:
            //
            //     DIM_NUM  Volume
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
            //    For the unit sphere, Volume(DIM_NUM) = 2 * PI * Volume(DIM_NUM-2)/ DIM_NUM
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    12 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the dimension of the space.
            //
            //    Output, double SPHERE_UNIT_VOLUME_ND, the volume of the sphere.
            //
        {
            int i;
            int m;
            double pi = 3.141592653589793;
            double volume;

            if ((dim_num % 2) == 0)
            {
                m = dim_num / 2;
                volume = Math.Pow(pi, m);
                for (i = 1; i <= m; i++)
                {
                    volume = volume / (double)(i);
                }
            }
            else
            {
                m = (dim_num - 1) / 2;
                volume = Math.Pow(pi, m) * Math.Pow(2.0, dim_num);
                for (i = m + 1; i <= 2 * m + 1; i++)
                {
                    volume = volume / (double)(i);
                }
            }

            return volume;
        }

        public static void sphere_unit_volume_values(ref int n_data, ref int n, ref double volume)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_UNIT_VOLUME_VALUES returns some volumes of the unit sphere in ND.
            //
            //  Discussion:
            //
            //    The formula for the volume of the unit sphere in N dimensions is
            //
            //      Volume(N) = 2 * PI**(N/2) / ( N * Gamma ( N / 2 ) )
            //
            //    This function satisfies the relationships:
            //
            //      Volume(N) = 2 * PI * Volume(N-2) / N
            //      Volume(N) = Area(N) / N
            //
            //    Some values of the function include:
            //
            //       N  Volume
            //
            //       1    1
            //       2    1        * PI
            //       3  ( 4 /   3) * PI
            //       4  ( 1 /   2) * PI^2
            //       5  ( 8 /  15) * PI^2
            //       6  ( 1 /   6) * PI^3
            //       7  (16 / 105) * PI^3
            //       8  ( 1 /  24) * PI^4
            //       9  (32 / 945) * PI^4
            //      10  ( 1 / 120) * PI^5
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      2 * Pi^(n/2) / ( n * Gamma[n/2] )
            //
            //  Modified:
            //
            //    21 August 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Stephen Wolfram,
            //    The Mathematica Book,
            //    Fourth Edition,
            //    Cambridge University Press, 1999,
            //    ISBN: 0-521-64314-7,
            //    LC: QA76.95.W65.
            //
            //  Parameters:
            //
            //    Input/output, int *N_DATA.
            //    On input, if N_DATA is 0, the first test data is returned, and
            //    N_DATA is set to the index of the test data.  On each subsequent
            //    call, N_DATA is incremented and that test data is returned.  When
            //    there is no more test data, N_DATA is set to 0.
            //
            //    Output, int *N, the spatial dimension.
            //
            //    Output, double *VOLUME, the volume of the unit 
            //    sphere in that dimension.
            //
        {
            int N_MAX = 20;

            int[] n_vec =
            {
                1, 2,
                3, 4,
                5, 6,
                7, 8,
                9, 10,
                11, 12,
                13, 14,
                15, 16,
                17, 18,
                19, 20
            };

            double[] volume_vec =
            {
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
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                n = 0;
                volume = 0.0;
            }
            else
            {
                n = n_vec[n_data - 1];
                volume = volume_vec[n_data - 1];
            }
        }

        public static double sphere_volume_2d(double r)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_VOLUME_2D computes the volume of an implicit sphere in 2D.
            //
            //  Discussion:
            //
            //    An implicit sphere in 2D satisfies the equation:
            //
            //      sum ( ( P(1:DIM_NUM) - CENTER(1:DIM_NUM) )^2 ) = R * R
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    12 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double R, the radius of the sphere.
            //
            //    Output, double SPHERE_VOLUME_2D, the volume of the sphere.
            //
        {
            double pi = 3.141592653589793;
            double volume;

            volume = pi * r * r;

            return volume;
        }

        public static double sphere_volume_3d(double r)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_VOLUME_3D computes the volume of an implicit sphere in 3D.
            //
            //  Discussion:
            //
            //    An implicit sphere in 3D satisfies the equation:
            //
            //      sum ( ( P(1:DIM_NUM) - CENTER(1:DIM_NUM) )^2 ) = R * R
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    15 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double R, the radius of the sphere.
            //
            //    Output, double SPHERE_VOLUME_3D, the volume of the sphere.
            //
        {
            double pi = 3.141592653589793;
            double volume;

            volume = (4.0 / 3.0) * pi * r * r * r;

            return volume;
        }

        public static double sphere_volume_nd(int dim_num, double r)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_VOLUME_ND computes the volume of an implicit sphere in ND.
            //
            //  Discussion:
            //
            //    An implicit sphere in ND satisfies the equation:
            //
            //      sum ( ( X(1:N) - CENTER(1:N) )^2 ) = R * R
            //
            //    where R is the radius and CENTER is the center.
            //
            //    Results for the first few values of N are:
            //
            //    DIM_NUM  Volume
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
            //    15 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the dimension of the space.
            //
            //    Input, double R, the radius of the sphere.
            //
            //    Output, double SPHERE_VOLUME_ND, the volume of the sphere.
            //
        {
            double volume;

            volume = Math.Pow(r, dim_num) * sphere_unit_volume_nd(dim_num);

            return volume;
        }


    }
}