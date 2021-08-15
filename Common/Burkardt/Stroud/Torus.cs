using System;
using Burkardt.Quadrature;

namespace Burkardt.Stroud
{
    public static class Torus
    {
        public static double torus_1(Func<double, double, double, double> func, double r1,
                double r2, int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TORUS_1 approximates an integral on the surface of a torus in 3D.
            //
            //  Integration region:
            //
            //    ( Math.Sqrt ( X*X + Y*Y ) - R1 )^2 + Z*Z = R2 * R2.
            //
            //  Discussion:
            //
            //    An (N+1)*(N+2) point N-th degree formula is used.
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
            //    Input, double R1, R2, the two radii that define the torus.
            //
            //    Input, int N, defines the degree of the formula
            //    used to approximate the integral.
            //
            //    Output, double TORUS_1, the approximate integral of the function.
            //
        {
            double angle;
            double ct1;
            int i;
            int j;
            double pi = 3.141592653589793;
            double quad;
            double result;
            double st1;
            double u;
            double volume;
            double w;
            double x;
            double y;
            double z;

            w = 1.0 / (r1 * (double)((n + 1) * (n + 2)));
            quad = 0.0;

            for (i = 0; i < n + 1; i++)
            {
                angle = 2.0 * pi * (double)(i + 1) / (double)(n + 1);
                ct1 = Math.Cos(angle);
                st1 = Math.Sin(angle);

                for (j = 0; j < n + 2; j++)
                {
                    angle = 2.0 * pi * (double)(j + 1) / (double)(n + 2);
                    u = r1 + r2 * Math.Cos(angle);
                    x = u * ct1;
                    y = u * st1;
                    z = r2 * Math.Sin(angle);

                    quad = quad + w * u * func(x, y, z);
                }
            }

            volume = torus_area_3d(r1, r2);
            result = quad * volume;

            return result;
        }

        public static double torus_14s(Func<double, double, double, double> func, double r1,
                double r2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TORUS_14S approximates an integral inside a torus in 3D.
            //
            //  Integration region:
            //
            //    ( Math.Sqrt ( X*X + Y*Y ) - R1 )^2 + Z*Z <= R2 * R2.
            //
            //  Discussion:
            //
            //    A 960 point 14-th degree formula is used.
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
            //    Input, Func< double, double, double, double > func, the name of the
            //    user function which is to be integrated.
            //
            //    Input, double R1, R2, the two radii that define the torus.
            //
            //    Output, double TORUS_14S, the approximate integral of the function.
            //
        {
            double angle;
            double ct;
            double cth;
            int i;
            int j;
            int n;
            int order = 4;
            double pi = 3.141592653589793;
            double quad;
            double[] r =
            {
                0.263499230, 0.574464514, 0.818529487, 0.964659606
            };
            double result;
            double st;
            double sth;
            double u;
            double volume;
            double[] weight =
            {
                0.086963711, 0.163036289, 0.163036289, 0.086963711
            };
            double x;
            double y;
            double z;

            quad = 0.0;

            for (n = 1; n <= 15; n++)
            {
                angle = 2.0 * pi * (double)(n) / 15.0;
                cth = Math.Cos(angle);
                sth = Math.Sin(angle);

                for (i = 1; i <= 16; i++)
                {
                    angle = 2.0 * pi * (double)(i) / 16.0;
                    ct = Math.Cos(angle);
                    st = Math.Sin(angle);

                    for (j = 0; j < order; j++)
                    {
                        u = r1 + r[j] * ct * r2;
                        x = u * cth;
                        y = u * sth;
                        z = r[j] * st * r2;
                        quad = quad + u * weight[j] * func(x, y, z) / (120.0 * r1);
                    }
                }
            }

            volume = torus_volume_3d(r1, r2);
            result = quad * volume;

            return result;
        }

        public static double torus_5s2(Func<double, double, double, double> func, double r1,
                double r2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TORUS_5S2 approximates an integral inside a torus in 3D.
            //
            //  Integration region:
            //
            //    ( Math.Sqrt ( X*X + Y*Y ) - R1 )^2 + Z*Z <= R2 * R2.
            //
            //  Discussion:
            //
            //    A 24 point, 5-th degree formula is used, Stroud number TOR3-S2:5-1.
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
            //    Input, Func< double, double, double, double > func, the name of 
            //    the user supplied function to be integrated.
            //
            //    Input, double R1, R2, the two radii that define the torus.
            //
            //    Output, double TORUS_5S2, the approximate integral of the function.
            //
        {
            double angle;
            double cs;
            int i;
            double pi = 3.141592653589793;
            double quad;
            double result;
            double sn;
            double u1;
            double u2;
            double u3;
            double volume;
            double w;
            double x;
            double y;
            double z;

            w = 1.0 / 24.0;

            quad = 0.0;

            u1 = Math.Sqrt(r1 * r1 + 0.5 * r2 * r2);
            u2 = Math.Sqrt(r1 * r1 + Math.Sqrt(2.0) * r1 * r2 + r2 * r2);
            u3 = Math.Sqrt(r1 * r1 - Math.Sqrt(2.0) * r1 * r2 + r2 * r2);

            for (i = 1; i <= 6; i++)
            {
                angle = 2.0 * pi * (double)(i) / 6.0;
                cs = Math.Cos(angle);
                sn = Math.Sin(angle);

                x = u1 * cs;
                y = u1 * sn;
                z = r2 / Math.Sqrt(2.0);
                quad = quad + w * func(x, y, z);

                x = u1 * cs;
                y = u1 * sn;
                z = -r2 / Math.Sqrt(2.0);
                quad = quad + w * func(x, y, z);

                x = u2 * cs;
                y = u2 * sn;
                z = 0.0;
                quad = quad + w * func(x, y, z);

                x = u3 * cs;
                y = u3 * sn;
                z = 0.0;
                quad = quad + w * func(x, y, z);
            }

            volume = torus_volume_3d(r1, r2);
            result = quad * volume;

            return result;
        }

        public static double torus_6s2(Func<double, double, double, double> func, double r1,
                double r2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TORUS_6S2 approximates an integral inside a torus in 3D.
            //
            //  Integration region:
            //
            //    ( Math.Sqrt ( X*X + Y*Y ) - R1 )^2 + Z*Z <= R2 * R2.
            //
            //  Discussion:
            //
            //    An 84 point 6-th degree formula is used, Stroud number TOR3-S2:6-1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    17 March 2007
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
            //    user defined function to be integrated.
            //
            //    Input, double R1, R2, the two radii that define the torus.
            //
            //    Output, double TORUS_6S2, the approximate integral of the function.
            //
        {
            double cth;
            int i;
            int j;
            int k;
            int n;
            int order = 2;
            double pi = 3.141592653589793;
            double quad;
            double result;
            double[] s = { 0.322914992, 0.644171310 };
            double sth;
            double u;
            double v;
            double volume;
            double w;
            double[] weight = { 0.387077796, 0.165609800 };
            double x;
            double y;
            double z;

            w = 1.0 / (7.0 * r1 * pi);

            quad = 0.0;

            for (n = 1; n <= 7; n++)
            {
                u = 0.5 * Math.Sqrt(3.0) * r2;
                cth = Math.Cos(2.0 * pi * (double)(n) / 7.0);
                sth = Math.Sin(2.0 * pi * (double)(n) / 7.0);

                for (i = 1; i <= 2; i++)
                {
                    u = -u;

                    x = (r1 + u) * cth;
                    y = (r1 + u) * sth;
                    z = 0.0;
                    quad = quad + 0.232710567 * w * (r1 + u) * func(x, y, z);

                    x = r1 * cth;
                    y = r1 * sth;
                    z = u;
                    quad = quad + 0.232710567 * w * r1 * func(x, y, z);

                }

                for (k = 0; k < order; k++)
                {
                    u = s[k] * r2;
                    v = u;

                    for (i = 1; i <= 2; i++)
                    {
                        u = -u;

                        for (j = 1; j <= 2; j++)
                        {
                            v = -v;

                            x = (r1 + u) * cth;
                            y = (r1 + u) * sth;
                            z = v;
                            quad = quad + weight[k] * w * (r1 + u) * func(x, y, z);
                        }
                    }
                }
            }

            volume = torus_volume_3d(r1, r2);
            result = quad * volume;

            return result;
        }

        public static double torus_area_3d(double r1, double r2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TORUS_AREA_3D returns the area of a torus in 3D.
            //
            //  Integration region:
            //
            //    ( Math.Sqrt ( X*X + Y*Y ) - R1 )^2 + Z*Z = R2*R2.
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
            //    Input, double R1, R2, the two radii that define the torus.
            //
            //    Output, double TORUS_AREA_3D, the area of the torus.
            //
        {
            double area;
            double pi = 3.141592653589793;

            area = 4.0 * pi * pi * r1 * r2;

            return area;
        }

        public static double torus_square_14c(Func<double, double, double, double> func,
                double r1, double r2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TORUS_SQUARE_14C approximates an integral in a "square" torus in 3D.
            //
            //  Discussion:
            //
            //    A 14-th degree 960 point formula is used.
            //
            //  Integration region:
            //
            //      R1 - R2 <= Math.Sqrt ( X*X + Y*Y ) <= R1 + R2,
            //    and
            //      -R2 <= Z <= R2.
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
            //    Input, double R1, R2, the radii that define the torus.
            //
            //    Output, double RESULT, the approximate integral of the function.
            //
        {
            double angle;
            double cth;
            int i;
            int j;
            int n;
            int order = 8;
            double pi = 3.141592653589793;
            double quad;
            double result;
            double[] rtab;
            double sth;
            double u;
            double volume;
            double w;
            double[] weight;
            double x;
            double y;
            double z;

            rtab = new double[order];
            weight = new double[order];

            LegendreQuadrature.legendre_set(order, ref rtab, ref weight);

            w = 1.0 / (60.0 * r1);
            quad = 0.0;

            for (n = 1; n <= 15; n++)
            {
                angle = 2.0 * pi * (double)(n) / 15.0;
                cth = Math.Cos(angle);
                sth = Math.Sin(angle);

                for (i = 0; i < order; i++)
                {
                    u = r1 + rtab[i] * r2;
                    x = u * cth;
                    y = u * sth;

                    for (j = 0; j < order; j++)
                    {
                        z = rtab[j] * r2;
                        quad = quad + u * w * weight[i] * weight[j] * func(x, y, z);
                    }
                }
            }

            volume = torus_square_volume_3d(r1, r2);
            result = quad * volume;

            return result;
        }

        public static double torus_square_5c2(Func<double, double, double, double> func,
                double r1, double r2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TORUS_SQUARE_5C2 approximates an integral in a "square" torus in 3D.
            //
            //  Integration region:
            //
            //      R1 - R2 <= Math.Sqrt ( X*X + Y*Y ) <= R1 + R2,
            //    and
            //      -R2 <= Z <= R2.
            //
            //  Discussion:
            //
            //    A 24 point 5-th degree formula is used, Stroud number TOR3-C2:5-1.
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
            //    Input, double R1, the primary radius of the torus.
            //
            //    Input, double R2, one-half the length of a side of the
            //    square cross-section.
            //
            //    Output, double TORUS_SQUARE_5C2, the approximate integral of the function.
            //
        {
            double b1 = 5.0 / 108.0;
            double b2 = 4.0 / 108.0;
            double cs;
            int i;
            double pi = 3.141592653589793;
            double quad;
            double result;
            double sn;
            double u1;
            double u2;
            double u3;
            double v;
            double volume;
            double x;
            double y;
            double z;

            quad = 0.0;

            u1 = Math.Sqrt(r1 * r1 + r2 * r2);

            v = r2 * Math.Sqrt(0.6);

            u2 = Math.Sqrt(r1 * r1 - Math.Sqrt(3.0) * r1 * r2 + r2 * r2);

            u3 = Math.Sqrt(r1 * r1 + Math.Sqrt(3.0) * r1 * r2 + r2 * r2);

            for (i = 1; i <= 6; i++)
            {
                cs = Math.Cos((double)(i) * pi / 3.0);
                sn = Math.Sin((double)(i) * pi / 3.0);

                x = u1 * cs;
                y = u1 * sn;
                z = v;
                quad = quad + b1 * func(x, y, z);

                z = -v;
                quad = quad + b1 * func(x, y, z);

                x = u2 * cs;
                y = u2 * sn;
                z = 0.0;
                quad = quad + b2 * func(x, y, z);

                x = u3 * cs;
                y = u3 * sn;
                z = 0.0;
                quad = quad + b2 * func(x, y, z);
            }

            volume = torus_square_volume_3d(r1, r2);
            result = quad * volume;

            return result;
        }

        public static double torus_square_area_3d(double r1, double r2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TORUS_SQUARE_AREA_3D returns the area of a square torus in 3D.
            //
            //  Integration region:
            //
            //      R1 - R2 <= Math.Sqrt ( X*X + Y*Y ) <= R1 + R2,
            //    and
            //      -R2 <= Z <= R2.
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
            //    Input, double R1, R2, the two radii that define the torus.
            //
            //    Output, double TORUS_SQUARE_AREA_3D, the area of the torus.
            //
        {
            double area;
            double pi = 3.141592653589793;

            area = 16.0 * pi * r1 * r2;

            return area;
        }

        public static double torus_square_volume_3d(double r1, double r2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TORUS_SQUARE_VOLUME_3D returns the volume of a square torus in 3D.
            //
            //  Integration region:
            //
            //      R1 - R2 <= Math.Sqrt ( X*X + Y*Y ) <= R1 + R2,
            //    and
            //      -R2 <= Z <= R2.
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
            //    Input, double R1, R2, the two radii that define the torus.
            //
            //    Output, double TORUS_SQUARE_VOLUME_3D, the volume of the torus.
            //
        {
            double pi = 3.141592653589793;
            double volume;

            volume = 8.0 * pi * r1 * r2 * r2;

            return volume;
        }

        public static double torus_volume_3d(double r1, double r2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TORUS_VOLUME_3D returns the volume of a torus in 3D.
            //
            //  Integration region:
            //
            //    ( Math.Sqrt ( X*X + Y*Y ) - R1 )^2 + Z * Z = R2 * R2.
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
            //    Input, double R1, R2, the two radii that define the torus.
            //
            //    Output, double TORUS_VOLUME_3D, the volume of the torus.
            //
        {
            double pi = 3.141592653589793;
            double volume;

            volume = 2.0 * pi * pi * r1 * r2 * r2;

            return volume;
        }

    }
}