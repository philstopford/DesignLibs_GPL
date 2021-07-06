using System;

namespace Burkardt.Edge
{
    public static class OneD
    {
        public static double fx1(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FX1 is the first 1D example, scalar version.
            //
            //  Discussion:
            //
            //    This function allows the user a more convenient interface when
            //    only a single input argument is supplied.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 February 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X, the argument.
            //
            //    Output, double FX1, the function value.
            //
        {
            double value;
            double[] value_vec;
            double[] x_vec = new double[1];

            x_vec[0] = x;
            value_vec = fx1_vec(1, x_vec);
            value = value_vec[0];

            return value;
        }

        public static double fx2(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FX2 is the second 1D example, scalar version.
            //
            //  Discussion:
            //
            //    This function allows the user a more convenient interface when
            //    only a single input argument is supplied.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 February 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X, the argument.
            //
            //    Output, double FX2, the function value.
            //
        {
            double value;
            double[] value_vec;
            double[] x_vec = new double[1];

            x_vec[0] = x;
            value_vec = fx2_vec(1, x_vec);
            value = value_vec[0];

            return value;
        }

        public static double fx3(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FX3 is the third 1D example, scalar version.
            //
            //  Discussion:
            //
            //    This function allows the user a more convenient interface when
            //    only a single input argument is supplied.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 February 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X, the argument.
            //
            //    Output, double FX3, the function value.
            //
        {
            double value;
            double[] value_vec;
            double[] x_vec = new double[1];

            x_vec[0] = x;
            value_vec = fx3_vec(1, x_vec);
            value = value_vec[0];

            return value;
        }

        public static double fx4(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FX4 is the fourth 1D example, scalar version.
            //
            //  Discussion:
            //
            //    This function allows the user a more convenient interface when
            //    only a single input argument is supplied.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 February 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X, the argument.
            //
            //    Output, double FX4, the function value.
            //
        {
            double value;
            double[] value_vec;
            double[] x_vec = new double[1];

            x_vec[0] = x;
            value_vec = fx4_vec(1, x_vec);
            value = value_vec[0];

            return value;
        }

        public static double fx5(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FX5 is 1D example #5, scalar version.
            //
            //  Discussion:
            //
            //    This function allows the user a more convenient interface when
            //    only a single input argument is supplied.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 February 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X, the argument.
            //
            //    Output, double FX5, the function value.
            //
        {
            double value;
            double[] value_vec;
            double[] x_vec = new double[1];

            x_vec[0] = x;
            value_vec = fx5_vec(1, x_vec);
            value = value_vec[0];

            return value;
        }

        public static double fx6(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FX6 is 1D example #6, scalar version.
            //
            //  Discussion:
            //
            //    This function allows the user a more convenient interface when
            //    only a single input argument is supplied.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 September 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X, the argument.
            //
            //    Output, double FX6, the function value.
            //
        {
            double value;
            double[] value_vec;
            double[] x_vec = new double[1];

            x_vec[0] = x;
            value_vec = fx6_vec(1, x_vec);
            value = value_vec[0];

            return value;
        }

        public static double fx7(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FX7 is 1D example #7, scalar version.
            //
            //  Discussion:
            //
            //    This function allows the user a more convenient interface when
            //    only a single input argument is supplied.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 September 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X, the argument.
            //
            //    Output, double FX7, the function value.
            //
        {
            double value;
            double[] value_vec;
            double[] x_vec = new double[1];

            x_vec[0] = x;
            value_vec = fx7_vec(1, x_vec);
            value = value_vec[0];

            return value;
        }

        public static double[] fx1_vec(int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FX1_VEC is the first 1D example, vector version.
            //
            //  Discussion:
            //
            //    This is example 3.1 in the reference.
            //
            //    The function should be plotted over [-1.0,+1.0].
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 February 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Rick Archibald, Anne Gelb, Jungho Yoon,
            //    Polynomial fitting for edge detection in irregularly sampled signals 
            //    and images,
            //    SIAM Journal on Numerical Analysis,
            //    Volume 43, Number 1, 2006, pages 259-279.
            //
            //  Parameters:
            //
            //    Input, int N, the number of points.
            //
            //    Input, double X[N], the arguments.
            //
            //    Output, double FX1_VEC[N], the function values.
            //
            //  Local parameters:
            //
            //    Local, real STEEP, controls the steepness of the slope.
            //    The default value is a moderate 5.  For a sharp rise, use 25 instead.  
            //
        {
            double[] f;
            int i;
            const double r8_pi = 3.141592653589793;
            const double steep = 5.0;

            f = new double[n];

            for (i = 0; i < n; i++)
            {
                if (x[i] < 0.0)
                {
                    f[i] = Math.Cos(3.0 * r8_pi * x[i]);
                }
                else if (0.0 <= x[i])
                {
                    f[i] = -1.0 + 2.0 / (1.0 + 3.0 * Math.Exp(-steep * (2.0 * x[i] - 1.0)));
                }
            }

            return f;
        }

        public static double[] fx2_vec(int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FX2_VEC is the second 1D example, vector version.
            //
            //  Discussion:
            //
            //    The function should be plotted over [-1,+1].
            //
            //    The "internal" coordinate range will be [-2.0,6.0*pi].
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 February 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Rick Archibald, Anne Gelb, Jungho Yoon,
            //    Polynomial fitting for edge detection in irregularly sampled signals 
            //    and images,
            //    SIAM Journal on Numerical Analysis,
            //    Volume 43, Number 1, 2006, pages 259-279.
            //
            //  Parameters:
            //
            //    Input, int N, the number of points.
            //
            //    Input, double X[N], the arguments.
            //
            //    Output, double FX2_VEC[N], the function values.
            //
        {
            double[] f;
            int i;
            const double r8_pi = 3.141592653589793;
            double x2;

            f = new double[n];
            //
            //  Map from the convenient range [-1,+1] to the physical range [-2,6pi].
            //
            for (i = 0; i < n; i++)
            {
                x2 = ((1.0 - x[i]) * (-2.0)
                      + (1.0 + x[i]) * 6.0 * r8_pi)
                     / 2.0;

                if (x2 < 0.0)
                {
                    f[i] = Math.Exp(x2);
                }
                else if (0.0 <= x2 && x2 < 3.0 * r8_pi / 2.0)
                {
                    f[i] = -Math.Exp(-x2);
                }
                else if (3.0 * r8_pi / 2.0 <= x2)
                {
                    f[i] = -1.5 * Math.Sin(x2);
                }
            }

            return f;
        }

        public static double[] fx3_vec(int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FX3_VEC is the third 1D example, vector version.
            //
            //  Discussion:
            //
            //    The function should be plotted over [-1.0,+1.0].
            //
            //    Internally, this range is mapped to [-3.0,+3.0].
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 February 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Rick Archibald, Anne Gelb, Jungho Yoon,
            //    Polynomial fitting for edge detection in irregularly sampled signals 
            //    and images,
            //    SIAM Journal on Numerical Analysis,
            //    Volume 43, Number 1, 2006, pages 259-279.
            //
            //  Parameters:
            //
            //    Input, int N, the number of points.
            //
            //    Input, double X[N], the arguments.
            //
            //    Output, double FX3_VEC[N], the function values.
            //
        {
            double[] f;
            int i;
            double x2;

            f = new double[n];
            //
            //  Map from the convenient range [-1,+1] to the physical range [-3,+3].
            //
            for (i = 0; i < n; i++)
            {
                x2 = ((1.0 - x[i]) * (-3.0)
                      + (1.0 + x[i]) * (+3.0))
                     / 2.0;

                if (-2.0 <= x2 && x2 <= -1.0)
                {
                    f[i] = 1.0;
                }
                else if (-0.5 <= x2 && x2 <= 0.5)
                {
                    f[i] = 0.5 + 4.0 * Math.Pow(x2 + 0.5, 2);
                }
                else if (1.25 <= x2 && 3.0 * x2 <= 7.0)
                {
                    f[i] = 3.0 * (2.0 - x2);
                }
                else
                {
                    f[i] = 0.0;
                }
            }

            return f;
        }

        public static double[] fx4_vec(int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FX4_VEC is the fourth 1D example, vector version.
            //
            //  Discussion:
            //
            //    The function should be plotted over [0.0,+1.0].
            //
            //    The function is continuous, but the derivative has a discontinuity at 0.5.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 February 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Rick Archibald, Anne Gelb, Jungho Yoon,
            //    Polynomial fitting for edge detection in irregularly sampled signals 
            //    and images,
            //    SIAM Journal on Numerical Analysis,
            //    Volume 43, Number 1, 2006, pages 259-279.
            //
            //  Parameters:
            //
            //    Input, int N, the number of points.
            //
            //    Input, double X[N], the arguments.
            //
            //    Output, double FX4_VEC[N], the function values.
            //
        {
            double[] f;
            int i;
            const double r8_pi = 3.141592653589793;
            double x2;

            f = new double[n];
            //
            //  Convert from -1 <= x <= 1 to 0 <= x <= 1:
            //
            for (i = 0; i < n; i++)
            {
                x2 = (x[i] + 1.0) / 2.0;

                if (x2 <= 0.5)
                {
                    f[i] = -(x2 - 0.5) + Math.Sin(4.0 * r8_pi * x2) / 6.0;
                }
                else if (0.5 < x2)
                {
                    f[i] = (x2 - 0.5) + Math.Sin(4.0 * r8_pi * x2) / 6.0;
                }
            }

            return f;
        }

        public static double[] fx5_vec(int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FX5_VEC is 1D example #5, vector version.
            //
            //  Discussion:
            //
            //    The function should be plotted over [-1.0,+1.0].
            //
            //    The function actually has no discontinuities, but does have a
            //    steep rise.  The local parameter S controls the steepness of the rise.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 February 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Rick Archibald, Anne Gelb, Jungho Yoon,
            //    Polynomial fitting for edge detection in irregularly sampled signals 
            //    and images,
            //    SIAM Journal on Numerical Analysis,
            //    Volume 43, Number 1, 2006, pages 259-279.
            //
            //  Parameters:
            //
            //    Input, int N, the number of points.
            //
            //    Input, double X[N], the arguments.
            //
            //    Output, double FX5_VEC[N], the function values.
            //
        {
            double[] f;
            int i;
            const double steep = 20.0;

            f = new double[n];

            for (i = 0; i < n; i++)
            {
                f[i] = Math.Tanh(steep * x[i]);
            }

            return f;
        }

        public static double[] fx6_vec(int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FX6_VEC is 1D example #6, vector version.
            //
            //  Discussion:
            //
            //    This is example 2.1 in the reference.
            //
            //    The function should be plotted over [0.0,+1.0].
            //
            //    The function has a discontinuous first derivative at 1/2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 September 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Rick Archibald, Anne Gelb, Jungho Yoon,
            //    Determining the location of discontinuities in the derivatives
            //    of functions,
            //    Applied Numerical Mathematics,
            //    Volume 58, 2008, pages 577-592.
            //
            //  Parameters:
            //
            //    Input, int N, the number of points.
            //
            //    Input, double X[N], the arguments.
            //
            //    Output, double FX6_VEC[N], the function values.
            //
        {
            double[] f;
            int i;
            const double r8_pi = 3.141592653589793;

            f = new double[n];

            for (i = 0; i < n; i++)
            {
                f[i] = Math.Sin(2.0 * r8_pi * x[i]) / 6.0;

                if (x[i] < 0.5)
                {
                    f[i] = f[i] - (x[i] - 0.5);
                }
                else
                {
                    f[i] = f[i] + (x[i] - 0.5);
                }
            }

            return f;
        }
        //****************************************************************************80

        public static double[] fx7_vec(int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FX7_VEC is 1D example #7, vector version.
            //
            //  Discussion:
            //
            //    This is example 2.1 in the reference.
            //
            //    The function should be plotted over [0.0,+1.0].
            //
            //    The function has a discontinuous second derivative at 1/2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 September 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Rick Archibald, Anne Gelb, Jungho Yoon,
            //    Determining the location of discontinuities in the derivatives
            //    of functions,
            //    Applied Numerical Mathematics,
            //    Volume 58, 2008, pages 577-592.
            //
            //  Parameters:
            //
            //    Input, int N, the number of points.
            //
            //    Input, double X[N], the arguments.
            //
            //    Output, double FX6_VEC[N], the function values.
            //
        {
            double[] f;
            int i;
            const double r8_pi = 3.141592653589793;

            f = new double[n];

            for (i = 0; i < n; i++)
            {
                f[i] = Math.Sin(2.0 * r8_pi * x[i]) / 6.0;

                if (x[i] < 0.5)
                {
                    f[i] = f[i] - 0.5 * Math.Pow(x[i] - 0.5, 2);
                }
                else
                {
                    f[i] = f[i] + 0.5 * Math.Pow(x[i] - 0.5, 2);
                }
            }

            return f;
        }

    }
}