using System;
using Burkardt.Types;

namespace Burkardt.Edge;

public static class TwoD
{
    public static double fxy1(double x, double y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FXY1 is the first 2D example, scalar version.
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
    {
        double value = 0;
        double[] value_vec;
        double[] x_vec = new double[1];
        double[] y_vec = new double[1];

        x_vec[0] = x;
        y_vec[0] = y;
        value_vec = fxy1_vec(1, x_vec, y_vec);
        value = value_vec[0];

        return value;
    }

    public static double fxy2(double x, double y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FXY2 is the second 2D example, scalar version.
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
    {
        double value = 0;
        double[] value_vec;
        double[] x_vec = new double[1];
        double[] y_vec = new double[1];

        x_vec[0] = x;
        y_vec[0] = y;
        value_vec = fxy2_vec(1, x_vec, y_vec);
        value = value_vec[0];

        return value;
    }

    public static double fxy3(double x, double y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FXY3 is the third 2D example, scalar version.
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
        //    19 September 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
    {
        double value = 0;
        double[] value_vec;
        double[] x_vec = new double[1];
        double[] y_vec = new double[1];

        x_vec[0] = x;
        y_vec[0] = y;
        value_vec = fxy3_vec(1, x_vec, y_vec);
        value = value_vec[0];

        return value;
    }

    public static double fxy4(double x, double y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FXY4 is the fourth 2D example, scalar version.
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
        //    21 September 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
    {
        double value = 0;
        double[] value_vec;
        double[] x_vec = new double[1];
        double[] y_vec = new double[1];

        x_vec[0] = x;
        y_vec[0] = y;
        value_vec = fxy4_vec(1, x_vec, y_vec);
        value = value_vec[0];

        return value;
    }

    public static double fxy5(double x, double y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FXY5 is the fifth 2D example, scalar version.
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
        //    21 September 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
    {
        double value = 0;
        double[] value_vec;
        double[] x_vec = new double[1];
        double[] y_vec = new double[1];

        x_vec[0] = x;
        y_vec[0] = y;
        value_vec = fxy5_vec(1, x_vec, y_vec);
        value = value_vec[0];

        return value;
    }

    public static double[] fxy1_vec(int n, double[] x, double[] y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FXY1_VEC is the first 2D example, vector version.
        //
        //  Discussion:
        //
        //    This is example 4.1 in the reference.
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
        //    Input, double X[N], Y[N], the arguments.
        //
        //    Output, double FXY1_VEC[N], the function values.
        //
    {
        double[] f;
        int i;
            

        f = new double[n];

        for (i = 0; i < n; i++)
        {
            f[i] = (x[i] * x[i] + y[i] * y[i]) switch
            {
                > 0.25 => f[i] + 10.0 * x[i] - 5.0,
                _ => x[i] * y[i] + Math.Cos(2.0 * Math.PI * x[i] * x[i]) - Math.Sin(2.0 * Math.PI * x[i] * x[i])
            };
        }

        return f;
    }

    public static double[] fxy2_vec(int n, double[] x, double[] y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FXY2_VEC is the second 2D example, vector version.
        //
        //  Discussion:
        //
        //    This is example 4.2 in the reference.
        //
        //    It is known as the Shepp-Logan phantom.
        //
        //    It should be plotted on [-1,+1] x [-1,+1].
        //
        //    Note that the Archibald reference incorrectly describes the divisor
        //    of x in the second ellipse as 0.06624, when it should be 0.6624.
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
        //    Larry Shepp, Ben Logan,
        //    The Fourier reconstruction of a head section,
        //    IEEE Transactions on Nuclear Science,
        //    Volume  NS-21, June 1974, pages 21-43.
        //
        //  Parameters:
        //
        //    Input, int N, the number of points.
        //
        //    Input, double X[N], Y[N], the arguments.
        //
        //    Output, double FXY2_VEC[N], the function values.
        //
        //  Local parameters:
        //
        //    Local, integer CHOICE:
        //    1, use Archibald's (and Shepp and Logan's) level values;
        //    2, use Matlab's level values;
        //    3, use Matlab's enhanced contrast level values.
        //
    {
        double[] c;
        double[] c1 =
        {
            2.0, -0.98, -0.02, +0.01
        };
        double[] c2 =
        {
            1.0, -0.98, -0.02, +0.01
        };
        double[] c3 =
        {
            1.0, -0.8, -0.2, +0.1
        };
        int choice;
        double eta1;
        double eta2;
        double[] f;
        int i;
            
        double xi1;
        double xi2;

        f = new double[n];

        choice = 3;

        c = choice switch
        {
            1 => typeMethods.r8vec_copy_new(4, c1),
            2 => typeMethods.r8vec_copy_new(4, c2),
            _ => typeMethods.r8vec_copy_new(4, c3)
        };

        for (i = 0; i < n; i++)
        {
            f[i] = 0.0;

            xi1 = (x[i] - 0.22) * Math.Cos(0.4 * Math.PI)
                  + y[i] * Math.Sin(0.4 * Math.PI);
            eta1 = -(x[i] - 0.22) * Math.Sin(0.4 * Math.PI)
                   + y[i] * Math.Cos(0.4 * Math.PI);

            xi2 = (x[i] + 0.22) * Math.Cos(0.6 * Math.PI)
                  + y[i] * Math.Sin(0.6 * Math.PI);
            eta2 = -(x[i] + 0.22) * Math.Sin(0.6 * Math.PI)
                   + y[i] * Math.Cos(0.6 * Math.PI);

            switch (Math.Pow(x[i] / 0.69, 2) + Math.Pow(y[i] / 0.92, 2))
            {
                case <= 1.0:
                    f[i] += c[0];
                    break;
            }

            switch (Math.Pow(x[i] / 0.6624, 2)
                    + Math.Pow((y[i] + 0.0184) / 0.874, 2))
            {
                case <= 1.0:
                    f[i] += c[1];
                    break;
            }

            if (Math.Pow(xi1 / 0.31, 2) + Math.Pow(eta1 / 0.11, 2) <= 1.0 ||
                Math.Pow(xi2 / 0.41, 2) + Math.Pow(eta2 / 0.16, 2) <= 1.0)
            {
                f[i] += c[2];
            }

            if (Math.Pow((x[i] - 0.35) / 0.3, 2)
                + Math.Pow(y[i] / 0.6, 2) <= 1.0 ||
                Math.Pow(x[i] / 0.21, 2)
                + Math.Pow((y[i] - 0.35) / 0.25, 2) <= 1.0 ||
                Math.Pow(x[i] / 0.046, 2)
                + Math.Pow((y[i] - 0.1) / 0.046, 2) <= 1.0 ||
                Math.Pow(x[i] / 0.046, 2)
                + Math.Pow((y[i] + 0.1) / 0.046, 2) <= 1.0 ||
                Math.Pow((x[i] + 0.08) / 0.046, 2)
                + Math.Pow((y[i] + 0.605) / 0.023, 2) <= 1.0 ||
                Math.Pow(x[i] / 0.023, 2)
                + Math.Pow((y[i] + 0.605) / 0.023, 2) <= 1.0 ||
                Math.Pow((x[i] - 0.06) / 0.023, 2)
                + Math.Pow((y[i] + 0.605) / 0.023, 2) <= 1.0)
            {
                f[i] += c[3];
            }
        }

        return f;
    }

    public static double[] fxy3_vec(int n, double[] x, double[] y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FXY3_VEC is the third 2D example, vector version.
        //
        //  Discussion:
        //
        //    This is example 3.2 in the reference.
        //
        //    It is known as the modified two-dimensional Harten function.
        //
        //    It should be plotted on [-1,+1] x [-1,+1].
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 September 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Rick Archibald, Anne Gelb, Jungho Yoon,
        //    Determining the locations and discontinuities in the derivatives
        //    of functions,
        //    Applied Numerical Mathematics,
        //    Volume 58, 2008, pages 577-592.
        //
        //  Parameters:
        //
        //    Input, int N, the number of points.
        //
        //    Input, double X[N], Y[N], the arguments.
        //
        //    Output, double FXY3_VEC[N], the function values.
        //
        //  Local parameters:
        //
        //    Local, integer CHOICE:
        //    1, use Archibald's (and Shepp and Logan's) level values;
        //    2, use Matlab's level values;
        //    3, use Matlab's enhanced contrast level values.
        //
    {
        double[] f;
        int i;
        double r;
            

        f = new double[n];

        for (i = 0; i < n; i++)
        {
            r = (4.0 * x[i] * x[i] + 4.0 * y[i] * y[i] - 1.0) / 6.0;

            f[i] = (3.0 * r) switch
            {
                <= -1.0 => -r * Math.Sin(0.5 * Math.PI * r * r),
                < 1.0 => Math.Abs(Math.Sin(2.0 * Math.PI * r)),
                _ => 2.0 * r - 1.0 - Math.Sin(3.0 * Math.PI * r) / 6.0
            };
        }

        return f;
    }

    public static double[] fxy4_vec(int n, double[] x, double[] y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FXY4_VEC is the fourth 2D example, vector version.
        //
        //  Discussion:
        //
        //    This is example 3.1 in the reference.
        //
        //    It is known as the discontinuous medium wave function.
        //
        //    Here, we are computing the first component of the solution, P(X,Y).
        //
        //    It should be plotted on (x,y) in [-1,0]x[0,0.1].
        //
        //    The second variable y actually represents time.
        //
        //    Note that in the reference, the formula reads:
        //     f(i) = 2.0D+00 * rhor * cr / ( rhol * cl + rhor * cr ) 
        //          *  Math.Sin( Math.PI * omega * ( y(i) - ( x(i) + 0.5D+00 ) / cr ) )
        //    but I believe this should be:
        //     f(i) = 2.0D+00 * rhor * cr / ( rhol * cl + rhor * cr ) 
        //          *  Math.Sin( Math.PI * omega * ( y(i) - ( x(i) + 0.5D+00 ) / cl ) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 September 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Rick Archibald, Anne Gelb, Jungho Yoon,
        //    Determining the locations and discontinuities in the derivatives
        //    of functions,
        //    Applied Numerical Mathematics,
        //    Volume 58, 2008, pages 577-592.
        //
        //  Parameters:
        //
        //    Input, int N, the number of points.
        //
        //    Input, double X[N], Y[N], the arguments.
        //
        //    Output, double FXY3_VEC[N], the function values.
        //
    {
        const double cl = 0.87879;
        const double cr = 1.0;
        double[] f;
        int i;
        const double omega = 12.0;
            
        const double rhol = 0.55556;
        const double rhor = 1.0;

        f = new double[n];

        for (i = 0; i < n; i++)
        {
            f[i] = x[i] switch
            {
                <= -0.5 => Math.Sin(Math.PI * omega * (y[i] - (x[i] + 0.5) / cl)) - (rhol * cl - rhor * cr) /
                    (rhol * cl + rhor * cr) * Math.Sin(Math.PI * omega * (y[i] + (x[i] + 0.5) / cl)),
                _ => 2.0 * rhor * cr / (rhol * cl + rhor * cr) * Math.Sin(Math.PI * omega * (y[i] - (x[i] + 0.5) / cl))
            };
        }

        return f;
    }

    public static double[] fxy5_vec(int n, double[] x, double[] y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FXY5_VEC is the fifth 2D example, vector version.
        //
        //  Discussion:
        //
        //    This is example 3.1 in the reference.
        //
        //    It is known as the discontinuous medium wave function.
        //
        //    Here, we are computing the second component of the solution, U(X,Y).
        //
        //    It should be plotted on (x,y) in [-1,0]x[0,0.1].
        //
        //    The second variable y actually represents time.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 September 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Rick Archibald, Anne Gelb, Jungho Yoon,
        //    Determining the locations and discontinuities in the derivatives
        //    of functions,
        //    Applied Numerical Mathematics,
        //    Volume 58, 2008, pages 577-592.
        //
        //  Parameters:
        //
        //    Input, int N, the number of points.
        //
        //    Input, double X[N], Y[N], the arguments.
        //
        //    Output, double FXY3_VEC[N], the function values.
        //
    {
        const double cl = 0.87879;
        const double cr = 1.0;
        double[] f;
        int i;
        const double omega = 12.0;
            
        const double rhol = 0.55556;
        const double rhor = 1.0;

        f = new double[n];

        for (i = 0; i < n; i++)
        {
            f[i] = x[i] switch
            {
                <= -0.5 => Math.Sin(Math.PI * omega * (y[i] - (x[i] + 0.5) / cl)) + (rhol * cl - rhor * cr) /
                    (rhol * cl + rhor * cr) / (rhol * cl) * Math.Sin(Math.PI * omega * (y[i] + (x[i] + 0.5) / cl)),
                _ => 2.0 / (rhol * cl + rhor * cr) * Math.Sin(Math.PI * omega * (y[i] - (x[i] + 0.5) / cl))
            };
        }

        return f;
    }
}