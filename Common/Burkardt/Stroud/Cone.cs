using System;

namespace Burkardt.Stroud;

public static class Cone
{
    public static double cone_unit_3d(int settings, Func<int, double, double, double, double> func)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CONE_UNIT_3D approximates an integral inside the unit cone in 3D.
        //
        //  Integration Region:
        //
        //      X^2 + Y^2 <= 1 - Z  
        //
        //    and
        //
        //      0 <= Z <= 1.
        //
        //  Discussion:
        //
        //    An 48 point degree 7 formula, Stroud CN:S2:7-1, is used.
        //
        //    (There is a typographical error in the S2:7-1 formula for B3.)
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
        //    Arthur Stroud,
        //    Approximate Calculation of Multiple Integrals,
        //    Prentice Hall, 1971,
        //    ISBN: 0130438936,
        //    LC: QA311.S85.
        //
        //  Parameters:
        //
        //    Input, Func< double, double, double, double > func, the name of the 
        //    user supplied function which evaluates F(X,Y,Z).
        //
        //    Output, double CONE_UNIT_3D, the approximate integral of the function.
        //
    {
        int i;
        double[] u =
        {
            0.04850054945, 0.2386007376,
            0.5170472951, 0.7958514179
        };
        double[] w1 =
        {
            0.1108884156, 0.1434587878,
            0.06863388717, 0.01035224075
        };
        double[] w2 = new double[3];
        double x;
        double y;
        double z;

        double a = Math.Sqrt(3.0) / 2.0;
        double b = Math.Sqrt((27.0 - 3.0 * Math.Sqrt(29.0)) / 104.0);
        double c = Math.Sqrt((27.0 + 3.0 * Math.Sqrt(29.0)) / 104.0);
        w2[0] = 2.0 / 9.0;
        w2[1] = 3.0 * (551.0 + 4.0 * Math.Sqrt(29.0)) / 6264.0;
        w2[2] = 3.0 * (551.0 - 4.0 * Math.Sqrt(29.0)) / 6264.0;

        double quad = 0.0;

        for (i = 0; i < 4; i++)
        {
            x = a * (1.0 - u[i]);
            y = 0.0;
            z = u[i];
            quad += w1[i] * w2[0] * func(settings, x, y, z);

            x = -a * (1.0 - u[i]);
            y = 0.0;
            z = u[i];
            quad += w1[i] * w2[0] * func(settings, x, y, z);

            x = 0.0;
            y = a * (1.0 - u[i]);
            z = u[i];
            quad += w1[i] * w2[0] * func(settings, x, y, z);

            x = 0.0;
            y = -a * (1.0 - u[i]);
            z = u[i];
            quad += w1[i] * w2[0] * func(settings, x, y, z);
        }

        for (i = 0; i < 4; i++)
        {
            x = b * (1.0 - u[i]);
            y = b * (1.0 - u[i]);
            z = u[i];
            quad += w1[i] * w2[1] * func(settings, x, y, z);

            x = -b * (1.0 - u[i]);
            y = b * (1.0 - u[i]);
            z = u[i];
            quad += w1[i] * w2[1] * func(settings, x, y, z);

            x = -b * (1.0 - u[i]);
            y = -b * (1.0 - u[i]);
            z = u[i];
            quad += w1[i] * w2[1] * func(settings, x, y, z);

            x = b * (1.0 - u[i]);
            y = -b * (1.0 - u[i]);
            z = u[i];
            quad += w1[i] * w2[1] * func(settings, x, y, z);

            x = c * (1.0 - u[i]);
            y = c * (1.0 - u[i]);
            z = u[i];
            quad += w1[i] * w2[2] * func(settings, x, y, z);

            x = -c * (1.0 - u[i]);
            y = c * (1.0 - u[i]);
            z = u[i];
            quad += w1[i] * w2[2] * func(settings, x, y, z);

            x = -c * (1.0 - u[i]);
            y = -c * (1.0 - u[i]);
            z = u[i];
            quad += w1[i] * w2[2] * func(settings, x, y, z);

            x = c * (1.0 - u[i]);
            y = -c * (1.0 - u[i]);
            z = u[i];
            quad += w1[i] * w2[2] * func(settings, x, y, z);

        }

        double r = 1.0;
        double h = 1.0;

        double volume = cone_volume_3d(r, h);
        double result = quad * volume;

        return result;
    }

    public static double cone_volume_3d(double r, double h)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CONE_VOLUME_3D returns the volume of a cone in 3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 March 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the base of the cone.
        //
        //    Input, double H, the height of the cone.
        //
        //    Output, double CONE_VOLUME_3D, the volume of the cone.
        //
    {
        double value = Math.PI / 3.0 * h * r * r;

        return value;
    }

}