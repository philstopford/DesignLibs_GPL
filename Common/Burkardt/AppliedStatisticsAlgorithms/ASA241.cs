using System;
using Burkardt.Types;

namespace Burkardt.AppliedStatistics;

public static partial class Algorithms
{
    public static double r8_normal_01_cdf_inverse(double p)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_NORMAL_01_CDF_INVERSE inverts the standard normal CDF.
        //
        //  Discussion:
        //
        //    The result is accurate to about 1 part in 10**16.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 March 2010
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Michael Wichura.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Michael Wichura,
        //    The Percentage Points of the Normal Distribution,
        //    Algorithm AS 241,
        //    Applied Statistics,
        //    Volume 37, Number 3, pages 477-484, 1988.
        //
        //  Parameters:
        //
        //    Input, double P, the value of the cumulative probability 
        //    densitity function.  0 < P < 1.  If P is outside this range, an "infinite"
        //    value is returned.
        //
        //    Output, double R8_NORMAL_01_CDF_INVERSE, the normal deviate value 
        //    with the property that the probability of a standard normal deviate being 
        //    less than or equal to this value is P.
        //
    {
        double[] a = {
                3.3871328727963666080, 1.3314166789178437745e+2,
                1.9715909503065514427e+3, 1.3731693765509461125e+4,
                4.5921953931549871457e+4, 6.7265770927008700853e+4,
                3.3430575583588128105e+4, 2.5090809287301226727e+3
            }
            ;
        double[] b = {
                1.0, 4.2313330701600911252e+1,
                6.8718700749205790830e+2, 5.3941960214247511077e+3,
                2.1213794301586595867e+4, 3.9307895800092710610e+4,
                2.8729085735721942674e+4, 5.2264952788528545610e+3
            }
            ;
        double[] c = {
                1.42343711074968357734, 4.63033784615654529590,
                5.76949722146069140550, 3.64784832476320460504,
                1.27045825245236838258, 2.41780725177450611770e-1,
                2.27238449892691845833e-2, 7.74545014278341407640e-4
            }
            ;
        const double const1 = 0.180625;
        const double const2 = 1.6;
        double[] d = {
                1.0, 2.05319162663775882187,
                1.67638483018380384940, 6.89767334985100004550e-1,
                1.48103976427480074590e-1, 1.51986665636164571966e-2,
                5.47593808499534494600e-4, 1.05075007164441684324e-9
            }
            ;
        double[] e = {
                6.65790464350110377720, 5.46378491116411436990,
                1.78482653991729133580, 2.96560571828504891230e-1,
                2.65321895265761230930e-2, 1.24266094738807843860e-3,
                2.71155556874348757815e-5, 2.01033439929228813265e-7
            }
            ;
        double[] f = {
                1.0, 5.99832206555887937690e-1,
                1.36929880922735805310e-1, 1.48753612908506148525e-2,
                7.86869131145613259100e-4, 1.84631831751005468180e-5,
                1.42151175831644588870e-7, 2.04426310338993978564e-15
            }
            ;
        double r;
        const double split1 = 0.425;
        const double split2 = 5.0;
        double value;

        switch (p)
        {
            case <= 0.0:
                value = -typeMethods.r8_huge();
                return value;
            case >= 1.0:
                value = typeMethods.r8_huge();
                return value;
        }

        double q = p - 0.5;

        if (Math.Abs(q) <= split1)
        {
            r = const1 - q * q;
            value = q * typeMethods.r8poly_value(8, a, r) / typeMethods.r8poly_value(8, b, r);
        }
        else
        {
            r = q switch
            {
                < 0.0 => p,
                _ => 1.0 - p
            };

            switch (r)
            {
                case <= 0.0:
                    value = -1.0;
                    return 1;
            }

            r = Math.Sqrt(-Math.Log(r));

            if (r <= split2)
            {
                r -= const2;
                value = typeMethods.r8poly_value(8, c, r) / typeMethods.r8poly_value(8, d, r);
            }
            else
            {
                r -= split2;
                value = typeMethods.r8poly_value(8, e, r) / typeMethods.r8poly_value(8, f, r);
            }

            value = q switch
            {
                < 0.0 => -value,
                _ => value
            };
        }

        return value;
    }



    public static float r4_normal_01_cdf_inverse(float p)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R4_NORMAL_01_CDF_INVERSE inverts the standard normal CDF.
        //
        //  Discussion:
        //
        //    The result is accurate to about 1 part in 10**7.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 March 2010
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Michael Wichura.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Michael Wichura,
        //    The Percentage Points of the Normal Distribution,
        //    Algorithm AS 241,
        //    Applied Statistics,
        //    Volume 37, Number 3, pages 477-484, 1988.
        //
        //  Parameters:
        //
        //    Input, float P, the value of the cumulative probability densitity function.
        //    0 < P < 1.  If P is outside this range, an "infinite" result is returned.
        //
        //    Output, float R4_NORMAL_01_CDF_INVERSE, the normal deviate value with the 
        //    property that the probability of a standard normal deviate being less than or
        //    equal to this value is P.
        //
    {
        float[] a = {
            3.3871327179f, 50.434271938f, 159.29113202f, 59.109374720f
        };
        float[] b = {
            1.0f, 17.895169469f, 78.757757664f, 67.187563600f
        };
        float[] c = {
            1.4234372777f, 2.7568153900f, 1.3067284816f, 0.17023821103f
        };
        const float const1 = 0.180625f;
        const float const2 = 1.6f;
        float[] d = {
            1.0f, 0.73700164250f, 0.12021132975f
        };
        float[] e = {
            6.6579051150f, 3.0812263860f, 0.42868294337f, 0.017337203997f
        };
        float[] f = {
            1.0f, 0.24197894225f, 0.012258202635f
        };
        float r;
        const float split1 = 0.425f;
        const float split2 = 5.0f;
        float value;

        if (p <= 0.0)
        {
            value = - typeMethods.r4_huge();
            return value;
        }

        if (1.0 <= p)
        {
            value = typeMethods.r4_huge();
            return value;
        }

        float q = p - 0.5f;

        if (Math.Abs(q) <= split1)
        {
            r = const1 - q * q;
            value = q * typeMethods.r4poly_value(4, a, r) / typeMethods.r4poly_value(4, b, r);
        }
        else
        {
            if (q < 0.0)
            {
                r = p;
            }
            else
            {
                r = 1.0f - p;
            }

            if (r <= 0.0)
            {
                value = -1.0f;
                return 1;
            }

            r = (float)Math.Sqrt(-Math.Log(r));

            if (r <= split2)
            {
                r -= const2;
                value = typeMethods.r4poly_value(4, c, r) / typeMethods.r4poly_value(3, d, r);
            }
            else
            {
                r -= split2;
                value = typeMethods.r4poly_value(4, e, r) / typeMethods.r4poly_value(3, f, r);
            }

            if (q < 0.0)
            {
                value = -value;
            }

        }

        return value;
    }

}