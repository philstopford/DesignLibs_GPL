using System;

namespace Burkardt.CDFLib
{
    public static partial class CDF
    {
        public static double error_f(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ERROR_F evaluates the error function ERF.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 February 2021
            //
            //  Author:
            //
            //    Barry Brown, James Lovato, Kathy Russell.
            //
            //  Parameters:
            //
            //    Input, double *X, the argument.
            //
            //    Output, double ERROR_F, the value of the error function at X.
            //
        {
            double[] a =  {
                .771058495001320e-04,-.133733772997339e-02,.323076579225834e-01,
                .479137145607681e-01,.128379167095513e+00
            }
            ;
            double ax;
            double[] b =  {
                .301048631703895e-02,.538971687740286e-01,.375795757275549e+00
            }
            ;
            double bot;
            double c = .564189583547756e0;
            double erf1;
            double[] p =  {
                -1.36864857382717e-07,5.64195517478974e-01,7.21175825088309e+00,
                4.31622272220567e+01,1.52989285046940e+02,3.39320816734344e+02,
                4.51918953711873e+02,3.00459261020162e+02
            }
            ;
            double[] q =  {
                1.00000000000000e+00,1.27827273196294e+01,7.70001529352295e+01,
                2.77585444743988e+02,6.38980264465631e+02,9.31354094850610e+02,
                7.90950925327898e+02,3.00459260956983e+02
            }
            ;
            double[] r =  {
                2.10144126479064e+00,2.62370141675169e+01,2.13688200555087e+01,
                4.65807828718470e+00,2.82094791773523e-01
            }
            ;
            double[] s =  {
                9.41537750555460e+01,1.87114811799590e+02,9.90191814623914e+01,
                1.80124575948747e+01
            }
            ;
            double t;
            double top;
            double x2;

            ax = Math.Abs(x);

            if (ax <= 0.5e0)
            {
                t = x * x;
                top = (((a[0] * t + a[1]) * t + a[2]) * t + a[3]) * t + a[4] + 1.0e0;
                bot = ((b[0] * t + b[1]) * t + b[2]) * t + 1.0e0;
                erf1 = x * (top / bot);
                return erf1;
            }

            if (ax > 4.0e0) goto S20;
            top = ((((((p[0] * ax + p[1]) * ax + p[2]) * ax + p[3]) * ax + p[4]) * ax + p[5]) * ax + p[6]) * ax + p[
                7];
            bot = ((((((q[0] * ax + q[1]) * ax + q[2]) * ax + q[3]) * ax + q[4]) * ax + q[5]) * ax + q[6]) * ax + q[
                7];
            erf1 = 0.5e0 + (0.5e0 - Math.Exp(-(x * x)) * top / bot);
            if (x < 0.0e0) erf1 = -erf1;
            return erf1;

            S20:
            if (ax >= 5.8e0) goto S30;
            x2 = x * x;
            t = 1.0e0 / x2;
            top = (((r[0] * t + r[1]) * t + r[2]) * t + r[3]) * t + r[4];
            bot = (((s[0] * t + s[1]) * t + s[2]) * t + s[3]) * t + 1.0e0;
            erf1 = (c - top / (x2 * bot)) / ax;
            erf1 = 0.5e0 + (0.5e0 - Math.Exp(-x2) * erf1);
            if (x < 0.0e0) erf1 = -erf1;
            return erf1;

            S30:
            erf1 = fifdsign(1.0e0, x);

            return erf1;
        }
    }
}