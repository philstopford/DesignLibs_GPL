using System;

namespace Burkardt
{
    public static class Crude
    {

        public class CrudeData
        {
            public double c13 = 0;
            public double em = 0;
            public double em2 = 0;
            public double em9 = 0;
            public double s2 = 0;
            public double s21 = 0;
            public double s22 = 0;
            public double s23 = 0;
            public int init = 0;

        }

        public static double crude(ref CrudeData data, double xx, int nb)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CRUDE returns a crude approximation for the W function.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 June 2014
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Andrew Barry, S. J. Barry, 
            //    Patricia Culligan-Hensley.
            //    This C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Andrew Barry, S. J. Barry, Patricia Culligan-Hensley,
            //    Algorithm 743: WAPR - A Fortran routine for calculating real 
            //    values of the W-function,
            //    ACM Transactions on Mathematical Software,
            //    Volume 21, Number 2, June 1995, pages 172-181.
            //
            //  Parameters:
            //
            //    Input, double XX, the argument.
            //
            //    Input, int NB, indicates the desired branch.
            //    * 0, the upper branch;
            //    * nonzero, the lower branch.
            //
            //    Output, double CRUDE, the crude approximation to W at XX.
            //
        {
            double an2;
            double eta;
            double reta;
            double t;
            double ts;
            double value;
            double zl;

            value = 0.0;
            //
            //  Various mathematical constants.
            //
            if (data.init == 0)
            {
                data.init = 1;
                data.em = -Math.Exp(-1.0);
                data.em9 = -Math.Exp(-9.0);
                data.c13 = 1.0 / 3.0;
                data.em2 = 2.0 / data.em;
                data.s2 = Math.Sqrt(2.0);
                data.s21 = 2.0 * data.s2 - 3.0;
                data.s22 = 4.0 - 3.0 * data.s2;
                data.s23 = data.s2 - 2.0;
            }

            //
            //  Crude Wp.
            //
            if (nb == 0)
            {
                if (xx <= 20.0)
                {
                    reta = data.s2 * Math.Sqrt(1.0 - xx / data.em);
                    an2 = 4.612634277343749 * Math.Sqrt(Math.Sqrt(reta +
                                                                  1.09556884765625));
                    value = reta / (1.0 + reta / (3.0
                                                  + (data.s21 * an2 + data.s22) * reta / (data.s23 * (an2 + reta)))) -
                            1.0;
                }
                else
                {
                    zl = Math.Log(xx);
                    value = Math.Log(xx / Math.Log(xx
                                                   / Math.Pow(zl, Math.Exp(-1.124491989777808 /
                                                                           (0.4225028202459761 + zl)))));
                }
            }
            //
            //  Crude Wm.
            //
            else
            {
                if (xx <= data.em9)
                {
                    zl = Math.Log(-xx);
                    t = -1.0 - zl;
                    ts = Math.Sqrt(t);
                    value = zl - (2.0 * ts) / (data.s2 + (data.c13 - t
                        / (270.0 + ts * 127.0471381349219)) * ts);
                }
                else
                {
                    zl = Math.Log(-xx);
                    eta = 2.0 - data.em2 * xx;
                    value = Math.Log(xx / Math.Log(-xx / ((1.0
                                                           - 0.5043921323068457 * (zl + 1.0))
                        * (Math.Sqrt(eta) + eta / 3.0) + 1.0)));
                }
            }

            return value;
        }
    }
}