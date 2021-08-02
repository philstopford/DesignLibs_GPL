using System;

namespace Burkardt
{
    public static class Lambert
    {
        public static double wew_a(double x, ref double en)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    WEW_A estimates Lambert's W function.
            //
            //  Discussion:
            //
            //    For a given X, this routine estimates the solution W of Lambert's 
            //    equation:
            //
            //      X = W * EXP ( W )
            //
            //    This routine has higher accuracy than WEW_B.
            //
            //  Modified:
            //
            //    11 June 2014
            //
            //  Reference:
            //
            //    Fred Fritsch, R Shafer, W Crowley,
            //    Algorithm 443: Solution of the transcendental equation w e^w = x,
            //    Communications of the ACM,
            //    October 1973, Volume 16, Number 2, pages 123-124.
            //
            //  Parameters:
            //
            //    Input, double X, the argument of W(X)
            //
            //    Output, double &EN, the last relative correction to W(X).
            //
            //    Output, double WEW_A, the estimated value of W(X).
            //
        {
            const double c1 = 4.0 / 3.0;
            const double c2 = 7.0 / 3.0;
            const double c3 = 5.0 / 6.0;
            const double c4 = 2.0 / 3.0;
            double f;
            double temp;
            double temp2;
            double wn;
            double y;
            double zn;
            //
            //  Initial guess.
            //
            f = Math.Log(x);

            if (x <= 6.46)
            {
                wn = x * (1.0 + c1 * x) / (1.0 + x * (c2 + c3 * x));
                zn = f - wn - Math.Log(wn);
            }
            else
            {
                wn = f;
                zn = -Math.Log(wn);
            }

            //
            //  Iteration 1.
            //
            temp = 1.0 + wn;
            y = 2.0 * temp * (temp + c4 * zn) - zn;
            wn = wn * (1.0 + zn * y / (temp * (y - zn)));
            //
            //  Iteration 2.
            //
            zn = f - wn - Math.Log(wn);
            temp = 1.0 + wn;
            temp2 = temp + c4 * zn;
            en = zn * temp2 / (temp * temp2 - 0.5 * zn);
            wn = wn * (1.0 + en);

            return wn;
        }

        public static double wew_b(double x, ref double en)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    WEW_B estimates Lambert's W function.
            //
            //  Discussion:
            //
            //    For a given X, this routine estimates the solution W of Lambert's 
            //    equation:
            //
            //      X = W * EXP ( W )
            //
            //    This routine has lower accuracy than WEW_A.
            //
            //  Modified:
            //
            //    11 June 2014
            //
            //  Reference:
            //
            //    Fred Fritsch, R Shafer, W Crowley,
            //    Algorithm 443: Solution of the transcendental equation w e^w = x,
            //    Communications of the ACM,
            //    October 1973, Volume 16, Number 2, pages 123-124.
            //
            //  Parameters:
            //
            //    Input, double X, the argument of W(X)
            //
            //    Output, double &EN, the last relative correction to W(X).
            //
            //    Output, double WEW_B, the estimated value of W(X).
            //
        {
            const double c1 = 4.0 / 3.0;
            const double c2 = 7.0 / 3.0;
            const double c3 = 5.0 / 6.0;
            const double c4 = 2.0 / 3.0;
            double f;
            double temp;
            double wn;
            double y;
            double zn;
            //
            //  Initial guess.
            //
            f = Math.Log(x);

            if (x <= 0.7385)
            {
                wn = x * (1.0 + c1 * x) / (1.0 + x * (c2 + c3 * x));
            }
            else
            {
                wn = f - 24.0 * ((f + 2.0) * f - 3.0)
                    / ((0.7 * f + 58.0) * f + 127.0);
            }

            //
            //  Iteration 1.
            //
            zn = f - wn - Math.Log(wn);
            temp = 1.0 + wn;
            y = 2.0 * temp * (temp + c4 * zn) - zn;
            en = zn * y / (temp * (y - zn));
            wn = wn * (1.0 + en);

            return wn;
        }
    }
}