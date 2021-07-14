using System;

namespace Burkardt.Laguerre
{
    public static partial class Integrands
    {
        public static double p05_alpha()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P05_ALPHA returns ALPHA for problem 5.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    28 July 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, double P05_ALPHA, the value of ALPHA.
            //
        {
            double alpha;

            alpha = 2.0;

            return alpha;
        }

        public static double p05_exact()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P05_EXACT returns the estimated integral for problem 5.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    06 October 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, double P05_EXACT, the estimated value of the integral.
            //
        {
            double exact;

            exact = 0.0015897286158592328774;

            return exact;
        }

        public static double[] p05_fun(int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P05_FUN evaluates the integrand for problem 5.
            //
            //  Discussion:
            //
            //    D&R gives "exact" value as  0.00158973...
            //    Mathematica returns         0.0015897286158592328774...
            //    D&R gives Laguerre(16) as  -0.067859545...
            //
            //  Integral:
            //
            //    exp ( -2 ) Integral ( 2 <= x < +oo ) cos ( pi * x^2 / 2 ) dx
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    28 July 2007
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
            //    Input, int N, the number of points.
            //
            //    Input, double X[N], the evaluation points.
            //
            //    Output, double P05_FUN[N], the function values.
            //
        {
            double[] f;
            int i;
            const double r8_pi = 3.1415926535897932385;

            f = new double[n];

            for (i = 0; i < n; i++)
            {
                f[i] = Math.Exp(-2.0) * Math.Cos(0.5 * r8_pi * x[i] * x[i]);
            }

            return f;
        }

        public static string p05_title()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P05_TITLE returns the title for problem 5.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    28 July 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, string P05_TITLE, the title of the problem.
            //
        {
            string title;

            title = "Fresnel integral";

            return title;
        }

    }
}