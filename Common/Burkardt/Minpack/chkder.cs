using System;
using Burkardt.Types;

namespace Burkardt.MinpackNS;

public static partial class Minpack
{
    public static void chkder(int m, int n, double[] x, double[] fvec, double[] fjac,
            int ldfjac, double[] xp, double[] fvecp, int mode, ref double[] err )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    chkder() checks the gradients of M functions in N variables.
        //
        //  Discussion:
        //
        //    This function checks the gradients of M nonlinear functions
        //    in N variables, evaluated at a point x, for consistency with
        //    the functions themselves.  The user must call chkder twice,
        //    first with mode = 1 and then with mode = 2.
        //
        //    mode = 1: on input, x must contain the point of evaluation.
        //    on output, xp is set to a neighboring point.
        //
        //    mode = 2. on input, fvec must contain the functions and the
        //    rows of fjac must contain the gradients
        //    of the respective functions each evaluated
        //    at x, and fvecp must contain the functions
        //    evaluated at xp.
        //    on output, err contains measures of correctness of
        //    the respective gradients.
        //
        //    The function does not perform reliably if cancellation or
        //    rounding errors cause a severe loss of significance in the
        //    evaluation of a function. therefore, none of the components
        //    of x should be unusually small (in particular, zero) or any
        //    other value which may cause loss of significance.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 April 2010
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Jorge More, Burt Garbow, Ken Hillstrom.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Jorge More, Burton Garbow, Kenneth Hillstrom,
        //    User Guide for MINPACK-1,
        //    Technical Report ANL-80-74,
        //    Argonne National Laboratory, 1980.
        //
        //  Parameters:
        //
        //    Input, int M, the number of functions.
        //
        //    Input, int N, the number of variables.
        //
        //    Input, double X[N], the point at which the jacobian is to be checked.
        //
        //    Input, double FVEC(M).  On input, when mode = 2,
        //    this must contain the functions evaluated at x.
        //
        //    Input, double FDJAC(LDFJAC,N).  On input when mode = 2,
        //    the rows of fjac must contain the gradients of
        //    the respective functions evaluated at x.
        //
        //    Input, int LDFJAC, the leading dimension of the array fjac,
        //    which must not be less than M.
        //
        //    Output, double XP[N], created when the input value of MODE is 1.
        //    This is a neighboring point of X.
        //
        //    Input, double FVECP[M].  On input when mode = 2,
        //    fvecp must contain the functions evaluated at xp.
        //
        //    Input, integer MODE.  This should be set to 1 on the first call
        //    and 2 thereafter.
        //
        //    Output, double ERR[M].  On output when mode = 2,
        //    err contains measures of correctness of the respective
        //    gradients.  If there is no severe loss of significance,
        //    then if err(i) is 1.0 the i-th gradient is correct,
        //    while if err(i) is 0.0 the i-th gradient is incorrect.
        //    For values of err between 0.0 and 1.0, the categorization
        //    is less certain.  In general, a value of err(i) greater
        //    than 0.5 indicates that the i-th gradient is probably
        //    correct, while a value of err(i) less than 0.5 indicates
        //    that the i-th gradient is probably incorrect.
        //
    {
        const double factor = 100.0;
        int j;
        double temp;
        //
        //  EPSMCH is the machine precision.
        //
        double epsmch = typeMethods.r8_epsilon();
        //
        double eps = Math.Sqrt(epsmch);
        switch (mode)
        {
            //
            //  MODE = 1.
            //
            case 1:
            {
                for (j = 0; j < n; j++)
                {
                    temp = x[j] switch
                    {
                        0.0 => eps,
                        _ => eps * Math.Abs(x[j])
                    };

                    xp[j] = x[j] + temp;
                }

                break;
            }
            //
            default:
            {
                double epsf = factor * epsmch;
                double epslog = Math.Log10(eps);
                int i;
                for (i = 0; i < m; i++)
                {
                    err[i] = 0.0;
                }

                for (j = 0; j < n; j++)
                {
                    temp = x[j] switch
                    {
                        0.0 => 1.0,
                        _ => Math.Abs(x[j])
                    };

                    for (i = 0; i < m; i++)
                    {
                        err[i] += temp * fjac[i + j * ldfjac];
                    }
                }

                for (i = 0; i < m; i++)
                {
                    if (fvec[i] == 0.0 || fvecp[i] == 0.0 ||
                        !(epsf * Math.Abs(fvec[i]) <= Math.Abs(fvecp[i] - fvec[i])))
                    {
                        continue;
                    }

                    temp = eps * Math.Abs((fvecp[i] - fvec[i]) / eps - err[i])
                           / (Math.Abs(fvec[i]) + Math.Abs(fvecp[i]));

                    if (temp <= epsmch)
                    {
                        err[i] = 1.0;
                    }
                    else if (temp < eps)
                    {
                        err[i] = (Math.Log10(temp) - epslog) / epslog;
                    }
                    else
                    {
                        err[i] = 0.0;
                    }
                }

                break;
            }
        }
    }
}