using System;
using Burkardt.Types;

namespace Burkardt.ODE.RungeKuttaFehlberg
{
    public static partial class RungeKuttaFehlberg
    {
        public static void r4_fehl(Func<float, float[], float[], float[]> f, int neqn,
                float[] y, float t, float h, float[] yp, ref float[] f1, ref float[] f2, ref float[] f3,
                ref float[] f4, ref float[] f5, ref float[] s)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R4_FEHL takes one Fehlberg fourth-fifth order step.
            //
            //  Discussion:
            //
            //    This version of the routine uses FLOAT real arithmetic.
            //
            //    This routine integrates a system of NEQN first order ordinary differential
            //    equations of the form
            //      dY(i)/dT = F(T,Y(1:NEQN))
            //    where the initial values Y and the initial derivatives
            //    YP are specified at the starting point T.
            //
            //    The routine advances the solution over the fixed step H and returns
            //    the fifth order (sixth order accurate locally) solution
            //    approximation at T+H in array S.
            //
            //    The formulas have been grouped to control loss of significance.
            //    The routine should be called with an H not smaller than 13 units of
            //    roundoff in T so that the various independent arguments can be
            //    distinguished.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    27 March 2004
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Herman Watts, Lawrence Shampine.
            //    C++ version by John Burkardt
            //
            //  Reference:
            //
            //    Erwin Fehlberg,
            //    Low-order Classical Runge-Kutta Formulas with Stepsize Control,
            //    NASA Technical Report R-315, 1969.
            //
            //    Lawrence Shampine, Herman Watts, S Davenport,
            //    Solving Non-stiff Ordinary Differential Equations - The State of the Art,
            //    SIAM Review,
            //    Volume 18, pages 376-411, 1976.
            //
            //  Parameters:
            //
            //    Input, external F, a user-supplied subroutine to evaluate the
            //    derivatives Y'(T), of the form:
            //
            //      void f ( double t, double y[], double yp[] )
            //
            //    Input, int NEQN, the number of equations to be integrated.
            //
            //    Input, float Y[NEQN], the current value of the dependent variable.
            //
            //    Input, float T, the current value of the independent variable.
            //
            //    Input, float H, the step size to take.
            //
            //    Input, float YP[NEQN], the current value of the derivative of the
            //    dependent variable.
            //
            //    Output, float F1[NEQN], F2[NEQN], F3[NEQN], F4[NEQN], F5[NEQN], derivative
            //    values needed for the computation.
            //
            //    Output, float S[NEQN], the estimate of the solution at T+H.
            //
        {
            float ch;
            int i;

            ch = h / 4.0f;

            for (i = 0; i < neqn; i++)
            {
                f5[i] = y[i] + ch * yp[i];
            }

            f1 = f(t + ch, f5, f1);

            ch = 3.0f * h / 32.0f;

            for (i = 0; i < neqn; i++)
            {
                f5[i] = y[i] + ch * (yp[i] + 3.0f * f1[i]);
            }

            f2 = f(t + 3.0f * h / 8.0f, f5, f2);

            ch = h / 2197.0f;

            for (i = 0; i < neqn; i++)
            {
                f5[i] = y[i] + ch *
                    (1932.0f * yp[i]
                     + (7296.0f * f2[i] - 7200.0f * f1[i])
                    );
            }

            f3 = f(t + 12.0f * h / 13.0f, f5, f3);

            ch = h / 4104.0f;

            for (i = 0; i < neqn; i++)
            {
                f5[i] = y[i] + ch *
                (
                    (8341.0f * yp[i] - 845.0f * f3[i])
                    + (29440.0f * f2[i] - 32832.0f * f1[i])
                );
            }

            f4 = f(t + h, f5, f4);

            ch = h / 20520.0f;

            for (i = 0; i < neqn; i++)
            {
                f1[i] = y[i] + ch *
                (
                    (-6080.0f * yp[i]
                     + (9295.0f * f3[i] - 5643.0f * f4[i])
                    )
                    + (41040.0f * f1[i] - 28352.0f * f2[i])
                );
            }

            f5 = f(t + h / 2.0f, f1, f5);
            //
            //  Ready to compute the approximate solution at T+H.
            //
            ch = h / 7618050.0f;

            for (i = 0; i < neqn; i++)
            {
                s[i] = y[i] + ch *
                (
                    (902880.0f * yp[i]
                     + (3855735.0f * f3[i] - 1371249.0f * f4[i]))
                    + (3953664.0f * f2[i] + 277020.0f * f5[i])
                );
            }

        }

        public class r4RKFData
        {
            public float abserr_save = -1.0f;
            public int flag_save = -1000;
            public float h = -1.0f;
            public int init = -1000;
            public int kflag = -1000;
            public int kop = -1;
            public float relerr_save = -1.0f;
            public float remin = 1.0E-12f;
            public int nfe = -1;
        }

        public static int r4_rkf45(ref r4RKFData data, Func<float, float[], float[], float[]> f, int neqn,
                ref float[] y, ref float[] yp, ref float t, float tout, ref float relerr, float abserr,
                int flag)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R4_RKF45 carries out the Runge-Kutta-Fehlberg method.
            //
            //  Discussion:
            //
            //    This version of the routine uses FLOAT real arithmetic.
            //
            //    This routine is primarily designed to solve non-stiff and mildly stiff
            //    differential equations when derivative evaluations are inexpensive.
            //    It should generally not be used when the user is demanding
            //    high accuracy.
            //
            //    This routine integrates a system of NEQN first-order ordinary differential
            //    equations of the form:
            //
            //      dY(i)/dT = F(T,Y(1),Y(2),...,Y(NEQN))
            //
            //    where the Y(1:NEQN) are given at T.
            //
            //    Typically the subroutine is used to integrate from T to TOUT but it
            //    can be used as a one-step integrator to advance the solution a
            //    single step in the direction of TOUT.  On return, the parameters in
            //    the call list are set for continuing the integration.  The user has
            //    only to call again (and perhaps define a new value for TOUT).
            //
            //    Before the first call, the user must 
            //
            //    * supply the subroutine F(T,Y,YP) to evaluate the right hand side;
            //      and declare F in an EXTERNAL statement;
            //
            //    * initialize the parameters:
            //      NEQN, Y(1:NEQN), T, TOUT, RELERR, ABSERR, FLAG.
            //      In particular, T should initially be the starting point for integration,
            //      Y should be the value of the initial conditions, and FLAG should 
            //      normally be +1.
            //
            //    Normally, the user only sets the value of FLAG before the first call, and
            //    thereafter, the program manages the value.  On the first call, FLAG should
            //    normally be +1 (or -1 for single step mode.)  On normal return, FLAG will
            //    have been reset by the program to the value of 2 (or -2 in single 
            //    step mode), and the user can continue to call the routine with that 
            //    value of FLAG.
            //
            //    (When the input magnitude of FLAG is 1, this indicates to the program 
            //    that it is necessary to do some initialization work.  An input magnitude
            //    of 2 lets the program know that that initialization can be skipped, 
            //    and that useful information was computed earlier.)
            //
            //    The routine returns with all the information needed to continue
            //    the integration.  If the integration reached TOUT, the user need only
            //    define a new TOUT and call again.  In the one-step integrator
            //    mode, returning with FLAG = -2, the user must keep in mind that 
            //    each step taken is in the direction of the current TOUT.  Upon 
            //    reaching TOUT, indicated by the output value of FLAG switching to 2,
            //    the user must define a new TOUT and reset FLAG to -2 to continue 
            //    in the one-step integrator mode.
            //
            //    In some cases, an error or difficulty occurs during a call.  In that case,
            //    the output value of FLAG is used to indicate that there is a problem
            //    that the user must address.  These values include:
            //
            //    * 3, integration was not completed because the input value of RELERR, the 
            //      relative error tolerance, was too small.  RELERR has been increased 
            //      appropriately for continuing.  If the user accepts the output value of
            //      RELERR, then simply reset FLAG to 2 and continue.
            //
            //    * 4, integration was not completed because more than MAXNFE derivative 
            //      evaluations were needed.  This is approximately (MAXNFE/6) steps.
            //      The user may continue by simply calling again.  The function counter 
            //      will be reset to 0, and another MAXNFE function evaluations are allowed.
            //
            //    * 5, integration was not completed because the solution vanished, 
            //      making a pure relative error test impossible.  The user must use 
            //      a non-zero ABSERR to continue.  Using the one-step integration mode 
            //      for one step is a good way to proceed.
            //
            //    * 6, integration was not completed because the requested accuracy 
            //      could not be achieved, even using the smallest allowable stepsize. 
            //      The user must increase the error tolerances ABSERR or RELERR before
            //      continuing.  It is also necessary to reset FLAG to 2 (or -2 when 
            //      the one-step integration mode is being used).  The occurrence of 
            //      FLAG = 6 indicates a trouble spot.  The solution is changing 
            //      rapidly, or a singularity may be present.  It often is inadvisable 
            //      to continue.
            //
            //    * 7, it is likely that this routine is inefficient for solving
            //      this problem.  Too much output is restricting the natural stepsize
            //      choice.  The user should use the one-step integration mode with 
            //      the stepsize determined by the code.  If the user insists upon 
            //      continuing the integration, reset FLAG to 2 before calling 
            //      again.  Otherwise, execution will be terminated.
            //
            //    * 8, invalid input parameters, indicates one of the following:
            //      NEQN <= 0;
            //      T = TOUT and |FLAG| /= 1;
            //      RELERR < 0 or ABSERR < 0;
            //      FLAG == 0  or FLAG < -2 or 8 < FLAG.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    05 April 2011
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Herman Watts, Lawrence Shampine.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Erwin Fehlberg,
            //    Low-order Classical Runge-Kutta Formulas with Stepsize Control,
            //    NASA Technical Report R-315, 1969.
            //
            //    Lawrence Shampine, Herman Watts, S Davenport,
            //    Solving Non-stiff Ordinary Differential Equations - The State of the Art,
            //    SIAM Review,
            //    Volume 18, pages 376-411, 1976.
            //
            //  Parameters:
            //
            //    Input, external F, a user-supplied subroutine to evaluate the
            //    derivatives Y'(T), of the form:
            //
            //      void f ( float t, float y[], float yp[] )
            //
            //    Input, int NEQN, the number of equations to be integrated.
            //
            //    Input/output, float Y[NEQN], the current solution vector at T.
            //
            //    Input/output, float YP[NEQN], the derivative of the current solution 
            //    vector at T.  The user should not set or alter this information!
            //
            //    Input/output, float[] T, the current value of the independent variable.
            //
            //    Input, float TOUT, the output point at which solution is desired.  
            //    TOUT = T is allowed on the first call only, in which case the routine
            //    returns with FLAG = 2 if continuation is possible.
            //
            //    Input, float[] RELERR, ABSERR, the relative and absolute error tolerances
            //    for the local error test.  At each step the code requires:
            //      Math.Abs ( local error ) <= RELERR * Math.Abs ( Y ) + ABSERR
            //    for each component of the local error and the solution vector Y.
            //    RELERR cannot be "too small".  If the routine believes RELERR has been
            //    set too small, it will reset RELERR to an acceptable value and return
            //    immediately for user action.
            //
            //    Input, int FLAG, indicator for status of integration. On the first call, 
            //    set FLAG to +1 for normal use, or to -1 for single step mode.  On 
            //    subsequent continuation steps, FLAG should be +2, or -2 for single 
            //    step mode.
            //
            //    Output, int RKF45_S, indicator for status of integration.  A value of 2 
            //    or -2 indicates normal progress, while any other value indicates a 
            //    problem that should be addressed.
            //
        {
            int MAXNFE = 3000;

            float ae;
            float dt;
            float ee;
            float eeoet;
            float eps;
            float esttol;
            float et;
            float[] f1;
            float[] f2;
            float[] f3;
            float[] f4;
            float[] f5;
            int flag_return;
            bool hfaild;
            float hmin;
            int i;
            int k;
            int mflag;
            bool output;
            float relerr_min;
            float s;
            float scale;
            float tol;
            float toln;
            float ypk;

            flag_return = flag;
            //
            //  Check the input parameters.
            //
            eps = typeMethods.r4_epsilon();

            if (neqn < 1)
            {
                flag_return = 8;
                return flag_return;
            }

            if ((relerr) < 0.0)
            {
                flag_return = 8;
                return flag_return;
            }

            if (abserr < 0.0)
            {
                flag_return = 8;
                return flag_return;
            }

            if (flag_return == 0 || 8 < flag_return || flag_return < -2)
            {
                flag_return = 8;
                return flag_return;
            }

            mflag = Math.Abs(flag_return);
            //
            //  Is this a continuation call?
            //
            if (mflag != 1)
            {
                if (t == tout && data.kflag != 3)
                {
                    flag_return = 8;
                    return flag_return;
                }

                //
                //  FLAG = -2 or +2:
                //
                if (mflag == 2)
                {
                    if (data.kflag == 3)
                    {
                        flag_return = data.flag_save;
                        mflag = Math.Abs(flag_return);
                    }
                    else if (data.init == 0)
                    {
                        flag_return = data.flag_save;
                    }
                    else if (data.kflag == 4)
                    {
                        data.nfe = 0;
                    }
                    else if (data.kflag == 5 && abserr == 0.0)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("R4_RKF45 - Fatal error!");
                        Console.WriteLine("  KFLAG = 5 and ABSERR = 0.0");
                        return 1;
                    }
                    else if (data.kflag == 6 && (relerr) <= data.relerr_save && abserr <= data.abserr_save)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("R4_RKF45 - Fatal error!");
                        Console.WriteLine("  KFLAG = 6 and");
                        Console.WriteLine("  RELERR <= RELERR_SAVE and");
                        Console.WriteLine("  ABSERR <= ABSERR_SAVE");
                        return (1);
                    }
                }
                //
                //  FLAG = 3, 4, 5, 6, 7 or 8.
                //
                else
                {
                    if (flag_return == 3)
                    {
                        flag_return = data.flag_save;
                        if (data.kflag == 3)
                        {
                            mflag = Math.Abs(flag_return);
                        }
                    }
                    else if (flag_return == 4)
                    {
                        data.nfe = 0;
                        flag_return = data.flag_save;
                        if (data.kflag == 3)
                        {
                            mflag = Math.Abs(flag_return);
                        }
                    }
                    else if (flag_return == 5 && 0.0 < abserr)
                    {
                        flag_return = data.flag_save;
                        if (data.kflag == 3)
                        {
                            mflag = Math.Abs(flag_return);
                        }
                    }
                    //
                    //  Integration cannot be continued because the user did not respond to
                    //  the instructions pertaining to FLAG = 5, 6, 7 or 8.
                    //
                    else
                    {
                        Console.WriteLine("");
                        Console.WriteLine("R4_RKF45 - Fatal error!");
                        Console.WriteLine("  Integration cannot be continued.");
                        Console.WriteLine("  The user did not respond to the output");
                        Console.WriteLine("  value FLAG = 5, 6, 7, or 8.");
                        return (1);
                    }
                }
            }

            //
            //  Save the input value of FLAG.  
            //  Set the continuation flag KFLAG for subsequent input checking.
            //
            data.flag_save = flag_return;
            data.kflag = 0;
            //
            //  Save RELERR and ABSERR for checking input on subsequent calls.
            //
            data.relerr_save = (relerr);
            data.abserr_save = abserr;
            //
            //  Restrict the relative error tolerance to be at least 
            //
            //    2*EPS+REMIN 
            //
            //  to avoid limiting precision difficulties arising from impossible 
            //  accuracy requests.
            //
            relerr_min = 2.0f * typeMethods.r4_epsilon() + data.remin;
            //
            //  Is the relative error tolerance too small?
            //
            if ((relerr) < relerr_min)
            {
                (relerr) = relerr_min;
                data.kflag = 3;
                flag_return = 3;
                return flag_return;
            }

            dt = tout - t;
            //
            //  Initialization:
            //
            //  Set the initialization completion indicator, INIT;
            //  set the indicator for too many output points, KOP;
            //  evaluate the initial derivatives
            //  set the counter for function evaluations, NFE;
            //  estimate the starting stepsize.
            //
            f1 = new float[neqn];
            f2 = new float[neqn];
            f3 = new float[neqn];
            f4 = new float[neqn];
            f5 = new float[neqn];

            if (mflag == 1)
            {
                data.init = 0;
                data.kop = 0;
                yp = f(t, y, yp);
                data.nfe = 1;

                if (t == tout)
                {
                    flag_return = 2;
                    return flag_return;
                }

            }

            if (data.init == 0)
            {
                data.init = 1;
                data.h = Math.Abs(dt);
                toln = 0.0f;

                for (k = 0; k < neqn; k++)
                {
                    tol = (relerr) * Math.Abs(y[k]) + abserr;
                    if (0.0 < tol)
                    {
                        toln = tol;
                        ypk = Math.Abs(yp[k]);
                        if (tol < ypk * Math.Pow(data.h, 5))
                        {
                            data.h = (float) Math.Pow((double) (tol / ypk), 0.2);
                        }
                    }
                }

                if (toln <= 0.0)
                {
                    data.h = 0.0f;
                }

                data.h = (float) Math.Max(data.h, 26.0 * eps * Math.Max(Math.Abs(t), Math.Abs(dt)));

                if (flag_return < 0)
                {
                    data.flag_save = -2;
                }
                else
                {
                    data.flag_save = 2;
                }
            }

            //
            //  Set stepsize for integration in the direction from T to TOUT.
            //
            data.h = typeMethods.r4_sign(dt) * Math.Abs(data.h);
            //
            //  Test to see if too may output points are being requested.
            //
            if (2.0 * Math.Abs(dt) <= Math.Abs(data.h))
            {
                data.kop = data.kop + 1;
            }

            //
            //  Unnecessary frequency of output.
            //
            if (data.kop == 100)
            {
                data.kop = 0;
                flag_return = 7;
                return flag_return;
            }

            //
            //  If we are too close to the output point, then simply extrapolate and return.
            //
            if (Math.Abs(dt) <= 26.0 * eps * Math.Abs(t))
            {
                t = tout;
                for (i = 0; i < neqn; i++)
                {
                    y[i] = y[i] + dt * yp[i];
                }

                yp = f(t, y, yp);
                data.nfe = data.nfe + 1;

                flag_return = 2;
                return flag_return;
            }

            //
            //  Initialize the output point indicator.
            //
            output = false;
            //
            //  To avoid premature underflow in the error tolerance function,
            //  scale the error tolerances.
            //
            scale = 2.0f / (relerr);
            ae = scale * abserr;
            //
            //  Step by step integration.
            //
            for (;;)
            {
                hfaild = false;
                //
                //  Set the smallest allowable stepsize.
                //
                hmin = 26.0f * eps * Math.Abs(t);
                //
                //  Adjust the stepsize if necessary to hit the output point.
                //
                //  Look ahead two steps to avoid drastic changes in the stepsize and
                //  thus lessen the impact of output points on the code.
                //
                dt = tout - t;

                if (2.0 * Math.Abs(data.h) <= Math.Abs(dt))
                {
                }
                else
                    //
                    //  Will the next successful step complete the integration to the output point?
                    //
                {
                    if (Math.Abs(dt) <= Math.Abs(data.h))
                    {
                        output = true;
                        data.h = dt;
                    }
                    else
                    {
                        data.h = 0.5f * dt;
                    }

                }

                //
                //  Here begins the core integrator for taking a single step.
                //
                //  The tolerances have been scaled to avoid premature underflow in
                //  computing the error tolerance function ET.
                //  To avoid problems with zero crossings, relative error is measured
                //  using the average of the magnitudes of the solution at the
                //  beginning and end of a step.
                //  The error estimate formula has been grouped to control loss of
                //  significance.
                //
                //  To distinguish the various arguments, H is not permitted
                //  to become smaller than 26 units of roundoff in T.
                //  Practical limits on the change in the stepsize are enforced to
                //  smooth the stepsize selection process and to avoid excessive
                //  chattering on problems having discontinuities.
                //  To prevent unnecessary failures, the code uses 9/10 the stepsize
                //  it estimates will succeed.
                //
                //  After a step failure, the stepsize is not allowed to increase for
                //  the next attempted step.  This makes the code more efficient on
                //  problems having discontinuities and more effective in general
                //  since local extrapolation is being used and extra caution seems
                //  warranted.
                //
                //  Test the number of derivative function evaluations.
                //  If okay, try to advance the integration from T to T+H.
                //
                for (;;)
                {
                    //
                    //  Have we done too much work?
                    //
                    if (MAXNFE < data.nfe)
                    {
                        data.kflag = 4;
                        flag_return = 4;
                        return flag_return;
                    }

                    //
                    //  Advance an approximate solution over one step of length H.
                    //
                    r4_fehl(f, neqn, y, t, data.h, yp, ref f1, ref f2, ref f3, ref f4, ref f5, ref f1);
                    data.nfe = data.nfe + 5;
                    //
                    //  Compute and test allowable tolerances versus local error estimates
                    //  and remove scaling of tolerances.  The relative error is
                    //  measured with respect to the average of the magnitudes of the
                    //  solution at the beginning and end of the step.
                    //
                    eeoet = 0.0f;

                    for (k = 0; k < neqn; k++)
                    {
                        et = Math.Abs(y[k]) + Math.Abs(f1[k]) + ae;

                        if (et <= 0.0)
                        {
                            flag_return = 5;
                            return flag_return;
                        }

                        ee = (float) Math.Abs
                        ((-2090.0 * yp[k]
                          + (21970.0 * f3[k] - 15048.0 * f4[k])
                         )
                         + (22528.0 * f2[k] - 27360.0 * f5[k])
                        );

                        eeoet = Math.Max(eeoet, ee / et);

                    }

                    esttol = Math.Abs(data.h) * eeoet * scale / 752400.0f;

                    if (esttol <= 1.0)
                    {
                        break;
                    }

                    //
                    //  Unsuccessful step.  Reduce the stepsize, try again.
                    //  The decrease is limited to a factor of 1/10.
                    //
                    hfaild = true;
                    output = false;

                    if (esttol < 59049.0)
                    {
                        s = 0.9f / (float) Math.Pow((double) esttol, 0.2);
                    }
                    else
                    {
                        s = 0.1f;
                    }

                    data.h = s * data.h;

                    if (Math.Abs(data.h) < hmin)
                    {
                        data.kflag = 6;
                        flag_return = 6;
                        return flag_return;
                    }

                }

                //
                //  We exited the loop because we took a successful step.  
                //  Store the solution for T+H, and evaluate the derivative there.
                //
                t = t + data.h;
                for (i = 0; i < neqn; i++)
                {
                    y[i] = f1[i];
                }

                yp = f(t, y, yp);
                data.nfe = data.nfe + 1;
                //
                //  Choose the next stepsize.  The increase is limited to a factor of 5.
                //  If the step failed, the next stepsize is not allowed to increase.
                //
                if (0.0001889568 < esttol)
                {
                    s = 0.9f / (float) Math.Pow((double) esttol, 0.2);
                }
                else
                {
                    s = 5.0f;
                }

                if (hfaild)
                {
                    s = Math.Min(s, 1.0f);
                }

                data.h = typeMethods.r4_sign(data.h) * Math.Max(s * Math.Abs(data.h), hmin);
                //
                //  End of core integrator
                //
                //  Should we take another step?
                //
                if (output)
                {
                    t = tout;
                    flag_return = 2;
                    return flag_return;
                }

                if (flag_return <= 0)
                {
                    flag_return = -2;
                    return flag_return;
                }

            }
        }

        public static void r8_fehl(Func<double, double[], double[], double[]> f, int neqn,
                double[] y, double t, double h, double[] yp, ref double[] f1, ref double[] f2,
                ref double[] f3, ref double[] f4, ref double[] f5, ref double[] s)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_FEHL takes one Fehlberg fourth-fifth order step.
            //
            //  Discussion:
            //
            //    This version of the routine uses DOUBLE real arithemtic.
            //
            //    This routine integrates a system of NEQN first order ordinary differential
            //    equations of the form
            //      dY(i)/dT = F(T,Y(1:NEQN))
            //    where the initial values Y and the initial derivatives
            //    YP are specified at the starting point T.
            //
            //    The routine advances the solution over the fixed step H and returns
            //    the fifth order (sixth order accurate locally) solution
            //    approximation at T+H in array S.
            //
            //    The formulas have been grouped to control loss of significance.
            //    The routine should be called with an H not smaller than 13 units of
            //    roundoff in T so that the various independent arguments can be
            //    distinguished.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    27 March 2004
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Herman Watts, Lawrence Shampine.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Erwin Fehlberg,
            //    Low-order Classical Runge-Kutta Formulas with Stepsize Control,
            //    NASA Technical Report R-315, 1969.
            //
            //    Lawrence Shampine, Herman Watts, S Davenport,
            //    Solving Non-stiff Ordinary Differential Equations - The State of the Art,
            //    SIAM Review,
            //    Volume 18, pages 376-411, 1976.
            //
            //  Parameters:
            //
            //    Input, external F, a user-supplied subroutine to evaluate the
            //    derivatives Y'(T), of the form:
            //
            //      void f ( double t, double y[], double yp[] )
            //
            //    Input, int NEQN, the number of equations to be integrated.
            //
            //    Input, double Y[NEQN], the current value of the dependent variable.
            //
            //    Input, double T, the current value of the independent variable.
            //
            //    Input, double H, the step size to take.
            //
            //    Input, double YP[NEQN], the current value of the derivative of the
            //    dependent variable.
            //
            //    Output, double F1[NEQN], F2[NEQN], F3[NEQN], F4[NEQN], F5[NEQN], derivative
            //    values needed for the computation.
            //
            //    Output, double S[NEQN], the estimate of the solution at T+H.
            //
        {
            double ch;
            int i;

            ch = h / 4.0;

            for (i = 0; i < neqn; i++)
            {
                f5[i] = y[i] + ch * yp[i];
            }

            f1 = f(t + ch, f5, f1);

            ch = 3.0 * h / 32.0;

            for (i = 0; i < neqn; i++)
            {
                f5[i] = y[i] + ch * (yp[i] + 3.0 * f1[i]);
            }

            f2 = f(t + 3.0 * h / 8.0, f5, f2);

            ch = h / 2197.0;

            for (i = 0; i < neqn; i++)
            {
                f5[i] = y[i] + ch *
                    (1932.0 * yp[i]
                     + (7296.0 * f2[i] - 7200.0 * f1[i])
                    );
            }

            f3 = f(t + 12.0 * h / 13.0, f5, f3);

            ch = h / 4104.0;

            for (i = 0; i < neqn; i++)
            {
                f5[i] = y[i] + ch *
                (
                    (8341.0 * yp[i] - 845.0 * f3[i])
                    + (29440.0 * f2[i] - 32832.0 * f1[i])
                );
            }

            f4 = f(t + h, f5, f4);

            ch = h / 20520.0;

            for (i = 0; i < neqn; i++)
            {
                f1[i] = y[i] + ch *
                (
                    (-6080.0 * yp[i]
                     + (9295.0 * f3[i] - 5643.0 * f4[i])
                    )
                    + (41040.0 * f1[i] - 28352.0 * f2[i])
                );
            }

            f5 = f(t + h / 2.0, f1, f5);
            //
            //  Ready to compute the approximate solution at T+H.
            //
            ch = h / 7618050.0;

            for (i = 0; i < neqn; i++)
            {
                s[i] = y[i] + ch *
                (
                    (902880.0 * yp[i]
                     + (3855735.0 * f3[i] - 1371249.0 * f4[i]))
                    + (3953664.0 * f2[i] + 277020.0 * f5[i])
                );
            }
        }

        public class r8RKFData
        {
            public double abserr_save = -1.0;
            public int flag_save = -1000;
            public double h = -1.0;
            public int init = -1000;
            public int kflag = -1000;
            public int kop = -1;
            public int nfe = -1;
            public double relerr_save = -1.0;
            public double remin = 1.0E-12;
        }

        public static int r8_rkf45(ref r8RKFData data, Func<double, double[], double[], double[]> f, int neqn,
                ref double[] y, ref double[] yp, ref double t, double tout, ref double relerr,
                double abserr, int flag)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_RKF45 carries out the Runge-Kutta-Fehlberg method.
            //
            //  Discussion:
            //
            //    This version of the routine uses DOUBLE real arithmetic.
            //
            //    This routine is primarily designed to solve non-stiff and mildly stiff
            //    differential equations when derivative evaluations are inexpensive.
            //    It should generally not be used when the user is demanding
            //    high accuracy.
            //
            //    This routine integrates a system of NEQN first-order ordinary differential
            //    equations of the form:
            //
            //      dY(i)/dT = F(T,Y(1),Y(2),...,Y(NEQN))
            //
            //    where the Y(1:NEQN) are given at T.
            //
            //    Typically the subroutine is used to integrate from T to TOUT but it
            //    can be used as a one-step integrator to advance the solution a
            //    single step in the direction of TOUT.  On return, the parameters in
            //    the call list are set for continuing the integration.  The user has
            //    only to call again (and perhaps define a new value for TOUT).
            //
            //    Before the first call, the user must 
            //
            //    * supply the subroutine F(T,Y,YP) to evaluate the right hand side;
            //      and declare F in an EXTERNAL statement;
            //
            //    * initialize the parameters:
            //      NEQN, Y(1:NEQN), T, TOUT, RELERR, ABSERR, FLAG.
            //      In particular, T should initially be the starting point for integration,
            //      Y should be the value of the initial conditions, and FLAG should 
            //      normally be +1.
            //
            //    Normally, the user only sets the value of FLAG before the first call, and
            //    thereafter, the program manages the value.  On the first call, FLAG should
            //    normally be +1 (or -1 for single step mode.)  On normal return, FLAG will
            //    have been reset by the program to the value of 2 (or -2 in single 
            //    step mode), and the user can continue to call the routine with that 
            //    value of FLAG.
            //
            //    (When the input magnitude of FLAG is 1, this indicates to the program 
            //    that it is necessary to do some initialization work.  An input magnitude
            //    of 2 lets the program know that that initialization can be skipped, 
            //    and that useful information was computed earlier.)
            //
            //    The routine returns with all the information needed to continue
            //    the integration.  If the integration reached TOUT, the user need only
            //    define a new TOUT and call again.  In the one-step integrator
            //    mode, returning with FLAG = -2, the user must keep in mind that 
            //    each step taken is in the direction of the current TOUT.  Upon 
            //    reaching TOUT, indicated by the output value of FLAG switching to 2,
            //    the user must define a new TOUT and reset FLAG to -2 to continue 
            //    in the one-step integrator mode.
            //
            //    In some cases, an error or difficulty occurs during a call.  In that case,
            //    the output value of FLAG is used to indicate that there is a problem
            //    that the user must address.  These values include:
            //
            //    * 3, integration was not completed because the input value of RELERR, the 
            //      relative error tolerance, was too small.  RELERR has been increased 
            //      appropriately for continuing.  If the user accepts the output value of
            //      RELERR, then simply reset FLAG to 2 and continue.
            //
            //    * 4, integration was not completed because more than MAXNFE derivative 
            //      evaluations were needed.  This is approximately (MAXNFE/6) steps.
            //      The user may continue by simply calling again.  The function counter 
            //      will be reset to 0, and another MAXNFE function evaluations are allowed.
            //
            //    * 5, integration was not completed because the solution vanished, 
            //      making a pure relative error test impossible.  The user must use 
            //      a non-zero ABSERR to continue.  Using the one-step integration mode 
            //      for one step is a good way to proceed.
            //
            //    * 6, integration was not completed because the requested accuracy 
            //      could not be achieved, even using the smallest allowable stepsize. 
            //      The user must increase the error tolerances ABSERR or RELERR before
            //      continuing.  It is also necessary to reset FLAG to 2 (or -2 when 
            //      the one-step integration mode is being used).  The occurrence of 
            //      FLAG = 6 indicates a trouble spot.  The solution is changing 
            //      rapidly, or a singularity may be present.  It often is inadvisable 
            //      to continue.
            //
            //    * 7, it is likely that this routine is inefficient for solving
            //      this problem.  Too much output is restricting the natural stepsize
            //      choice.  The user should use the one-step integration mode with 
            //      the stepsize determined by the code.  If the user insists upon 
            //      continuing the integration, reset FLAG to 2 before calling 
            //      again.  Otherwise, execution will be terminated.
            //
            //    * 8, invalid input parameters, indicates one of the following:
            //      NEQN <= 0;
            //      T = TOUT and |FLAG| /= 1;
            //      RELERR < 0 or ABSERR < 0;
            //      FLAG == 0  or FLAG < -2 or 8 < FLAG.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    13 October 2012
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Herman Watts, Lawrence Shampine.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Erwin Fehlberg,
            //    Low-order Classical Runge-Kutta Formulas with Stepsize Control,
            //    NASA Technical Report R-315, 1969.
            //
            //    Lawrence Shampine, Herman Watts, S Davenport,
            //    Solving Non-stiff Ordinary Differential Equations - The State of the Art,
            //    SIAM Review,
            //    Volume 18, pages 376-411, 1976.
            //
            //  Parameters:
            //
            //    Input, external F, a user-supplied subroutine to evaluate the
            //    derivatives Y'(T), of the form:
            //
            //      void f ( double t, double y[], double yp[] )
            //
            //    Input, int NEQN, the number of equations to be integrated.
            //
            //    Input/output, double Y[NEQN], the current solution vector at T.
            //
            //    Input/output, double YP[NEQN], the derivative of the current solution 
            //    vector at T.  The user should not set or alter this information!
            //
            //    Input/output, double[] T, the current value of the independent variable.
            //
            //    Input, double TOUT, the output point at which solution is desired.  
            //    TOUT = T is allowed on the first call only, in which case the routine
            //    returns with FLAG = 2 if continuation is possible.
            //
            //    Input, double[] RELERR, ABSERR, the relative and absolute error tolerances
            //    for the local error test.  At each step the code requires:
            //      Math.Abs ( local error ) <= RELERR * Math.Abs ( Y ) + ABSERR
            //    for each component of the local error and the solution vector Y.
            //    RELERR cannot be "too small".  If the routine believes RELERR has been
            //    set too small, it will reset RELERR to an acceptable value and return
            //    immediately for user action.
            //
            //    Input, int FLAG, indicator for status of integration. On the first call, 
            //    set FLAG to +1 for normal use, or to -1 for single step mode.  On 
            //    subsequent continuation steps, FLAG should be +2, or -2 for single 
            //    step mode.
            //
            //    Output, int RKF45_D, indicator for status of integration.  A value of 2 
            //    or -2 indicates normal progress, while any other value indicates a 
            //    problem that should be addressed.
            //
        {
            int MAXNFE = 3000;

            double ae;
            double dt;
            double ee;
            double eeoet;
            double eps;
            double esttol;
            double et;
            double[] f1;
            double[] f2;
            double[] f3;
            double[] f4;
            double[] f5;
            int flag_return;
            bool hfaild;
            double hmin;
            int i;
            int k;
            int mflag;
            bool output;
            double relerr_min;
            double s;
            double scale;
            double tol;
            double toln;
            double ypk;

            flag_return = flag;
            //
            //  Check the input parameters.
            //
            eps = typeMethods.r8_epsilon();

            if (neqn < 1)
            {
                flag_return = 8;
                Console.WriteLine("");
                Console.WriteLine("R8_RKF45 - Fatal error!");
                Console.WriteLine("  Invalid input value of NEQN.");
                return flag_return;
            }

            if ((relerr) < 0.0)
            {
                flag_return = 8;
                Console.WriteLine("");
                Console.WriteLine("R8_RKF45 - Fatal error!");
                Console.WriteLine("  Invalid input value of RELERR.");
                return flag_return;
            }

            if (abserr < 0.0)
            {
                flag_return = 8;
                Console.WriteLine("");
                Console.WriteLine("R8_RKF45 - Fatal error!");
                Console.WriteLine("  Invalid input value of ABSERR.");
                return flag_return;
            }

            if (flag_return == 0 || 8 < flag_return || flag_return < -2)
            {
                flag_return = 8;
                Console.WriteLine("");
                Console.WriteLine("R8_RKF45 - Fatal error!");
                Console.WriteLine("  Invalid input.");
                return flag_return;
            }

            mflag = Math.Abs(flag_return);
            //
            //  Is this a continuation call?
            //
            if (mflag != 1)
            {
                if (t == tout && data.kflag != 3)
                {
                    flag_return = 8;
                    return flag_return;
                }

                //
                //  FLAG = -2 or +2:
                //
                if (mflag == 2)
                {
                    if (data.kflag == 3)
                    {
                        flag_return = data.flag_save;
                        mflag = Math.Abs(flag_return);
                    }
                    else if (data.init == 0)
                    {
                        flag_return = data.flag_save;
                    }
                    else if (data.kflag == 4)
                    {
                        data.nfe = 0;
                    }
                    else if (data.kflag == 5 && abserr == 0.0)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("R8_RKF45 - Fatal error!");
                        Console.WriteLine("  KFLAG = 5 and ABSERR = 0.0");
                        return (1);
                    }
                    else if (data.kflag == 6 && (relerr) <= data.relerr_save && abserr <= data.abserr_save)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("R8_RKF45 - Fatal error!");
                        Console.WriteLine("  KFLAG = 6 and");
                        Console.WriteLine("  RELERR <= RELERR_SAVE and");
                        Console.WriteLine("  ABSERR <= ABSERR_SAVE");
                        return (1);
                    }
                }
                //
                //  FLAG = 3, 4, 5, 6, 7 or 8.
                //
                else
                {
                    if (flag_return == 3)
                    {
                        flag_return = data.flag_save;
                        if (data.kflag == 3)
                        {
                            mflag = Math.Abs(flag_return);
                        }
                    }
                    else if (flag_return == 4)
                    {
                        data.nfe = 0;
                        flag_return = data.flag_save;
                        if (data.kflag == 3)
                        {
                            mflag = Math.Abs(flag_return);
                        }
                    }
                    else if (flag_return == 5 && 0.0 < abserr)
                    {
                        flag_return = data.flag_save;
                        if (data.kflag == 3)
                        {
                            mflag = Math.Abs(flag_return);
                        }
                    }
                    //
                    //  Integration cannot be continued because the user did not respond to
                    //  the instructions pertaining to FLAG = 5, 6, 7 or 8.
                    //
                    else
                    {
                        Console.WriteLine("");
                        Console.WriteLine("R8_RKF45 - Fatal error!");
                        Console.WriteLine("  Integration cannot be continued.");
                        Console.WriteLine("  The user did not respond to the output");
                        Console.WriteLine("  value FLAG = 5, 6, 7, or 8.");
                        return (1);
                    }
                }
            }

            //
            //  Save the input value of FLAG.  
            //  Set the continuation flag KFLAG for subsequent input checking.
            //
            data.flag_save = flag_return;
            data.kflag = 0;
            //
            //  Save RELERR and ABSERR for checking input on subsequent calls.
            //
            data.relerr_save = (relerr);
            data.abserr_save = abserr;
            //
            //  Restrict the relative error tolerance to be at least 
            //
            //    2*EPS+REMIN 
            //
            //  to avoid limiting precision difficulties arising from impossible 
            //  accuracy requests.
            //
            relerr_min = 2.0 * typeMethods.r8_epsilon() + data.remin;
            //
            //  Is the relative error tolerance too small?
            //
            if ((relerr) < relerr_min)
            {
                (relerr) = relerr_min;
                data.kflag = 3;
                flag_return = 3;
                return flag_return;
            }

            dt = tout - t;
            //
            //  Initialization:
            //
            //  Set the initialization completion indicator, INIT;
            //  set the indicator for too many output points, KOP;
            //  evaluate the initial derivatives
            //  set the counter for function evaluations, NFE;
            //  estimate the starting stepsize.
            //
            f1 = new double[neqn];
            f2 = new double[neqn];
            f3 = new double[neqn];
            f4 = new double[neqn];
            f5 = new double[neqn];

            if (mflag == 1)
            {
                data.init = 0;
                data.kop = 0;
                yp = f(t, y, yp);
                data.nfe = 1;

                if (t == tout)
                {
                    flag_return = 2;
                    return flag_return;
                }

            }

            if (data.init == 0)
            {
                data.init = 1;
                data.h = Math.Abs(dt);
                toln = 0.0;

                for (k = 0; k < neqn; k++)
                {
                    tol = (relerr) * Math.Abs(y[k]) + abserr;
                    if (0.0 < tol)
                    {
                        toln = tol;
                        ypk = Math.Abs(yp[k]);
                        if (tol < ypk * Math.Pow(data.h, 5))
                        {
                            data.h = Math.Pow((tol / ypk), 0.2);
                        }
                    }
                }

                if (toln <= 0.0)
                {
                    data.h = 0.0;
                }

                data.h = Math.Max(data.h, 26.0 * eps * Math.Max(Math.Abs(t), Math.Abs(dt)));

                if (flag_return < 0)
                {
                    data.flag_save = -2;
                }
                else
                {
                    data.flag_save = 2;
                }
            }

            //
            //  Set stepsize for integration in the direction from T to TOUT.
            //
            data.h = typeMethods.r8_sign(dt) * Math.Abs(data.h);
            //
            //  Test to see if too may output points are being requested.
            //
            if (2.0 * Math.Abs(dt) <= Math.Abs(data.h))
            {
                data.kop = data.kop + 1;
            }

            //
            //  Unnecessary frequency of output.
            //
            if (data.kop == 100)
            {
                data.kop = 0;
                flag_return = 7;
                return flag_return;
            }

            //
            //  If we are too close to the output point, then simply extrapolate and return.
            //
            if (Math.Abs(dt) <= 26.0 * eps * Math.Abs(t))
            {
                t = tout;
                for (i = 0; i < neqn; i++)
                {
                    y[i] = y[i] + dt * yp[i];
                }

                yp = f(t, y, yp);
                data.nfe = data.nfe + 1;

                flag_return = 2;
                return flag_return;
            }

            //
            //  Initialize the output point indicator.
            //
            output = false;
            //
            //  To avoid premature underflow in the error tolerance function,
            //  scale the error tolerances.
            //
            scale = 2.0 / (relerr);
            ae = scale * abserr;
            //
            //  Step by step integration.
            //
            for (;;)
            {
                hfaild = false;
                //
                //  Set the smallest allowable stepsize.
                //
                hmin = 26.0 * eps * Math.Abs(t);
                //
                //  Adjust the stepsize if necessary to hit the output point.
                //
                //  Look ahead two steps to avoid drastic changes in the stepsize and
                //  thus lessen the impact of output points on the code.
                //
                dt = tout - t;

                if (2.0 * Math.Abs(data.h) <= Math.Abs(dt))
                {
                }
                else
                    //
                    //  Will the next successful step complete the integration to the output point?
                    //
                {
                    if (Math.Abs(dt) <= Math.Abs(data.h))
                    {
                        output = true;
                        data.h = dt;
                    }
                    else
                    {
                        data.h = 0.5 * dt;
                    }

                }

                //
                //  Here begins the core integrator for taking a single step.
                //
                //  The tolerances have been scaled to avoid premature underflow in
                //  computing the error tolerance function ET.
                //  To avoid problems with zero crossings, relative error is measured
                //  using the average of the magnitudes of the solution at the
                //  beginning and end of a step.
                //  The error estimate formula has been grouped to control loss of
                //  significance.
                //
                //  To distinguish the various arguments, H is not permitted
                //  to become smaller than 26 units of roundoff in T.
                //  Practical limits on the change in the stepsize are enforced to
                //  smooth the stepsize selection process and to avoid excessive
                //  chattering on problems having discontinuities.
                //  To prevent unnecessary failures, the code uses 9/10 the stepsize
                //  it estimates will succeed.
                //
                //  After a step failure, the stepsize is not allowed to increase for
                //  the next attempted step.  This makes the code more efficient on
                //  problems having discontinuities and more effective in general
                //  since local extrapolation is being used and extra caution seems
                //  warranted.
                //
                //  Test the number of derivative function evaluations.
                //  If okay, try to advance the integration from T to T+H.
                //
                for (;;)
                {
                    //
                    //  Have we done too much work?
                    //
                    if (MAXNFE < data.nfe)
                    {
                        data.kflag = 4;
                        flag_return = 4;
                        return flag_return;
                    }

                    //
                    //  Advance an approximate solution over one step of length H.
                    //
                    r8_fehl(f, neqn, y, t, data.h, yp, ref f1, ref f2, ref f3, ref f4, ref f5, ref f1);
                    data.nfe = data.nfe + 5;
                    //
                    //  Compute and test allowable tolerances versus local error estimates
                    //  and remove scaling of tolerances.  The relative error is
                    //  measured with respect to the average of the magnitudes of the
                    //  solution at the beginning and end of the step.
                    //
                    eeoet = 0.0;

                    for (k = 0; k < neqn; k++)
                    {
                        et = Math.Abs(y[k]) + Math.Abs(f1[k]) + ae;

                        if (et <= 0.0)
                        {
                            flag_return = 5;
                            return flag_return;
                        }

                        ee = Math.Abs
                        ((-2090.0 * yp[k]
                          + (21970.0 * f3[k] - 15048.0 * f4[k])
                         )
                         + (22528.0 * f2[k] - 27360.0 * f5[k])
                        );

                        eeoet = Math.Max(eeoet, ee / et);

                    }

                    esttol = Math.Abs(data.h) * eeoet * scale / 752400.0;

                    if (esttol <= 1.0)
                    {
                        break;
                    }

                    //
                    //  Unsuccessful step.  Reduce the stepsize, try again.
                    //  The decrease is limited to a factor of 1/10.
                    //
                    hfaild = true;
                    output = false;

                    if (esttol < 59049.0)
                    {
                        s = 0.9 / Math.Pow(esttol, 0.2);
                    }
                    else
                    {
                        s = 0.1;
                    }

                    data.h = s * data.h;

                    if (Math.Abs(data.h) < hmin)
                    {
                        data.kflag = 6;
                        flag_return = 6;
                        return flag_return;
                    }

                }

                //
                //  We exited the loop because we took a successful step.  
                //  Store the solution for T+H, and evaluate the derivative there.
                //
                t = t + data.h;
                for (i = 0; i < neqn; i++)
                {
                    y[i] = f1[i];
                }

                yp = f(t, y, yp);
                data.nfe = data.nfe + 1;
                //
                //  Choose the next stepsize.  The increase is limited to a factor of 5.
                //  If the step failed, the next stepsize is not allowed to increase.
                //
                if (0.0001889568 < esttol)
                {
                    s = 0.9 / Math.Pow(esttol, 0.2);
                }
                else
                {
                    s = 5.0;
                }

                if (hfaild)
                {
                    s = Math.Min(s, 1.0);
                }

                data.h = typeMethods.r8_sign(data.h) * Math.Max(s * Math.Abs(data.h), hmin);
                //
                //  End of core integrator
                //
                //  Should we take another step?
                //
                if (output)
                {
                    t = tout;
                    flag_return = 2;
                    return flag_return;
                }

                if (flag_return <= 0)
                {
                    flag_return = -2;
                    return flag_return;
                }

            }
        }

    }
}