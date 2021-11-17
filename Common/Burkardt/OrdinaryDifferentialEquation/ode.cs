using System;
using Burkardt.Types;

namespace Burkardt.ODENS;

public static class ODE
{
    public static void de
        (
            Func<double, double[], int, double[], int, double[]> f,
            int neqn,
            ref double[] y,
            ref int yIndex,
            ref double t,
            double tout,
            double relerr,
            double abserr,
            ref int iflag,
            ref double[] yy,
            ref int yyIndex,
            double[] wt,
            ref int wtIndex,
            ref double[] p,
            ref int pIndex,
            ref double[] yp,
            ref int ypIndex,
            ref double[] ypout,
            ref int ypoutIndex,
            ref double[] phi,
            ref int phiIndex,
            ref double[] alpha,
            ref int alphaIndex,
            ref double[] beta,
            ref int betaIndex,
            ref double[] sig,
            ref int sigIndex,
            ref double[] v,
            ref int vIndex,
            ref double[] w,
            ref int wIndex,
            ref double[] g,
            ref int gIndex,
            ref bool phase1,
            ref double[] psi,
            ref int psiIndex,
            ref double x,
            ref double h,
            ref double hold,
            ref bool start,
            ref double told,
            ref double delsgn,
            ref int ns,
            ref bool nornd,
            ref int k,
            ref int kold,
            ref int isnold
        )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DE carries out the ODE solution algorithm.
        //
        //  Discussion:
        //
        //    ODE merely allocates storage for DE, to relieve the user of the
        //    inconvenience of a long call list.  Consequently, DE is used as
        //    described in the comments for ODE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 February 2012
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Lawrence Shampine, Marilyn Gordon.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Lawrence Shampine, Marilyn Gordon,
        //    Computer Solution of Ordinary Differential Equations:
        //    The Initial Value Problem,
        //    Freeman, 1975,
        //    ISBN: 0716704617,
        //    LC: QA372.S416.
        //
        //  Parameters:
        //
        //    Input, void F ( double T, double Y[], double YP[] ), the user-supplied function
        //    which accepts input values T and Y[], evaluates the right hand
        //    sides of the ODE, and stores the result in YP[].
        //
        //    Input, int NEQN, the number of equations.
        //
        //    Input/output, double Y[NEQN], the current solution.
        //
        //    Input/output, double &T, the current value of the independent
        //    variable.
        //
        //    Input, double TOUT, the desired value of T on output.
        //
        //    Input, double RELERR, ABSERR, the relative and absolute error
        //    tolerances.  At each step, the code requires
        //      abs ( local error ) <= abs ( Y ) * RELERR + ABSERR
        //    for each component of the local error and solution vectors.
        //
        //    Input/output, int &IFLAG, indicates the status of integration.
        //    On input, IFLAG is normally 1 (or -1 in the special case where TOUT is
        //    not to be exceeded.)  On normal output, IFLAG is 2.  Other output values
        //    are:
        //    * 3, integration did not reach TOUT because the error tolerances were
        //         too small.
        //         But RELERR and ABSERR were increased appropriately for continuing;
        //    * 4, integration did not reach TOUT because more than 500 steps were taken;
        //    * 5, integration did not reach TOUT because the equations appear to be
        //         stiff;
        //    * 6, invalid input parameters (fatal error).
        //    The value of IFLAG is returned negative when the input value is negative
        //    and the integration does not reach TOUT.
        //
        //    Workspace, double YY[NEQN], used to hold old solution data.
        //
        //    Input, double WT[NEQN], the error weight vector.
        //
        //    Workspace, double P[NEQN].
        //
        //    Workspace, double YP[NEQN], used to hold values of the
        //    solution derivative.
        //
        //    Workspace, double YPOUT[NEQN], used to hold values of the
        //    solution derivative.
        //
        //    Workspace, double PHI[NEQN*16], contains divided difference
        //    information about the polynomial interpolant to the solution.
        //
        //    Workspace, double ALPHA[12], BETA[12], SIG[13].
        //
        //    Workspace, double V[12], W[12], G[13].
        //
        //    Input/output, bool &PHASE1, indicates whether the program is in the
        //    first phase, when it always wants to increase the ODE method order.
        //
        //    Workspace, double PSI[12], contains information about
        //    the polynomial interpolant to the solution.
        //
        //    Input/output, double &X, a "working copy" of T, the current value
        //    of the independent variable, which is adjusted as the code attempts
        //    to take a step.
        //
        //    Input/output, double &H, the current stepsize.
        //
        //    Input/output, double &HOLD, the last successful stepsize.
        //
        //    Input/output, bool &START, is TRUE on input for the first step.
        //    The program initializes data, and sets START to FALSE.
        //
        //    Input/output, double &TOLD, the previous value of T.
        //
        //    Input/output, double &DELSGN, the sign (+1 or -1) of TOUT - T.
        //
        //    Input/output, int &NS, the number of steps taken with stepsize H.
        //
        //    Input/output, bool &NORND, ?
        //
        //    Input/output, int &K, the order of the current ODE method.
        //
        //    Input/output, int &KOLD, the order of the ODE method on the previous step.
        //
        //    Input/output, int &ISNOLD, the previous value of ISN, the sign
        //    of IFLAG.
        //
        //  Local parameters:
        //
        //    Local, integer MAXNUM, the maximum number of steps allowed in one
        //    call to DE.
        //
    {
        double absdel;
        double abseps;
        bool crash = false;
        double del;
        double eps;
        double fouru;
        int isn;
        int kle4;
        int l;
        const int maxnum = 500;
        int nostep;
        double releps;
        bool stiff;
        double tend;
        //
        //  Test for improper parameters.
        //
        fouru = 4.0 * typeMethods.r8_epsilon();

        switch (neqn)
        {
            case < 1:
                iflag = 6;
                Console.WriteLine("");
                Console.WriteLine("DE - Fatal error!");
                Console.WriteLine("  NEQN < 1.");
                return;
        }

        if (t == tout)
        {
            iflag = 6;
            Console.WriteLine("");
            Console.WriteLine("DE - Fatal error!");
            Console.WriteLine("  T = TOUT.");
            return;
        }

        if (relerr < 0.0 || abserr < 0.0)
        {
            iflag = 6;
            Console.WriteLine("");
            Console.WriteLine("DE - Fatal error!");
            Console.WriteLine("  RELERR < 0 or ABSERR < 0.");
            return;
        }

        eps = Math.Max(relerr, abserr);

        switch (eps)
        {
            case <= 0.0:
                iflag = 6;
                Console.WriteLine("");
                Console.WriteLine("DE - Fatal error!");
                Console.WriteLine("  max ( RELERR, ABSERR ) <= 0.");
                return;
        }

        switch (iflag)
        {
            case 0:
                iflag = 6;
                Console.WriteLine("");
                Console.WriteLine("DE - Fatal error!");
                Console.WriteLine("  IFLAG = 0 on input.");
                return;
        }

        isn = typeMethods.i4_sign(iflag);
        iflag = Math.Abs(iflag);

        if (iflag != 1)
        {
            if (t != told)
            {
                iflag = 6;
                Console.WriteLine("");
                Console.WriteLine("DE - Fatal error!");
                Console.WriteLine("  IFLAG is not 1, and T is not equal to TOLD.");
                return;
            }

            switch (iflag)
            {
                case < 2:
                case > 5:
                    iflag = 6;
                    return;
            }
        }

        //
        //  On each call set interval of integration and counter for number of
        //  steps.  Adjust input error tolerances to define weight vector for
        //  subroutine STEP.
        //
        del = tout - t;
        absdel = Math.Abs(del);

        tend = isn switch
        {
            < 0 => tout,
            _ => t + 10.0 * del
        };

        nostep = 0;
        kle4 = 0;
        stiff = false;
        releps = relerr / eps;
        abseps = abserr / eps;
        //
        //  On start and restart, also set work variables X and YY(*), store the
        //  direction of integration, and initialize the step size.
        //
        if (iflag == 1 || isnold < 0 || delsgn * del <= 0.0)
        {
            start = true;
            x = t;
            for (l = 1; l <= neqn; l++)
            {
                yy[l - 1] = y[l - 1];
            }

            delsgn = typeMethods.r8_sign(del);

            h = Math.Max(Math.Abs(tout - x), fouru * Math.Abs(x)) * typeMethods.r8_sign(tout - x);
        }

        //
        //  If already past the output point, then interpolate and return.
        //
        for (;;)
        {
            if (absdel <= Math.Abs(x - t))
            {
                intrp(x, yy, yyIndex, tout, ref y, yIndex, ref ypout, ypoutIndex, neqn, kold, phi, phiIndex, psi,
                    psiIndex);
                iflag = 2;
                t = tout;
                told = t;
                isnold = isn;
                break;
            }

            //
            //  If we cannot go past the output point, and we are sufficiently
            //  close to it, then extrapolate and return.
            //
            if (isn <= 0 && Math.Abs(tout - x) < fouru * Math.Abs(x))
            {
                h = tout - x;
                yp = f(x, yy, yyIndex, yp, ypIndex);
                for (l = 1; l <= neqn; l++)
                {
                    y[l - 1] = yy[(l - 1 + yyIndex) % yy.Length] + h * yp[(l - 1 + ypIndex) % yp.Length];
                }

                iflag = 2;
                t = tout;
                told = t;
                isnold = isn;
                break;
            }

            //
            //  Test for too many steps.
            //
            if (maxnum <= nostep)
            {
                iflag = stiff switch
                {
                    true => isn * 5,
                    _ => isn * 4
                };

                for (l = 1; l <= neqn; l++)
                {
                    y[l - 1] = yy[(l - 1 + yyIndex) % yy.Length];
                }

                t = x;
                told = t;
                isnold = 1;
                break;
            }

            //
            //  Limit the step size, set the weight vector and take a step.
            //
            h = Math.Min(Math.Abs(h), Math.Abs(tend - x)) * typeMethods.r8_sign(h);

            for (l = 1; l <= neqn; l++)
            {
                wt[(l - 1 + wtIndex) % wt.Length] =
                    releps * Math.Abs(yy[(l - 1 + yyIndex) % yy.Length]) + abseps;
            }

            step(ref x, ref yy, yyIndex, f, neqn, ref h, ref eps, wt, wtIndex, ref start,
                ref hold, ref k, ref kold, ref crash, ref phi, phiIndex, ref p, pIndex, ref yp, ypIndex, ref psi,
                psiIndex,
                ref alpha, alphaIndex, ref beta, betaIndex, ref sig, sigIndex, ref v, vIndex, ref w, wIndex, ref g,
                gIndex, ref phase1, ref ns, ref nornd);
            //
            //  Test for tolerances too small.
            //
            if (crash)
            {
                iflag = isn * 3;
                relerr = eps * releps;
                abserr = eps * abseps;
                for (l = 1; l <= neqn; l++)
                {
                    y[l - 1] = yy[(l - 1 + yyIndex) % yy.Length];
                }

                t = x;
                told = t;
                isnold = 1;
                break;
            }

            //
            //  Augment the step counter and test for stiffness.
            //
            nostep += 1;
            kle4 += 1;

            kle4 = kold switch
            {
                > 4 => 0,
                _ => kle4
            };

            stiff = kle4 switch
            {
                >= 50 => true,
                _ => stiff
            };
        }
    }

    public static void intrp
        (
            double x,
            double[] y,
            int yIndex,
            double xout,
            ref double[] yout,
            int youtIndex,
            ref double[] ypout,
            int ypoutIndex,
            int neqn,
            int kold,
            double[] phi,
            int phiIndex,
            double[] psi,
            int psiIndex
        )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INTRP approximates the solution at XOUT by polynomial interpolation.
        //
        //  Discussion:
        //
        //    The methods in STEP approximate the solution near X by a polynomial.
        //    This routine approximates the solution at XOUT by evaluating the
        //    polynomial there.  Information defining this polynomial is passed
        //    from STEP, so INTRP cannot be used alone.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 January 2012
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Lawrence Shampine, Marilyn Gordon.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Lawrence Shampine, Marilyn Gordon,
        //    Computer Solution of Ordinary Differential Equations:
        //    The Initial Value Problem,
        //    Freeman, 1975,
        //    ISBN: 0716704617,
        //    LC: QA372.S416.
        //
        //  Parameters:
        //
        //    Input, double X, the point where the solution has been computed.
        //
        //    Input, double Y[NEQN], the computed solution at X.
        //
        //    Input, double XOUT, the point at which the solution is desired.
        //
        //    Output, double YOUT[NEQN], the solution at XOUT.
        //
        //    Output, double YPOUT[NEQN], the derivative of the solution
        //    at XOUT.
        //
        //    Input, int NEQN, the number of equations.
        //
        //    Input, int KOLD, the order used for the last
        //    successful step.
        //
        //    Input, double PHI[NEQN*16], contains information about the
        //    interpolating polynomial.
        //
        //    Input, double PSI[12], contains information about the
        //    interpolating polynomial.
        //
    {
        double eta;
        double[] g = new double[13];
        double gamma;
        double hi;
        int i;
        int j;
        int k;
        int ki;
        double psijm1;
        double[] rho = new double[13];
        double term;
        double[] w = new double[13];

        hi = xout - x;
        ki = kold + 1;
        //
        //  Initialize W for computing G.
        //
        for (i = 1; i <= ki; i++)
        {
            w[i - 1] = 1.0 / i;
        }

        //
        //  Compute G.
        //
        g[0] = 1.0;
        rho[0] = 1.0;
        term = 0.0;

        for (j = 2; j <= ki; j++)
        {
            psijm1 = psi[(j - 2 + psiIndex) % psi.Length];
            gamma = (hi + term) / psijm1;
            eta = hi / psijm1;
            for (i = 1; i <= ki + 1 - j; i++)
            {
                w[i - 1] = gamma * w[i - 1] - eta * w[i];
            }

            g[j - 1] = w[0];
            rho[j - 1] = gamma * rho[j - 2];
            term = psijm1;
        }

        //
        //  Interpolate.
        //
        for (k = 0; k < neqn; k++)
        {
            ypout[(k + ypoutIndex) % ypout.Length] = 0.0;
            yout[(k + youtIndex) % yout.Length] = 0.0;
        }

        for (j = 1; j <= ki; j++)
        {
            i = ki + 1 - j;
            for (k = 0; k < neqn; k++)
            {
                yout[(k + youtIndex) % yout.Length] += g[i - 1] * phi[
                    (k + (i - 1) * neqn + phiIndex) % phi.Length];
                ypout[(k + ypoutIndex) % ypout.Length] += rho[i - 1] *
                                                          phi[(k + (i - 1) * neqn + phiIndex) % phi.Length];
            }
        }

        for (k = 0; k < neqn; k++)
        {
            yout[(k + youtIndex) % yout.Length] =
                y[(k + yIndex) % y.Length] + hi * yout[(k + youtIndex) % yout.Length];
        }
    }


    public static void ode
        (
            Func<double, double[], int, double[], int, double[]> f,
            int neqn,
            double[] y,
            ref double t,
            double tout,
            double relerr,
            double abserr,
            ref int iflag,
            ref double[] work,
            ref int[] iwork
        )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ODE is the user interface to an ordinary differential equation solver.
        //
        //  Discussion:
        //
        //    ODE integrates a system of NEQN first order ordinary differential
        //    equations of the form:
        //      dY(i)/dT = F(T,Y(1),Y(2),...,Y(NEQN))
        //      Y(i) given at T.
        //    The subroutine integrates from T to TOUT.  On return, the
        //    parameters in the call list are set for continuing the integration.
        //    The user has only to define a new value TOUT and call ODE again.
        //
        //    The differential equations are actually solved by a suite of codes
        //    DE, STEP, and INTRP.  ODE allocates virtual storage in the
        //    arrays WORK and IWORK and calls DE.  DE is a supervisor which
        //    directs the solution.  It calls the routines STEP and INTRP
        //    to advance the integration and to interpolate at output points.
        //
        //    STEP uses a modified divided difference form of the Adams PECE
        //    formulas and local extrapolation.  It adjusts the order and step
        //    size to control the local error per unit step in a generalized
        //    sense.  Normally each call to STEP advances the solution one step
        //    in the direction of TOUT.  For reasons of efficiency, DE integrates
        //    beyond TOUT internally, though never beyond T+10*(TOUT-T), and
        //    calls INTRP to interpolate the solution at TOUT.  An option is
        //    provided to stop the integration at TOUT but it should be used
        //    only if it is impossible to continue the integration beyond TOUT.
        //
        //    On the first call to ODE, the user must provide storage in the calling
        //    program for the arrays in the call list,
        //      Y(NEQN), WORK(100+21*NEQN), IWORK(5),
        //    declare F in an external statement, supply the double precision
        //      SUBROUTINE F ( T, Y, YP )
        //    to evaluate dy(i)/dt = yp(i) = f(t,y(1),y(2),...,y(neqn))
        //    and initialize the parameters:
        //    * NEQN, the number of equations to be integrated;
        //    * Y(1:NEQN), the vector of initial conditions;
        //    * T, the starting point of integration;
        //    * TOUT, the point at which a solution is desired;
        //    * RELERR, ABSERR, the relative and absolute local error tolerances;
        //    * IFLAG, an indicator to initialize the code.  Normal input
        //      is +1.  The user should set IFLAG = -1 only if it is
        //      impossible to continue the integration beyond TOUT.
        //    All parameters except F, NEQN and TOUT may be altered by the
        //    code on output, and so must be variables in the calling program.
        //
        //    On normal return from ODE, IFLAG is 2, indicating that T has been
        //    set to TOUT, and Y has been set to the approximate solution at TOUT.
        //
        //    If IFLAG is 3, then the program noticed that RELERR or ABSERR was
        //    too small; the output values of these variables are more appropriate,
        //    and integration can be resumed by setting IFLAG to 1.
        //
        //    IFLAG is -2 if the user input IFLAG = -1, and the code was able to
        //    reach TOUT exactly.  In that case, the output value of T is TOUT,
        //    and the output value of Y is the solution at TOUT, which was computed
        //    directly, and not by interpolation.
        //
        //    Other values of IFLAG generally indicate an error.
        //
        //    Normally, it is desirable to compute many points along the solution
        //    curve.  After the first successful step, more steps may be taken
        //    simply by updating the value of TOUT and calling ODE again.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 January 2012
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Lawrence Shampine, Marilyn Gordon.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Lawrence Shampine, Marilyn Gordon,
        //    Computer Solution of Ordinary Differential Equations:
        //    The Initial Value Problem,
        //    Freeman, 1975,
        //    ISBN: 0716704617,
        //    LC: QA372.S416.
        //
        //  Parameters:
        //
        //    Input, void F ( double T, double Y[], double YP[] ), the user-supplied function
        //    which accepts input values T and Y[], evaluates the right hand
        //    sides of the ODE, and stores the result in YP[].
        //
        //    Input, int NEQN, the number of equations.
        //
        //    Input/output, double Y[NEQN], the current solution.
        //
        //    Input/output, double &T, the current value of the independent
        //    variable.
        //
        //    Input, double TOUT, the desired value of T on output.
        //
        //    Input, double RELERR, ABSERR, the relative and absolute error
        //    tolerances.  At each step, the code requires
        //      abs ( local error ) <= abs ( y ) * relerr + abserr
        //    for each component of the local error and solution vectors.
        //
        //    Input/output, int &IFLAG, indicates the status of integration.
        //    On input, IFLAG is normally 1 (or -1 in the special case where TOUT is
        //    not to be exceeded.)  On normal output, IFLAG is 2.  Other output values
        //    are:
        //    * 3, integration did not reach TOUT because the error tolerances
        //      were too small.  But RELERR and ABSERR were increased appropriately
        //      for continuing;
        //    * 4, integration did not reach TOUT because more than 500 steps were taken;
        //    * 5, integration did not reach TOUT because the equations appear to
        //      be stiff;
        //    * 6, invalid input parameters (fatal error).
        //    The value of IFLAG is returned negative when the input value is
        //    negative and the integration does not reach TOUT.
        //
        //    Input/output, double WORK[100+21*NEQN], workspace.
        //
        //    Input/output, int IWORK[5], workspace.
        //
    {
        const int ialpha = 1;
        const int ibeta = 13;
        const int idelsn = 93;
        const int ig = 62;
        const int ih = 89;
        const int ihold = 90;
        int ip;
        const int iphase = 75;
        int iphi;
        const int ipsi = 76;
        const int isig = 25;
        const int istart = 91;
        const int itold = 92;
        const int iv = 38;
        const int iw = 50;
        int iwt;
        const int ix = 88;
        int iyp;
        int iypout;
        const int iyy = 100;
        bool nornd = false;
        bool phase1 = false;
        bool start = false;

        iwt = iyy + neqn;
        ip = iwt + neqn;
        iyp = ip + neqn;
        iypout = iyp + neqn;
        iphi = iypout + neqn;

        if (Math.Abs(iflag) != 1)
        {
            start = 0.0 < work[istart - 1];
            phase1 = 0.0 < work[iphase - 1];
            nornd = iwork[1] != -1;
        }

        int yIndex = 0;
        int yyIndex = iyy - 1;
        int wtIndex = +iwt - 1;
        int pIndex = +ip - 1;
        int ypIndex = +iyp - 1;
        int ypoutIndex = +iypout - 1;
        int phiIndex = +iphi - 1;
        int alphaIndex = +ialpha - 1;
        int betaIndex = +ibeta - 1;
        int sigIndex = +isig - 1;
        int vIndex = +iv - 1;
        int wIndex = +iw - 1;
        int gIndex = +ig - 1;
        int psiIndex = +ipsi - 1;

        de(f, neqn, ref y, ref yIndex, ref t, tout, relerr, abserr, ref iflag, ref work, ref yyIndex,
            work, ref wtIndex, ref work, ref pIndex, ref work, ref ypIndex, ref work, ref ypoutIndex, ref work,
            ref phiIndex,
            ref work, ref alphaIndex, ref work, ref betaIndex, ref work, ref sigIndex, ref work, ref vIndex,
            ref work, ref wIndex, ref work, ref gIndex,
            ref phase1, ref work, ref psiIndex, ref work[ix - 1], ref work[ih - 1], ref work[ihold - 1], ref start,
            ref work[itold - 1], ref work[idelsn - 1], ref iwork[0], ref nornd, ref iwork[2], ref iwork[3],
            ref iwork[4]);

        work[istart - 1] = start switch
        {
            true => 1.0,
            _ => -1.0
        };

        work[iphase - 1] = phase1 switch
        {
            true => 1.0,
            _ => -1.0
        };

        iwork[1] = nornd switch
        {
            true => 1,
            _ => -1
        };
    }

    public static void step
        (
            ref double x,
            ref double[] y,
            int yIndex,
            Func<double, double[], int, double[], int, double[]> f,
            int neqn,
            ref double h,
            ref double eps,
            double[] wt,
            int wtIndex,
            ref bool start,
            ref double hold,
            ref int k,
            ref int kold,
            ref bool crash,
            ref double[] phi,
            int phiIndex,
            ref double[] p,
            int pIndex,
            ref double[] yp,
            int ypIndex,
            ref double[] psi,
            int psiIndex,
            ref double[] alpha,
            int alphaIndex,
            ref double[] beta,
            int betaIndex,
            ref double[] sig,
            int sigIndex,
            ref double[] v,
            int vIndex,
            ref double[] w,
            int wIndex,
            ref double[] g,
            int gIndex,
            ref bool phase1,
            ref int ns,
            ref bool nornd
        )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    STEP integrates the system of ODEs one step, from X to X+H.
        //
        //  Discussion:
        //
        //    This routine integrates a system of first order ordinary differential
        //    equations one step, normally from x to x+h, using a modified divided
        //    difference form of the Adams PECE formulas.  Local extrapolation is
        //    used to improve absolute stability and accuracy.  The code adjusts its
        //    order and step size to control the local error per unit step in a
        //    generalized sense.  Special devices are included to control roundoff
        //    error and to detect when the user is requesting too much accuracy.
        //
        //    STEP is normally not called directly by the user.  However, it is
        //    possible to do so.
        //
        //    On the first call to STEP, the user must pass in values for:
        //    * X, the initial value of the independent variable;
        //    * Y, the vector of initial values of dependent variables;
        //    * NEQN, the number of equations to be integrated;
        //    * H, the nominal step size indicating direction of integration
        //      and maximum size of step.  H must be a variable, not a constant;
        //    * EPS, the local error tolerance per step.  EPS must be variable;
        //    * WT, the vector of non-zero weights for error criterion;
        //    * START, set to TRUE.
        //
        //    STEP requires the L2 norm of the vector with components
        //      local error(1:NEQN) / WT(1:NEQN)
        //    to be less than EPS for a successful step.  The array WT allows the user
        //    to specify an error test appropriate for the problem.  For example,
        //    if WT(L):
        //    = 1.0, specifies absolute error,
        //    = abs(Y(L)), specifies error relative to the most recent value of
        //      the L-th component of the solution,
        //    = abs(YP(L)), specifies error relative to the most recent value of
        //      the L-th component of the derivative,
        //    = max (WT(L),abs(Y(L))), specifies error relative to the largest
        //      magnitude of L-th component obtained so far,
        //    = abs(Y(L))*RELERR/EPS + ABSERR/EPS, specifies a mixed
        //      relative-absolute test where EPS = max ( RELERR, ABSERR ).
        //
        //    On subsequent calls to STEP, the routine is designed so that all
        //    information needed to continue the integration, including the next step
        //    size H and the next order K, is returned with each step.  With the
        //    exception of the step size, the error tolerance, and the weights, none
        //    of the parameters should be altered.  The array WT must be updated after
        //    each step to maintain relative error tests like those above.
        //
        //    Normally the integration is continued just beyond the desired endpoint
        //    and the solution interpolated there with subroutine INTRP.  If it is
        //    impossible to integrate beyond the endpoint, the step size may be
        //    reduced to hit the endpoint since the code will not take a step
        //    larger than the H input.
        //
        //    Changing the direction of integration, that is, the sign of H, requires
        //    the user to set START = TRUE before calling STEP again.  This is the
        //    only situation in which START should be altered.
        //
        //    A successful step: the subroutine returns after each successful step with
        //    START and CRASH both set to FALSE.  X represents the independent variable
        //    advanced by one step of length HOLD from its input value; Y has been
        //    updated to the solution vector at the new value of X.  All other parameters
        //    represent information corresponding to the new X needed to continue the
        //    integration.
        //
        //    Unsuccessful steps: when the error tolerance is too small, the subroutine
        //    returns without taking a step and sets CRASH to TRUE. An appropriate step
        //    size and error tolerance for continuing are estimated and all other
        //    information is restored as upon input before returning.  To continue
        //    with the larger tolerance, the user just calls the code again.  A
        //    restart is neither required nor desirable.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 January 2012
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Lawrence Shampine, Marilyn Gordon.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Lawrence Shampine, Marilyn Gordon,
        //    Computer Solution of Ordinary Differential Equations:
        //    The Initial Value Problem,
        //    Freeman, 1975,
        //    ISBN: 0716704617,
        //    LC: QA372.S416.
        //
        //  Parameters:
        //
        //    Input/output, double &X, the value of the independent variable.
        //
        //    Input/output, double Y[NEQN], the approximate solution at the
        //    current value of the independent variable.
        //
        //    Input, void F ( double T, double Y[], double YP[] ), the user-supplied function
        //    which accepts input values T and Y[], evaluates the right hand
        //    sides of the ODE, and stores the result in YP[].
        //
        //    Input, int NEQN, the number of equations.
        //
        //    Input/output, double &H, the suggested stepsize.
        //
        //    Input/output, double &EPS, the local error tolerance.
        //
        //    Input, double WT[NEQN], the vector of error weights.
        //
        //    Input/output, bool &START, is set to TRUE before the first step.
        //    The program initializes data, and resets START to FALSE.
        //
        //    Input/output, double &HOLD, the step size used on the last
        //    successful step.
        //
        //    Input/output, int &K, the appropriate order for the next step.
        //
        //    Input/output, int &KOLD, the order used on the last
        //    successful step.
        //
        //    Output, bool &CRASH, is set to TRUE if no step can be taken.
        //
        //    Workspace, double PHI[NEQN*16], contains divided difference
        //    information about the polynomial interpolant to the solution.
        //
        //    Workspace, double P[NEQN].
        //
        //    Workspace, double YP[NEQN], used to hold values of the
        //    solution derivative.
        //
        //    Workspace, double PSI[12], contains information about
        //    the polynomial interpolant to the solution.
        //
        //    Workspace, double ALPHA[12], BETA[12], SIG[13].
        //
        //    Workspace, double V[12], W[12], G[13].
        //
        //    Input/output, bool &PHASE1, indicates whether the program is in the
        //    first phase, when it always wants to increase the ODE method order.
        //
        //    Input/output, int &NS, the number of steps taken with
        //    stepsize H.
        //
        //    Input/output, bool &NORND, ?
        //
    {
        double absh;
        double erk;
        double erkm1;
        double erkm2;
        double erkp1;
        double err;
        double fouru;
        double[] gstr =
        {
            0.50, 0.0833, 0.0417, 0.0264, 0.0188,
            0.0143, 0.0114, 0.00936, 0.00789, 0.00679,
            0.00592, 0.00524, 0.00468
        };
        double hnew;
        int i;
        int ifail;
        int iq;
        int j;
        int km1;
        int km2;
        int knew;
        int kp1;
        int kp2;
        int l;
        int nsp1;
        double p5eps;
        double r;
        double rho;
        double round;
        double total;
        double tau;
        double temp1;
        double temp2;
        double[] two =
        {
            2.0, 4.0, 8.0, 16.0, 32.0,
            64.0, 128.0, 256.0, 512.0, 1024.0,
            2048.0, 4096.0, 8192.0
        };
        double twou;
        double xold;

        twou = 2.0 * typeMethods.r8_epsilon();
        fouru = 2.0 * twou;
        //
        //  Check if the step size or error tolerance is too small.  If this is the
        //  first step, initialize the PHI array and estimate a starting step size.
        //
        //  If the step size is too small, determine an acceptable one.
        //
        crash = true;

        if (Math.Abs(h) < fouru * Math.Abs(x))
        {
            h = fouru * Math.Abs(x) * typeMethods.r8_sign(h);
            return;
        }

        p5eps = 0.5 * eps;
        //
        //  If the error tolerance is too small, increase it to an acceptable value.
        //
        round = 0.0;
        for (i = 0; i < neqn; i++)
        {
            round += Math.Pow(y[(i + yIndex) % y.Length] / wt[(i + wtIndex) % wt.Length], 2);
        }

        round = twou * Math.Sqrt(round);

        if (p5eps < round)
        {
            eps = 2.0 * round * (1.0 + fouru);
            return;
        }

        crash = false;
        g[(0 + gIndex) % g.Length] = 1.0;
        g[(1 + gIndex) % g.Length] = 0.5;
        sig[(0 + sigIndex) % sig.Length] = 1.0;
        switch (start)
        {
            //
            //  Initialize.  Compute an appropriate step size for the first step.
            //
            case true:
            {
                yp = f(x, y, yIndex, yp, ypIndex);
                for (l = 1; l <= neqn; l++)
                {
                    phi[(l - 1 + 0 * neqn + phiIndex) % phi.Length] = yp[(l - 1 + ypIndex) % yp.Length];
                    phi[(l - 1 + 1 * neqn + phiIndex) % phi.Length] = 0.0;
                }

                total = 0.0;
                for (l = 1; l <= neqn; l++)
                {
                    total += Math.Pow(yp[(l - 1 + ypIndex) % yp.Length] / wt[(l - 1 + wtIndex) % wt.Length],
                        2);
                }

                total = Math.Sqrt(total);
                absh = Math.Abs(h);
                if (eps < 16.0 * total * h * h)
                {
                    absh = 0.25 * Math.Sqrt(eps / total);
                }

                h = Math.Max(absh, fouru * Math.Abs(x)) * typeMethods.r8_sign(h);
                hold = 0.0;
                k = 1;
                kold = 0;
                start = false;
                phase1 = true;
                nornd = true;

                if (p5eps <= 100.0 * round)
                {
                    nornd = false;
                    for (l = 1; l <= neqn; l++)
                    {
                        phi[(l - 1 + 14 * neqn + phiIndex) % phi.Length] = 0.0;
                    }
                }

                break;
            }
        }

        ifail = 0;
        //
        //  Compute coefficients of formulas for this step.  Avoid computing
        //  those quantities not changed when step size is not changed.
        //
        for (;;)
        {
            kp1 = k + 1;
            kp2 = k + 2;
            km1 = k - 1;
            km2 = k - 2;
            //
            //  NS is the number of steps taken with size H, including the current
            //  one.  When K < NS, no coefficients change.
            //
            if (h != hold)
            {
                ns = 0;
            }

            if (ns <= kold)
            {
                ns += 1;
            }

            nsp1 = ns + 1;
            //
            //  Compute those components of ALPHA, BETA, PSI and SIG which change.
            //
            if (ns <= k)
            {
                beta[(ns - 1 + betaIndex) % beta.Length] = 1.0;
                alpha[(ns - 1 + alphaIndex) % alpha.Length] = 1.0 / ns;
                temp1 = h * ns;
                sig[(nsp1 - 1 + sigIndex) % sig.Length] = 1.0;

                for (i = nsp1; i <= k; i++)
                {
                    temp2 = psi[(i - 2 + psiIndex) % psi.Length];
                    psi[(i - 2 + psiIndex) % psi.Length] = temp1;
                    beta[(i - 1 + betaIndex) % beta.Length] = beta[(i - 2 + betaIndex) % beta.Length] *
                        psi[(i - 2 + psiIndex) % psi.Length] / temp2;
                    temp1 = temp2 + h;
                    alpha[(i - 1 + alphaIndex) % alpha.Length] = h / temp1;
                    sig[(i + sigIndex) % sig.Length] = i *
                                                       alpha[(i - 1 + alphaIndex) % alpha.Length] *
                                                       sig[(i - 1 + sigIndex) % sig.Length];
                }

                psi[(k - 1 + psiIndex) % psi.Length] = temp1;
                switch (ns)
                {
                    //
                    //  Compute coefficients G.
                    //
                    //  Initialize V and set W.
                    //
                    case <= 1:
                    {
                        for (iq = 1; iq <= k; iq++)
                        {
                            v[(iq - 1 + vIndex) % v.Length] = 1.0 / (iq * (iq + 1));
                            w[(iq - 1 + wIndex) % w.Length] = v[(iq - 1 + vIndex) % v.Length];
                        }

                        break;
                    }
                    //
                    default:
                    {
                        if (kold < k)
                        {
                            v[(k - 1 + vIndex) % v.Length] = 1.0 / (k * kp1);

                            for (j = 1; j <= ns - 2; j++)
                            {
                                i = k - j;
                                v[(i - 1 + vIndex) % v.Length] -= alpha[(j + alphaIndex) % alpha.Length] *
                                                                  v[(i + vIndex) % v.Length];
                            }
                        }

                        //
                        //  Update V and set W.
                        //
                        for (iq = 1; iq <= kp1 - ns; iq++)
                        {
                            v[(iq - 1 + vIndex) % v.Length] -= alpha[(ns - 1 + alphaIndex) % alpha.Length] *
                                                               v[(iq + vIndex) % v.Length];
                            w[(iq - 1 + wIndex) % w.Length] = v[(iq - 1 + vIndex) % v.Length];
                        }

                        g[(nsp1 - 1 + gIndex) % g.Length] = w[(0 + wIndex) % w.Length];
                        break;
                    }
                }

                //
                //  Compute the G in the work vector W.
                //
                for (i = ns + 2; i <= kp1; i++)
                {
                    for (iq = 1; iq <= kp2 - i; iq++)
                    {
                        w[(iq - 1 + wIndex) % w.Length] -= alpha[(i - 2 + alphaIndex) % alpha.Length] *
                                                           w[(iq + wIndex) % w.Length];
                    }

                    g[(i - 1 + gIndex) % g.Length] = w[(0 + wIndex) % w.Length];
                }
            }

            //
            //  Predict a solution P, evaluate derivatives using predicted
            //  solution, estimate local error at order K and errors at orders K,
            //  K-1, K-2 as if a constant step size were used.
            //
            //  Change PHI to PHI star.
            //
            for (i = nsp1; i <= k; i++)
            {
                for (l = 1; l <= neqn; l++)
                {
                    phi[(l - 1 + (i - 1) * neqn + phiIndex) % phi.Length] =
                        beta[(i - 1 + betaIndex) % beta.Length] *
                        phi[(l - 1 + (i - 1) * neqn + phiIndex) % phi.Length];
                }
            }

            //
            //  Predict solution and differences.
            //
            for (l = 1; l <= neqn; l++)
            {
                phi[(l - 1 + (kp2 - 1) * neqn + phiIndex) % phi.Length] =
                    phi[(l - 1 + (kp1 - 1) * neqn + phiIndex) % phi.Length];
                phi[(l - 1 + (kp1 - 1) * neqn + phiIndex) % phi.Length] = 0.0;
                p[(l - 1 + pIndex) % p.Length] = 0.0;
            }

            for (j = 1; j <= k; j++)
            {
                i = kp1 - j;
                for (l = 1; l <= neqn; l++)
                {
                    p[(l - 1 + pIndex) % p.Length] += g[(i - 1 + gIndex) % g.Length] *
                                                      phi[(l - 1 + (i - 1) * neqn + phiIndex) % phi.Length];
                    phi[(l - 1 + (i - 1) * neqn + phiIndex) % phi.Length] += phi[(l - 1 + i * neqn + phiIndex) % phi.Length];
                }
            }

            switch (nornd)
            {
                case false:
                {
                    for (l = 1; l <= neqn; l++)
                    {
                        tau = h * p[(l - 1 + pIndex) % p.Length] -
                              phi[(l - 1 + (15 - 1) * neqn + phiIndex) % phi.Length];
                        p[(l - 1 + pIndex) % p.Length] = y[(l - 1 + yIndex) % y.Length] + tau;
                        phi[(l - 1 + (16 - 1) * neqn + phiIndex) % phi.Length] =
                            p[(l - 1 + pIndex) % p.Length] - y[(l - 1 + yIndex) % y.Length] - tau;
                    }

                    break;
                }
                default:
                {
                    for (l = 1; l <= neqn; l++)
                    {
                        p[(l - 1 + pIndex) % p.Length] =
                            y[(l - 1 + yIndex) % y.Length] + h * p[(l - 1 + pIndex) % p.Length];
                    }

                    break;
                }
            }

            xold = x;
            x += h;
            absh = Math.Abs(h);
            yp = f(x, p, pIndex, yp, ypIndex);
            //
            //  Estimate the errors at orders K, K-1 and K-2.
            //
            erkm2 = 0.0;
            erkm1 = 0.0;
            erk = 0.0;

            for (l = 1; l <= neqn; l++)
            {
                switch (km2)
                {
                    case > 0:
                        erkm2 += Math.Pow(
                            (phi[(l - 1 + (km1 - 1) * neqn + phiIndex) % phi.Length] +
                             yp[(l - 1 + ypIndex) % yp.Length] -
                             phi[(l - 1 + 0 * neqn + phiIndex) % phi.Length]) /
                            wt[(l - 1 + wtIndex) % wt.Length], 2);
                        break;
                }

                switch (km2)
                {
                    case >= 0:
                        erkm1 += Math.Pow(
                            (phi[(l - 1 + (k - 1) * neqn + phiIndex) % phi.Length] +
                             yp[(l - 1 + ypIndex) % yp.Length] -
                             phi[(l - 1 + 0 * neqn + phiIndex) % phi.Length]) /
                            wt[(l - 1 + wtIndex) % wt.Length], 2);
                        break;
                }

                erk += Math.Pow(
                    (yp[(l - 1 + ypIndex) % yp.Length] - phi[(l - 1 + 0 * neqn + phiIndex) % phi.Length]) /
                    wt[(l - 1 + wtIndex) % wt.Length], 2);
            }

            erkm2 = km2 switch
            {
                > 0 => absh * sig[(km1 - 1 + sigIndex) % sig.Length] * gstr[km2 - 1] * Math.Sqrt(erkm2),
                _ => erkm2
            };

            erkm1 = km2 switch
            {
                >= 0 => absh * sig[(k - 1 + sigIndex) % sig.Length] * gstr[km1 - 1] * Math.Sqrt(erkm1),
                _ => erkm1
            };

            err = absh * Math.Sqrt(erk) * (g[(k - 1 + gIndex) % g.Length] - g[(kp1 - 1 + gIndex) % g.Length]);
            erk = absh * Math.Sqrt(erk) * sig[(kp1 - 1 + sigIndex) % sig.Length] * gstr[k - 1];
            knew = k;
            switch (km2)
            {
                //
                //  Test if the order should be lowered.
                //
                case > 0:
                {
                    if (Math.Max(erkm1, erkm2) <= erk)
                    {
                        knew = km1;
                    }

                    break;
                }
                case 0:
                {
                    if (erkm1 <= 0.5 * erk)
                    {
                        knew = km1;
                    }

                    break;
                }
            }

            //
            //  Test if the step was successful.
            //
            if (err <= eps)
            {
                break;
            }

            //
            //  The step is unsuccessful.  Restore X, PHI and PSI.
            //  If third consecutive failure, set order to one.  If the step fails more
            //  than three times, consider an optimal step size.  Double the error
            //  tolerance and return if the estimated step size is too small for machine
            //  precision.
            //
            //  Restore X, PHI and PSI.
            //
            phase1 = false;
            x = xold;
            for (i = 1; i <= k; i++)
            {
                for (l = 1; l <= neqn; l++)
                {
                    phi[(l - 1 + (i - 1) * neqn + phiIndex) % phi.Length] =
                        (phi[(l - 1 + (i - 1) * neqn + phiIndex) % phi.Length] -
                         phi[(l - 1 + i * neqn + phiIndex) % phi.Length]) /
                        beta[(i - 1 + betaIndex) % beta.Length];
                }
            }

            for (i = 2; i <= k; i++)
            {
                psi[(i - 2 + psiIndex) % psi.Length] = psi[(i - 1 + psiIndex) % psi.Length] - h;
            }

            //
            //  On third failure, set the order to one.  Thereafter, use optimal step size.
            //
            ifail += 1;
            temp2 = 0.5;

            switch (ifail)
            {
                case > 3:
                {
                    if (p5eps < 0.25 * erk)
                    {
                        temp2 = Math.Sqrt(p5eps / erk);
                    }

                    break;
                }
            }

            knew = ifail switch
            {
                >= 3 => 1,
                _ => knew
            };

            h = temp2 * h;
            k = knew;

            if (Math.Abs(h) < fouru * Math.Abs(x))
            {
                crash = true;
                h = Math.Abs(fouru * Math.Abs(x)) * typeMethods.r8_sign(h);
                eps += eps;
                return;
            }
        }

        //
        //  The step is successful.  Correct the predicted solution, evaluate
        //  the derivatives using the corrected solution and update the
        //  differences.  Determine best order and step size for next step.
        //
        kold = k;
        hold = h;
        switch (nornd)
        {
            //
            //  Correct and evaluate.
            //
            case false:
            {
                for (l = 1; l <= neqn; l++)
                {
                    rho = h * g[(kp1 - 1 + gIndex) % g.Length] * (yp[(l - 1 + ypIndex) % yp.Length] -
                                                                  phi[(l - 1 + 0 * neqn + phiIndex) % phi.Length]) -
                          phi[(l - 1 + (16 - 1) * neqn + phiIndex) % phi.Length];
                    y[(l - 1 + yIndex) % y.Length] = p[(l - 1 + pIndex) % p.Length] + rho;
                    phi[(l - 1 + (15 - 1) * neqn + phiIndex) % phi.Length] =
                        y[(l - 1 + phiIndex) % phi.Length] - p[(l - 1 + pIndex) % p.Length] - rho;
                }

                break;
            }
            default:
            {
                for (l = 1; l <= neqn; l++)
                {
                    y[(l - 1 + phiIndex) % phi.Length] = p[(l - 1 + pIndex) % p.Length] +
                                                         h * g[(kp1 - 1 + gIndex) % g.Length] *
                                                         (yp[(l - 1 + ypIndex) % yp.Length] -
                                                          phi[(l - 1 + 0 * neqn + phiIndex) % phi.Length]);
                }

                break;
            }
        }

        yp = f(x, y, yIndex, yp, ypIndex);
        //
        //  Update differences for the next step.
        //
        for (l = 1; l <= neqn; l++)
        {
            phi[(l - 1 + (kp1 - 1) * neqn + phiIndex) % phi.Length] = yp[(l - 1 + ypIndex) % yp.Length] -
                                                                      phi[
                                                                          (l - 1 + 0 * neqn + phiIndex) %
                                                                          phi.Length];
            phi[(l - 1 + (kp2 - 1) * neqn + phiIndex) % phi.Length] =
                phi[(l - 1 + (kp1 - 1) * neqn + phiIndex) % phi.Length] -
                phi[(l - 1 + (kp2 - 1) * neqn + phiIndex) % phi.Length];
        }

        for (i = 1; i <= k; i++)
        {
            for (l = 1; l <= neqn; l++)
            {
                phi[(l - 1 + (i - 1) * neqn + phiIndex) % phi.Length] += phi[(l - 1 + (kp1 - 1) * neqn + phiIndex) % phi.Length];
            }
        }

        //
        //  Estimate error at order K+1 unless:
        //  * in first phase when always raise order, or,
        //  * already decided to lower order, or,
        //  * step size not constant so estimate unreliable.
        //
        erkp1 = 0.0;

        if (knew == km1 || k == 12)
        {
            phase1 = false;
        }

        switch (phase1)
        {
            case true:
                k = kp1;
                erk = erkp1;
                break;
            default:
            {
                if (knew == km1)
                {
                    k = km1;
                    erk = erkm1;
                }
                else if (kp1 <= ns)
                {
                    for (l = 1; l <= neqn; l++)
                    {
                        erkp1 += Math.Pow(
                            phi[(l - 1 + (kp2 - 1) * neqn + phiIndex) % phi.Length] /
                            wt[(l - 1 + wtIndex) % wt.Length], 2);
                    }

                    erkp1 = absh * gstr[kp1 - 1] * Math.Sqrt(erkp1);
                    switch (k)
                    {
                        //
                        //  Using estimated error at order K+1, determine appropriate order
                        //  for next step.
                        //
                        case 1:
                        {
                            if (erkp1 < 0.5 * erk)
                            {
                                k = kp1;
                                erk = erkp1;
                            }

                            break;
                        }
                        default:
                        {
                            if (erkm1 <= Math.Min(erk, erkp1))
                            {
                                k = km1;
                                erk = erkm1;
                            }
                            else if (erkp1 < erk && k < 12)
                            {
                                k = kp1;
                                erk = erkp1;
                            }

                            break;
                        }
                    }
                }

                break;
            }
        }

        //
        //  With the new order, determine appropriate step size for next step.
        //
        hnew = h + h;

        switch (phase1)
        {
            case false:
            {
                if (p5eps < erk * two[k])
                {
                    hnew = h;

                    if (p5eps < erk)
                    {
                        temp2 = k + 1;
                        r = Math.Pow(p5eps / erk, 1.0 / temp2);
                        hnew = absh * Math.Max(0.5, Math.Min(0.9, r));
                        hnew = Math.Abs(Math.Max(hnew, fouru * Math.Abs(x))) * typeMethods.r8_sign(h);
                    }
                }

                break;
            }
        }

        h = hnew;
    }

}