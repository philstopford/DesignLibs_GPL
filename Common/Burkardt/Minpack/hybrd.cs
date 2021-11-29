using System;
using Burkardt.SolveNS;
using Burkardt.Types;

namespace Burkardt.MinpackNS;

public static partial class Minpack
{
    // Last two ints before fcnData are for the array index baselines.
    public static int hybrd(Func < int, double[], double[], int, int, int, fcnData > fcn,
            int n, double[] x,
            double[] fvec, double xtol, int maxfev, int ml, int mu, double epsfcn,
            double[] diag, int mode, double factor, int nprint, int nfev,
            double[] fjac, int ldfjac, double[] r, int lr, double[] qtf, double[] wa1,
            double[] wa2, double[] wa3, double[] wa4, int fjacIndex, int rIndex, int qtfIndex, int wa1Index,  int wa2Index,  int wa3Index,  int wa4Index  )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    hybrd() finds a zero of a system of N nonlinear equations.
        //
        //  Discussion:
        //
        //    The purpose of HYBRD is to find a zero of a system of
        //    N nonlinear functions in N variables by a modification
        //    of the Powell hybrid method. 
        //
        //    The user must provide FCN, which calculates the functions. 
        //
        //    The jacobian is calculated by a forward-difference approximation.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 April 2010
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
        //       fcn is the name of the user-supplied subroutine which
        //         calculates the functions. fcn must be declared
        //         in an external statement in the user calling
        //         program, and should be written as follows.
        //
        //         subroutine fcn(n,x,fvec,iflag)
        //         integer n,iflag
        //         double precision x(n),fvec(n)
        //         ----------
        //         calculate the functions at x and
        //         return this vector in fvec.
        //         ---------
        //         return
        //         end
        //
        //         the value of iflag should not be changed by fcn unless
        //         the user wants to terminate execution of hybrd.
        //         in this case set iflag to a negative integer.
        //
        //    Input, int N, the number of functions and variables.
        //
        //    Input/output, double X[N].  On input an initial estimate of the solution.
        //    On output, the final estimate of the solution.
        //
        //    Output, double FVEC[N], the functions evaluated at the output value of X.
        //
        //    Input, double XTOL, a nonnegative value.  Termination occurs when the 
        //    relative error between two consecutive iterates is at most XTOL.
        //
        //    Input, int MAXFEV.  Termination occurs when the number of calls to FCN 
        //    is at least MAXFEV by the end of an iteration.
        //
        //    Input, int ML, specifies the number of subdiagonals within the band of 
        //    the jacobian matrix.  If the jacobian is not banded, set
        //    ml to at least n - 1.
        //
        //       mu is a nonnegative integer input variable which specifies
        //         the number of superdiagonals within the band of the
        //         jacobian matrix. if the jacobian is not banded, set
        //         mu to at least n - 1.
        //
        //       epsfcn is an input variable used in determining a suitable
        //         step length for the forward-difference approximation. this
        //         approximation assumes that the relative errors in the
        //         functions are of the order of epsfcn. if epsfcn is less
        //         than the machine precision, it is assumed that the relative
        //         errors in the functions are of the order of the machine
        //         precision.
        //
        //       diag is an array of length n. if mode = 1 (see
        //         below), diag is internally set. if mode = 2, diag
        //         must contain positive entries that serve as
        //         multiplicative scale factors for the variables.
        //
        //       mode is an integer input variable. if mode = 1, the
        //         variables will be scaled internally. if mode = 2,
        //         the scaling is specified by the input diag. other
        //         values of mode are equivalent to mode = 1.
        //
        //       factor is a positive input variable used in determining the
        //         initial step bound. this bound is set to the product of
        //         factor and the euclidean norm of diag*x if nonzero, or else
        //         to factor itself. in most cases factor should lie in the
        //         interval (.1,100.). 100. is a generally recommended value.
        //
        //       nprint is an integer input variable that enables controlled
        //         printing of iterates if it is positive. in this case,
        //         fcn is called with iflag = 0 at the beginning of the first
        //         iteration and every nprint iterations thereafter and
        //         immediately prior to return, with x and fvec available
        //         for printing. if nprint is not positive, no special calls
        //         of fcn with iflag = 0 are made.
        //
        //       info is an integer output variable. if the user has
        //         terminated execution, info is set to the (negative)
        //         value of iflag. see description of fcn. otherwise,
        //         info is set as follows.
        //
        //         info = 0   improper input parameters.
        //
        //         info = 1   relative error between two consecutive iterates
        //                    is at most xtol.
        //
        //         info = 2   number of calls to fcn has reached or exceeded
        //                    maxfev.
        //
        //         info = 3   xtol is too small. no further improvement in
        //                    the approximate solution x is possible.
        //
        //         info = 4   iteration is not making good progress, as
        //                    measured by the improvement from the last
        //                    five jacobian evaluations.
        //
        //         info = 5   iteration is not making good progress, as
        //                    measured by the improvement from the last
        //                    ten iterations.
        //
        //       nfev is an integer output variable set to the number of
        //         calls to fcn.
        //
        //       fjac is an output n by n array which contains the
        //         orthogonal matrix q produced by the qr factorization
        //         of the final approximate jacobian.
        //
        //       ldfjac is a positive integer input variable not less than n
        //         which specifies the leading dimension of the array fjac.
        //
        //       r is an output array of length lr which contains the
        //         upper triangular matrix produced by the qr factorization
        //         of the final approximate jacobian, stored rowwise.
        //
        //       lr is a positive integer input variable not less than
        //         (n*(n+1))/2.
        //
        //       qtf is an output array of length n which contains
        //         the vector (q transpose)*fvec.
        //
        //       wa1, wa2, wa3, and wa4 are work arrays of length n.
        //
    {
        double delta = 0;
        int[] iwa = new int[1];
        int j;
        const double p001 = 0.001;
        const double p0001 = 0.0001;
        const double p1 = 0.1;
        const double p5 = 0.5;
        double xnorm = 0;
        //
        //  Certain loops in this function were kept closer to their original FORTRAN77
        //  format, to avoid confusing issues with the array index L.  These loops are
        //  marked "DO NOT ADJUST", although they certainly could be adjusted (carefully)
        //  once the initial translated code is tested.
        //

        //
        //  EPSMCH is the machine precision.
        //
        double epsmch = typeMethods.r8_epsilon();

        int info = 0;
        int iflag = 0;
        nfev = 0;
        switch (n)
        {
            //
            //  Check the input parameters.
            //
            case <= 0:
                info = 0;
                return info;
        }

        switch (xtol)
        {
            case < 0.0:
                info = 0;
                return info;
        }

        switch (maxfev)
        {
            case <= 0:
                info = 0;
                return info;
        }

        switch (ml)
        {
            case < 0:
                info = 0;
                return info;
        }

        switch (mu)
        {
            case < 0:
                info = 0;
                return info;
        }

        switch (factor)
        {
            case <= 0.0:
                info = 0;
                return info;
        }

        if (ldfjac < n)
        {
            info = 0;
            return info;
        }

        if (lr < n * (n + 1) / 2)
        {
            info = 0;
            return info;
        }

        switch (mode)
        {
            case 2:
            {
                for (j = 0; j < n; j++)
                {
                    switch (diag[j])
                    {
                        case <= 0.0:
                            info = 0;
                            return info;
                    }
                }

                break;
            }
        }

        //
        //  Evaluate the function at the starting point and calculate its norm.
        //
        iflag = 1;
        fcnData res = fcn(n, x, fvec, iflag, 0, 0);
        fvec = res.fvec;
        iflag = res.iflag;
        nfev = 1;
        switch (iflag)
        {
            case < 0:
                info = iflag;
                return info;
        }

        double fnorm = Helpers.enorm(n, fvec);
        //
        //  Determine the number of calls to FCN needed to compute the jacobian matrix.
        //
        int msum = Math.Min(ml + mu + 1, n);
        //
        //  Initialize iteration counter and monitors.
        //
        int iter = 1;
        int ncsuc = 0;
        int ncfail = 0;
        int nslow1 = 0;
        int nslow2 = 0;
        //
        //  Beginning of the outer loop.
        //
        for (;;)
        {
            bool jeval = true;
            //
            //  Calculate the jacobian matrix.
            //
            iflag = 2;
            fdjac1(fcn, n, x, fvec, fjac, ldfjac, ref iflag, ml, mu, epsfcn, wa1, wa2, fjacIndex: fjacIndex, wa1Index: wa1Index, wa2Index: wa2Index);

            nfev += msum;
            switch (iflag)
            {
                case < 0:
                    info = iflag;
                    return info;
            }

            //
            //  Compute the QR factorization of the jacobian.
            //
            int tmpi = 1;
            QRSolve.qrfac(n, n, ref fjac, ldfjac, false, ref iwa, ref tmpi, ref wa1, ref wa2, rdiagIndex:wa1Index, acnormIndex:wa2Index);
            switch (iter)
            {
                //
                //  On the first iteration and if MODE is 1, scale according
                //  to the norms of the columns of the initial jacobian.
                //
                case 1:
                {
                    switch (mode)
                    {
                        case 1:
                        {
                            for (j = 0; j < n; j++)
                            {
                                if (wa2[wa2Index + j] != 0.0)
                                {
                                    diag[j] = wa2[wa2Index + j];
                                }
                                else
                                {
                                    diag[j] = 1.0;
                                }
                            }

                            break;
                        }
                    }

                    //
                    //  On the first iteration, calculate the norm of the scaled X
                    //  and initialize the step bound DELTA.
                    //
                    for (j = 0; j < n; j++)
                    {
                        wa3[wa3Index + j] = diag[j] * x[j];
                    }

                    xnorm = Helpers.enorm(n, wa3);

                    delta = xnorm switch
                    {
                        0.0 => factor,
                        _ => factor * xnorm
                    };

                    break;
                }
            }

            //
            //  Form Q' * FVEC and store in QTF.
            //
            int i;
            for (i = 0; i < n; i++)
            {
                qtf[qtfIndex + i] = fvec[i];
            }

            double temp;
            double sum;
            for (j = 0; j < n; j++)
            {
                if (fjac[fjacIndex + j + j * ldfjac] == 0.0)
                {
                    continue;
                }

                sum = 0.0;
                for (i = j; i < n; i++)
                {
                    sum += fjac[fjacIndex + i + j * ldfjac] * qtf[qtfIndex + i];
                }

                temp = -sum / fjac[fjacIndex + j + j * ldfjac];
                for (i = j; i < n; i++)
                {
                    qtf[qtfIndex + i] += fjac[fjacIndex + i + j * ldfjac] * temp;
                }
            }

            //
            //  Copy the triangular factor of the QR factorization into R.
            //
            //  DO NOT ADJUST THIS LOOP, BECAUSE OF L.
            //
            int l;
            for (j = 1; j <= n; j++)
            {
                l = j;
                for (i = 1; i <= j - 1; i++)
                {
                    r[rIndex + (l - 1)] = fjac[fjacIndex + (i - 1) + (j - 1) * ldfjac];
                    l = l + n - i;
                }

                r[rIndex + (l - 1)] = wa1[wa1Index + (j - 1)];
                switch (wa1[wa1Index + (j - 1)])
                {
                    case 0.0:
                        Console.WriteLine("  Matrix is singular.");
                        break;
                }
            }

            //
            //  Accumulate the orthogonal factor in FJAC.
            //
            QRSolve.qform(n, n, ref fjac, ldfjac, qIndex:fjacIndex);
            switch (mode)
            {
                //
                //  Rescale if necessary.
                //
                case 1:
                {
                    for (j = 0; j < n; j++)
                    {
                        diag[j] = Math.Max(diag[j], wa2[wa2Index + j]);
                    }

                    break;
                }
            }

            //
            //  Beginning of the inner loop.
            //
            for (;;)
            {
                switch (nprint)
                {
                    //
                    //  If requested, call FCN to enable printing of iterates.
                    //
                    case > 0:
                    {
                        switch ((iter - 1) % nprint)
                        {
                            case 0:
                            {
                                iflag = 0;
                                fcn(n, x, fvec, iflag, 0, 0);
                                switch (iflag)
                                {
                                    case < 0:
                                        info = iflag;
                                        return info;
                                }

                                break;
                            }
                        }

                        break;
                    }
                }

                //
                //  Determine the direction P.
                //
                dogleg(n, r, lr, diag, qtf, delta, ref wa1, wa2, wa3, rIndex: rIndex, qtbIndex:qtfIndex, xIndex:wa1Index, wa1Index:wa2Index, wa2Index: wa3Index);
                //
                //  Store the direction P and X + P.  Calculate the norm of P.
                //
                for (j = 0; j < n; j++)
                {
                    wa1[wa1Index + j] = -wa1[j];
                    wa2[wa2Index + j] = x[j] + wa1[j];
                    wa3[wa3Index + j] = diag[j] * wa1[wa1Index + j];
                }

                double pnorm = Helpers.enorm(n, wa3);
                delta = iter switch
                {
                    //
                    //  On the first iteration, adjust the initial step bound.
                    //
                    1 => Math.Min(delta, pnorm),
                    _ => delta
                };

                //
                //  Evaluate the function at X + P and calculate its norm.
                //
                iflag = 1;
                fcn(n, wa2, wa4, iflag, wa2Index, wa4Index);
                nfev += 1;
                switch (iflag)
                {
                    case < 0:
                        info = iflag;
                        return info;
                }

                double fnorm1 = Helpers.enorm(n, wa4);
                //
                //  Compute the scaled actual reduction.
                //
                double actred;
                if (fnorm1 < fnorm)
                {
                    actred = 1.0 - fnorm1 / fnorm * (fnorm1 / fnorm);
                }
                else
                {
                    actred = -1.0;
                }

                //
                //  Compute the scaled predicted reduction.
                //
                //  DO NOT ADJUST THIS LOOP, BECAUSE OF L.
                //
                l = 1;
                for (i = 1; i <= n; i++)
                {
                    sum = 0.0;
                    for (j = i; j <= n; j++)
                    {
                        sum += r[rIndex + (l - 1)] * wa1[wa1Index + (j - 1)];
                        l += 1;
                    }

                    wa3[i - 1] = qtf[i - 1] + sum;
                }

                temp = Helpers.enorm(n, wa3);

                double prered;
                if (temp < fnorm)
                {
                    prered = 1.0 - temp / fnorm * (temp / fnorm);
                }
                else
                {
                    prered = 0.0;
                }

                double ratio = prered switch
                {
                    //
                    //  Compute the ratio of the actual to the predicted reduction.
                    //
                    > 0.0 => actred / prered,
                    _ => 0.0
                };

                switch (ratio)
                {
                    //
                    //  Update the step bound.
                    //
                    case < p1:
                        ncsuc = 0;
                        ncfail += 1;
                        delta = p5 * delta;
                        break;
                    default:
                    {
                        ncfail = 0;
                        ncsuc += 1;
                        if (p5 <= ratio || 1 < ncsuc)
                        {
                            delta = Math.Max(delta, pnorm / p5);
                        }

                        delta = Math.Abs(ratio - 1.0) switch
                        {
                            <= p1 => pnorm / p5,
                            _ => delta
                        };

                        break;
                    }
                }

                switch (ratio)
                {
                    //
                    //  On successful iteration, update X, FVEC, and their norms.
                    //
                    case >= p0001:
                    {
                        for (j = 0; j < n; j++)
                        {
                            x[j] = wa2[j];
                            wa2[j] = diag[j] * x[j];
                            fvec[j] = wa4[j];
                        }

                        xnorm = Helpers.enorm(n, wa2);
                        fnorm = fnorm1;
                        iter += 1;
                        break;
                    }
                }

                //
                //  Determine the progress of the iteration.
                //
                nslow1 += 1;
                nslow1 = actred switch
                {
                    >= p001 => 0,
                    _ => nslow1
                };

                switch (jeval)
                {
                    case true:
                        nslow2 += 1;
                        break;
                }

                nslow2 = actred switch
                {
                    >= p1 => 0,
                    _ => nslow2
                };

                //
                //  Test for convergence.
                //
                if (delta <= xtol * xnorm || fnorm == 0.0)
                {
                    info = 1;
                    return info;
                }

                //
                //  Tests for termination and stringent tolerances.
                //
                if (maxfev <= nfev)
                {
                    info = 2;
                    return info;
                }

                if (p1 * Math.Max(p1 * delta, pnorm) <= epsmch * xnorm)
                {
                    info = 3;
                    return info;
                }

                switch (nslow2)
                {
                    case 5:
                        info = 4;
                        return info;
                }

                switch (nslow1)
                {
                    case 10:
                        info = 5;
                        return info;
                }

                //
                //  Criterion for recalculating jacobian approximation by forward differences.
                //
                if (ncfail == 2)
                {
                    break;
                }

                //
                //  Calculate the rank one modification to the jacobian
                //  and update QTF if necessary.
                //
                for (j = 0; j < n; j++)
                {
                    sum = 0.0;
                    for (i = 0; i < n; i++)
                    {
                        sum += fjac[fjacIndex + i + j * ldfjac] * wa4[i];
                    }

                    wa2[j] = (sum - wa3[j]) / pnorm;
                    wa1[wa1Index + j] = diag[j] * (diag[j] * wa1[j] / pnorm);
                    qtf[j] = ratio switch
                    {
                        >= p0001 => sum,
                        _ => qtf[j]
                    };
                }

                //
                //  Compute the QR factorization of the updated jacobian.
                //
                r1updt(n, n, ref r, lr, wa1, ref wa2, wa3, sIndex:rIndex, uIndex:wa1Index);
                r1mpyq(n, n, ref fjac, ldfjac, wa2, wa3, aIndex:fjacIndex);
                r1mpyq(1, n, ref qtf, 1, wa2, wa3, aIndex:qtfIndex);

                jeval = false;
            }
            //
            //  End of the inner loop.
            //
        }
        //
        //  End of the outer loop.
        //
    }

        
    // Last two ints before fcnData are the array baseline indices.
    public static int hybrd1(Func< int, double[], double[], int, int, int, fcnData > fcn, int n,
            double[] x, double[] fvec, double tol, double[] wa, int lwa )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    hybrd1() finds a zero of a system of N nonlinear equations. 
        //
        //  Discussion:
        //
        //    The purpose of HYBRD1 is to find a zero of a system of
        //    N nonlinear functions in N variables by a modification
        //    of the Powell hybrid method.  
        //
        //    This is done by using the more general nonlinear equation solver HYBRD.  
        //
        //    The user must provide FCN, which calculates the functions.
        //
        //    The jacobian is calculated by a forward-difference approximation.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 July 2014
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
        //       fcn is the name of the user-supplied subroutine which
        //         calculates the functions. fcn must be declared
        //         in an external statement in the user calling
        //         program, and should be written as follows.
        //
        //         subroutine fcn(n,x,fvec,iflag)
        //         integer n,iflag
        //         double precision x(n),fvec(n)
        //         ----------
        //         calculate the functions at x and
        //         return this vector in fvec.
        //         ---------
        //         return
        //         end
        //
        //         the value of iflag should not be changed by fcn unless
        //         the user wants to terminate execution of hybrd1.
        //         in this case set iflag to a negative integer.
        //
        //    Input, int N, the number of functions and variables.
        //
        //    Input/output, double X[N].  On input, an initial estimate of the solution.
        //    On output, the final estimate of the solution.
        //
        //    Output, double FVEC[N], the functions evaluated at the output X.
        //
        //    Input, double TOL, a nonnegative variable. tTermination occurs when the 
        //    algorithm estimates that the relative error between X and the solution 
        //    is at most TOL.
        //
        //       info is an integer output variable. if the user has
        //         terminated execution, info is set to the (negative)
        //         value of iflag. see description of fcn. otherwise,
        //         info is set as follows.
        //         info = 0   improper input parameters.
        //         info = 1   algorithm estimates that the relative error
        //                    between x and the solution is at most tol.
        //         info = 2   number of calls to fcn has reached or exceeded
        //                    200*(n+1).
        //         info = 3   tol is too small. no further improvement in
        //                    the approximate solution x is possible.
        //         info = 4   iteration is not making good progress.
        //
        //    Workspace, double WA[LWA].
        //
        //       lwa is a positive integer input variable not less than
        //         (n*(3*n+13))/2.
        //
    {
        int j;

        int info = 0;
        switch (n)
        {
            //
            //  Check the input.
            //
            case <= 0:
                return info;
        }

        switch (tol)
        {
            case <= 0.0:
                return info;
        }

        if (lwa < n * (3 * n + 13) / 2)
        {
            return info;
        }

        //
        //  Call HYBRD.
        //
        int maxfev = 200 * (n + 1);
        int ml = n - 1;
        int mu = n - 1;
        const double epsfcn = 0.0;
        for (j = 0; j < n; j++)
        {
            wa[j] = 1.0;
        }

        const int mode = 2;
        const double factor = 100.0;
        const int nprint = 0;
        const int nfev = 0;
        int lr = n * (n + 1) / 2;
        int index = 6 * n + lr;

        info = hybrd(fcn, n, x, fvec, tol, maxfev, ml, mu, epsfcn, wa, mode,
            factor, nprint, nfev, wa, n, wa, lr,
            wa, wa, wa, wa, wa, fjacIndex:index, rIndex:  + 6 * n, qtfIndex : n, wa1Index:  + 2 * n,
            wa2Index: + 3 * n, wa3Index:  + 4 * n, wa4Index:  + 5 * n);

        switch (info)
        {
            case 5:
                info = 4;
                break;
        }

        return info;
    }
}