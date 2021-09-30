using System;
using Burkardt.Types;

namespace Burkardt.MinpackNS
{
    public class fcnData
    {
        public double[] fvec { get; set; }
        public int iflag;
    }

    public static partial class Minpack
    {
        
        // Last two ints before fcnData are the array index baseline for each array
        public static void fdjac1(Func < int, double[], double[], int, int, int, fcnData > fcn,
        int n, double[] x, double[] fvec, double[] fjac, int ldfjac, ref int iflag,
        int ml, int mu, double epsfcn, double[] wa1, double[] wa2, int fjacIndex = 0, int wa1Index = 0, int wa2Index = 0 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    fdjac1() estimates an N by N Jacobian matrix using forward differences.
        //
        //  Discussion:
        //
        //    This function computes a forward-difference approximation
        //    to the N by N jacobian matrix associated with a specified
        //    problem of N functions in N variables. 
        //
        //    If the jacobian has a banded form, then function evaluations are saved 
        //    by only approximating the nonzero terms.
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
        //         ----------
        //         return
        //         end
        //
        //         the value of iflag should not be changed by fcn unless
        //         the user wants to terminate execution of fdjac1.
        //         in this case set iflag to a negative integer.
        //
        //    Input, int N, the number of functions and variables.
        //
        //    Input, double X[N], the evaluation point.
        //
        //    Input, double FVEC[N], the functions evaluated at X.
        //
        //    Output, double FJAC[N*N], the approximate jacobian matrix at X.
        //
        //       ldfjac is a positive integer input variable not less than n
        //         which specifies the leading dimension of the array fjac.
        //
        //       iflag is an integer variable which can be used to terminate
        //         the execution of fdjac1. see description of fcn.
        //
        //       ml is a nonnegative integer input variable which specifies
        //         the number of subdiagonals within the band of the
        //         jacobian matrix. if the jacobian is not banded, set
        //         ml to at least n - 1.
        //
        //       epsfcn is an input variable used in determining a suitable
        //         step length for the forward-difference approximation. this
        //         approximation assumes that the relative errors in the
        //         functions are of the order of epsfcn. if epsfcn is less
        //         than the machine precision, it is assumed that the relative
        //         errors in the functions are of the order of the machine
        //         precision.
        //
        //       mu is a nonnegative integer input variable which specifies
        //         the number of superdiagonals within the band of the
        //         jacobian matrix. if the jacobian is not banded, set
        //         mu to at least n - 1.
        //
        //       wa1 and wa2 are work arrays of length n. if ml + mu + 1 is at
        //         least n, then the jacobian is considered dense, and wa2 is
        //         not referenced.
        {
            double eps;
            double epsmch;
            double h;
            int i;
            int j;
            int k;
            int msum;
            double temp;
            //
            //  EPSMCH is the machine precision.
            //
            epsmch = typeMethods.r8_epsilon();

            eps = Math.Sqrt(Math.Max(epsfcn, epsmch));
            msum = ml + mu + 1;
            //
            //  Computation of dense approximate jacobian.
            //
            if (n <= msum)
            {
                for (j = 0; j < n; j++)
                {
                    temp = x[j];
                    h = eps * Math.Abs(temp);
                    if (h == 0.0)
                    {
                        h = eps;
                    }

                    x[j] = temp + h;
                    fcnData r = fcn(n, x, wa1, iflag, 0, wa2Index);
                    iflag = r.iflag;
                    fvec = r.fvec;
                    if (iflag < 0)
                    {
                        break;
                    }

                    x[j] = temp;
                    for (i = 0; i < n; i++)
                    {
                        fjac[fjacIndex + (i + j * ldfjac)] = (wa1[wa1Index + (i)] - fvec[i]) / h;
                    }
                }
            }
            //
            //  Computation of a banded approximate jacobian.
            //
            else
            {
                for (k = 0; k < msum; k++)
                {
                    for (j = k; j < n; j = j + msum)
                    {
                        wa2[wa2Index + (j)] = x[j];
                        h = eps * Math.Abs(wa2[wa2Index + (j)]);
                        if (h == 0.0)
                        {
                            h = eps;
                        }

                        x[j] = wa2[wa2Index + (j)] + h;
                    }

                    fcn(n, x, wa1, iflag, 0, wa1Index);
                    if (iflag < 0)
                    {
                        break;
                    }

                    for (j = k; j < n; j = j + msum)
                    {
                        x[j] = wa2[wa2Index + (j)];
                        h = eps * Math.Abs(wa2[wa2Index + (j)]);
                        if (h == 0.0)
                        {
                            h = eps;
                        }

                        for (i = 0; i < n; i++)
                        {
                            if (j - mu <= i && i <= j + ml)
                            {
                                fjac[fjacIndex + (i + j * ldfjac)] = (wa1[wa1Index + (i)] - fvec[i]) / h;
                            }
                            else
                            {
                                fjac[fjacIndex + (i + j * ldfjac)] = 0.0;
                            }
                        }
                    }
                }
            }
        }
        
        public static void fdjac2(Func < int, int, double[], double[], int, fcnData > fcn,
        int m, int n, double[] x, double[] fvec, double[] fjac, int ldfjac,
        ref int iflag, double epsfcn, double[] wa )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    fdjac2() estimates an M by N Jacobian matrix using forward differences.
        //
        //  Discussion:
        //
        //    This function computes a forward-difference approximation
        //    to the M by N jacobian matrix associated with a specified
        //    problem of M functions in N variables.
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
        //       fcn is the name of the user-supplied subroutine which
        //         calculates the functions. fcn must be declared
        //         in an external statement in the user calling
        //         program, and should be written as follows.
        //
        //         subroutine fcn(m,n,x,fvec,iflag)
        //         integer m,n,iflag
        //         double precision x(n),fvec(m)
        //         ----------
        //         calculate the functions at x and
        //         return this vector in fvec.
        //         ----------
        //         return
        //         end
        //
        //         the value of iflag should not be changed by fcn unless
        //         the user wants to terminate execution of fdjac2.
        //         in this case set iflag to a negative integer.
        //
        //    Input, int M, the number of functions.
        //
        //    Input, int N, the number of variables.  N must not exceed M.
        //
        //    Input, double X[N], the point at which the jacobian is to be estimated.
        //
        //       fvec is an input array of length m which must contain the
        //         functions evaluated at x.
        //
        //       fjac is an output m by n array which contains the
        //         approximation to the jacobian matrix evaluated at x.
        //
        //       ldfjac is a positive integer input variable not less than m
        //         which specifies the leading dimension of the array fjac.
        //
        //       iflag is an integer variable which can be used to terminate
        //         the execution of fdjac2. see description of fcn.
        //
        //       epsfcn is an input variable used in determining a suitable
        //         step length for the forward-difference approximation. this
        //         approximation assumes that the relative errors in the
        //         functions are of the order of epsfcn. if epsfcn is less
        //         than the machine precision, it is assumed that the relative
        //         errors in the functions are of the order of the machine
        //         precision.
        //
        //    Workspace, double WA[M].
        //
        {
            double eps;
            double epsmch;
            double h;
            int i;
            int j;
            double temp;
            //
            //  EPSMCH is the machine precision.
            //
            epsmch = typeMethods.r8_epsilon();
            eps = Math.Sqrt(Math.Max(epsfcn, epsmch));

            for (j = 0; j < n; j++)
            {
                temp = x[j];
                if (temp == 0.0)
                {
                    h = eps;
                }
                else
                {
                    h = eps * Math.Abs(temp);
                }

                x[j] = temp + h;
                fcnData r = fcn(m, n, x, wa, iflag);
                iflag = r.iflag;
                fvec = r.fvec;
                if (iflag < 0)
                {
                    break;
                }

                x[j] = temp;
                for (i = 0; i < m; i++)
                {
                    fjac[i + j * ldfjac] = (wa[i] - fvec[i]) / h;
                }
            }
        }
    }
}