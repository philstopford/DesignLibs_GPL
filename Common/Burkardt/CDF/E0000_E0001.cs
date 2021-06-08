using System;

namespace Burkardt.CDFLib
{
    public class E0000_E0001_Data
    {
        public int status;
        public double x;
        public double fx;
        public bool qleft;
        public bool qhi;

        public double e0000_zabsst;
        public double e0000_zabsto;
        public double e0000_zbig;
        public double e0000_zrelst;
        public double e0000_zrelto;
        public double e0000_zsmall;
        public double e0000_zstpmu;

        public double e0001_xlo;
        public double e0001_xhi;
        public double e0001_zabstl;
        public double e0001_zreltl;
        public double e0001_zxhi;
        public double e0001_zxlo;
    }

    public static partial class E0000E0001
    {
        public static void dinvr(ref E0000_E0001_Data data, double x_, double fx_, bool qleft_, bool qhi_)
        {
            data.x = x_;
            data.fx = fx_;
            data.qleft = qleft_;
            data.qhi = qhi_;
            dinvr(ref data);
        }
        
        public static void dinvr(ref E0000_E0001_Data data)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DINVR bounds the zero of the function and invokes DZROR.
            //
            //  Discussion:
            //
            //    This routine seeks to find bounds on a root of the function and
            //    invokes ZROR to perform the zero finding.  STINVR must have been
            //    called before this routine in order to set its parameters.
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
            //  Reference:
            //
            //    J C P Bus and T J Dekker,
            //    Two Efficient Algorithms with Guaranteed Convergence for
            //    Finding a Zero of a Function,
            //    ACM Transactions on Mathematical Software,
            //    Volume 1, Number 4, pages 330-345, 1975.
            //
            //  Parameters:
            //
            //    Input/output, integer STATUS.  At the beginning of a zero finding
            //    problem, STATUS should be set to 0 and INVR invoked.  The value
            //    of parameters other than X will be ignored on this call.
            //    If INVR needs the function to be evaluated, it will set STATUS to 1
            //    and return.  The value of the function should be set in FX and INVR
            //    again called without changing any of its other parameters.
            //    If INVR finishes without error, it returns with STATUS 0, and X an
            //    approximate root of F(X).
            //    If INVR cannot bound the function, it returns a negative STATUS and
            //    sets QLEFT and QHI.
            //
            //    Output, double precision X, the value at which F(X) is to be evaluated.
            //
            //    Input, double precision FX, the value of F(X) calculated by the user
            //    on the previous call, when INVR returned with STATUS = 1.
            //
            //    Output, logical QLEFT, is defined only if QMFINV returns FALSE.  In that
            //    case, QLEFT is TRUE if the stepping search terminated unsucessfully
            //    at SMALL, and FALSE if the search terminated unsucessfully at BIG.
            //
            //    Output, logical QHI, is defined only if QMFINV returns FALSE.  In that
            //    case, it is TRUE if Y < F(X) at the termination of the search and FALSE
            //    if F(X) < Y.
            //
        {
            E0000(0, ref data);
        }

        public static void dstinv(ref E0000_E0001_Data data, double zsmall_, double zbig_, double zabsst_,
            double zrelst_, double zstpmu_, double zabsto_, double zrelto_)
        {
            data.e0000_zsmall = zsmall_;
            data.e0000_zbig = zbig_;
            data.e0000_zabsst = zabsst_;
            data.e0000_zrelst = zrelst_;
            data.e0000_zstpmu = zstpmu_;
            data.e0000_zabsto = zabsto_;
            data.e0000_zrelto = zrelto_;
            
            dstinv(ref data);
        }
        
        public static void dstinv(ref E0000_E0001_Data data)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DSTINV seeks a value X such that F(X) = Y.
            //
            //  Discussion:
            //
            //    DSTINV is the double precision set inverse finder.
            //    It uses reverse communication.
            //
            //    Given a monotone function F and a value Y, it finds X
            //    such that F(X) = Y.
            //
            //    This routine sets quantities needed by INVR.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 February 2021
            //
            //  More Precise Description of INVR -
            //
            //     F must be a monotone function, the results of QMFINV are
            //     otherwise undefined.  QINCR must be .TRUE. if F is non-
            //     decreasing and .FALSE. if F is non-increasing.
            //     QMFINV will return .TRUE. if and only if F(SMALL) and
            //     F(BIG) bracket Y, i. e.,
            //          QINCR is .TRUE. and F(SMALL)<=Y<=F(BIG) or
            //          QINCR is .FALSE. and F(BIG)<=Y<=F(SMALL)
            //     if QMFINV returns .TRUE., then the X returned satisfies
            //     the following condition.  let
            //               TOL(X) = MAX(ABSTOL,RELTOL*ABS(X))
            //     then if QINCR is .TRUE.,
            //          F(X-TOL(X)) <= Y <= F(X+TOL(X))
            //     and if QINCR is .FALSE.
            //          F(X-TOL(X)) .GE. Y .GE. F(X+TOL(X))
            //
            //                              Method
            //     Compares F(X) with Y for the input value of X then uses QINCR
            //     to determine whether to step left or right to bound the
            //     desired x.  the initial step size is
            //          MAX(ABSSTP,RELSTP*ABS(S)) for the input value of X.
            //     Iteratively steps right or left until it bounds X.
            //     At each step which doesn't bound X, the step size is doubled.
            //     The routine is careful never to step beyond SMALL or BIG.  If
            //     it hasn't bounded X at SMALL or BIG, QMFINV returns .FALSE.
            //     after setting QLEFT and QHI.
            //     If X is successfully bounded then Algorithm R of the paper
            //     'Two Efficient Algorithms with Guaranteed Convergence for
            //     Finding a Zero of a Function' by J. C. P. Bus and
            //     T. J. Dekker in ACM Transactions on Mathematical
            //     Software, Volume 1, No. 4 page 330 (DEC. '75) is employed
            //     to find the zero of the function F(X)-Y. This is routine
            //     QRZERO.
            //
            //  Parameters:
            //
            //    double SMALL --> The left endpoint of the interval to be
            //          searched for a solution.
            //
            //    double BIG --> The right endpoint of the interval to be
            //          searched for a solution.
            //
            //    double ABSSTP, RELSTP --> The initial step size in the search
            //          is MAX(ABSSTP,RELSTP*ABS(X)). See algorithm.
            //
            //    double STPMUL --> When a step doesn't bound the zero, the step
            //                size is multiplied by STPMUL and another step
            //                taken.  A popular value is 2.0
            //
            //    double ABSTOL, RELTOL --> Two numbers that determine the accuracy
            //          of the solution.  See function for a precise definition.
            //
        {
            E0000(1, ref data);
        }

        public static void dstzr(ref E0000_E0001_Data data, double zxlo, double zxhi, double zabstl, double zreltl)
        {
            data.e0001_zxlo = zxlo;
            data.e0001_zxhi = zxhi;
            data.e0001_zabstl = zabstl;
            data.e0001_zreltl = zreltl;
            dstzr(ref data);
        }

        public static void dstzr(ref E0000_E0001_Data data)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DSTXR sets quantities needed by the zero finder.
            //
            //  Discussion:
            //
            //     Double precision SeT ZeRo finder - Reverse communication version
            //                              Function
            //     Sets quantities needed by ZROR.  The function of ZROR
            //     and the quantities set is given here.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 February 2021
            //
            //  Concise Description
            //
            //    Given a function F
            //     find XLO such that F(XLO) = 0.
            //          More Precise Description -
            //     Input condition. F is a double function of a single
            //     double argument and XLO and XHI are such that
            //          F(XLO)*F(XHI)  <=  0.0
            //     If the input condition is met, QRZERO returns .TRUE.
            //     and output values of XLO and XHI satisfy the following
            //          F(XLO)*F(XHI)  <= 0.
            //          ABS(F(XLO)  <= ABS(F(XHI)
            //          ABS(XLO-XHI)  <= TOL(X)
            //     where
            //          TOL(X) = MAX(ABSTOL,RELTOL*ABS(X))
            //     If this algorithm does not find XLO and XHI satisfying
            //     these conditions then QRZERO returns .FALSE.  This
            //     implies that the input condition was not met.
            //
            //  Parameters:
            //
            //     XLO --> The left endpoint of the interval to be
            //           searched for a solution.
            //                    XLO is DOUBLE PRECISION
            //     XHI --> The right endpoint of the interval to be
            //           for a solution.
            //                    XHI is DOUBLE PRECISION
            //     ABSTOL, RELTOL --> Two numbers that determine the accuracy
            //                      of the solution.  See function for a
            //                      precise definition.
            //                    ABSTOL is DOUBLE PRECISION
            //                    RELTOL is DOUBLE PRECISION
            //
            //                              Method
            //     Algorithm R of the paper 'Two Efficient Algorithms with
            //     Guaranteed Convergence for Finding a Zero of a Function'
            //     by J. C. P. Bus and T. J. Dekker in ACM Transactions on
            //     Mathematical Software, Volume 1, no. 4 page 330
            //     (Dec. '75) is employed to find the zero of F(X)-Y.
            //
        {
            E0001(1, ref data);
        }

        public static void dzror(ref E0000_E0001_Data data, double x_)
        {
            
        }
        
        public static void dzror(ref E0000_E0001_Data data)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DZROR seeks the zero of a function using reverse communication.
            //
            //  Discussion:
            //
            //    Performs the zero finding.  STZROR must have been called before
            //    this routine in order to set its parameters.
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
            //    int STATUS <--> At the beginning of a zero finding problem, STATUS
            //                 should be set to 0 and ZROR invoked.  (The value
            //                 of other parameters will be ignored on this call.)
            //
            //                 When ZROR needs the function evaluated, it will set
            //                 STATUS to 1 and return.  The value of the function
            //                 should be set in FX and ZROR again called without
            //                 changing any of its other parameters.
            //
            //                 When ZROR has finished without error, it will return
            //                 with STATUS 0.  In that case (XLO,XHI) bound the answe
            //
            //                 If ZROR finds an error (which implies that F(XLO)-Y an
            //                 F(XHI)-Y have the same sign, it returns STATUS -1.  In
            //                 this case, XLO and XHI are undefined.
            //
            //    double X <-- The value of X at which F(X) is to be evaluated.
            //
            //    double FX --> The value of F(X) calculated when ZROR returns with
            //            STATUS = 1.
            //
            //    double XLO <-- When ZROR returns with STATUS = 0, XLO bounds the
            //             inverval in X containing the solution below.
            //
            //    double XHI <-- When ZROR returns with STATUS = 0, XHI bounds the
            //             inverval in X containing the solution above.
            //
            //    bool QLEFT <-- .TRUE. if the stepping search terminated unsucessfully
            //                at XLO.  If it is .FALSE. the search terminated
            //                unsucessfully at XHI.
            //
            //    bool QHI <-- .TRUE. if F(X) .GT. Y at the termination of the
            //              search and .FALSE. if F(X) < Y at the
            //              termination of the search.
            //
        {
            E0001(0, ref data);
        }

        static bool qxmon(double zx, double zy, double zz)
        {
            return ((zx) <= (zy) && (zy) <= (zz));
        }


        public static void E0000(int IENTRY, ref E0000_E0001_Data data)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    E0000 is a reverse-communication zero bounder.
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
        {

            double absstp = 0;
            double abstol = 0;
            double big = data.x;
            double fbig = 0,
                fsmall = 0,
                relstp = 0,
                reltol = 0,
                small = data.x,
                step = 0,
                stpmul = 0,
                xlb = 0,
                xsave = 0,
                xub = 0,
                yy = 0;
            int i99999 = 0;
            bool qbdd = false, qcond = false, qdum1 = false, qdum2 = false, qincr = false, qlim = false;
            bool qup = false;

            switch (IENTRY)
            {
                case 0: goto DINVR;
                case 1: goto DSTINV;
            }

            DINVR:
            if (data.status > 0) goto S310;
            qcond = !qxmon(small, data.x, big);
            if (qcond)
            {
                throw new Exception(" SMALL, X, BIG not monotone in INVR");
            }

            xsave = data.x;
            //
            //     See that SMALL and BIG bound the zero and set QINCR
            //
            data.x = small;
            //
            //     GET-FUNCTION-VALUE
            //
            i99999 = 1;
            goto S300;
            S10:
            fsmall = data.fx;
            data.x = big;
            //
            //     GET-FUNCTION-VALUE
            //
            i99999 = 2;
            goto S300;
            S20:
            fbig = data.fx;
            qincr = fbig > fsmall;
            if (!qincr) goto S50;
            if (fsmall <= 0.0e0) goto S30;
            data.status = -1;
            data.qleft = data.qhi = true;
            return;
            S30:
            if (fbig >= 0.0e0) goto S40;
            data.status = -1;
            data.qleft = data.qhi = false;
            return;
            S40:
            goto S80;
            S50:
            if (fsmall >= 0.0e0) goto S60;
            data.status = -1;
            data.qleft = true;
            data.qhi = false;
            return;
            S60:
            if (fbig <= 0.0e0) goto S70;
            data.status = -1;
            data.qleft = false;
            data.qhi = true;
            return;
            S80:
            S70:
            data.x = xsave;
            step = Math.Max(absstp, relstp * Math.Abs(data.x));
            //
            //      YY = F(X) - Y
            //     GET-FUNCTION-VALUE
            //
            i99999 = 3;
            goto S300;
            S90:
            yy = data.fx;
            if (!(yy == 0.0e0)) goto S100;
            data.status = 0;
            //  qok = 1;
            return;
            S100:
            qup = (qincr && yy < 0.0e0) || (!qincr && yy > 0.0e0);
            //
            //     HANDLE CASE IN WHICH WE MUST STEP HIGHER
            //
            if (!qup) goto S170;
            xlb = xsave;
            xub = Math.Min(xlb + step, big);
            goto S120;
            S110:
            if (qcond) goto S150;
            S120:
            //
            //      YY = F(XUB) - Y
            //
            data.x = xub;
            //
            //     GET-FUNCTION-VALUE
            //
            i99999 = 4;
            goto S300;
            S130:
            yy = data.fx;
            qbdd = (qincr && yy >= 0.0e0) || (!qincr && yy <= 0.0e0);
            qlim = xub >= big;
            qcond = qbdd || qlim;
            if (qcond) goto S140;
            step = stpmul * step;
            xlb = xub;
            xub = Math.Min(xlb + step, big);
            S140:
            goto S110;
            S150:
            if (!(qlim && !qbdd)) goto S160;
            data.status = -1;
            data.qleft = false;
            data.qhi = !qincr;
            data.x = big;
            return;
            S160:
            goto S240;
            S170:
            //
            //     HANDLE CASE IN WHICH WE MUST STEP LOWER
            //
            xub = xsave;
            xlb = Math.Max(xub - step, small);
            goto S190;
            S180:
            if (qcond) goto S220;
            S190:
            //
            //      YY = F(XLB) - Y
            //
            data.x = xlb;
            //
            //     GET-FUNCTION-VALUE
            //
            i99999 = 5;
            goto S300;
            S200:
            yy = data.fx;
            qbdd = (qincr && yy <= 0.0e0) || (!qincr && yy >= 0.0e0);
            qlim = xlb <= small;
            qcond = qbdd || qlim;
            if (qcond) goto S210;
            step = stpmul * step;
            xub = xlb;
            xlb = Math.Max(xub - step, small);
            S210:
            goto S180;
            S220:
            if (!(qlim && !qbdd)) goto S230;
            data.status = -1;
            data.qleft = true;
            data.qhi = qincr;
            data.x = small;
            return;
            S240:
            S230:
            dstzr(ref data, xlb, xub, abstol, reltol);
            //
            //  IF WE REACH HERE, XLB AND XUB BOUND THE ZERO OF F.
            //
            data.status = 0;
            goto S260;
            S250:
            if (!(data.status == 1)) goto S290;
            S260:
            data.qleft = qdum1;
            data.qhi = qdum2;
            dzror(ref data);
            if (!(data.status == 1)) goto S280;
            //
            //     GET-FUNCTION-VALUE
            //
            i99999 = 6;
            goto S300;
            S280:
            S270:
            goto S250;
            S290:
            data.x = data.e0001_xlo;
            data.status = 0;
            return;
            DSTINV:
            small = data.e0000_zsmall;
            big = data.e0000_zbig;
            absstp = data.e0000_zabsst;
            relstp = data.e0000_zrelst;
            stpmul = data.e0000_zstpmu;
            abstol = data.e0000_zabsto;
            reltol = data.e0000_zrelto;
            return;
            S300:
            //
            //     TO GET-FUNCTION-VALUE
            //
            data.status = 1;
            return;
            S310:
            switch ((int) i99999)
            {
                case 1: goto S10;
                case 2: goto S20;
                case 3: goto S90;
                case
                    4: goto S130;
                case 5: goto S200;
                case 6: goto S270;
                default: break;
            }
        }

        public static double fifdsign ( double mag, double sign )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FIFDSIGN transfers the sign of the variable "sign" to the variable "mag"
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
            //  mag     -     magnitude
            //  sign    -     sign to be transfered
            //
        {
            if (mag < 0) mag = -mag;
            if (sign < 0) mag = -mag;
            return mag;

        }

        public static void E0001(int IENTRY, ref E0000_E0001_Data data)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    E00001 is a reverse-communication zero finder.
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
        {

            double a = 0;
            double abstol = 0, b = 0, c = 0, d = 0, fa = 0, fb = 0, fc = 0, fd = 0, fda = 0;
            double fdb = 0, m = 0, mb = 0, p = 0, q = 0, reltol = 0, tol = 0, w = 0, xxhi = 0, xxlo = 0;
            int ext = 0, i99999 = 0;
            bool first = false, qrzero = false;

            double ftol(double zx)
            {
                return (0.5e0 * Math.Max(abstol, reltol * Math.Abs((zx))));
            }

            switch (IENTRY)
            {
                case 0:
                    goto DZROR;
                case 1:
                    goto DSTZR;
            }

            DZROR:
            if (data.status > 0) goto S280;
            data.e0001_xlo = xxlo;
            data.e0001_xhi = xxhi;
            b = data.x = data.e0001_xlo;
            //
            //     GET-FUNCTION-VALUE
            //
            i99999 = 1;
            goto S270;
            S10:
            fb = data.fx;
            data.e0001_xlo = data.e0001_xhi;
            a = data.x = data.e0001_xlo;
            //
            //     GET-FUNCTION-VALUE
            //
            i99999 = 2;
            goto S270;
            S20:
            //
            //     Check that F(ZXLO) < 0 < F(ZXHI)  or
            //                F(ZXLO) > 0 > F(ZXHI)
            //
            if (!(fb < 0.0e0)) goto S40;
            if (!(data.fx < 0.0e0)) goto S30;
            data.status = -1;
            data.qleft = data.fx < fb;
            data.qhi = false;
            return;
            S40:
            S30:
            if (!(fb > 0.0e0)) goto S60;
            if (!(data.fx > 0.0e0)) goto S50;
            data.status = -1;
            data.qleft = data.fx > fb;
            data.qhi = true;
            return;
            S60:
            S50:
            fa = data.fx;
            first = true;
            S70:
            c = a;
            fc = fa;
            ext = 0;
            S80:
            if (!(Math.Abs(fc) < Math.Abs(fb))) goto S100;
            if (!(c != a)) goto S90;
            d = a;
            fd = fa;
            S90:
            a = b;
            fa = fb;
            data.e0001_xlo = c;
            b = data.e0001_xlo;
            fb = fc;
            c = a;
            fc = fa;
            S100:
            tol = ftol(data.e0001_xlo);
            m = (c + b) * .5e0;
            mb = m - b;
            if (!(Math.Abs(mb) > tol)) goto S240;
            if (!(ext > 3)) goto S110;
            w = mb;
            goto S190;
            S110:
            tol = fifdsign(tol, mb);
            p = (b - a) * fb;
            if (!first) goto S120;
            q = fa - fb;
            first = false;
            goto S130;
            S120:
            fdb = (fd - fb) / (d - b);
            fda = (fd - fa) / (d - a);
            p = fda * p;
            q = fdb * fa - fda * fb;
            S130:
            if (!(p < 0.0e0)) goto S140;
            p = -p;
            q = -q;
            S140:
            if (ext == 3) p = p * 2.0e0;
            if (!(p * 1.0e0 == 0.0e0 || p <= q * tol)) goto S150;
            w = tol;
            goto S180;
            S150:
            if (!(p < mb * q)) goto S160;
            w = p / q;
            goto S170;
            S160:
            w = mb;
            S190:
            S180:
            S170:
            d = a;
            fd = fa;
            a = b;
            fa = fb;
            b = b + w;
            data.e0001_xlo = b;
            data.x = data.e0001_xlo;
            //
            //  GET-FUNCTION-VALUE
            //
            i99999 = 3;
            goto S270;
            S200:
            fb = data.fx;
            if (!(fc * fb >= 0.0e0)) goto S210;
            goto S70;
            S210:
            if (!(w == mb)) goto S220;
            ext = 0;
            goto S230;
            S220:
            ext = ext + 1;
            S230:
            goto S80;
            S240:
            data.e0001_xhi = c;
            qrzero = (fc >= 0.0e0 && fb <= 0.0e0) || (fc < 0.0e0 && fb >= 0.0e0);
            if (!qrzero) goto S250;
            data.status = 0;
            goto S260;
            S250:
            data.status = -1;
            S260:
            return;
            DSTZR:
            xxlo = data.e0001_zxlo;
            xxhi = data.e0001_zxhi;
            abstol = data.e0001_zabstl;
            reltol = data.e0001_zreltl;
            return;
            S270:
            //
            //     TO GET-FUNCTION-VALUE
            //
            data.status = 1;
            return;
            S280:
            switch ((int) i99999)
            {
                case 1: goto S10;
                case 2: goto S20;
                case 3: goto S200;
                default: break;
            }
        }

    }
}