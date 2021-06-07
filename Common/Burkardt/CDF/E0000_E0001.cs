using System;

namespace Burkardt.CDFLib
{
    public class E0000_E0001_Data
    {
        public int e0000_status;
        public double e0000_x;
        public double e0000_fx;
        public bool e0000_qleft;
        public bool e0000_qhi;

        public int e0001_status;
        public double e0001_x;
        public double e0001_fx;
        public bool e0001_qleft;
        public bool e0001_qhi;

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
            E000(1, ref data);
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

            double absstp;
            double abstol;
            double big;
            double fbig,
                fsmall,
                relstp,
                reltol,
                small,
                step,
                stpmul,
                xhi,
                xlb,
                xlo,
                xsave,
                xub,
                yy;
            int i99999
                ;
            bool qbdd, qcond, qdum1, qdum2, qincr, qlim;
            bool qup;

            switch (IENTRY)
            {
                case 0: goto DINVR;
                case 1: goto DSTINV;
            }

            DINVR:
            if (data.e0000_status > 0) goto S310;
            qcond = !qxmon(small, data.e0000_x, big);
            if (qcond)
            {
                throw new Exception(" SMALL, X, BIG not monotone in INVR");
            }

            xsave = data.e0000_x;
            //
            //     See that SMALL and BIG bound the zero and set QINCR
            //
            data.e0000_x = small;
            //
            //     GET-FUNCTION-VALUE
            //
            i99999 = 1;
            goto S300;
            S10:
            fsmall = data.e0000_fx;
            data.e0000_x = big;
            //
            //     GET-FUNCTION-VALUE
            //
            i99999 = 2;
            goto S300;
            S20:
            fbig = data.e0000_fx;
            qincr = fbig > fsmall;
            if (!qincr) goto S50;
            if (fsmall <= 0.0e0) goto S30;
            data.e0000_status = -1;
            data.e0000_qleft = data.e0000_qhi = true;
            return;
            S30:
            if (fbig >= 0.0e0) goto S40;
            data.e0000_status = -1;
            data.e0000_qleft = data.e0000_qhi = false
            return;
            S40:
            goto S80;
            S50:
            if (fsmall >= 0.0e0) goto S60;
            data.e0000_status = -1;
            data.e0000_qleft = true;
            data.e0000_qhi = false;
            return;
            S60:
            if (fbig <= 0.0e0) goto S70;
            data.e0000_status = -1;
            data.e0000_qleft = false;
            data.e0000_qhi = true;
            return;
            S80:
            S70:
            data.e0000_x = xsave;
            step = Math.Max(absstp, relstp * Math.Abs(data.e0000_x));
            //
            //      YY = F(X) - Y
            //     GET-FUNCTION-VALUE
            //
            i99999 = 3;
            goto S300;
            S90:
            yy = data.e0000_fx;
            if (!(yy == 0.0e0)) goto S100;
            data.e0000_status = 0;
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
            data.e0000_x = xub;
            //
            //     GET-FUNCTION-VALUE
            //
            i99999 = 4;
            goto S300;
            S130:
            yy = data.e0000_fx;
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
            data.e0000_status = -1;
            data.e0000_qleft = 0;
            data.e0000_qhi = !qincr;
            data.e0000_x = big;
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
            data.e0000_x = xlb;
            //
            //     GET-FUNCTION-VALUE
            //
            i99999 = 5;
            goto S300;
            S200:
            yy = data.e0000_fx;
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
            data.e0000_status = -1;
            data.e0000_qleft = 1;
            data.e0000_qhi = qincr;
            data.e0000_x = small;
            return;
            S240:
            S230:
            dstzr(xlb, xub, abstol, reltol);
            //
            //  IF WE REACH HERE, XLB AND XUB BOUND THE ZERO OF F.
            //
            data.e0000_status = 0;
            goto S260;
            S250:
            if (!(data.e0000_status == 1)) goto S290;
            S260:
            dzror(data.e0000_status, data.e0000_x, data.e0000_fx, xlo, xhi, qdum1, qdum2);
            if (!(data.e0000_status == 1)) goto S280;
            //
            //     GET-FUNCTION-VALUE
            //
            i99999 = 6;
            goto S300;
            S280:
            S270:
            goto S250;
            S290:
            data.e0000_x = xlo;
            data.e0000_status = 0;
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
            data.e0000_status = 1;
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

            double a;
            double abstol, b, c, d, fa, fb, fc, fd, fda;
            double fdb, m, mb, p, q, reltol, tol, w, xxhi, xxlo;
            int ext, i99999;
            bool first, qrzero;

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
            if (data.e0001_status > 0) goto S280;
            data.e0001_xlo = xxlo;
            data.e0001_xhi = xxhi;
            b = data.e0001_x = data.e0001_xlo;
            //
            //     GET-FUNCTION-VALUE
            //
            i99999 = 1;
            goto S270;
            S10:
            fb = data.e0001_fx;
            data.e0001_xlo = data.e0001_xhi;
            a = data.e0001_x = data.e0001_xlo;
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
            if (!(data.e0001_fx < 0.0e0)) goto S30;
            data.e0001_status = -1;
            data.e0001_qleft = data.e0001_fx < fb;
            data.e0001_qhi = 0;
            return;
            S40:
            S30:
            if (!(fb > 0.0e0)) goto S60;
            if (!(data.e0001_fx > 0.0e0)) goto S50;
            data.e0001_status = -1;
            data.e0001_qleft = data.e0001_fx > fb;
            data.e0001_qhi = 1;
            return;
            S60:
            S50:
            fa = data.e0001_fx;
            first = 1;
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
            first = 0;
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
            data.e0001_x = data.e0001_xlo;
            //
            //  GET-FUNCTION-VALUE
            //
            i99999 = 3;
            goto S270;
            S200:
            fb = data.e0001_fx;
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
            data.e0001_status = 0;
            goto S260;
            S250:
            data.e0001_status = -1;
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
            data.e0001_status = 1;
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