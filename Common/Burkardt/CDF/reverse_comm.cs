﻿using System;

namespace Burkardt.CDFLib
{
    public static partial class CDF
    {
        static void dzror(ref int status, ref double x, ref double fx, ref double xlo,
                ref double xhi, ref bool qleft, ref bool qhi)

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
            double zabstl = 0;
            double zreltl = 0;
            double zxhi = 0;
            double zxlo = 0;
            E0001(0, ref status, ref x, ref fx, ref xlo, ref xhi, ref qleft, ref qhi, ref zabstl, ref zreltl, ref zxhi,
                ref zxlo);
        }

        static void E0000(int IENTRY, ref int status, ref double x, ref double fx,
                ref bool qleft, ref bool qhi, ref double zabsst,
                ref double zabsto, ref double zbig, ref double zrelst,
                ref double zrelto, ref double zsmall, ref double zstpmu)

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
            bool qxmon(double zx, double zy, double zz)
            {
                return ((zx) <= (zy) && (zy) <= (zz));
            }

            double absstp = 0;
            double abstol = 0;
            double big = 0;
            double fbig = 0,
                fsmall = 0,
                relstp = 0,
                reltol = 0,
                small = 0,
                step = 0,
                stpmul = 0,
                xhi = 0,
                xlb = 0,
                xlo = 0,
                xsave = 0,
                xub = 0,
                yy;
            int i99999 = 0;
            bool qbdd, qcond, qdum1 = false, qdum2 = false, qincr = false, qlim;
            bool qup = false;

            switch (IENTRY)
            {
                case 0: goto DINVR;
                case 1: goto DSTINV;
            }

            DINVR:
            if (status > 0) goto S310;
            qcond = !qxmon(small, x, big);
            if (qcond)
            {
                throw new Exception(" SMALL, X, BIG not monotone in INVR");
            }

            xsave = x;
            //
            //     See that SMALL and BIG bound the zero and set QINCR
            //
            x = small;
            //
            //     GET-FUNCTION-VALUE
            //
            i99999 = 1;
            goto S300;
            S10:
            fsmall = fx;
            x = big;
            //
            //     GET-FUNCTION-VALUE
            //
            i99999 = 2;
            goto S300;
            S20:
            fbig = fx;
            qincr = fbig > fsmall;
            if (!qincr) goto S50;
            if (fsmall <= 0.0e0) goto S30;
            status = -1;
            qleft = qhi = true;
            return;
            S30:
            if (fbig >= 0.0e0) goto S40;
            status = -1;
            qleft = qhi = false;
            return;
            S40:
            goto S80;
            S50:
            if (fsmall >= 0.0e0) goto S60;
            status = -1;
            qleft = true;
            qhi = false;
            return;
            S60:
            if (fbig <= 0.0e0) goto S70;
            status = -1;
            qleft = false;
            qhi = true;
            return;
            S80:
            S70:
            x = xsave;
            step = Math.Max(absstp, relstp * Math.Abs(x));
            //
            //      YY = F(X) - Y
            //     GET-FUNCTION-VALUE
            //
            i99999 = 3;
            goto S300;
            S90:
            yy = fx;
            if (!(yy == 0.0e0)) goto S100;
            status = 0;
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
            x = xub;
            //
            //     GET-FUNCTION-VALUE
            //
            i99999 = 4;
            goto S300;
            S130:
            yy = fx;
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
            status = -1;
            qleft = false;
            qhi = !qincr;
            x = big;
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
            x = xlb;
            //
            //     GET-FUNCTION-VALUE
            //
            i99999 = 5;
            goto S300;
            S200:
            yy = fx;
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
            status = -1;
            qleft = true;
            qhi = qincr;
            x = small;
            return;
            S240:
            S230:
            dstzr(xlb, xub, abstol, reltol);
            //
            //  IF WE REACH HERE, XLB AND XUB BOUND THE ZERO OF F.
            //
            status = 0;
            goto S260;
            S250:
            if (!(status == 1)) goto S290;
            S260:
            dzror(ref status, ref x, ref fx, ref xlo, ref xhi, ref qdum1, ref qdum2);
            if (!(status == 1)) goto S280;
            //
            //     GET-FUNCTION-VALUE
            //
            i99999 = 6;
            goto S300;
            S280:
            S270:
            goto S250;
            S290:
            x = xlo;
            status = 0;
            return;
            DSTINV:
            small = zsmall;
            big = zbig;
            absstp = zabsst;
            relstp = zrelst;
            stpmul = zstpmu;
            abstol = zabsto;
            reltol = zrelto;
            return;
            S300:
            //
            //     TO GET-FUNCTION-VALUE
            //
            status = 1;
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

        static void E0001(int IENTRY, ref int status, ref double x, ref double fx,
                ref double xlo, ref double xhi, ref bool qleft,
                ref bool qhi, ref double zabstl, ref double zreltl,
                ref double zxhi, ref double zxlo)

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
            double abstol = 0,
                b = 0,
                c = 0,
                d = 0,
                fa = 0,
                fb = 0,
                fc = 0,
                fd = 0,
                fda = 0;
            double fdb = 0,
                m = 0,
                mb = 0,
                p = 0,
                q = 0,
                reltol = 0,
                tol = 0,
                w = 0,
                xxhi = 0,
                xxlo = 0;
            int ext = 0,
                i99999 = 0;
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
            if (status > 0) goto S280;
            xlo = xxlo;
            xhi = xxhi;
            b = x = xlo;
            //
            //     GET-FUNCTION-VALUE
            //
            i99999 = 1;
            goto S270;
            S10:
            fb = fx;
            xlo = xhi;
            a = x = xlo;
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
            if (!(fx < 0.0e0)) goto S30;
            status = -1;
            qleft = fx < fb;
            qhi = false;
            return;
            S40:
            S30:
            if (!(fb > 0.0e0)) goto S60;
            if (!(fx > 0.0e0)) goto S50;
            status = -1;
            qleft = fx > fb;
            qhi = true;
            return;
            S60:
            S50:
            fa = fx;
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
            xlo = c;
            b = xlo;
            fb = fc;
            c = a;
            fc = fa;
            S100:
            tol = ftol(xlo);
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
            xlo = b;
            x = xlo;
            //
            //  GET-FUNCTION-VALUE
            //
            i99999 = 3;
            goto S270;
            S200:
            fb = fx;
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
            xhi = c;
            qrzero = (fc >= 0.0e0 && fb <= 0.0e0) || (fc < 0.0e0 && fb >= 0.0e0);
            if (!qrzero) goto S250;
            status = 0;
            goto S260;
            S250:
            status = -1;
            S260:
            return;
            DSTZR:
            xxlo = zxlo;
            xxhi = zxhi;
            abstol = zabstl;
            reltol = zreltl;
            return;
            S270:
            //
            //     TO GET-FUNCTION-VALUE
            //
            status = 1;
            return;
            S280:
            switch ((int) i99999)
            {
                case 1: goto S10;
                case 2: goto S20;
                case 3: goto S200;
                default: break;
            }
# undef ftol
        }
    }
}