using System;

namespace Burkardt.CDFLib
{
    public static partial class CDF
    {
        static void E0000(int IENTRY, ref E0000Data e0000Data )

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

            double absstp = double.NaN;
            double abstol = double.NaN;
            double big = double.NaN;
            double fbig = double.NaN,
            fsmall = double.NaN,relstp = double.NaN,reltol = double.NaN,small = double.NaN,step = double.NaN,stpmul = double.NaN,xhi = double.NaN,
            xlb = double.NaN,xlo = double.NaN,xsave = double.NaN,xub = double.NaN,yy = double.NaN;
            int i99999 = 0;
            bool qbdd = false, qcond = false, qdum1 = false, qdum2 = false, qincr = false, qlim = false;
            bool qup = false;

            E0001Data e0001Data = new E0001Data() { };

            bool qxmon(double zx, double zy, double zz)
            {
                if (zx == Double.NaN)
                    zx = 0;
                if (zy == Double.NaN)
                    zy = 0;
                if (zz == Double.NaN)
                    zz = 0;
                return !((zx) <= (zy) && (zy) <= (zz));
            }

            switch (IENTRY)
            {
                case 0: goto DINVR;
                case 1: goto DSTINV;
            }

            DINVR:
            if (e0000Data.status > 0) goto S310;
            qcond = !qxmon(small, e0000Data.x, big);
            if (qcond)
            {
                throw new Exception(" SMALL, X, BIG not monotone in INVR");
            }

            xsave = e0000Data.x;
            //
            //     See that SMALL and BIG bound the zero and set QINCR
            //
            e0000Data.x = small;
            //
            //     GET-FUNCTION-VALUE
            //
            i99999 = 1;
            goto S300;
            S10:
            fsmall = e0000Data.fx;
            e0000Data.x = big;
            //
            //     GET-FUNCTION-VALUE
            //
            i99999 = 2;
            goto S300;
            S20:
            fbig = e0000Data.fx;
            qincr = fbig > fsmall;
            if (!qincr) goto S50;
            if (fsmall <= 0.0e0) goto S30;
            e0000Data.status = -1;
            e0000Data.qleft = e0000Data.qhi = true;
            return;
            S30:
            if (fbig >= 0.0e0) goto S40;
            e0000Data.status = -1;
            e0000Data.qleft = e0000Data.qhi = false;
            return;
            S40:
            goto S80;
            S50:
            if (fsmall >= 0.0e0) goto S60;
            e0000Data.status = -1;
            e0000Data.qleft = true;
            e0000Data.qhi = false;
            return;
            S60:
            if (fbig <= 0.0e0) goto S70;
            e0000Data.status = -1;
            e0000Data.qleft = false;
            e0000Data.qhi = true;
            return;
            S80:
            S70:
            e0000Data.x = xsave;
            step = Math.Max(absstp, relstp * Math.Abs(e0000Data.x));
            //
            //      YY = F(X) - Y
            //     GET-FUNCTION-VALUE
            //
            i99999 = 3;
            goto S300;
            S90:
            yy = e0000Data.fx;
            if (!(yy == 0.0e0)) goto S100;
            e0000Data.status = 0;
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
            e0000Data.x = xub;
            //
            //     GET-FUNCTION-VALUE
            //
            i99999 = 4;
            goto S300;
            S130:
            yy = e0000Data.fx;
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
            e0000Data.status = -1;
            e0000Data.qleft = false;
            e0000Data.qhi = !qincr;
            e0000Data.x = big;
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
            e0000Data.x = xlb;
            //
            //     GET-FUNCTION-VALUE
            //
            i99999 = 5;
            goto S300;
            S200:
            yy = e0000Data.fx;
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
            e0000Data.status = -1;
            e0000Data.qleft = true;
            e0000Data.qhi = qincr;
            e0000Data.x = small;
            return;
            S240:
            S230:
            e0001Data.zxlo = xlb;
            e0001Data.zxhi = xub;
            e0001Data.zabstl = abstol;
            e0001Data.zreltl = reltol;
            dstzr(ref e0001Data);
            //
            //  IF WE REACH HERE, XLB AND XUB BOUND THE ZERO OF F.
            //
            e0000Data.status = 0;
            goto S260;
            S250:
            if (!(e0000Data.status == 1)) goto S290;
            S260:
            dzror(ref e0000Data, ref e0001Data);
            if (!(e0000Data.status == 1)) goto S280;
            //
            //     GET-FUNCTION-VALUE
            //
            i99999 = 6;
            goto S300;
            S280:
            S270:
            goto S250;
            S290:
            e0000Data.x = xlo;
            e0000Data.status = 0;
            return;
            DSTINV:
            small = e0000Data.zsmall;
            big = e0000Data.zbig;
            absstp = e0000Data.zabsst;
            relstp = e0000Data.zrelst;
            stpmul = e0000Data.zstpmu;
            abstol = e0000Data.zabsto;
            reltol = e0000Data.zrelto;
            return;
            S300:
            //
            //     TO GET-FUNCTION-VALUE
            //
            e0000Data.status = 1;
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

        static void E0001(int IENTRY, ref E0001Data edata )

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

            double a = double.NaN;
            double abstol = double.NaN,
            b = double.NaN,c = double.NaN,d = double.NaN,fa = double.NaN,fb = double.NaN,fc = double.NaN,fd = double.NaN,fda = double.NaN;
            double fdb = double.NaN,
            m = double.NaN,mb = 0,p = double.NaN,q = double.NaN,reltol = double.NaN,tol = double.NaN,w = double.NaN,xxhi = double.NaN,xxlo = double.NaN;
            int ext = 0, i99999 = 0;
            bool first = false, qrzero = false;

            switch (IENTRY)
            {
                case 0:
                    goto DZROR;
                case 1:
                    goto DSTZR;
            }

            double ftol(double zx)
            {
                return (0.5e0 * Math.Max(abstol,reltol*Math.Abs((zx))));
            }

            DZROR:
            if (edata.status > 0) goto S280;
            edata.xlo = xxlo;
            edata.xhi = xxhi;
            b = edata.x = edata.xlo;
            //
            //     GET-FUNCTION-VALUE
            //
            i99999 = 1;
            goto S270;
            S10:
            fb = edata.fx;
            edata.xlo = edata.xhi;
            a = edata.x = edata.xlo;
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
            if (!(edata.fx < 0.0e0)) goto S30;
            edata.status = -1;
            edata.qleft = edata.fx < fb;
            edata.qhi = false;
            return;
            S40:
            S30:
            if (!(fb > 0.0e0)) goto S60;
            if (!(edata.fx > 0.0e0)) goto S50;
            edata.status = -1;
            edata.qleft = edata.fx > fb;
            edata.qhi = true;
            return;
            S60:
            S50:
            fa = edata.fx;
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
            edata.xlo = c;
            b = edata.xlo;
            fb = fc;
            c = a;
            fc = fa;
            S100:
            tol = ftol(edata.xlo);
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
            edata.xlo = b;
            edata.x = edata.xlo;
            //
            //  GET-FUNCTION-VALUE
            //
            i99999 = 3;
            goto S270;
            S200:
            fb = edata.fx;
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
            edata.xhi = c;
            qrzero = (fc >= 0.0e0 && fb <= 0.0e0) || (fc < 0.0e0 && fb >= 0.0e0);
            if (!qrzero) goto S250;
            edata.status = 0;
            goto S260;
            S250:
            edata.status = -1;
            S260:
            return;
            DSTZR:
            xxlo = edata.zxlo;
            xxhi = edata.zxhi;
            abstol = edata.zabstl;
            reltol = edata.zreltl;
            return;
            S270:
            //
            //     TO GET-FUNCTION-VALUE
            //
            edata.status = 1;
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