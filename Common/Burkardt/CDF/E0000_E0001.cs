using System;

namespace Burkardt.CDFLib;

public class E0000_E0001_Data
{
    public int status;
    public double x;
    public double fx;
    public bool qleft;
    public bool qhi;

    public double zabsst;
    public double zabsto;
    public double zbig;
    public double zrelst;
    public double zrelto;
    public double zsmall;
    public double zstpmu;

    public double xlo;
    public double xhi;
    public double zabstl;
    public double zreltl;
        
    public class E0000Variables
    {
        public double absstp;
        public double abstol;
        public double big = double.MaxValue;
        public double fbig,
            fsmall,
            relstp,
            reltol,
            small = double.MinValue,
            step,
            stpmul,
            xlb,
            xsave,
            xub,
            yy;
        public int i99999;
        public bool qbdd, qcond, qdum1 = false, qdum2 = false, qincr, qlim;
        public bool qup;
        
    }
        
    public class E0001Variables
    {
        public double a;
        public double abstol, b, c, d, fa, fb, fc, fd, fda;
        public double fdb, m, mb, p, q, reltol, tol, w, xxhi, xxlo;
        public int ext, i99999;
        public bool first, qrzero;        
    }

    public E0000Variables e0000vars;
    public E0001Variables e0001vars;
        
    public E0000_E0001_Data()
    {
        e0000vars = new E0000Variables();
        e0001vars = new E0001Variables();
    }
        
}


    
public static class E0000E0001
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
        data.zsmall = zsmall_;
        data.zbig = zbig_;
        data.zabsst = zabsst_;
        data.zrelst = zrelst_;
        data.zstpmu = zstpmu_;
        data.zabsto = zabsto_;
        data.zrelto = zrelto_;
            
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
        data.xlo = zxlo;
        data.xhi = zxhi;
        data.zabstl = zabstl;
        data.zreltl = zreltl;
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

    public static bool qxmon(double zx, double zy, double zz)
    {
        return zx <= zy && zy <= zz;
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


        switch (IENTRY)
        {
            case 0: goto DINVR;
            case 1: goto DSTINV;
        }

        DINVR:
        switch (data.status)
        {
            case > 0:
                goto S310;
        }

        data.e0000vars.qcond = !qxmon(data.e0000vars.small, data.x, data.e0000vars.big);
        switch (data.e0000vars.qcond)
        {
            case true:
                throw new Exception(" SMALL, X, BIG not monotone in INVR");
        }

        data.e0000vars.xsave = data.x;
        //
        //     See that SMALL and BIG bound the zero and set QINCR
        //
        data.x = data.e0000vars.small;
        //
        //     GET-FUNCTION-VALUE
        //
        data.e0000vars.i99999 = 1;
        goto S300;
        S10:
        data.e0000vars.fsmall = data.fx;
        data.x = data.e0000vars.big;
        //
        //     GET-FUNCTION-VALUE
        //
        data.e0000vars.i99999 = 2;
        goto S300;
        S20:
        data.e0000vars.fbig = data.fx;
        data.e0000vars.qincr = data.e0000vars.fbig > data.e0000vars.fsmall;
        switch (data.e0000vars.qincr)
        {
            case false:
                goto S50;
        }

        switch (data.e0000vars.fsmall)
        {
            case <= 0.0e0:
                goto S30;
        }

        data.status = -1;
        data.qleft = data.qhi = true;
        return;
        S30:
        switch (data.e0000vars.fbig)
        {
            case >= 0.0e0:
                goto S40;
        }

        data.status = -1;
        data.qleft = data.qhi = false;
        return;
        S40:
        goto S80;
        S50:
        switch (data.e0000vars.fsmall)
        {
            case >= 0.0e0:
                goto S60;
        }

        data.status = -1;
        data.qleft = true;
        data.qhi = false;
        return;
        S60:
        switch (data.e0000vars.fbig)
        {
            case <= 0.0e0:
                goto S70;
        }

        data.status = -1;
        data.qleft = false;
        data.qhi = true;
        return;
        S80:
        S70:
        data.x = data.e0000vars.xsave;
        data.e0000vars.step = Math.Max(data.e0000vars.absstp, data.e0000vars.relstp * Math.Abs(data.x));
        //
        //      YY = F(X) - Y
        //     GET-FUNCTION-VALUE
        //
        data.e0000vars.i99999 = 3;
        goto S300;
        S90:
        data.e0000vars.yy = data.fx;
        if (data.e0000vars.yy != 0.0e0)
        {
            goto S100;
        }

        data.status = 0;
        //  qok = 1;
        return;
        S100:
        data.e0000vars.qup = data.e0000vars.qincr && data.e0000vars.yy < 0.0e0 || !data.e0000vars.qincr && data.e0000vars.yy > 0.0e0;
        switch (data.e0000vars.qup)
        {
            //
            //     HANDLE CASE IN WHICH WE MUST STEP HIGHER
            //
            case false:
                goto S170;
        }

        data.e0000vars.xlb = data.e0000vars.xsave;
        data.e0000vars.xub = Math.Min(data.e0000vars.xlb + data.e0000vars.step, data.e0000vars.big);
        goto S120;
        S110:
        switch (data.e0000vars.qcond)
        {
            case true:
                goto S150;
        }

        S120:
        //
        //      YY = F(XUB) - Y
        //
        data.x = data.e0000vars.xub;
        //
        //     GET-FUNCTION-VALUE
        //
        data.e0000vars.i99999 = 4;
        goto S300;
        S130:
        data.e0000vars.yy = data.fx;
        data.e0000vars.qbdd = data.e0000vars.qincr && data.e0000vars.yy >= 0.0e0 || !data.e0000vars.qincr && data.e0000vars.yy <= 0.0e0;
        data.e0000vars.qlim = data.e0000vars.xub >= data.e0000vars.big;
        data.e0000vars.qcond = data.e0000vars.qbdd || data.e0000vars.qlim;
        switch (data.e0000vars.qcond)
        {
            case true:
                goto S140;
        }

        data.e0000vars.step = data.e0000vars.stpmul * data.e0000vars.step;
        data.e0000vars.xlb = data.e0000vars.xub;
        data.e0000vars.xub = Math.Min(data.e0000vars.xlb + data.e0000vars.step, data.e0000vars.big);
        S140:
        goto S110;
        S150:
        switch ((data.e0000vars.qlim && !data.e0000vars.qbdd))
        {
            case false:
                goto S160;
        }

        data.status = -1;
        data.qleft = false;
        data.qhi = !data.e0000vars.qincr;
        data.x = data.e0000vars.big;
        return;
        S160:
        goto S240;
        S170:
        //
        //     HANDLE CASE IN WHICH WE MUST STEP LOWER
        //
        data.e0000vars.xub = data.e0000vars.xsave;
        data.e0000vars.xlb = Math.Max(data.e0000vars.xub - data.e0000vars.step, data.e0000vars.small);
        goto S190;
        S180:
        switch (data.e0000vars.qcond)
        {
            case true:
                goto S220;
        }

        S190:
        //
        //      YY = F(XLB) - Y
        //
        data.x = data.e0000vars.xlb;
        //
        //     GET-FUNCTION-VALUE
        //
        data.e0000vars.i99999 = 5;
        goto S300;
        S200:
        data.e0000vars.yy = data.fx;
        data.e0000vars.qbdd = data.e0000vars.qincr && data.e0000vars.yy <= 0.0e0 || !data.e0000vars.qincr && data.e0000vars.yy >= 0.0e0;
        data.e0000vars.qlim = data.e0000vars.xlb <= data.e0000vars.small;
        data.e0000vars.qcond = data.e0000vars.qbdd || data.e0000vars.qlim;
        switch (data.e0000vars.qcond)
        {
            case true:
                goto S210;
        }

        data.e0000vars.step = data.e0000vars.stpmul * data.e0000vars.step;
        data.e0000vars.xub = data.e0000vars.xlb;
        data.e0000vars.xlb = Math.Max(data.e0000vars.xub - data.e0000vars.step, data.e0000vars.small);
        S210:
        goto S180;
        S220:
        switch ((data.e0000vars.qlim && !data.e0000vars.qbdd))
        {
            case false:
                goto S230;
        }

        data.status = -1;
        data.qleft = true;
        data.qhi = data.e0000vars.qincr;
        data.x = data.e0000vars.small;
        return;
        S240:
        S230:
        dstzr(ref data, data.e0000vars.xlb, data.e0000vars.xub, data.e0000vars.abstol, data.e0000vars.reltol);
        //
        //  IF WE REACH HERE, XLB AND XUB BOUND THE ZERO OF F.
        //
        data.status = 0;
        goto S260;
        S250:
        if (data.status != 1)
        {
            goto S290;
        }

        S260:
        data.qleft = data.e0000vars.qdum1;
        data.qhi = data.e0000vars.qdum2;
        dzror(ref data);
        if (data.status != 1)
        {
            goto S280;
        }

        //
        //     GET-FUNCTION-VALUE
        //
        data.e0000vars.i99999 = 6;
        goto S300;
        S280:
        S270:
        goto S250;
        S290:
        data.x = data.xlo;
        data.status = 0;
        return;
        DSTINV:
        data.e0000vars.small = data.zsmall;
        data.e0000vars.big = data.zbig;
        data.e0000vars.absstp = data.zabsst;
        data.e0000vars.relstp = data.zrelst;
        data.e0000vars.stpmul = data.zstpmu;
        data.e0000vars.abstol = data.zabsto;
        data.e0000vars.reltol = data.zrelto;
        return;
        S300:
        //
        //     TO GET-FUNCTION-VALUE
        //
        data.status = 1;
        return;
        S310:
        switch (data.e0000vars.i99999)
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
        mag = sign switch
        {
            < 0 => -mag,
            _ => mag switch
            {
                < 0 => -mag,
                _ => mag
            }
        };

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

        double ftol(double zx, ref E0000_E0001_Data data)
        {
            return 0.5e0 * Math.Max(data.e0001vars.abstol, data.e0001vars.reltol * Math.Abs(zx));
        }

        switch (IENTRY)
        {
            case 0:
                goto DZROR;
            case 1:
                goto DSTZR;
        }

        DZROR:
        switch (data.status)
        {
            case > 0:
                goto S280;
        }

        data.xlo = data.e0001vars.xxlo;
        data.xhi = data.e0001vars.xxhi;
        data.e0001vars.b = data.x = data.xlo;
        //
        //     GET-FUNCTION-VALUE
        //
        data.e0001vars.i99999 = 1;
        goto S270;
        S10:
        data.e0001vars.fb = data.fx;
        data.xlo = data.xhi;
        data.e0001vars.a = data.x = data.xlo;
        //
        //     GET-FUNCTION-VALUE
        //
        data.e0001vars.i99999 = 2;
        goto S270;
        S20:
        switch ((data.e0001vars.fb < 0.0e0))
        {
            //
            //     Check that F(ZXLO) < 0 < F(ZXHI)  or
            //                F(ZXLO) > 0 > F(ZXHI)
            //
            case false:
                goto S40;
        }

        switch ((data.fx < 0.0e0))
        {
            case false:
                goto S30;
        }

        data.status = -1;
        data.qleft = data.fx < data.e0001vars.fb;
        data.qhi = false;
        return;
        S40:
        S30:
        switch ((data.e0001vars.fb > 0.0e0))
        {
            case false:
                goto S60;
        }

        switch ((data.fx > 0.0e0))
        {
            case false:
                goto S50;
        }

        data.status = -1;
        data.qleft = data.fx > data.e0001vars.fb;
        data.qhi = true;
        return;
        S60:
        S50:
        data.e0001vars.fa = data.fx;
        data.e0001vars.first = true;
        S70:
        data.e0001vars.c = data.e0001vars.a;
        data.e0001vars.fc = data.e0001vars.fa;
        data.e0001vars.ext = 0;
        S80:
        switch ((Math.Abs(data.e0001vars.fc) < Math.Abs(data.e0001vars.fb)))
        {
            case false:
                goto S100;
        }

        switch ((data.e0001vars.c != data.e0001vars.a))
        {
            case false:
                goto S90;
        }

        data.e0001vars.d = data.e0001vars.a;
        data.e0001vars.fd = data.e0001vars.fa;
        S90:
        data.e0001vars.a = data.e0001vars.b;
        data.e0001vars.fa = data.e0001vars.fb;
        data.xlo = data.e0001vars.c;
        data.e0001vars.b = data.xlo;
        data.e0001vars.fb = data.e0001vars.fc;
        data.e0001vars.c = data.e0001vars.a;
        data.e0001vars.fc = data.e0001vars.fa;
        S100:
        data.e0001vars.tol = ftol(data.xlo, ref data);
        data.e0001vars.m = (data.e0001vars.c + data.e0001vars.b) * .5e0;
        data.e0001vars.mb = data.e0001vars.m - data.e0001vars.b;
        switch ((Math.Abs(data.e0001vars.mb) > data.e0001vars.tol))
        {
            case false:
                goto S240;
        }

        switch ((data.e0001vars.ext > 3))
        {
            case false:
                goto S110;
        }

        data.e0001vars.w = data.e0001vars.mb;
        goto S190;
        S110:
        data.e0001vars.tol = fifdsign(data.e0001vars.tol, data.e0001vars.mb);
        data.e0001vars.p = (data.e0001vars.b - data.e0001vars.a) * data.e0001vars.fb;
        switch (data.e0001vars.first)
        {
            case false:
                goto S120;
        }

        data.e0001vars.q = data.e0001vars.fa - data.e0001vars.fb;
        data.e0001vars.first = false;
        goto S130;
        S120:
        data.e0001vars.fdb = (data.e0001vars.fd - data.e0001vars.fb) / (data.e0001vars.d - data.e0001vars.b);
        data.e0001vars.fda = (data.e0001vars.fd - data.e0001vars.fa) / (data.e0001vars.d - data.e0001vars.a);
        data.e0001vars.p = data.e0001vars.fda * data.e0001vars.p;
        data.e0001vars.q = data.e0001vars.fdb * data.e0001vars.fa - data.e0001vars.fda * data.e0001vars.fb;
        S130:
        switch ((data.e0001vars.p < 0.0e0))
        {
            case false:
                goto S140;
        }

        data.e0001vars.p = -data.e0001vars.p;
        data.e0001vars.q = -data.e0001vars.q;
        S140:
        switch (data.e0001vars.ext)
        {
            case 3:
                data.e0001vars.p *= 2.0e0;
                break;
        }

        switch ((data.e0001vars.p * 1.0e0 == 0.0e0 || data.e0001vars.p <= data.e0001vars.q * data.e0001vars.tol))
        {
            case false:
                goto S150;
        }

        data.e0001vars.w = data.e0001vars.tol;
        goto S180;
        S150:
        switch ((data.e0001vars.p < data.e0001vars.mb * data.e0001vars.q))
        {
            case false:
                goto S160;
        }

        data.e0001vars.w = data.e0001vars.p / data.e0001vars.q;
        goto S170;
        S160:
        data.e0001vars.w = data.e0001vars.mb;
        S190:
        S180:
        S170:
        data.e0001vars.d = data.e0001vars.a;
        data.e0001vars.fd = data.e0001vars.fa;
        data.e0001vars.a = data.e0001vars.b;
        data.e0001vars.fa = data.e0001vars.fb;
        data.e0001vars.b += data.e0001vars.w;
        data.xlo = data.e0001vars.b;
        data.x = data.xlo;
        //
        //  GET-FUNCTION-VALUE
        //
        data.e0001vars.i99999 = 3;
        goto S270;
        S200:
        data.e0001vars.fb = data.fx;
        switch ((data.e0001vars.fc * data.e0001vars.fb >= 0.0e0))
        {
            case false:
                goto S210;
        }

        goto S70;
        S210:
        switch ((data.e0001vars.w == data.e0001vars.mb))
        {
            case false:
                goto S220;
        }

        data.e0001vars.ext = 0;
        goto S230;
        S220:
        data.e0001vars.ext += 1;
        S230:
        goto S80;
        S240:
        data.xhi = data.e0001vars.c;
        data.e0001vars.qrzero = data.e0001vars.fc >= 0.0e0 && data.e0001vars.fb <= 0.0e0 || data.e0001vars.fc < 0.0e0 && data.e0001vars.fb >= 0.0e0;
        switch (data.e0001vars.qrzero)
        {
            case false:
                goto S250;
        }

        data.status = 0;
        goto S260;
        S250:
        data.status = -1;
        S260:
        return;
        DSTZR:
        data.e0001vars.xxlo = data.xlo;
        data.e0001vars.xxhi = data.xhi;
        data.e0001vars.abstol = data.zabstl;
        data.e0001vars.reltol = data.zreltl;
        return;
        S270:
        //
        //     TO GET-FUNCTION-VALUE
        //
        data.status = 1;
        return;
        S280:
        switch (data.e0001vars.i99999)
        {
            case 1: goto S10;
            case 2: goto S20;
            case 3: goto S200;
            default: break;
        }
    }

}