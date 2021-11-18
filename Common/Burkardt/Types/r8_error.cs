﻿using System;
using Burkardt.Probability;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public class r8ErrorData
    {
        public int nterf;
        public double sqeps;
        public const double sqrtpi = 1.77245385090551602729816748334115;
        public double xbig;

    }
        
    public static double r8_error( ref r8ErrorData data, ref r8ErrorcData cdata, double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_ERROR evaluates the error function of an R8 argument.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 September 2011
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Wayne Fullerton.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Wayne Fullerton,
        //    Portable Special Function Routines,
        //    in Portability of Numerical Software,
        //    edited by Wayne Cowell,
        //    Lecture Notes in Computer Science, Volume 57,
        //    Springer 1977,
        //    ISBN: 978-3-540-08446-4,
        //    LC: QA297.W65.
        //
        //  Parameters:
        //
        //    Input, double X, the argument.
        //
        //    Output, double R8_ERROR, the error function of X.
        //
    {
        double[] erfcs =
            {
                -0.49046121234691808039984544033376E-01,
                -0.14226120510371364237824741899631,
                +0.10035582187599795575754676712933E-01,
                -0.57687646997674847650827025509167E-03,
                +0.27419931252196061034422160791471E-04,
                -0.11043175507344507604135381295905E-05,
                +0.38488755420345036949961311498174E-07,
                -0.11808582533875466969631751801581E-08,
                +0.32334215826050909646402930953354E-10,
                -0.79910159470045487581607374708595E-12,
                +0.17990725113961455611967245486634E-13,
                -0.37186354878186926382316828209493E-15,
                +0.71035990037142529711689908394666E-17,
                -0.12612455119155225832495424853333E-18,
                +0.20916406941769294369170500266666E-20,
                -0.32539731029314072982364160000000E-22,
                +0.47668672097976748332373333333333E-24,
                -0.65980120782851343155199999999999E-26,
                +0.86550114699637626197333333333333E-28,
                -0.10788925177498064213333333333333E-29,
                +0.12811883993017002666666666666666E-31
            }
            ;
        double value = 0;

        switch (data.nterf)
        {
            case 0:
                data.nterf = inits(erfcs, 21, 0.1 * r8_mach(3));
                data.xbig = Math.Sqrt(-Math.Log(r8ErrorData.sqrtpi * r8_mach(3)));
                data.sqeps = Math.Sqrt(2.0 * r8_mach(3));
                break;
        }

        double y = Math.Abs(x);

        if (y <= data.sqeps)
        {
            value = 2.0 * x / r8ErrorData.sqrtpi;
        }
        else
        {
            switch (y)
            {
                case <= 1.0:
                    value = x * (1.0 + csevl(2.0 * x * x - 1.0, erfcs, data.nterf));
                    break;
                default:
                {
                    if (y <= data.xbig)
                    {
                        value = 1.0 - r8_errorc(ref cdata, y);
                        switch (x)
                        {
                            case < 0.0:
                                value = -value;
                                break;
                        }
                    }
                    else
                    {
                        value = x switch
                        {
                            < 0.0 => -value,
                            _ => 1.0
                        };
                    }

                    break;
                }
            }
        }

        return value;
    }

    public class r8ErrorcData
    {
        public int nterc2;
        public int nterf;
        public int nterfc;
        public double sqeps;
        public const double sqrtpi = 1.77245385090551602729816748334115;
        public double xmax;
        public double xsml;

    }

    public static double r8_errorc( ref r8ErrorcData data, double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_ERRORC evaluates the co-error function of an R8 argument.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 September 2011
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Wayne Fullerton.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Wayne Fullerton,
        //    Portable Special Function Routines,
        //    in Portability of Numerical Software,
        //    edited by Wayne Cowell,
        //    Lecture Notes in Computer Science, Volume 57,
        //    Springer 1977,
        //    ISBN: 978-3-540-08446-4,
        //    LC: QA297.W65.
        //
        //  Parameters:
        //
        //    Input, double X, the argument.
        //
        //    Output, double R8_ERRORC, the co-error function of X.
        //
    {
        double[] erc2cs =
            {
                -0.6960134660230950112739150826197E-01,
                -0.4110133936262089348982212084666E-01,
                +0.3914495866689626881561143705244E-02,
                -0.4906395650548979161280935450774E-03,
                +0.7157479001377036380760894141825E-04,
                -0.1153071634131232833808232847912E-04,
                +0.1994670590201997635052314867709E-05,
                -0.3642666471599222873936118430711E-06,
                +0.6944372610005012589931277214633E-07,
                -0.1371220902104366019534605141210E-07,
                +0.2788389661007137131963860348087E-08,
                -0.5814164724331161551864791050316E-09,
                +0.1238920491752753181180168817950E-09,
                -0.2690639145306743432390424937889E-10,
                +0.5942614350847910982444709683840E-11,
                -0.1332386735758119579287754420570E-11,
                +0.3028046806177132017173697243304E-12,
                -0.6966648814941032588795867588954E-13,
                +0.1620854541053922969812893227628E-13,
                -0.3809934465250491999876913057729E-14,
                +0.9040487815978831149368971012975E-15,
                -0.2164006195089607347809812047003E-15,
                +0.5222102233995854984607980244172E-16,
                -0.1269729602364555336372415527780E-16,
                +0.3109145504276197583836227412951E-17,
                -0.7663762920320385524009566714811E-18,
                +0.1900819251362745202536929733290E-18,
                -0.4742207279069039545225655999965E-19,
                +0.1189649200076528382880683078451E-19,
                -0.3000035590325780256845271313066E-20,
                +0.7602993453043246173019385277098E-21,
                -0.1935909447606872881569811049130E-21,
                +0.4951399124773337881000042386773E-22,
                -0.1271807481336371879608621989888E-22,
                +0.3280049600469513043315841652053E-23,
                -0.8492320176822896568924792422399E-24,
                +0.2206917892807560223519879987199E-24,
                -0.5755617245696528498312819507199E-25,
                +0.1506191533639234250354144051199E-25,
                -0.3954502959018796953104285695999E-26,
                +0.1041529704151500979984645051733E-26,
                -0.2751487795278765079450178901333E-27,
                +0.7290058205497557408997703680000E-28,
                -0.1936939645915947804077501098666E-28,
                +0.5160357112051487298370054826666E-29,
                -0.1378419322193094099389644800000E-29,
                +0.3691326793107069042251093333333E-30,
                -0.9909389590624365420653226666666E-31,
                +0.2666491705195388413323946666666E-31
            }
            ;
        double[] erfccs =
            {
                +0.715179310202924774503697709496E-01,
                -0.265324343376067157558893386681E-01,
                +0.171115397792085588332699194606E-02,
                -0.163751663458517884163746404749E-03,
                +0.198712935005520364995974806758E-04,
                -0.284371241276655508750175183152E-05,
                +0.460616130896313036969379968464E-06,
                -0.822775302587920842057766536366E-07,
                +0.159214187277090112989358340826E-07,
                -0.329507136225284321486631665072E-08,
                +0.722343976040055546581261153890E-09,
                -0.166485581339872959344695966886E-09,
                +0.401039258823766482077671768814E-10,
                -0.100481621442573113272170176283E-10,
                +0.260827591330033380859341009439E-11,
                -0.699111056040402486557697812476E-12,
                +0.192949233326170708624205749803E-12,
                -0.547013118875433106490125085271E-13,
                +0.158966330976269744839084032762E-13,
                -0.472689398019755483920369584290E-14,
                +0.143587337678498478672873997840E-14,
                -0.444951056181735839417250062829E-15,
                +0.140481088476823343737305537466E-15,
                -0.451381838776421089625963281623E-16,
                +0.147452154104513307787018713262E-16,
                -0.489262140694577615436841552532E-17,
                +0.164761214141064673895301522827E-17,
                -0.562681717632940809299928521323E-18,
                +0.194744338223207851429197867821E-18,
                -0.682630564294842072956664144723E-19,
                +0.242198888729864924018301125438E-19,
                -0.869341413350307042563800861857E-20,
                +0.315518034622808557122363401262E-20,
                -0.115737232404960874261239486742E-20,
                +0.428894716160565394623737097442E-21,
                -0.160503074205761685005737770964E-21,
                +0.606329875745380264495069923027E-22,
                -0.231140425169795849098840801367E-22,
                +0.888877854066188552554702955697E-23,
                -0.344726057665137652230718495566E-23,
                +0.134786546020696506827582774181E-23,
                -0.531179407112502173645873201807E-24,
                +0.210934105861978316828954734537E-24,
                -0.843836558792378911598133256738E-25,
                +0.339998252494520890627359576337E-25,
                -0.137945238807324209002238377110E-25,
                +0.563449031183325261513392634811E-26,
                -0.231649043447706544823427752700E-26,
                +0.958446284460181015263158381226E-27,
                -0.399072288033010972624224850193E-27,
                +0.167212922594447736017228709669E-27,
                -0.704599152276601385638803782587E-28,
                +0.297976840286420635412357989444E-28,
                -0.126252246646061929722422632994E-28,
                +0.539543870454248793985299653154E-29,
                -0.238099288253145918675346190062E-29,
                +0.109905283010276157359726683750E-29,
                -0.486771374164496572732518677435E-30,
                +0.152587726411035756763200828211E-30
            }
            ;
        double[] erfcs =
            {
                -0.49046121234691808039984544033376E-01,
                -0.14226120510371364237824741899631,
                +0.10035582187599795575754676712933E-01,
                -0.57687646997674847650827025509167E-03,
                +0.27419931252196061034422160791471E-04,
                -0.11043175507344507604135381295905E-05,
                +0.38488755420345036949961311498174E-07,
                -0.11808582533875466969631751801581E-08,
                +0.32334215826050909646402930953354E-10,
                -0.79910159470045487581607374708595E-12,
                +0.17990725113961455611967245486634E-13,
                -0.37186354878186926382316828209493E-15,
                +0.71035990037142529711689908394666E-17,
                -0.12612455119155225832495424853333E-18,
                +0.20916406941769294369170500266666E-20,
                -0.32539731029314072982364160000000E-22,
                +0.47668672097976748332373333333333E-24,
                -0.65980120782851343155199999999999E-26,
                +0.86550114699637626197333333333333E-28,
                -0.10788925177498064213333333333333E-29,
                +0.12811883993017002666666666666666E-31
            }
            ;
        double value = 0;

        switch (data.nterf)
        {
            case 0:
                double eta = 0.1 * r8_mach(3);
                data.nterf = inits(erfcs, 21, eta);
                data.nterfc = inits(erfccs, 59, eta);
                data.nterc2 = inits(erc2cs, 49, eta);

                data.xsml = -Math.Sqrt(-Math.Log(r8ErrorcData.sqrtpi * r8_mach(3)));
                data.xmax = Math.Sqrt(-Math.Log(r8ErrorcData.sqrtpi * r8_mach(1)));
                data.xmax = data.xmax - 0.5 * Math.Log(data.xmax) / data.xmax - 0.01;
                data.sqeps = Math.Sqrt(2.0 * r8_mach(3));
                break;
        }

        if (x <= data.xsml)
        {
            value = 2.0;
            return value;
        }

        if (data.xmax < x)
        {
            Console.WriteLine("");
            Console.WriteLine("R8_ERRORC - Warning!");
            Console.WriteLine("  X so big that ERFC underflows.");
            value = 0.0;
            return value;
        }

        double y = Math.Abs(x);

        if (y < data.sqeps)
        {
            value = 1.0 - 2.0 * x / r8ErrorcData.sqrtpi;
            return value;
        }

        switch (y)
        {
            case <= 1.0:
                value = 1.0 - x * (1.0
                                   + csevl(2.0 * x * x - 1.0, erfcs, data.nterf));
                return value;
        }

        y *= y;

        value = y switch
        {
            <= 4.0 => Math.Exp(-y) / Math.Abs(x) * (0.5 + csevl((8.0 / y - 5.0) / 3.0, erc2cs, data.nterc2)),
            _ => Math.Exp(-y) / Math.Abs(x) * (0.5 + csevl(8.0 / y - 1.0, erfccs, data.nterfc))
        };

        switch (x)
        {
            case < 0.0:
                value = 2.0 - value;
                break;
        }

        return value;
    }

    public static double r8_error_f(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_ERROR_F evaluates the error function ERF.
        //
        //  Discussion:
        //
        //    Since some compilers already supply a routine named ERF which evaluates
        //    the error function, this routine has been given a distinct, if
        //    somewhat unnatural, name.
        //
        //    The function is defined by:
        //
        //      ERF(X) = ( 2 / sqrt ( PI ) ) * Integral ( 0 <= T <= X ) EXP ( - T^2 ) dT.
        //
        //    Properties of the function include:
        //
        //      Limit ( X -> -Infinity ) ERF(X) =          -1.0;
        //                               ERF(0) =           0.0;
        //                               ERF(0.476936...) = 0.5;
        //      Limit ( X -> +Infinity ) ERF(X) =          +1.0.
        //
        //      0.5 * ( ERF(X/sqrt(2)) + 1 ) = Normal_01_CDF(X)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 November 2006
        //
        //  Author:
        //
        //    Original FORTRAN77 versino by William Cody.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    William Cody,
        //    "Rational Chebyshev Approximations for the Error Function",
        //    Mathematics of Computation,
        //    1969, pages 631-638.
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the error function.
        //
        //    Output, double R8_ERROR_F, the value of the error function.
        //
    {
        double[] a =
            {
                3.16112374387056560,
                1.13864154151050156E+02,
                3.77485237685302021E+02,
                3.20937758913846947E+03,
                1.85777706184603153E-01
            }
            ;
        double[] b =
            {
                2.36012909523441209E+01,
                2.44024637934444173E+02,
                1.28261652607737228E+03,
                2.84423683343917062E+03
            }
            ;
        double[] c =
            {
                5.64188496988670089E-01,
                8.88314979438837594,
                6.61191906371416295E+01,
                2.98635138197400131E+02,
                8.81952221241769090E+02,
                1.71204761263407058E+03,
                2.05107837782607147E+03,
                1.23033935479799725E+03,
                2.15311535474403846E-08
            }
            ;
        double[] d =
            {
                1.57449261107098347E+01,
                1.17693950891312499E+02,
                5.37181101862009858E+02,
                1.62138957456669019E+03,
                3.29079923573345963E+03,
                4.36261909014324716E+03,
                3.43936767414372164E+03,
                1.23033935480374942E+03
            }
            ;
        double erfxd;
        int i;
        double[] p =
            {
                3.05326634961232344E-01,
                3.60344899949804439E-01,
                1.25781726111229246E-01,
                1.60837851487422766E-02,
                6.58749161529837803E-04,
                1.63153871373020978E-02
            }
            ;
        double[] q =
            {
                2.56852019228982242,
                1.87295284992346047,
                5.27905102951428412E-01,
                6.05183413124413191E-02,
                2.33520497626869185E-03
            }
            ;
        const double sqrpi = 0.56418958354775628695;
        const double thresh = 0.46875;
        const double xbig = 26.543;
        double xden;
        double xnum;
        const double xsmall = 1.11E-16;
        double xsq;

        double xabs = Math.Abs(x);
        //
        //  Evaluate ERF(X) for |X| <= 0.46875.
        //
        if (xabs <= thresh)
        {
            if (xsmall < xabs)
            {
                xsq = xabs * xabs;
            }
            else
            {
                xsq = 0.0;
            }

            xnum = a[4] * xsq;
            xden = xsq;

            for (i = 0; i < 3; i++)
            {
                xnum = (xnum + a[i]) * xsq;
                xden = (xden + b[i]) * xsq;
            }

            erfxd = x * (xnum + a[3]) / (xden + b[3]);
        }
        else
        {
            double del;
            switch (xabs)
            {
                //
                //  Evaluate ERFC(X) for 0.46875 <= |X| <= 4.0.
                //
                case <= 4.0:
                {
                    xnum = c[8] * xabs;
                    xden = xabs;
                    for (i = 0; i < 7; i++)
                    {
                        xnum = (xnum + c[i]) * xabs;
                        xden = (xden + d[i]) * xabs;
                    }

                    erfxd = (xnum + c[7]) / (xden + d[7]);
                    xsq = (int) (xabs * 16.0) / 16.0;
                    del = (xabs - xsq) * (xabs + xsq);
                    erfxd = Math.Exp(-xsq * xsq) * Math.Exp(-del) * erfxd;

                    erfxd = x switch
                    {
                        < 0.0 => -erfxd,
                        _ => 0.5 - erfxd + 0.5
                    };

                    break;
                }
                //
                default:
                {
                    if (xbig <= xabs)
                    {
                        erfxd = x switch
                        {
                            > 0.0 => 1.0,
                            _ => -1.0
                        };
                    }
                    else
                    {
                        xsq = 1.0 / (xabs * xabs);

                        xnum = p[5] * xsq;
                        xden = xsq;

                        for (i = 0; i < 4; i++)
                        {
                            xnum = (xnum + p[i]) * xsq;
                            xden = (xden + q[i]) * xsq;
                        }

                        erfxd = xsq * (xnum + p[4]) / (xden + q[4]);
                        erfxd = (sqrpi - erfxd) / xabs;
                        xsq = (int) (xabs * 16.0) / 16.0;
                        del = (xabs - xsq) * (xabs + xsq);
                        erfxd = Math.Exp(-xsq * xsq) * Math.Exp(-del) * erfxd;

                        erfxd = x switch
                        {
                            < 0.0 => -erfxd,
                            _ => 0.5 - erfxd + 0.5
                        };
                    }

                    break;
                }
            }
        }

        return erfxd;
    }

    public static double r8_error_f_inverse(double y)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_ERROR_F_INVERSE inverts the error function ERF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 August 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double Y, the value of the error function.
        //
        //    Output, double R8_ERROR_F_INVERSE, the value X such that
        //    ERF(X) = Y.
        //
    {
        double z = (y + 1.0) / 2.0;

        double x = Normal.normal_01_cdf_inv(z);

        double value = x / Math.Sqrt(2.0);

        return value;
    }
        
}