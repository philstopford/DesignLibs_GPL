using System;

namespace Burkardt.FullertonFnLib
{
    public static partial class FullertonLib
    {
        public static double r8_tan(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_TAN evaluates the tangent of an R8 argument.
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
            //    Output, double R8_TAN, the tangent of X.
            //
        {
            double ainty;
            double ainty2;
            int ifn;
            int nterms = 0;
            const double pi2rec = 0.011619772367581343075535053490057;
            double prodbg;
            double sqeps = 0.0;
            double[] tancs = {
                +0.22627932763129357846578636531752,
                +0.43017913146548961775583410748067E-01,
                +0.68544610682565088756929473623461E-03,
                +0.11045326947597098383578849369696E-04,
                +0.17817477903926312943238512588940E-06,
                +0.28744968582365265947529646832471E-08,
                +0.46374854195902995494137478234363E-10,
                +0.74817609041556138502341633308215E-12,
                +0.12070497002957544801644516947824E-13,
                +0.19473610812823019305513858584533E-15,
                +0.31417224874732446504614586026666E-17,
                +0.50686132555800153941904891733333E-19,
                +0.81773105159836540043979946666666E-21,
                +0.13192643412147384408951466666666E-22,
                +0.21283995497042377309866666666666E-24,
                +0.34337960192345945292800000000000E-26,
                +0.55398222121173811200000000000000E-28,
                +0.89375227794352810666666666666666E-30,
                +0.14419111371369130666666666666666E-31
            }
            ;
            double value;
            double xmax = 0.0;
            double xsml = 0.0;
            double y;
            double yrem;

            if (nterms == 0)
            {
                nterms = r8_inits(tancs, 19, 0.1 * r8_mach(3));
                xmax = 1.0 / r8_mach(4);
                xsml = Math.Sqrt(3.0 * r8_mach(3));
                sqeps = Math.Sqrt(r8_mach(4));
            }

            y = Math.Abs(x);

            if (xmax < y)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_TAN - Warning");
                Console.WriteLine("  No precision because |X| is big.");
                value = 0.0;
                return value;
            }

            //
            //  Carefully compute y * (2/pi) = (aint(y) + rem(y)) * (.625 + pi2rec)
            //  = aint(.625*y) + rem(.625*y) + y*pi2rec  =  aint(.625*y) + z
            //  = aint(.625*y) + aint(z) + rem(z)
            //
            ainty = r8_aint(y);
            yrem = y - ainty;
            prodbg = 0.625 * ainty;
            ainty = r8_aint(prodbg);
            y = (prodbg - ainty) + 0.625 * yrem + pi2rec * y;
            ainty2 = r8_aint(y);
            ainty = ainty + ainty2;
            y = y - ainty2;

            ifn = (int) (ainty % 2.0);

            if (ifn == 1)
            {
                y = 1.0 - y;
            }

            if (1.0 - y < Math.Abs(x) * sqeps)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_TAN - Warning!");
                Console.WriteLine("  Answer < half precision.");
                Console.WriteLine("  |X| big or X near pi/2 or 3*pi/2.");
            }

            if (y == 1.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_TAN - Fatal error!");
                Console.WriteLine("  X is pi/2 or 3*pi/2.");
                return (1);
            }

            if (y <= 0.25)
            {
                value = y;
                if (xsml < y)
                {
                    value = y * (1.5 + r8_csevl(32.0 * y * y - 1.0, tancs, nterms));
                }
            }
            else if (y <= 0.5)
            {
                value = 0.5 * y * (1.5 + r8_csevl(
                    8.0 * y * y - 1.0, tancs, nterms));
                value = 2.0 * value / (1.0 - value * value);
            }
            else
            {
                value = 0.25 * y * (1.5 + r8_csevl(
                    2.0 * y * y - 1.0, tancs, nterms));
                value = 2.0 * value / (1.0 - value * value);
                value = 2.0 * value / (1.0 - value * value);
            }

            if (x < 0.0)
            {
                value = -Math.Abs(value);
            }
            else if (0.0 < x)
            {
                value = +Math.Abs(value);
            }

            if (ifn == 1)
            {
                value = -value;
            }

            return value;
        }

        public static double r8_tanh(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_TANH evaluates the hyperbolic tangent of an R8 argument.
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
            //    Output, double R8_TANH, the hyperbolic tangent of X.
            //
        {
            int nterms = 0;
            double sqeps = 0.0;
            double[] tanhcs = {
                -0.25828756643634710438338151450605,
                -0.11836106330053496535383671940204,
                +0.98694426480063988762827307999681E-02,
                -0.83579866234458257836163690398638E-03,
                +0.70904321198943582626778034363413E-04,
                -0.60164243181207040390743479001010E-05,
                +0.51052419080064402965136297723411E-06,
                -0.43320729077584087216545467387192E-07,
                +0.36759990553445306144930076233714E-08,
                -0.31192849612492011117215651480953E-09,
                +0.26468828199718962579377758445381E-10,
                -0.22460239307504140621870997006196E-11,
                +0.19058733768288196054319468396139E-12,
                -0.16172371446432292391330769279701E-13,
                +0.13723136142294289632897761289386E-14,
                -0.11644826870554194634439647293781E-15,
                +0.98812684971669738285540514338133E-17,
                -0.83847933677744865122269229055999E-18,
                +0.71149528869124351310723506176000E-19,
                -0.60374242229442045413288837119999E-20,
                +0.51230825877768084883404663466666E-21,
                -0.43472140157782110106047829333333E-22,
                +0.36888473639031328479423146666666E-23,
                -0.31301874774939399883325439999999E-24,
                +0.26561342006551994468488533333333E-25,
                -0.22538742304145029883494399999999E-26,
                +0.19125347827973995102208000000000E-27,
                -0.16228897096543663117653333333333E-28,
                +0.13771101229854738786986666666666E-29,
                -0.11685527840188950118399999999999E-30,
                +0.99158055384640389120000000000000E-32
            }
            ;
            double value;
            double xmax = 0.0;
            double y;
            double yrec;

            if (nterms == 0)
            {
                nterms = r8_inits(tanhcs, 31, 0.1 * r8_mach(3));
                sqeps = Math.Sqrt(3.0 * r8_mach(3));
                xmax = -0.5 * Math.Log(r8_mach(3));
            }

            y = Math.Abs(x);

            if (y <= sqeps)
            {
                value = x;
            }
            else if (y <= 1.0)
            {
                value = x * (1.0 + r8_csevl(2.0 * x * x - 1.0, tanhcs, nterms));
            }
            else if (y <= xmax)
            {
                y = Math.Exp(y);
                yrec = 1.0 / y;
                value = (y - yrec) / (y + yrec);

                if (x < 0.0)
                {
                    value = -value;
                }
            }
            else
            {
                if (x < 0.0)
                {
                    value = -1.0;
                }
                else
                {
                    value = +1.0;
                }
            }

            return value;
        }
    }
}