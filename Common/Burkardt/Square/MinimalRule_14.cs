﻿using Burkardt.Types;

namespace Burkardt.Square
{
    public static partial class MinimalRule
    {
        public static double[] smr14()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SMR14 returns the SMR rule of degree 14.
            //
            //  Discussion:
            // 
            //    DEGREE: 14
            //    POINTS CARDINALITY: 40
            //    NORM INF MOMS. RESIDUAL: 9.57567e-16
            //    SUM NEGATIVE WEIGHTS: 0.00000e+00,
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    21 February 2018
            //
            //  Author:
            //
            //    Original MATLAB version by Mattia Festa, Alvise Sommariva,
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Mattia Festa, Alvise Sommariva,
            //    Computing almost minimal formulas on the square,
            //    Journal of Computational and Applied Mathematics,
            //    Volume 17, Number 236, November 2012, pages 4296-4302.
            //
            //  Parameters:
            //
            //    Output, double *SMR14[3*40], the requested rule.
            //
        {
            int degree = 14;
            int order;
            double[] xw =
            {
                9.784247248091626e-01, 9.849379037345679e-01, 7.042609290974809e-03,
                9.696033367128390e-01, -3.331633352358324e-01, 5.127000401619231e-02,
                9.635907900433699e-01, -8.235135827696192e-01, 2.835904116834926e-02,
                9.609680853689336e-01, -9.665848054209388e-01, 5.802552341290735e-03,
                8.759175244824323e-01, 8.440519016975554e-01, 6.406839799781101e-02,
                9.850622977861553e-01, 5.927006882272775e-01, 3.059827453183588e-02,
                9.132378164262995e-01, 1.686644107446046e-01, 1.015293304177829e-01,
                8.210255575251280e-01, -6.312138484114799e-01, 1.070566792196576e-01,
                7.253802697646771e-01, 5.271299424251080e-01, 1.491850201682654e-01,
                7.129972234617307e-01, -1.746725364206650e-01, 1.672674729613243e-01,
                5.551781984544302e-01, -8.617677569495097e-01, 1.035166248728787e-01,
                8.047627946516351e-01, -9.674432786379303e-01, 2.787008138269169e-02,
                6.424373183439024e-01, 9.660008737161317e-01, 4.265704272492368e-02,
                4.202620588036903e-01, 7.809596116887254e-01, 1.425075466703280e-01,
                4.547166615470812e-01, 1.927401717576898e-01, 2.105516162791665e-01,
                4.259990143599163e-01, -4.962208715768625e-01, 1.867185172918136e-01,
                9.710021027954509e-02, 9.536977082453567e-01, 6.972182442589203e-02,
                1.114008919498024e-01, 4.923616289896904e-01, 2.043948618313033e-01,
                1.209667130758268e-01, -1.503469100876254e-01, 2.335930539280487e-01,
                7.524097400070873e-02, -7.452253206810735e-01, 1.601269026668343e-01,
                -4.518656196438389e-01, 9.427795716582286e-01, 6.465270834286081e-02,
                -2.040271575380492e-01, 7.591224031418166e-01, 1.466929221365241e-01,
                -2.427862218874625e-01, 1.736708825032474e-01, 2.280281871073875e-01,
                -2.459782585741061e-01, -4.449955936922590e-01, 2.055272577034085e-01,
                2.124852450722984e-01, -9.810732833242932e-01, 4.316769231929027e-02,
                -7.568108239381779e-01, 7.901122315570046e-01, 9.576721669303179e-02,
                -5.448202693114659e-01, 5.034324548840123e-01, 1.739315171493671e-01,
                -5.908385723762482e-01, -1.268935867110808e-01, 1.941490024297708e-01,
                -2.895046302710989e-01, -9.116133911949098e-01, 9.663345199603172e-02,
                -6.417036777526361e-01, -9.886715561256082e-01, 2.266331268661072e-02,
                -7.921595150813847e-01, 9.955322518565568e-01, 1.393382316893550e-02,
                -9.527313016559716e-01, 5.688243284550245e-01, 5.657603379663307e-02,
                -8.300912488610223e-01, 2.262029630794567e-01, 1.327223448878936e-01,
                -5.758839055540124e-01, -6.756366015563920e-01, 1.369481741372994e-01,
                -7.954880187581953e-01, -8.506684869408652e-01, 6.033321912862107e-02,
                -9.562007514872694e-01, 9.139906777934969e-01, 2.419546581425405e-02,
                -9.839876523259704e-01, -5.381597268716472e-02, 4.116110545021220e-02,
                -8.552518042944189e-01, -4.002028380944249e-01, 1.164178174263354e-01,
                -9.694591723778574e-01, -6.823556992833991e-01, 3.472756253802061e-02,
                -9.468282975367817e-01, -9.506255113908446e-01, 1.793373090014833e-02
            };
            double[] xw_copy;

            order = square_minimal_rule_order(degree);
            xw_copy = typeMethods.r8mat_copy_new(3, order, xw);

            return xw_copy;
        }
    }
}