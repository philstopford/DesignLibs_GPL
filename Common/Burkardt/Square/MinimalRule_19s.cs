﻿using Burkardt.Types;

namespace Burkardt.Square
{
    public static partial class MinimalRule
    {
        public static double[] smr19s()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SMR19S returns the rotationally invariant SMR rule of degree 19.
            //
            //  Discussion:
            //
            //    DEGREE: 19
            //    ROTATIONALLY INVARIANT: (X,Y),(-Y,X),(-X,-Y),(Y,-X).
            //    POINTS CARDINALITY: 68
            //    NORM INF MOMS. RESIDUAL: 2.63678e-16
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
            //    Output, double *SMR19S[3*68], the requested rule.
            //
        {
            int degree = 19;
            int order;
            double[] xw =
            {
                6.772712303303220e-02, 1.091623307294788e-01, 5.494509440087280e-02,
                1.911320239476111e-01, 7.684428790753622e-01, 9.924039371581232e-02,
                1.293584417049455e-01, 9.843401047067886e-01, 2.298164089064930e-02,
                2.653188588974263e-01, 3.187377804093571e-01, 1.356277329376223e-01,
                4.805884103058745e-01, 5.883672104526556e-01, 1.092830463202312e-01,
                4.353820421144143e-01, 9.253041863610407e-01, 5.332719030254342e-02,
                5.327122733511525e-01, 5.307888414194004e-02, 1.342811802398570e-01,
                6.817257899596626e-01, 8.006754071171294e-01, 6.668259225153295e-02,
                6.897274550636672e-01, 9.904094449990978e-01, 1.399432148751916e-02,
                8.960307328244307e-01, 1.256777467976840e-01, 6.752580403681208e-02,
                7.189012994269330e-01, 3.612251785232622e-01, 1.014916530440626e-01,
                8.636764064784566e-01, 6.273567750914366e-01, 5.639122590523545e-02,
                8.689071260438135e-01, 9.214958682077944e-01, 2.840189695206082e-02,
                9.851138463638244e-01, 3.506711929883037e-01, 1.440758011223270e-02,
                9.604323781459383e-01, 4.868528161195464e-01, 1.812287875530189e-02,
                9.757200071559576e-01, 7.955215584636620e-01, 1.826958037918748e-02,
                9.734329880957604e-01, 9.798840213938864e-01, 5.026188268466660e-03,
                -1.091623307294788e-01, 6.772712303303220e-02, 5.494509440087280e-02,
                -7.684428790753622e-01, 1.911320239476111e-01, 9.924039371581232e-02,
                -9.843401047067886e-01, 1.293584417049455e-01, 2.298164089064930e-02,
                -3.187377804093571e-01, 2.653188588974263e-01, 1.356277329376223e-01,
                -5.883672104526556e-01, 4.805884103058745e-01, 1.092830463202312e-01,
                -9.253041863610407e-01, 4.353820421144143e-01, 5.332719030254342e-02,
                -5.307888414194004e-02, 5.327122733511525e-01, 1.342811802398570e-01,
                -8.006754071171294e-01, 6.817257899596626e-01, 6.668259225153295e-02,
                -9.904094449990978e-01, 6.897274550636672e-01, 1.399432148751916e-02,
                -1.256777467976840e-01, 8.960307328244307e-01, 6.752580403681208e-02,
                -3.612251785232622e-01, 7.189012994269330e-01, 1.014916530440626e-01,
                -6.273567750914366e-01, 8.636764064784566e-01, 5.639122590523545e-02,
                -9.214958682077944e-01, 8.689071260438135e-01, 2.840189695206082e-02,
                -3.506711929883037e-01, 9.851138463638244e-01, 1.440758011223270e-02,
                -4.868528161195464e-01, 9.604323781459383e-01, 1.812287875530189e-02,
                -7.955215584636620e-01, 9.757200071559576e-01, 1.826958037918748e-02,
                -9.798840213938864e-01, 9.734329880957604e-01, 5.026188268466660e-03,
                -6.772712303303220e-02, -1.091623307294788e-01, 5.494509440087280e-02,
                -1.911320239476111e-01, -7.684428790753622e-01, 9.924039371581232e-02,
                -1.293584417049455e-01, -9.843401047067886e-01, 2.298164089064930e-02,
                -2.653188588974263e-01, -3.187377804093571e-01, 1.356277329376223e-01,
                -4.805884103058745e-01, -5.883672104526556e-01, 1.092830463202312e-01,
                -4.353820421144143e-01, -9.253041863610407e-01, 5.332719030254342e-02,
                -5.327122733511525e-01, -5.307888414194004e-02, 1.342811802398570e-01,
                -6.817257899596626e-01, -8.006754071171294e-01, 6.668259225153295e-02,
                -6.897274550636672e-01, -9.904094449990978e-01, 1.399432148751916e-02,
                -8.960307328244307e-01, -1.256777467976840e-01, 6.752580403681208e-02,
                -7.189012994269330e-01, -3.612251785232622e-01, 1.014916530440626e-01,
                -8.636764064784566e-01, -6.273567750914366e-01, 5.639122590523545e-02,
                -8.689071260438135e-01, -9.214958682077944e-01, 2.840189695206082e-02,
                -9.851138463638244e-01, -3.506711929883037e-01, 1.440758011223270e-02,
                -9.604323781459383e-01, -4.868528161195464e-01, 1.812287875530189e-02,
                -9.757200071559576e-01, -7.955215584636620e-01, 1.826958037918748e-02,
                -9.734329880957604e-01, -9.798840213938864e-01, 5.026188268466660e-03,
                1.091623307294788e-01, -6.772712303303220e-02, 5.494509440087280e-02,
                7.684428790753622e-01, -1.911320239476111e-01, 9.924039371581232e-02,
                9.843401047067886e-01, -1.293584417049455e-01, 2.298164089064930e-02,
                3.187377804093571e-01, -2.653188588974263e-01, 1.356277329376223e-01,
                5.883672104526556e-01, -4.805884103058745e-01, 1.092830463202312e-01,
                9.253041863610407e-01, -4.353820421144143e-01, 5.332719030254342e-02,
                5.307888414194004e-02, -5.327122733511525e-01, 1.342811802398570e-01,
                8.006754071171294e-01, -6.817257899596626e-01, 6.668259225153295e-02,
                9.904094449990978e-01, -6.897274550636672e-01, 1.399432148751916e-02,
                1.256777467976840e-01, -8.960307328244307e-01, 6.752580403681208e-02,
                3.612251785232622e-01, -7.189012994269330e-01, 1.014916530440626e-01,
                6.273567750914366e-01, -8.636764064784566e-01, 5.639122590523545e-02,
                9.214958682077944e-01, -8.689071260438135e-01, 2.840189695206082e-02,
                3.506711929883037e-01, -9.851138463638244e-01, 1.440758011223270e-02,
                4.868528161195464e-01, -9.604323781459383e-01, 1.812287875530189e-02,
                7.955215584636620e-01, -9.757200071559576e-01, 1.826958037918748e-02,
                9.798840213938864e-01, -9.734329880957604e-01, 5.026188268466660e-03
            };
            double[] xw_copy;

            order = square_minimal_rule_order(degree);
            xw_copy = typeMethods.r8mat_copy_new(3, order, xw);

            return xw_copy;
        }
    }
}