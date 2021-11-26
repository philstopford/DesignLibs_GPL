﻿using System;
using Burkardt.ExactnessNS;
using Burkardt.Types;

namespace ExactnessTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for EXACTNESS_TEST.
        //
        //  Discussion:
        //
        //    EXACTNESS_TEST tests the EXACTNESS library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 January 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("EXACTNESS_TEST");

        Console.WriteLine("  Test the EXACTNESS library.");

        chebyshev1_exactness_test();
        chebyshev2_exactness_test();
        chebyshev3_exactness_test();
        clenshaw_curtis_exactness_test();
        fejer1_exactness_test();
        fejer2_exactness_test();
        gegenbauer_exactness_test();
        hermite_exactness_test();
        hermite_1_exactness_test();
        laguerre_exactness_test();
        laguerre_1_exactness_test();
        legendre_exactness_test();

        Console.WriteLine("");
        Console.WriteLine("EXACTNESS_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void chebyshev1_exactness_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHEBYSHEV1_EXACTNESS_TEST tests rules for the Chebyshev1 integral.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n;

        Console.WriteLine("");
        Console.WriteLine("CHEBYSHEV1_EXACTNESS_TEST");
        Console.WriteLine("  Gauss-Chebyshev1 rules for the Chebyshev1 integral.");
        Console.WriteLine("  Density function rho(x) = 1/sqrt(1-x^2).");
        Console.WriteLine("  Region: -1 <= x <= +1.");
        Console.WriteLine("  Exactness: 2N-1.");

        for (n = 1; n <= 5; n++)
        {
            double[] x = new double[n];
            double[] w = new double[n];
            chebyshev1_set(n, ref x, ref w);
            int p_max = 2 * n;
            Exactness.chebyshev1_exactness(n, x, w, p_max);
        }
    }

    private static void chebyshev1_set(int n, ref double[] x, ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHEBYSHEV1_SET sets a Chebyshev Type 1 quadrature rule.
        //
        //  Discussion:
        //
        //    The integral:
        //
        //      integral ( -1 <= x <= 1 ) f(x) / sqrt ( 1 - x * x ) dx
        //
        //    The quadrature rule:
        //
        //      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order.
        //    N must be between 1 and 10.
        //
        //    Output, double X[N], the abscissas.
        //
        //    Output, double W[N], the weights.
        //
    {
        switch (n)
        {
            case 1:
                x[0] = 0.0;
                w[0] = 3.141592653589793;
                break;
            case 2:
                x[0] = -0.7071067811865475;
                x[1] = 0.7071067811865476;
                w[0] = 1.570796326794897;
                w[1] = 1.570796326794897;
                break;
            case 3:
                x[0] = -0.8660254037844387;
                x[1] = 0.0;
                x[2] = 0.8660254037844387;
                w[0] = 1.047197551196598;
                w[1] = 1.047197551196598;
                w[2] = 1.047197551196598;
                break;
            case 4:
                x[0] = -0.9238795325112867;
                x[1] = -0.3826834323650897;
                x[2] = 0.3826834323650898;
                x[3] = 0.9238795325112867;
                w[0] = 0.7853981633974483;
                w[1] = 0.7853981633974483;
                w[2] = 0.7853981633974483;
                w[3] = 0.7853981633974483;
                break;
            case 5:
                x[0] = -0.9510565162951535;
                x[1] = -0.5877852522924730;
                x[2] = 0.0;
                x[3] = 0.5877852522924731;
                x[4] = 0.9510565162951535;
                w[0] = 0.6283185307179586;
                w[1] = 0.6283185307179586;
                w[2] = 0.6283185307179586;
                w[3] = 0.6283185307179586;
                w[4] = 0.6283185307179586;
                break;
            case 6:
                x[0] = -0.9659258262890682;
                x[1] = -0.7071067811865475;
                x[2] = -0.2588190451025206;
                x[3] = 0.2588190451025207;
                x[4] = 0.7071067811865476;
                x[5] = 0.9659258262890683;
                w[0] = 0.5235987755982988;
                w[1] = 0.5235987755982988;
                w[2] = 0.5235987755982988;
                w[3] = 0.5235987755982988;
                w[4] = 0.5235987755982988;
                w[5] = 0.5235987755982988;
                break;
            case 7:
                x[0] = -0.9749279121818237;
                x[1] = -0.7818314824680295;
                x[2] = -0.4338837391175581;
                x[3] = 0.0;
                x[4] = 0.4338837391175582;
                x[5] = 0.7818314824680298;
                x[6] = 0.9749279121818236;
                w[0] = 0.4487989505128276;
                w[1] = 0.4487989505128276;
                w[2] = 0.4487989505128276;
                w[3] = 0.4487989505128276;
                w[4] = 0.4487989505128276;
                w[5] = 0.4487989505128276;
                w[6] = 0.4487989505128276;
                break;
            case 8:
                x[0] = -0.9807852804032304;
                x[1] = -0.8314696123025453;
                x[2] = -0.5555702330196020;
                x[3] = -0.1950903220161282;
                x[4] = 0.1950903220161283;
                x[5] = 0.5555702330196023;
                x[6] = 0.8314696123025452;
                x[7] = 0.9807852804032304;
                w[0] = 0.3926990816987241;
                w[1] = 0.3926990816987241;
                w[2] = 0.3926990816987241;
                w[3] = 0.3926990816987241;
                w[4] = 0.3926990816987241;
                w[5] = 0.3926990816987241;
                w[6] = 0.3926990816987241;
                w[7] = 0.3926990816987241;
                break;
            case 9:
                x[0] = -0.9848077530122080;
                x[1] = -0.8660254037844385;
                x[2] = -0.6427876096865394;
                x[3] = -0.3420201433256685;
                x[4] = 0.0;
                x[5] = 0.3420201433256688;
                x[6] = 0.6427876096865394;
                x[7] = 0.8660254037844387;
                x[8] = 0.9848077530122080;
                w[0] = 0.3490658503988659;
                w[1] = 0.3490658503988659;
                w[2] = 0.3490658503988659;
                w[3] = 0.3490658503988659;
                w[4] = 0.3490658503988659;
                w[5] = 0.3490658503988659;
                w[6] = 0.3490658503988659;
                w[7] = 0.3490658503988659;
                w[8] = 0.3490658503988659;
                break;
            case 10:
                x[0] = -0.9876883405951377;
                x[1] = -0.8910065241883678;
                x[2] = -0.7071067811865475;
                x[3] = -0.4539904997395467;
                x[4] = -0.1564344650402306;
                x[5] = 0.1564344650402309;
                x[6] = 0.4539904997395468;
                x[7] = 0.7071067811865476;
                x[8] = 0.8910065241883679;
                x[9] = 0.9876883405951378;
                w[0] = 0.3141592653589793;
                w[1] = 0.3141592653589793;
                w[2] = 0.3141592653589793;
                w[3] = 0.3141592653589793;
                w[4] = 0.3141592653589793;
                w[5] = 0.3141592653589793;
                w[6] = 0.3141592653589793;
                w[7] = 0.3141592653589793;
                w[8] = 0.3141592653589793;
                w[9] = 0.3141592653589793;
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("CHEBYSHEV1_SET - Fatal error!");
                Console.WriteLine("  Illegal value of N = " + n + "");
                Console.WriteLine("  Legal values are 1 through 10.");
                break;
        }
    }

    private static void chebyshev2_exactness_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHEBYSHEV2_EXACTNESS_TEST tests rules for the Chebyshev2 integral.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 January 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("CHEBYSHEV2_EXACTNESS_TEST");
        Console.WriteLine("  Gauss-Chebyshev2 rules for the Chebyshev2 integral.");
        Console.WriteLine("  Density function rho(x) = sqrt(1-x^2).");
        Console.WriteLine("  Region: -1 <= x <= +1.");
        Console.WriteLine("  Exactness: 2N-1.");

        for (int n = 1; n <= 5; n++)
        {
            double[] x = new double[n];
            double[] w = new double[n];
            chebyshev2_set(n, ref x, ref w);
            int p_max = 2 * n;
            Exactness.chebyshev2_exactness(n, x, w, p_max);
        }
    }

    private static void chebyshev2_set(int n, ref double[] x, ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHEBYSHEV2_SET sets a Chebyshev Type 2 quadrature rule.
        //
        //  Discussion:
        //
        //    The integral:
        //
        //      integral ( -1 <= x <= 1 ) f(x) * sqrt ( 1 - x * x ) dx
        //
        //    The quadrature rule:
        //
        //      sum ( 1 <= i <= n ) w[i) * f ( x[i) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order.
        //    N must be between 1 and 10.
        //
        //    Output, double X[N], the abscissas.
        //
        //    Output, double W[N], the weights.
        //
    {
        switch (n)
        {
            case 1:
                x[0] = 0.0;
                w[0] = 1.570796326794897;
                break;
            case 2:
                x[0] = -0.5000000000000000;
                x[1] = 0.5000000000000000;
                w[0] = 0.7853981633974484;
                w[1] = 0.7853981633974481;
                break;
            case 3:
                x[0] = -0.7071067811865475;
                x[1] = 0.0;
                x[2] = 0.7071067811865476;
                w[0] = 0.3926990816987243;
                w[1] = 0.7853981633974483;
                w[2] = 0.3926990816987240;
                break;
            case 4:
                x[0] = -0.8090169943749473;
                x[1] = -0.3090169943749473;
                x[2] = 0.3090169943749475;
                x[3] = 0.8090169943749475;
                w[0] = 0.2170787134227061;
                w[1] = 0.5683194499747424;
                w[2] = 0.5683194499747423;
                w[3] = 0.2170787134227060;
                break;
            case 5:
                x[0] = -0.8660254037844387;
                x[1] = -0.5000000000000000;
                x[2] = 0.0;
                x[3] = 0.5000000000000000;
                x[4] = 0.8660254037844387;
                w[0] = 0.1308996938995747;
                w[1] = 0.3926990816987242;
                w[2] = 0.5235987755982988;
                w[3] = 0.3926990816987240;
                w[4] = 0.1308996938995747;
                break;
            case 6:
                x[0] = -0.9009688679024190;
                x[1] = -0.6234898018587335;
                x[2] = -0.2225209339563143;
                x[3] = 0.2225209339563144;
                x[4] = 0.6234898018587336;
                x[5] = 0.9009688679024191;
                w[0] = 0.08448869089158863;
                w[1] = 0.2743330560697779;
                w[2] = 0.4265764164360819;
                w[3] = 0.4265764164360819;
                w[4] = 0.2743330560697778;
                w[5] = 0.08448869089158857;
                break;
            case 7:
                x[0] = -0.9238795325112867;
                x[1] = -0.7071067811865475;
                x[2] = -0.3826834323650897;
                x[3] = 0.0;
                x[4] = 0.3826834323650898;
                x[5] = 0.7071067811865476;
                x[6] = 0.9238795325112867;
                w[0] = 0.05750944903191316;
                w[1] = 0.1963495408493621;
                w[2] = 0.3351896326668110;
                w[3] = 0.3926990816987241;
                w[4] = 0.3351896326668110;
                w[5] = 0.1963495408493620;
                w[6] = 0.05750944903191313;
                break;
            case 8:
                x[0] = -0.9396926207859083;
                x[1] = -0.7660444431189779;
                x[2] = -0.5000000000000000;
                x[3] = -0.1736481776669303;
                x[4] = 0.1736481776669304;
                x[5] = 0.5000000000000000;
                x[6] = 0.7660444431189780;
                x[7] = 0.9396926207859084;
                w[0] = 0.04083294770910712;
                w[1] = 0.1442256007956728;
                w[2] = 0.2617993877991495;
                w[3] = 0.3385402270935190;
                w[4] = 0.3385402270935190;
                w[5] = 0.2617993877991494;
                w[6] = 0.1442256007956727;
                w[7] = 0.04083294770910708;
                break;
            case 9:
                x[0] = -0.9510565162951535;
                x[1] = -0.8090169943749473;
                x[2] = -0.5877852522924730;
                x[3] = -0.3090169943749473;
                x[4] = 0.0;
                x[5] = 0.3090169943749475;
                x[6] = 0.5877852522924731;
                x[7] = 0.8090169943749475;
                x[8] = 0.9510565162951535;
                w[0] = 0.02999954037160818;
                w[1] = 0.1085393567113530;
                w[2] = 0.2056199086476263;
                w[3] = 0.2841597249873712;
                w[4] = 0.3141592653589793;
                w[5] = 0.2841597249873711;
                w[6] = 0.2056199086476263;
                w[7] = 0.1085393567113530;
                w[8] = 0.02999954037160816;
                break;
            case 10:
                x[0] = -0.9594929736144974;
                x[1] = -0.8412535328311811;
                x[2] = -0.6548607339452850;
                x[3] = -0.4154150130018863;
                x[4] = -0.1423148382732850;
                x[5] = 0.1423148382732851;
                x[6] = 0.4154150130018864;
                x[7] = 0.6548607339452851;
                x[8] = 0.8412535328311812;
                x[9] = 0.9594929736144974;
                w[0] = 0.02266894250185884;
                w[1] = 0.08347854093418908;
                w[2] = 0.1631221774548166;
                w[3] = 0.2363135602034873;
                w[4] = 0.2798149423030966;
                w[5] = 0.2798149423030965;
                w[6] = 0.2363135602034873;
                w[7] = 0.1631221774548166;
                w[8] = 0.08347854093418902;
                w[9] = 0.02266894250185884;
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("CHEBYSHEV2_SET - Fatal error!");
                Console.WriteLine("  Illegal value of N = " + n + "");
                Console.WriteLine("  Legal values are 1 through 10.");
                break;
        }
    }

    private static void chebyshev3_exactness_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHEBYSHEV3_EXACTNESS_TEST tests rules for the Chebyshev1 integral.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 January 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n;

        Console.WriteLine("");
        Console.WriteLine("CHEBYSHEV3_EXACTNESS_TEST");
        Console.WriteLine("  Gauss-Chebyshev3 rules for the Chebyshev1 integral.");
        Console.WriteLine("  Density function rho(x) = 1/sqrt(1-x^2).");
        Console.WriteLine("  Region: -1 <= x <= +1.");
        Console.WriteLine("  Exactness: 2N-3.");

        for (n = 1; n <= 5; n++)
        {
            double[] x = new double[n];
            double[] w = new double[n];
            chebyshev3_set(n, x, w);
            int p_max = n switch
            {
                1 => 2 * n,
                _ => 2 * n - 2
            };

            Exactness.chebyshev1_exactness(n, x, w, p_max);
        }
    }

    private static void chebyshev3_set(int n, double[] x, double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHEBYSHEV3_SET sets a Chebyshev Type 3 quadrature rule.
        //
        //  Discussion:
        //
        //    The integral:
        //
        //      integral ( -1 <= x <= 1 ) f(x) / sqrt ( 1 - x * x ) dx
        //
        //    The quadrature rule:
        //
        //      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, integer N, the order.
        //    N must be between 1 and 10.
        //
        //    Output, double X[N], the abscissas.
        //
        //    Output, double W[N], the weights.
        //
    {
        switch (n)
        {
            case 1:
                x[0] = 0.000000000000000;
                w[0] = 3.141592653589793;
                break;
            case 2:
                x[0] = -1.000000000000000;
                x[1] = 1.000000000000000;
                w[0] = 1.570796326794897;
                w[1] = 1.570796326794897;
                break;
            case 3:
                x[0] = -1.000000000000000;
                x[1] = 0.0;
                x[2] = 1.000000000000000;
                w[0] = 0.7853981633974483;
                w[1] = 1.570796326794897;
                w[2] = 0.7853981633974483;
                break;
            case 4:
                x[0] = -1.000000000000000;
                x[1] = -0.5000000000000000;
                x[2] = 0.5000000000000000;
                x[3] = 1.000000000000000;
                w[0] = 0.5235987755982988;
                w[1] = 1.047197551196598;
                w[2] = 1.047197551196598;
                w[3] = 0.5235987755982988;
                break;
            case 5:
                x[0] = -1.000000000000000;
                x[1] = -0.7071067811865475;
                x[2] = 0.0;
                x[3] = 0.7071067811865476;
                x[4] = 1.000000000000000;
                w[0] = 0.3926990816987241;
                w[1] = 0.7853981633974483;
                w[2] = 0.7853981633974483;
                w[3] = 0.7853981633974483;
                w[4] = 0.3926990816987241;
                break;
            case 6:
                x[0] = -1.000000000000000;
                x[1] = -0.8090169943749473;
                x[2] = -0.3090169943749473;
                x[3] = 0.3090169943749475;
                x[4] = 0.8090169943749475;
                x[5] = 1.000000000000000;
                w[0] = 0.3141592653589793;
                w[1] = 0.6283185307179586;
                w[2] = 0.6283185307179586;
                w[3] = 0.6283185307179586;
                w[4] = 0.6283185307179586;
                w[5] = 0.3141592653589793;
                break;
            case 7:
                x[0] = -1.000000000000000;
                x[1] = -0.8660254037844387;
                x[2] = -0.5000000000000000;
                x[3] = 0.0;
                x[4] = 0.5000000000000001;
                x[5] = 0.8660254037844387;
                x[6] = 1.000000000000000;
                w[0] = 0.2617993877991494;
                w[1] = 0.5235987755982988;
                w[2] = 0.5235987755982988;
                w[3] = 0.5235987755982988;
                w[4] = 0.5235987755982988;
                w[5] = 0.5235987755982988;
                w[6] = 0.2617993877991494;
                break;
            case 8:
                x[0] = -1.000000000000000;
                x[1] = -0.9009688679024190;
                x[2] = -0.6234898018587335;
                x[3] = -0.2225209339563143;
                x[4] = 0.2225209339563144;
                x[5] = 0.6234898018587336;
                x[6] = 0.9009688679024191;
                x[7] = 1.000000000000000;
                w[0] = 0.2243994752564138;
                w[1] = 0.4487989505128276;
                w[2] = 0.4487989505128276;
                w[3] = 0.4487989505128276;
                w[4] = 0.4487989505128276;
                w[5] = 0.4487989505128276;
                w[6] = 0.4487989505128276;
                w[7] = 0.2243994752564138;
                break;
            case 9:
                x[0] = -1.000000000000000;
                x[1] = -0.9238795325112867;
                x[2] = -0.7071067811865475;
                x[3] = -0.3826834323650897;
                x[4] = 0.0;
                x[5] = 0.3826834323650898;
                x[6] = 0.7071067811865476;
                x[7] = 0.9238795325112867;
                x[8] = 1.000000000000000;
                w[0] = 0.1963495408493621;
                w[1] = 0.3926990816987241;
                w[2] = 0.3926990816987241;
                w[3] = 0.3926990816987241;
                w[4] = 0.3926990816987241;
                w[5] = 0.3926990816987241;
                w[6] = 0.3926990816987241;
                w[7] = 0.3926990816987241;
                w[8] = 0.1963495408493621;
                break;
            case 10:
                x[0] = -1.000000000000000;
                x[1] = -0.9396926207859083;
                x[2] = -0.7660444431189779;
                x[3] = -0.5000000000000000;
                x[4] = -0.1736481776669303;
                x[5] = 0.1736481776669304;
                x[6] = 0.5000000000000001;
                x[7] = 0.7660444431189780;
                x[8] = 0.9396926207859084;
                x[9] = 1.000000000000000;
                w[0] = 0.1745329251994329;
                w[1] = 0.3490658503988659;
                w[2] = 0.3490658503988659;
                w[3] = 0.3490658503988659;
                w[4] = 0.3490658503988659;
                w[5] = 0.3490658503988659;
                w[6] = 0.3490658503988659;
                w[7] = 0.3490658503988659;
                w[8] = 0.3490658503988659;
                w[9] = 0.1745329251994329;
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("CHEBYSHEV3_SET - Fatal error!");
                Console.WriteLine("  Illegal value of N = " + n + "");
                Console.WriteLine("  Legal values are 1 through 10.");
                break;
        }
    }

    private static void clenshaw_curtis_exactness_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CLENSHAW_CURTIS_EXACTNESS_TEST tests Clenshaw-Curtis rules for the Legendre integral.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 January 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n;

        Console.WriteLine("");
        Console.WriteLine("CLENSHAW_CURTIS_EXACTNESS_TEST");
        Console.WriteLine("  Clenshaw-Curtis rules for the Legendre integral.");
        Console.WriteLine("  Density function rho(x) = 1.");
        Console.WriteLine("  Region: -1 <= x <= +1.");
        Console.WriteLine("  Exactness: N   for N odd,");
        Console.WriteLine("             N-1 for N even.");

        for (n = 1; n <= 5; n++)
        {
            double[] x = new double[n];
            double[] w = new double[n];
            clenshaw_curtis_set(n, x, w);
            int p_max = (n % 2) switch
            {
                1 => n + 1,
                _ => n
            };

            Exactness.legendre_exactness(n, x, w, p_max);
        }
    }

    private static void clenshaw_curtis_set(int n, double[] x, double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CLENSHAW_CURTIS_SET sets a Clenshaw-Curtis quadrature rule.
        //
        //  Discussion:
        //
        //    The integral:
        //
        //      Integral ( -1 <= X <= 1 ) F(X) dX
        //
        //    The quadrature rule:
        //
        //      Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
        //
        //    The abscissas for the rule of order N can be regarded 
        //    as the cosines of equally spaced angles between 180 and 0 degrees:
        //
        //      X(I) = cos ( ( I - 1 ) * PI / ( N - 1 ) )
        //
        //    except for the basic case N = 1, when
        //
        //      X(1) = 0.
        //
        //    A Clenshaw-Curtis rule that uses N points will integrate
        //    exactly all polynomials of degrees 0 through N-1.  If N
        //    is odd, then by symmetry the polynomial of degree N will
        //    also be integrated exactly.
        //
        //    If the value of N is increased in a sensible way, then
        //    the new set of abscissas will include the old ones.  One such
        //    sequence would be N(K) = 2*K+1 for K = 0, 1, 2, ...
        //    Thus, in the table below, the abscissas for order 9 include
        //    those for order 5.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Charles Clenshaw, Alan Curtis,
        //    A Method for Numerical Integration on an Automatic Computer,
        //    Numerische Mathematik,
        //    Volume 2, Number 1, December 1960, pages 197-205.
        //
        //  Parameters:
        //
        //    Input, int N, the order.
        //    N must be between 1 and 17, 33, 65 or 129.
        //
        //    Output, double X[N], the abscissas.
        //
        //    Output, double W[N], the weights.
        //
    {
        switch (n)
        {
            case 1:
                x[0] = 0.00000000000000000000;
                w[0] = 2.00000000000000000000;
                break;
            case 2:
                x[0] = -1.00000000000000000000;
                x[1] = 1.00000000000000000000;

                w[0] = 1.00000000000000000000;
                w[1] = 1.00000000000000000000;
                break;
            case 3:
                x[0] = -1.00000000000000000000;
                x[1] = 0.00000000000000000000;
                x[2] = 1.00000000000000000000;

                w[0] = 0.33333333333333333333;
                w[1] = 1.33333333333333333333;
                w[2] = 0.33333333333333333333;
                break;
            case 4:
                x[0] = -1.00000000000000000000;
                x[1] = -0.50000000000000000000;
                x[2] = 0.50000000000000000000;
                x[3] = 1.00000000000000000000;

                w[0] = 0.11111111111111111111;
                w[1] = 0.88888888888888888889;
                w[2] = 0.88888888888888888889;
                w[3] = 0.11111111111111111111;
                break;
            case 5:
                x[0] = -1.00000000000000000000;
                x[1] = -0.70710678118654752440;
                x[2] = 0.00000000000000000000;
                x[3] = 0.70710678118654752440;
                x[4] = 1.00000000000000000000;

                w[0] = 0.06666666666666666667;
                w[1] = 0.53333333333333333333;
                w[2] = 0.80000000000000000000;
                w[3] = 0.53333333333333333333;
                w[4] = 0.06666666666666666667;
                break;
            case 6:
                x[0] = -1.00000000000000000000;
                x[1] = -0.80901699437494742410;
                x[2] = -0.30901699437494742410;
                x[3] = 0.30901699437494742410;
                x[4] = 0.80901699437493732410;
                x[5] = 1.00000000000000000000;

                w[0] = 0.04000000000000000000;
                w[1] = 0.36074304120001121619;
                w[2] = 0.59925695879998878381;
                w[3] = 0.59925695879998878381;
                w[4] = 0.36074304120001121619;
                w[5] = 0.04000000000000000000;
                break;
            case 7:
                x[0] = -1.00000000000000000000;
                x[1] = -0.86602540378443864676;
                x[2] = -0.50000000000000000000;
                x[3] = 0.00000000000000000000;
                x[4] = 0.50000000000000000000;
                x[5] = 0.86602540378443864676;
                x[6] = 1.00000000000000000000;

                w[0] = 0.02857142857142857143;
                w[1] = 0.25396825396825396825;
                w[2] = 0.45714285714285714286;
                w[3] = 0.52063492063492063492;
                w[4] = 0.45714285714285714286;
                w[5] = 0.25396825396825396825;
                w[6] = 0.02857142857142857143;
                break;
            case 8:
                x[0] = -1.00000000000000000000;
                x[1] = -0.90096886790241912624;
                x[2] = -0.62348980185873353053;
                x[3] = -0.22252093395631440429;
                x[4] = 0.22252093395631440429;
                x[5] = 0.62348980185873353053;
                x[6] = 0.90096886790241910624;
                x[7] = 1.00000000000000000000;

                w[0] = 0.02040816326530612245;
                w[1] = 0.19014100721820835178;
                w[2] = 0.35224242371815911533;
                w[3] = 0.43720840579832641044;
                w[4] = 0.43720840579832641044;
                w[5] = 0.35224242371815911533;
                w[6] = 0.19014100721820835178;
                w[7] = 0.02040816326530612245;
                break;
            case 9:
                x[0] = -1.00000000000000000000;
                x[1] = -0.92387953251128675613;
                x[2] = -0.70710678118654752440;
                x[3] = -0.38268343236508977173;
                x[4] = 0.00000000000000000000;
                x[5] = 0.38268343236508977173;
                x[6] = 0.70710678118654752440;
                x[7] = 0.92387953251128675613;
                x[8] = 1.00000000000000000000;

                w[0] = 0.01587301587301587302;
                w[1] = 0.14621864921601815501;
                w[2] = 0.27936507936507936508;
                w[3] = 0.36171785872048978150;
                w[4] = 0.39365079365079365079;
                w[5] = 0.36171785872048978150;
                w[6] = 0.27936507936507936508;
                w[7] = 0.14621864921601815501;
                w[8] = 0.01587301587301587302;
                break;
            case 10:
                x[0] = -1.00000000000000000000;
                x[1] = -0.93969262078590838405;
                x[2] = -0.76604444311897903520;
                x[3] = -0.50000000000000000000;
                x[4] = -0.17364817766693034885;
                x[5] = 0.17364817766693034885;
                x[6] = 0.50000000000000000000;
                x[7] = 0.76604444311897903520;
                x[8] = 0.93969262078590838405;
                x[9] = 1.00000000000000000000;

                w[0] = 0.01234567901234567901;
                w[1] = 0.11656745657203712296;
                w[2] = 0.22528432333810440813;
                w[3] = 0.30194003527336860670;
                w[4] = 0.34386250580414418320;
                w[5] = 0.34386250580414418320;
                w[6] = 0.30194003527336860670;
                w[7] = 0.22528432333810440813;
                w[8] = 0.11656745657203712296;
                w[9] = 0.01234567901234567901;
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("CLENSHAW_CURTIS_SET - Fatal error!");
                Console.WriteLine("  Illegal value of N = " + n + "");
                Console.WriteLine("  Legal values are 1 to 10.");
                break;
        }
    }

    private static void fejer1_exactness_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FEJER1_EXACTNESS_TEST tests Fejer Type 1 rules for the Legendre integral.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n;

        Console.WriteLine("");
        Console.WriteLine("FEJER1_EXACTNESS_TEST");
        Console.WriteLine("  Fejer Type 1 rules for the Legendre integral.");
        Console.WriteLine("  Density function rho(x) = 1.");
        Console.WriteLine("  Region: -1 <= x <= +1.");
        Console.WriteLine("  Exactness: N   for N odd,");
        Console.WriteLine("             N-1 for N even.");

        for (n = 1; n <= 5; n++)
        {
            double[] x = new double[n];
            double[] w = new double[n];
            fejer1_set(n, x, w);
            int p_max = (n % 2) switch
            {
                1 => n + 1,
                _ => n
            };

            Exactness.legendre_exactness(n, x, w, p_max);
        }
    }

    private static void fejer1_set(int n, double[] xtab, double[] weight )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FEJER1_SET sets abscissas and weights for Fejer type 1 quadrature.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 March 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Philip Davis, Philip Rabinowitz,
        //    Methods of Numerical Integration,
        //    Second Edition,
        //    Dover, 2007,
        //    ISBN: 0486453391,
        //    LC: QA299.3.D28.
        //
        //    Walter Gautschi,
        //    Numerical Quadrature in the Presence of a Singularity,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 4, Number 3, 1967, pages 357-362.
        //
        //    Joerg Waldvogel,
        //    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
        //    BIT Numerical Mathematics,
        //    Volume 43, Number 1, 2003, pages 1-18.
        //
        //  Parameters:
        //
        //    Input, int N, the order.  N should be
        //    between 1 and 10.
        //
        //    Output, double XTAB[N], the abscissas.
        //
        //    Output, double WEIGHT[N], the weights.
        //
    {
        switch (n)
        {
            case 1:
                xtab[0] = 0.000000000000000;
                weight[0] = 2.000000000000000;
                break;
            case 2:
                xtab[0] = -0.7071067811865475;
                xtab[1] = 0.7071067811865475;

                weight[0] = 1.000000000000000;
                weight[1] = 1.000000000000000;
                break;
            case 3:
                xtab[0] = -0.8660254037844387;
                xtab[1] = 0.0000000000000000;
                xtab[2] = 0.8660254037844387;

                weight[0] = 0.4444444444444444;
                weight[1] = 1.111111111111111;
                weight[2] = 0.4444444444444444;
                break;
            case 4:
                xtab[0] = -0.9238795325112867;
                xtab[1] = -0.3826834323650897;
                xtab[2] = 0.3826834323650898;
                xtab[3] = 0.9238795325112867;

                weight[0] = 0.2642977396044841;
                weight[1] = 0.7357022603955158;
                weight[2] = 0.7357022603955158;
                weight[3] = 0.2642977396044841;
                break;
            case 5:
                xtab[0] = -0.9510565162951535;
                xtab[1] = -0.5877852522924730;
                xtab[2] = 0.0000000000000000;
                xtab[3] = 0.5877852522924731;
                xtab[4] = 0.9510565162951535;

                weight[0] = 0.1677812284666835;
                weight[1] = 0.5255521048666498;
                weight[2] = 0.6133333333333333;
                weight[3] = 0.5255521048666498;
                weight[4] = 0.1677812284666835;
                break;
            case 6:
                xtab[0] = -0.9659258262890682;
                xtab[1] = -0.7071067811865475;
                xtab[2] = -0.2588190451025206;
                xtab[3] = 0.2588190451025207;
                xtab[4] = 0.7071067811865476;
                xtab[5] = 0.9659258262890683;

                weight[0] = 0.1186610213812358;
                weight[1] = 0.3777777777777778;
                weight[2] = 0.5035612008409863;
                weight[3] = 0.5035612008409863;
                weight[4] = 0.3777777777777778;
                weight[5] = 0.1186610213812358;
                break;
            case 7:
                xtab[0] = -0.9749279121818237;
                xtab[1] = -0.7818314824680295;
                xtab[2] = -0.4338837391175581;
                xtab[3] = 0.0000000000000000;
                xtab[4] = 0.4338837391175582;
                xtab[5] = 0.7818314824680298;
                xtab[6] = 0.9749279121818236;

                weight[0] = 0.08671618072672234;
                weight[1] = 0.2878313947886919;
                weight[2] = 0.3982415401308441;
                weight[3] = 0.4544217687074830;
                weight[4] = 0.3982415401308441;
                weight[5] = 0.2878313947886919;
                weight[6] = 0.08671618072672234;
                break;
            case 8:
                xtab[0] = -0.9807852804032304;
                xtab[1] = -0.8314696123025453;
                xtab[2] = -0.5555702330196020;
                xtab[3] = -0.1950903220161282;
                xtab[4] = 0.1950903220161283;
                xtab[5] = 0.5555702330196023;
                xtab[6] = 0.8314696123025452;
                xtab[7] = 0.9807852804032304;

                weight[0] = 0.06698294569858981;
                weight[1] = 0.2229879330145788;
                weight[2] = 0.3241525190645244;
                weight[3] = 0.3858766022223071;
                weight[4] = 0.3858766022223071;
                weight[5] = 0.3241525190645244;
                weight[6] = 0.2229879330145788;
                weight[7] = 0.06698294569858981;
                break;
            case 9:
                xtab[0] = -0.9848077530122080;
                xtab[1] = -0.8660254037844385;
                xtab[2] = -0.6427876096865394;
                xtab[3] = -0.3420201433256685;
                xtab[4] = 0.0000000000000000;
                xtab[5] = 0.3420201433256688;
                xtab[6] = 0.6427876096865394;
                xtab[7] = 0.8660254037844387;
                xtab[8] = 0.9848077530122080;

                weight[0] = 0.05273664990990676;
                weight[1] = 0.1791887125220458;
                weight[2] = 0.2640372225410044;
                weight[3] = 0.3308451751681364;
                weight[4] = 0.3463844797178130;
                weight[5] = 0.3308451751681364;
                weight[6] = 0.2640372225410044;
                weight[7] = 0.1791887125220458;
                weight[8] = 0.05273664990990676;
                break;
            case 10:
                xtab[0] = -0.9876883405951377;
                xtab[1] = -0.8910065241883678;
                xtab[2] = -0.7071067811865475;
                xtab[3] = -0.4539904997395467;
                xtab[4] = -0.1564344650402306;
                xtab[5] = 0.1564344650402309;
                xtab[6] = 0.4539904997395468;
                xtab[7] = 0.7071067811865476;
                xtab[8] = 0.8910065241883679;
                xtab[9] = 0.9876883405951378;

                weight[0] = 0.04293911957413078;
                weight[1] = 0.1458749193773909;
                weight[2] = 0.2203174603174603;
                weight[3] = 0.2808792186638755;
                weight[4] = 0.3099892820671425;
                weight[5] = 0.3099892820671425;
                weight[6] = 0.2808792186638755;
                weight[7] = 0.2203174603174603;
                weight[8] = 0.1458749193773909;
                weight[9] = 0.04293911957413078;
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("FEJER1_SET - Fatal error!");
                Console.WriteLine("  Illegal value of N = " + n + "");
                Console.WriteLine("  Legal values are 1 through 10.");
                break;
        }
    }

    private static void fejer2_exactness_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FEJER2_EXACTNESS_TEST tests Fejer Type 2 rules for the Legendre integral.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n;

        Console.WriteLine("");
        Console.WriteLine("FEJER2_EXACTNESS_TEST");
        Console.WriteLine("  Fejer Type 2 rules for the Legendre integral.");
        Console.WriteLine("  Density function rho(x) = 1.");
        Console.WriteLine("  Region: -1 <= x <= +1.");
        Console.WriteLine("  Exactness: N   for N odd,");
        Console.WriteLine("             N-1 for N even.");

        for (n = 1; n <= 5; n++)
        {
            double[] x = new double[n];
            double[] w = new double[n];
            fejer2_set(n, x, w);
            int p_max = (n % 2) switch
            {
                1 => n + 1,
                _ => n
            };

            Exactness.legendre_exactness(n, x, w, p_max);
        }
    }

    private static void fejer2_set(int n, double[] xtab, double[] weight )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FEJER2_SET sets abscissas and weights for Fejer type 2 quadrature.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 March 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Philip Davis, Philip Rabinowitz,
        //    Methods of Numerical Integration,
        //    Second Edition,
        //    Dover, 2007,
        //    ISBN: 0486453391,
        //    LC: QA299.3.D28.
        //
        //    Walter Gautschi,
        //    Numerical Quadrature in the Presence of a Singularity,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 4, Number 3, 1967, pages 357-362.
        //
        //    Joerg Waldvogel,
        //    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
        //    BIT Numerical Mathematics,
        //    Volume 43, Number 1, 2003, pages 1-18.
        //
        //  Parameters:
        //
        //    Input, int N, the order.  
        //    N should be between 1 and 10.
        //
        //    Output, double XTAB[N], the abscissas.
        //
        //    Output, double WEIGHT[N], the weights.
        //
    {
        switch (n)
        {
            case 1:
                xtab[0] = 0.000000000000000;
                weight[0] = 2.000000000000000;
                break;
            case 2:
                xtab[0] = -0.5000000000000000;
                xtab[1] = 0.5000000000000000;

                weight[0] = 1.0000000000000000;
                weight[1] = 1.0000000000000000;
                break;
            case 3:
                xtab[0] = -0.7071067811865476;
                xtab[1] = 0.0000000000000000;
                xtab[2] = 0.7071067811865476;

                weight[0] = 0.6666666666666666;
                weight[1] = 0.6666666666666666;
                weight[2] = 0.6666666666666666;
                break;
            case 4:
                xtab[0] = -0.8090169943749475;
                xtab[1] = -0.3090169943749475;
                xtab[2] = 0.3090169943749475;
                xtab[3] = 0.8090169943749475;

                weight[0] = 0.4254644007500070;
                weight[1] = 0.5745355992499930;
                weight[2] = 0.5745355992499930;
                weight[3] = 0.4254644007500070;
                break;
            case 5:
                xtab[0] = -0.8660254037844387;
                xtab[1] = -0.5000000000000000;
                xtab[2] = 0.0000000000000000;
                xtab[3] = 0.5000000000000000;
                xtab[4] = 0.8660254037844387;

                weight[0] = 0.3111111111111111;
                weight[1] = 0.4000000000000000;
                weight[2] = 0.5777777777777777;
                weight[3] = 0.4000000000000000;
                weight[4] = 0.3111111111111111;
                break;
            case 6:
                xtab[0] = -0.9009688679024191;
                xtab[1] = -0.6234898018587336;
                xtab[2] = -0.2225209339563144;
                xtab[3] = 0.2225209339563144;
                xtab[4] = 0.6234898018587336;
                xtab[5] = 0.9009688679024191;

                weight[0] = 0.2269152467244296;
                weight[1] = 0.3267938603769863;
                weight[2] = 0.4462908928985841;
                weight[3] = 0.4462908928985841;
                weight[4] = 0.3267938603769863;
                weight[5] = 0.2269152467244296;
                break;
            case 7:
                xtab[0] = -0.9238795325112867;
                xtab[1] = -0.7071067811865476;
                xtab[2] = -0.3826834323650898;
                xtab[3] = 0.0000000000000000;
                xtab[4] = 0.3826834323650898;
                xtab[5] = 0.7071067811865476;
                xtab[6] = 0.9238795325112867;

                weight[0] = 0.1779646809620499;
                weight[1] = 0.2476190476190476;
                weight[2] = 0.3934638904665215;
                weight[3] = 0.3619047619047619;
                weight[4] = 0.3934638904665215;
                weight[5] = 0.2476190476190476;
                weight[6] = 0.1779646809620499;
                break;
            case 8:
                xtab[0] = -0.9396926207859084;
                xtab[1] = -0.7660444431189780;
                xtab[2] = -0.5000000000000000;
                xtab[3] = -0.1736481776669304;
                xtab[4] = 0.1736481776669304;
                xtab[5] = 0.5000000000000000;
                xtab[6] = 0.7660444431189780;
                xtab[7] = 0.9396926207859084;

                weight[0] = 0.1397697435050225;
                weight[1] = 0.2063696457302284;
                weight[2] = 0.3142857142857143;
                weight[3] = 0.3395748964790348;
                weight[4] = 0.3395748964790348;
                weight[5] = 0.3142857142857143;
                weight[6] = 0.2063696457302284;
                weight[7] = 0.1397697435050225;
                break;
            case 9:
                xtab[0] = -0.9510565162951535;
                xtab[1] = -0.8090169943749475;
                xtab[2] = -0.5877852522924731;
                xtab[3] = -0.3090169943749475;
                xtab[4] = 0.0000000000000000;
                xtab[5] = 0.3090169943749475;
                xtab[6] = 0.5877852522924731;
                xtab[7] = 0.8090169943749475;
                xtab[8] = 0.9510565162951535;

                weight[0] = 0.1147810750857217;
                weight[1] = 0.1654331942222276;
                weight[2] = 0.2737903534857068;
                weight[3] = 0.2790112502222170;
                weight[4] = 0.3339682539682539;
                weight[5] = 0.2790112502222170;
                weight[6] = 0.2737903534857068;
                weight[7] = 0.1654331942222276;
                weight[8] = 0.1147810750857217;
                break;
            case 10:
                xtab[0] = -0.9594929736144974;
                xtab[1] = -0.8412535328311812;
                xtab[2] = -0.6548607339452851;
                xtab[3] = -0.4154150130018864;
                xtab[4] = -0.1423148382732851;
                xtab[5] = 0.1423148382732851;
                xtab[6] = 0.4154150130018864;
                xtab[7] = 0.6548607339452851;
                xtab[8] = 0.8412535328311812;
                xtab[9] = 0.9594929736144974;

                weight[0] = 0.09441954173982806;
                weight[1] = 0.1411354380109716;
                weight[2] = 0.2263866903636005;
                weight[3] = 0.2530509772156453;
                weight[4] = 0.2850073526699544;
                weight[5] = 0.2850073526699544;
                weight[6] = 0.2530509772156453;
                weight[7] = 0.2263866903636005;
                weight[8] = 0.1411354380109716;
                weight[9] = 0.09441954173982806;
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("FEJER2_SET - Fatal error!");
                Console.WriteLine("  Illegal value of N = " + n + "");
                Console.WriteLine("  Legal values are 1 through 10.");
                break;
        }
    }

    private static void gegenbauer_exactness_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GEGENBAUER_EXACTNESS_TEST tests rules for the Gegenbauer integral.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] w1 =  {
                1.2485988353771993
            }
            ;
        double[] w2 =  {
                0.6242994176885995,
                0.6242994176885995
            }
            ;
        double[] w3 =  {
                0.2837724625857273,
                0.6810539102057455,
                0.2837724625857273
            }
            ;
        double[] w4 =  {
                0.1379302690657785,
                0.4863691486228214,
                0.4863691486228208,
                0.1379302690657786
            }
            ;
        double[] w5 =  {
                0.0725955752894624,
                0.3156051535278124,
                0.4721973777426502,
                0.3156051535278118,
                0.0725955752894624
            }
            ;
        double[] x1 =  {
                0.0000000000000000
            }
            ;
        double[] x2 =  {
                -0.4264014327112208,
                0.4264014327112208
            }
            ;
        double[] x3 =  {
                -0.6324555320336757,
                0.0000000000000000,
                0.6324555320336757
            }
            ;
        double[] x4 =  {
                -0.7455376618816977,
                -0.2752317970082527,
                0.2752317970082527,
                0.7455376618816980
            }
            ;
        double[] x5 =  {
                -0.8137803260309515,
                -0.4553315257658559,
                0.0000000000000001,
                0.4553315257658557,
                0.8137803260309517
            }
            ;

        double lambda = 1.75;

        Console.WriteLine("");
        Console.WriteLine("GEGENBAUER_EXACTNESS_TEST");
        Console.WriteLine("  Gauss-Gegenbauer rules for the Gegenbauer integral.");
        Console.WriteLine("  Density function rho(x) = 1.");
        Console.WriteLine("  Region: -1 <= x <= +1.");
        Console.WriteLine("  Exactness: 2*N-1.");
        Console.WriteLine("  Lambda = " + lambda + "");

        for (int n = 1; n <= 5; n++)
        {
            double[] w;
            double[] x;
            switch (n)
            {
                case 1:
                    x = typeMethods.r8vec_copy_new(n, x1);
                    w = typeMethods.r8vec_copy_new(n, w1);
                    break;
                case 2:
                    x = typeMethods.r8vec_copy_new(n, x2);
                    w = typeMethods.r8vec_copy_new(n, w2);
                    break;
                case 3:
                    x = typeMethods.r8vec_copy_new(n, x3);
                    w = typeMethods.r8vec_copy_new(n, w3);
                    break;
                case 4:
                    x = typeMethods.r8vec_copy_new(n, x4);
                    w = typeMethods.r8vec_copy_new(n, w4);
                    break;
                case 5:
                    x = typeMethods.r8vec_copy_new(n, x5);
                    w = typeMethods.r8vec_copy_new(n, w5);
                    break;
                default:
                    x = new double[1];
                    w = new double[1];
                    break;
            }

            legendre_set(n, ref x, ref w);
            int p_max = 2 * n;
            Exactness.legendre_exactness(n, x, w, p_max);
        }
    }

    private static void hermite_exactness_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HERMITE_EXACTNESS_TEST tests Gauss-Hermite rules for the Hermite integral.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("HERMITE_EXACTNESS_TEST");
        Console.WriteLine("  Gauss-Hermite rules for the Hermite integral.");
        Console.WriteLine("  Density function rho(x) = exp(-x^2).");
        Console.WriteLine("  Region: -oo < x < +oo.");
        Console.WriteLine("  Exactness: 2N-1.");

        for (int n = 1; n <= 5; n++)
        {
            double[] x = new double[n];
            double[] w = new double[n];
            hermite_set(n, ref x, ref w);
            int p_max = 2 * n;
            Exactness.hermite_exactness(n, x, w, p_max);
        }
    }

    private static void hermite_set(int n, ref double[] x, ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HERMITE_SET sets abscissas and weights for Hermite quadrature.
        //
        //  Discussion:
        //
        //    The integral:
        //
        //      integral ( -oo < x < +oo ) f(x) * rho(x) dx
        //
        //    The weight:
        //
        //      rho(x) = exp ( - x * x )
        //
        //    The quadrature rule:
        //
        //      sum ( 1 <= i <= n ) w(i) * f ( x(i) ).
        //
        //    Mathematica can numerically estimate the abscissas 
        //    of order N to P digits by the command:
        //
        //      NSolve [ HermiteH [ n, x ] == 0, x, p ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 April 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz, Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    National Bureau of Standards, 1964,
        //    ISBN: 0-486-61272-4,
        //    LC: QA47.A34.
        //
        //    Vladimir Krylov,
        //    Approximate Calculation of Integrals,
        //    Dover, 2006,
        //    ISBN: 0486445798,
        //    LC: QA311.K713.
        //
        //    Arthur Stroud, Don Secrest,
        //    Gaussian Quadrature Formulas,
        //    Prentice Hall, 1966,
        //    LC: QA299.4G3S7.
        //
        //    Stephen Wolfram,
        //    The Mathematica Book,
        //    Fourth Edition,
        //    Cambridge University Press, 1999,
        //    ISBN: 0-521-64314-7,
        //    LC: QA76.95.W65.
        //
        //    Daniel Zwillinger, editor,
        //    CRC Standard Mathematical Tables and Formulae,
        //    30th Edition,
        //    CRC Press, 1996,
        //    ISBN: 0-8493-2479-3,
        //    LC: QA47.M315.
        //
        //  Parameters:
        //
        //    Input, int N, the order.
        //    N must be between 1 and 20, or 31/32/33, 63/64/65, 127/128/129.
        //
        //    Output, double X[N], the abscissas.
        //
        //    Output, double W[N], the weights.
        //
    {
        switch (n)
        {
            case 1:
                x[0] = 0.0;

                w[0] = 1.77245385090551602729816748334;
                break;
            case 2:
                x[0] = -0.707106781186547524400844362105;
                x[1] = 0.707106781186547524400844362105;

                w[0] = 0.886226925452758013649083741671;
                w[1] = 0.886226925452758013649083741671;
                break;
            case 3:
                x[0] = -0.122474487139158904909864203735E+01;
                x[1] = 0.0;
                x[2] = 0.122474487139158904909864203735E+01;

                w[0] = 0.295408975150919337883027913890;
                w[1] = 0.118163590060367735153211165556E+01;
                w[2] = 0.295408975150919337883027913890;
                break;
            case 4:
                x[0] = -0.165068012388578455588334111112E+01;
                x[1] = -0.524647623275290317884060253835;
                x[2] = 0.524647623275290317884060253835;
                x[3] = 0.165068012388578455588334111112E+01;

                w[0] = 0.813128354472451771430345571899E-01;
                w[1] = 0.804914090005512836506049184481;
                w[2] = 0.804914090005512836506049184481;
                w[3] = 0.813128354472451771430345571899E-01;
                break;
            case 5:
                x[0] = -0.202018287045608563292872408814E+01;
                x[1] = -0.958572464613818507112770593893;
                x[2] = 0.0;
                x[3] = 0.958572464613818507112770593893;
                x[4] = 0.202018287045608563292872408814E+01;

                w[0] = 0.199532420590459132077434585942E-01;
                w[1] = 0.393619323152241159828495620852;
                w[2] = 0.945308720482941881225689324449;
                w[3] = 0.393619323152241159828495620852;
                w[4] = 0.199532420590459132077434585942E-01;
                break;
            case 6:
                x[0] = -0.235060497367449222283392198706E+01;
                x[1] = -0.133584907401369694971489528297E+01;
                x[2] = -0.436077411927616508679215948251;
                x[3] = 0.436077411927616508679215948251;
                x[4] = 0.133584907401369694971489528297E+01;
                x[5] = 0.235060497367449222283392198706E+01;

                w[0] = 0.453000990550884564085747256463E-02;
                w[1] = 0.157067320322856643916311563508;
                w[2] = 0.724629595224392524091914705598;
                w[3] = 0.724629595224392524091914705598;
                w[4] = 0.157067320322856643916311563508;
                w[5] = 0.453000990550884564085747256463E-02;
                break;
            case 7:
                x[0] = -0.265196135683523349244708200652E+01;
                x[1] = -0.167355162876747144503180139830E+01;
                x[2] = -0.816287882858964663038710959027;
                x[3] = 0.0;
                x[4] = 0.816287882858964663038710959027;
                x[5] = 0.167355162876747144503180139830E+01;
                x[6] = 0.265196135683523349244708200652E+01;

                w[0] = 0.971781245099519154149424255939E-03;
                w[1] = 0.545155828191270305921785688417E-01;
                w[2] = 0.425607252610127800520317466666;
                w[3] = 0.810264617556807326764876563813;
                w[4] = 0.425607252610127800520317466666;
                w[5] = 0.545155828191270305921785688417E-01;
                w[6] = 0.971781245099519154149424255939E-03;
                break;
            case 8:
                x[0] = -0.293063742025724401922350270524E+01;
                x[1] = -0.198165675669584292585463063977E+01;
                x[2] = -0.115719371244678019472076577906E+01;
                x[3] = -0.381186990207322116854718885584;
                x[4] = 0.381186990207322116854718885584;
                x[5] = 0.115719371244678019472076577906E+01;
                x[6] = 0.198165675669584292585463063977E+01;
                x[7] = 0.293063742025724401922350270524E+01;

                w[0] = 0.199604072211367619206090452544E-03;
                w[1] = 0.170779830074134754562030564364E-01;
                w[2] = 0.207802325814891879543258620286;
                w[3] = 0.661147012558241291030415974496;
                w[4] = 0.661147012558241291030415974496;
                w[5] = 0.207802325814891879543258620286;
                w[6] = 0.170779830074134754562030564364E-01;
                w[7] = 0.199604072211367619206090452544E-03;
                break;
            case 9:
                x[0] = -0.319099320178152760723004779538E+01;
                x[1] = -0.226658058453184311180209693284E+01;
                x[2] = -0.146855328921666793166701573925E+01;
                x[3] = -0.723551018752837573322639864579;
                x[4] = 0.0;
                x[5] = 0.723551018752837573322639864579;
                x[6] = 0.146855328921666793166701573925E+01;
                x[7] = 0.226658058453184311180209693284E+01;
                x[8] = 0.319099320178152760723004779538E+01;

                w[0] = 0.396069772632643819045862946425E-04;
                w[1] = 0.494362427553694721722456597763E-02;
                w[2] = 0.884745273943765732879751147476E-01;
                w[3] = 0.432651559002555750199812112956;
                w[4] = 0.720235215606050957124334723389;
                w[5] = 0.432651559002555750199812112956;
                w[6] = 0.884745273943765732879751147476E-01;
                w[7] = 0.494362427553694721722456597763E-02;
                w[8] = 0.396069772632643819045862946425E-04;
                break;
            case 10:
                x[0] = -0.343615911883773760332672549432E+01;
                x[1] = -0.253273167423278979640896079775E+01;
                x[2] = -0.175668364929988177345140122011E+01;
                x[3] = -0.103661082978951365417749191676E+01;
                x[4] = -0.342901327223704608789165025557;
                x[5] = 0.342901327223704608789165025557;
                x[6] = 0.103661082978951365417749191676E+01;
                x[7] = 0.175668364929988177345140122011E+01;
                x[8] = 0.253273167423278979640896079775E+01;
                x[9] = 0.343615911883773760332672549432E+01;

                w[0] = 0.764043285523262062915936785960E-05;
                w[1] = 0.134364574678123269220156558585E-02;
                w[2] = 0.338743944554810631361647312776E-01;
                w[3] = 0.240138611082314686416523295006;
                w[4] = 0.610862633735325798783564990433;
                w[5] = 0.610862633735325798783564990433;
                w[6] = 0.240138611082314686416523295006;
                w[7] = 0.338743944554810631361647312776E-01;
                w[8] = 0.134364574678123269220156558585E-02;
                w[9] = 0.764043285523262062915936785960E-05;
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("HERMITE_SET - Fatal error!");
                Console.WriteLine("  Illegal value of N = " + n + "");
                Console.WriteLine("  Legal values are 1 to 10.");
                break;
        }
    }

    private static void hermite_1_exactness_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HERMITE_1_EXACTNESS_TEST tests Gauss-Hermite rules for the Hermite integral.
        //
        //  Discussion:
        //
        //    Instead of the usual density rho(x)=exp(-x*x), these rules apply to
        //    approximating the integral:
        //      I(f) = integral ( -oo < x < +oo ) f(x) dx
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {

        Console.WriteLine("");
        Console.WriteLine("HERMITE_1_EXACTNESS_TEST");
        Console.WriteLine("  Gauss-Hermite rules for the Hermite integral.");
        Console.WriteLine("  Density function rho(x) = 1.");
        Console.WriteLine("  Region: -oo < x < +oo.");
        Console.WriteLine("  Exactness: 2N-1.");

        for (int n = 1; n <= 5; n++)
        {
            double[] x = new double[n];
            double[] w = new double[n];
            hermite_1_set(n, ref x, ref w);
            //
            //  Standardize the rule by multiplying every weight w(i) by exp(-x(i)^2).
            //
            for (int i = 0; i < n; i++)
            {
                w[i] = Math.Exp(-x[i] * x[i]) * w[i];
            }

            //
            //  Now test the rule in standard form.
            //
            int p_max = 2 * n;
            Exactness.hermite_exactness(n, x, w, p_max);
        }
    }

    private static void hermite_1_set(int n, ref double[] x, ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HERMITE_1_SET sets abscissas and weights for Hermite quadrature.
        //
        //  Discussion:
        //  
        //    This routine is for the case with unit density:
        //      integral ( -oo < x < +oo ) f(x) dx
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 May 2014
        //  
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order.
        //    N must be between 1 and 10.
        //
        //    Output, double X[N], the abscissas.
        //
        //    Output, double W[N], the weights.
        //
    {
        switch (n)
        {
            case 1:
                x[0] = 0.0;

                w[0] = 1.7724538509055161;
                break;
            case 2:
                x[0] = -0.707106781186547524400844362105;
                x[1] = 0.707106781186547524400844362105;

                w[0] = 1.4611411826611391;
                w[1] = 1.4611411826611391;
                break;
            case 3:
                x[0] = -0.122474487139158904909864203735E+01;
                x[1] = 0.0;
                x[2] = 0.122474487139158904909864203735E+01;

                w[0] = 1.3239311752136438;
                w[1] = 1.1816359006036774;
                w[2] = 1.3239311752136438;
                break;
            case 4:
                x[0] = -0.165068012388578455588334111112E+01;
                x[1] = -0.524647623275290317884060253835;
                x[2] = 0.524647623275290317884060253835;
                x[3] = 0.165068012388578455588334111112E+01;

                w[0] = 1.2402258176958150;
                w[1] = 1.0599644828949693;
                w[2] = 1.0599644828949693;
                w[3] = 1.2402258176958150;
                break;
            case 5:
                x[0] = -0.202018287045608563292872408814E+01;
                x[1] = -0.958572464613818507112770593893;
                x[2] = 0.0;
                x[3] = 0.958572464613818507112770593893;
                x[4] = 0.202018287045608563292872408814E+01;

                w[0] = 1.1814886255359869;
                w[1] = 0.98658099675142830;
                w[2] = 0.94530872048294190;
                w[3] = 0.98658099675142830;
                w[4] = 1.1814886255359869;
                break;
            case 6:
                x[0] = -0.235060497367449222283392198706E+01;
                x[1] = -0.133584907401369694971489528297E+01;
                x[2] = -0.436077411927616508679215948251;
                x[3] = 0.436077411927616508679215948251;
                x[4] = 0.133584907401369694971489528297E+01;
                x[5] = 0.235060497367449222283392198706E+01;

                w[0] = 1.1369083326745253;
                w[1] = 0.93558055763118075;
                w[2] = 0.87640133443623058;
                w[3] = 0.87640133443623058;
                w[4] = 0.93558055763118075;
                w[5] = 1.1369083326745253;
                break;
            case 7:
                x[0] = -0.265196135683523349244708200652E+01;
                x[1] = -0.167355162876747144503180139830E+01;
                x[2] = -0.816287882858964663038710959027;
                x[3] = 0.0;
                x[4] = 0.816287882858964663038710959027;
                x[5] = 0.167355162876747144503180139830E+01;
                x[6] = 0.265196135683523349244708200652E+01;

                w[0] = 1.1013307296103216;
                w[1] = 0.89718460022518409;
                w[2] = 0.82868730328363926;
                w[3] = 0.81026461755680734;
                w[4] = 0.82868730328363926;
                w[5] = 0.89718460022518409;
                w[6] = 1.1013307296103216;
                break;
            case 8:
                x[0] = -0.293063742025724401922350270524E+01;
                x[1] = -0.198165675669584292585463063977E+01;
                x[2] = -0.115719371244678019472076577906E+01;
                x[3] = -0.381186990207322116854718885584;
                x[4] = 0.381186990207322116854718885584;
                x[5] = 0.115719371244678019472076577906E+01;
                x[6] = 0.198165675669584292585463063977E+01;
                x[7] = 0.293063742025724401922350270524E+01;

                w[0] = 1.0719301442479805;
                w[1] = 0.86675260656338138;
                w[2] = 0.79289004838640131;
                w[3] = 0.76454412865172916;
                w[4] = 0.76454412865172916;
                w[5] = 0.79289004838640131;
                w[6] = 0.86675260656338138;
                w[7] = 1.0719301442479805;
                break;
            case 9:
                x[0] = -0.319099320178152760723004779538E+01;
                x[1] = -0.226658058453184311180209693284E+01;
                x[2] = -0.146855328921666793166701573925E+01;
                x[3] = -0.723551018752837573322639864579;
                x[4] = 0.0;
                x[5] = 0.723551018752837573322639864579;
                x[6] = 0.146855328921666793166701573925E+01;
                x[7] = 0.226658058453184311180209693284E+01;
                x[8] = 0.319099320178152760723004779538E+01;

                w[0] = 1.0470035809766838;
                w[1] = 0.84175270147867043;
                w[2] = 0.76460812509455023;
                w[3] = 0.73030245274509220;
                w[4] = 0.72023521560605097;
                w[5] = 0.73030245274509220;
                w[6] = 0.76460812509455023;
                w[7] = 0.84175270147867043;
                w[8] = 1.0470035809766838;
                break;
            case 10:
                x[0] = -0.343615911883773760332672549432E+01;
                x[1] = -0.253273167423278979640896079775E+01;
                x[2] = -0.175668364929988177345140122011E+01;
                x[3] = -0.103661082978951365417749191676E+01;
                x[4] = -0.342901327223704608789165025557;
                x[5] = 0.342901327223704608789165025557;
                x[6] = 0.103661082978951365417749191676E+01;
                x[7] = 0.175668364929988177345140122011E+01;
                x[8] = 0.253273167423278979640896079775E+01;
                x[9] = 0.343615911883773760332672549432E+01;

                w[0] = 1.0254516913657352;
                w[1] = 0.82066612640481640;
                w[2] = 0.74144193194356511;
                w[3] = 0.70329632310490608;
                w[4] = 0.68708185395127341;
                w[5] = 0.68708185395127341;
                w[6] = 0.70329632310490608;
                w[7] = 0.74144193194356511;
                w[8] = 0.82066612640481640;
                w[9] = 1.0254516913657352;
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("HERMITE_1_SET - Fatal error!");
                Console.WriteLine("  Illegal value of N = " + n + "");
                Console.WriteLine("  Legal values are 1 to 10.");
                break;
        }
    }

    private static void laguerre_exactness_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGUERRE_EXACTNESS_TEST tests Gauss-Laguerre rules for the Laguerre integral.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 January 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("LAGUERRE_EXACTNESS_TEST");
        Console.WriteLine("  Gauss-Laguerre rules for the Laguerre integral.");
        Console.WriteLine("  Density function rho(x) = exp(-x).");
        Console.WriteLine("  Region: 0 <= x < +oo.");
        Console.WriteLine("  Exactness: 2N-1.");

        for (int n = 1; n <= 5; n++)
        {
            double[] x = new double[n];
            double[] w = new double[n];
            laguerre_set(n, ref x, ref w);
            int p_max = 2 * n;
            Exactness.laguerre_exactness(n, x, w, p_max);
        }
    }

    private static void laguerre_set(int n, ref double[] x, ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGUERRE_SET sets abscissas and weights for Laguerre quadrature.
        //
        //  Discussion:
        //
        //    The abscissas are the zeroes of the Laguerre polynomial L(N)(X).
        //
        //    The integral:
        //
        //      Integral ( 0 <= X < +oo ) exp ( -X ) * F(X) dX
        //
        //    The quadrature rule:
        //
        //      Sum ( 1 <= I <= N ) W(I) * f ( X(I) )
        //
        //    The integral:
        //
        //      Integral ( 0 <= X < +oo ) F(X) dX
        //
        //    The quadrature rule:
        //
        //      Sum ( 1 <= I <= N ) W(I) * exp ( X(I) ) * f ( X(I) )
        //
        //    Mathematica can numerically estimate the abscissas for the
        //    n-th order polynomial to p digits of precision by the command:
        //
        //      NSolve [ LaguerreL[n,x] == 0, x, p ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 April 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz, Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    National Bureau of Standards, 1964,
        //    ISBN: 0-486-61272-4,
        //    LC: QA47.A34.
        //
        //    Vladimir Krylov,
        //    Approximate Calculation of Integrals,
        //    Dover, 2006,
        //    ISBN: 0486445798,
        //    LC: QA311.K713.
        //
        //    Arthur Stroud, Don Secrest,
        //    Gaussian Quadrature Formulas,
        //    Prentice Hall, 1966,
        //    LC: QA299.4G3S7.
        //
        //    Stephen Wolfram,
        //    The Mathematica Book,
        //    Fourth Edition,
        //    Cambridge University Press, 1999,
        //    ISBN: 0-521-64314-7,
        //    LC: QA76.95.W65.
        //
        //    Daniel Zwillinger, editor,
        //    CRC Standard Mathematical Tables and Formulae,
        //    30th Edition,
        //    CRC Press, 1996,
        //    ISBN: 0-8493-2479-3.
        //
        //  Parameters:
        //
        //    Input, int N, the order.
        //    N must be between 1 and 20, 31/32/33, 63/64/65, 127/128/129.
        //
        //    Output, double X[N], the abscissas.
        //
        //    Output, double W[N], the weights.
        //
    {
        switch (n)
        {
            case 1:
                x[0] = 1.00000000000000000000000000000E+00;

                w[0] = 1.00000000000000000000000000000E+00;
                break;
            case 2:
                x[0] = 0.585786437626904951198311275790E+00;
                x[1] = 3.41421356237309504880168872421E+00;

                w[0] = 0.85355339059327376220042218105E+00;
                w[1] = 0.146446609406726237799577818948E+00;
                break;
            case 3:
                x[0] = 0.415774556783479083311533873128E+00;
                x[1] = 2.29428036027904171982205036136E+00;
                x[2] = 6.28994508293747919686641576551E+00;

                w[0] = 0.71109300992917301544959019114E+00;
                w[1] = 0.27851773356924084880144488846E+00;
                w[2] = 0.010389256501586135748964920401E+00;
                break;
            case 4:
                x[0] = 0.322547689619392311800361459104E+00;
                x[1] = 1.74576110115834657568681671252E+00;
                x[2] = 4.53662029692112798327928538496E+00;
                x[3] = 9.39507091230113312923353644342E+00;

                w[0] = 0.60315410434163360163596602382E+00;
                w[1] = 0.35741869243779968664149201746E+00;
                w[2] = 0.03888790851500538427243816816E+00;
                w[3] = 0.0005392947055613274501037905676E+00;
                break;
            case 5:
                x[0] = 0.263560319718140910203061943361E+00;
                x[1] = 1.41340305910651679221840798019E+00;
                x[2] = 3.59642577104072208122318658878E+00;
                x[3] = 7.08581000585883755692212418111E+00;
                x[4] = 12.6408008442757826594332193066E+00;

                w[0] = 0.52175561058280865247586092879E+00;
                w[1] = 0.3986668110831759274541333481E+00;
                w[2] = 0.0759424496817075953876533114E+00;
                w[3] = 0.00361175867992204845446126257E+00;
                w[4] = 0.00002336997238577622789114908455E+00;
                break;
            case 6:
                x[0] = 0.222846604179260689464354826787E+00;
                x[1] = 1.18893210167262303074315092194E+00;
                x[2] = 2.99273632605931407769132528451E+00;
                x[3] = 5.77514356910451050183983036943E+00;
                x[4] = 9.83746741838258991771554702994E+00;
                x[5] = 15.9828739806017017825457915674E+00;

                w[0] = 0.45896467394996359356828487771E+00;
                w[1] = 0.4170008307721209941133775662E+00;
                w[2] = 0.1133733820740449757387061851E+00;
                w[3] = 0.01039919745314907489891330285E+00;
                w[4] = 0.000261017202814932059479242860E+00;
                w[5] = 8.98547906429621238825292053E-07;
                break;
            case 7:
                x[0] = 0.193043676560362413838247885004E+00;
                x[1] = 1.02666489533919195034519944317E+00;
                x[2] = 2.56787674495074620690778622666E+00;
                x[3] = 4.90035308452648456810171437810E+00;
                x[4] = 8.18215344456286079108182755123E+00;
                x[5] = 12.7341802917978137580126424582E+00;
                x[6] = 19.3957278622625403117125820576E+00;

                w[0] = 0.40931895170127390213043288002E+00;
                w[1] = 0.4218312778617197799292810054E+00;
                w[2] = 0.1471263486575052783953741846E+00;
                w[3] = 0.0206335144687169398657056150E+00;
                w[4] = 0.00107401014328074552213195963E+00;
                w[5] = 0.0000158654643485642012687326223E+00;
                w[6] = 3.17031547899558056227132215E-08;
                break;
            case 8:
                x[0] = 0.170279632305100999788861856608E+00;
                x[1] = 0.903701776799379912186020223555E+00;
                x[2] = 2.25108662986613068930711836697E+00;
                x[3] = 4.26670017028765879364942182690E+00;
                x[4] = 7.04590540239346569727932548212E+00;
                x[5] = 10.7585160101809952240599567880E+00;
                x[6] = 15.7406786412780045780287611584E+00;
                x[7] = 22.8631317368892641057005342974E+00;

                w[0] = 0.36918858934163752992058283938E+00;
                w[1] = 0.4187867808143429560769785813E+00;
                w[2] = 0.175794986637171805699659867E+00;
                w[3] = 0.033343492261215651522132535E+00;
                w[4] = 0.0027945362352256725249389241E+00;
                w[5] = 0.00009076508773358213104238501E+00;
                w[6] = 8.4857467162725315448680183E-07;
                w[7] = 1.04800117487151038161508854E-09;
                break;
            case 9:
                x[0] = 0.152322227731808247428107073127E+00;
                x[1] = 0.807220022742255847741419210952E+00;
                x[2] = 2.00513515561934712298303324701E+00;
                x[3] = 3.78347397333123299167540609364E+00;
                x[4] = 6.20495677787661260697353521006E+00;
                x[5] = 9.37298525168757620180971073215E+00;
                x[6] = 13.4662369110920935710978818397E+00;
                x[7] = 18.8335977889916966141498992996E+00;
                x[8] = 26.3740718909273767961410072937E+00;

                w[0] = 0.336126421797962519673467717606E+00;
                w[1] = 0.411213980423984387309146942793E+00;
                w[2] = 0.199287525370885580860575607212E+00;
                w[3] = 0.0474605627656515992621163600479E+00;
                w[4] = 0.00559962661079458317700419900556E+00;
                w[5] = 0.000305249767093210566305412824291E+00;
                w[6] = 6.59212302607535239225572284875E-06;
                w[7] = 4.1107693303495484429024104033E-08;
                w[8] = 3.29087403035070757646681380323E-11;
                break;
            case 10:
                x[0] = 0.137793470540492430830772505653E+00;
                x[1] = 0.729454549503170498160373121676E+00;
                x[2] = 1.80834290174031604823292007575E+00;
                x[3] = 3.40143369785489951448253222141E+00;
                x[4] = 5.55249614006380363241755848687E+00;
                x[5] = 8.33015274676449670023876719727E+00;
                x[6] = 11.8437858379000655649185389191E+00;
                x[7] = 16.2792578313781020995326539358E+00;
                x[8] = 21.9965858119807619512770901956E+00;
                x[9] = 29.9206970122738915599087933408E+00;

                w[0] = 0.30844111576502014154747083468E+00;
                w[1] = 0.4011199291552735515157803099E+00;
                w[2] = 0.218068287611809421588648523E+00;
                w[3] = 0.062087456098677747392902129E+00;
                w[4] = 0.009501516975181100553839072E+00;
                w[5] = 0.0007530083885875387754559644E+00;
                w[6] = 0.00002825923349599565567422564E+00;
                w[7] = 4.249313984962686372586577E-07;
                w[8] = 1.839564823979630780921535E-09;
                w[9] = 9.911827219609008558377547E-13;
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("LAGUERRE_SET - Fatal error!");
                Console.WriteLine("  Illegal value of N = " + n + "");
                Console.WriteLine("  Legal values are 1 to 10");
                break;
        }
    }

    private static void laguerre_1_exactness_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGUERRE_1_EXACTNESS_TEST tests Gauss-Laguerre rules for the Laguerre integral with rho=1.
        //
        //  Discussion:
        //
        //    Instead of the usual density rho(x)=exp(-x), these rules apply to
        //    approximating the integral:
        //      I(f) = integral ( 0 <= x < +oo ) f(x) dx
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n;

        Console.WriteLine("");
        Console.WriteLine("LAGUERRE_1_EXACTNESS_TEST");
        Console.WriteLine("  Gauss-Laguerre rules for the Laguerre integral.");
        Console.WriteLine("  Density function rho(x) = 1.");
        Console.WriteLine("  Region: 0 <= x < +oo.");
        Console.WriteLine("  Exactness: 2N-1.");

        for (n = 1; n <= 5; n++)
        {
            double[] x = new double[n];
            double[] w = new double[n];
            laguerre_1_set(n, ref x, ref w);
            //
            //  Standardize the rule by multiplying every weight w(i) by exp(-x(i)).
            //
            int i;
            for (i = 0; i < n; i++)
            {
                w[i] = Math.Exp(-x[i]) * w[i];
            }

            //
            //  Now test the rule in standard form.
            //
            int p_max = 2 * n;
            Exactness.laguerre_exactness(n, x, w, p_max);
        }
    }

    private static void laguerre_1_set(int n, ref double[] x, ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGUERRE_1_SET sets abscissas and weights for Laguerre quadrature.
        //
        //  Discussion:
        //
        //    This routine is specialized for the case where the density function is 1.
        //
        //    The integral:
        //      I(f) = integral ( 0 <= x < +oo ) f(x) dx
        //    The quadrature rule:
        //      Q(f) = sum ( 1 <= i <= n ) w(i) * f ( x(i) )
        //
        //  Licensing:
        //
        //   This code is distributed under the GNU LGPL license.
        // 
        //  Modified:
        //  
        //    16 May 2014
        //  
        //  Author:
        //  
        //    John Burkardt
        //  
        //  Parameters:
        //  
        //    Input, int N, the order.
        //    N must be between 1 and 10.
        //  
        //    Output, double X[N], the abscissas.
        //  
        //    Output, double W[N], the weights.
        //
    {
        switch (n)
        {
            case 1:
                x[0] = 1.00000000000000000000000000000E+00;

                w[0] = 2.7182818284590451;
                break;
            case 2:
                x[0] = 0.585786437626904951198311275790E+00;
                x[1] = 3.41421356237309504880168872421E+00;

                w[0] = 1.5333260331194167;
                w[1] = 4.4509573350545928;
                break;
            case 3:
                x[0] = 0.415774556783479083311533873128E+00;
                x[1] = 2.29428036027904171982205036136E+00;
                x[2] = 6.28994508293747919686641576551E+00;

                w[0] = 1.0776928592709207;
                w[1] = 2.7621429619015876;
                w[2] = 5.6010946254344267;
                break;
            case 4:
                x[0] = 0.322547689619392311800361459104E+00;
                x[1] = 1.74576110115834657568681671252E+00;
                x[2] = 4.53662029692112798327928538496E+00;
                x[3] = 9.39507091230113312923353644342E+00;

                w[0] = 0.83273912383788917;
                w[1] = 2.0481024384542965;
                w[2] = 3.6311463058215168;
                w[3] = 6.4871450844076604;
                break;
            case 5:
                x[0] = 0.263560319718140910203061943361E+00;
                x[1] = 1.41340305910651679221840798019E+00;
                x[2] = 3.59642577104072208122318658878E+00;
                x[3] = 7.08581000585883755692212418111E+00;
                x[4] = 12.6408008442757826594332193066E+00;

                w[0] = 0.67909404220775038;
                w[1] = 1.6384878736027471;
                w[2] = 2.7694432423708375;
                w[3] = 4.3156569009208940;
                w[4] = 7.2191863543544450;
                break;
            case 6:
                x[0] = 0.222846604179260689464354826787E+00;
                x[1] = 1.18893210167262303074315092194E+00;
                x[2] = 2.99273632605931407769132528451E+00;
                x[3] = 5.77514356910451050183983036943E+00;
                x[4] = 9.83746741838258991771554702994E+00;
                x[5] = 15.9828739806017017825457915674E+00;

                w[0] = 0.57353550742273818;
                w[1] = 1.3692525907123045;
                w[2] = 2.2606845933826722;
                w[3] = 3.3505245823554555;
                w[4] = 4.8868268002108213;
                w[5] = 7.8490159455958279;
                break;
            case 7:
                x[0] = 0.193043676560362413838247885004E+00;
                x[1] = 1.02666489533919195034519944317E+00;
                x[2] = 2.56787674495074620690778622666E+00;
                x[3] = 4.90035308452648456810171437810E+00;
                x[4] = 8.18215344456286079108182755123E+00;
                x[5] = 12.7341802917978137580126424582E+00;
                x[6] = 19.3957278622625403117125820576E+00;

                w[0] = 0.49647759753997234;
                w[1] = 1.1776430608611976;
                w[2] = 1.9182497816598063;
                w[3] = 2.7718486362321113;
                w[4] = 3.8412491224885148;
                w[5] = 5.3806782079215330;
                w[6] = 8.4054324868283103;
                break;
            case 8:
                x[0] = 0.170279632305100999788861856608E+00;
                x[1] = 0.903701776799379912186020223555E+00;
                x[2] = 2.25108662986613068930711836697E+00;
                x[3] = 4.26670017028765879364942182690E+00;
                x[4] = 7.04590540239346569727932548212E+00;
                x[5] = 10.7585160101809952240599567880E+00;
                x[6] = 15.7406786412780045780287611584E+00;
                x[7] = 22.8631317368892641057005342974E+00;

                w[0] = 0.43772341049291136;
                w[1] = 1.0338693476655976;
                w[2] = 1.6697097656587756;
                w[3] = 2.3769247017585995;
                w[4] = 3.2085409133479259;
                w[5] = 4.2685755108251344;
                w[6] = 5.8180833686719184;
                w[7] = 8.9062262152922216;
                break;
            case 9:
                x[0] = 0.152322227731808247428107073127E+00;
                x[1] = 0.807220022742255847741419210952E+00;
                x[2] = 2.00513515561934712298303324701E+00;
                x[3] = 3.78347397333123299167540609364E+00;
                x[4] = 6.20495677787661260697353521006E+00;
                x[5] = 9.37298525168757620180971073215E+00;
                x[6] = 13.4662369110920935710978818397E+00;
                x[7] = 18.8335977889916966141498992996E+00;
                x[8] = 26.3740718909273767961410072937E+00;

                w[0] = 0.39143112431563987;
                w[1] = 0.92180502852896307;
                w[2] = 1.4801279099429154;
                w[3] = 2.0867708075492613;
                w[4] = 2.7729213897119713;
                w[5] = 3.5916260680922663;
                w[6] = 4.6487660021402037;
                w[7] = 6.2122754197471348;
                w[8] = 9.3632182377057980;
                break;
            case 10:
                x[0] = 0.137793470540492430830772505653E+00;
                x[1] = 0.729454549503170498160373121676E+00;
                x[2] = 1.80834290174031604823292007575E+00;
                x[3] = 3.40143369785489951448253222141E+00;
                x[4] = 5.55249614006380363241755848687E+00;
                x[5] = 8.33015274676449670023876719727E+00;
                x[6] = 11.8437858379000655649185389191E+00;
                x[7] = 16.2792578313781020995326539358E+00;
                x[8] = 21.9965858119807619512770901956E+00;
                x[9] = 29.9206970122738915599087933408E+00;

                w[0] = 0.35400973860699630;
                w[1] = 0.83190230104358065;
                w[2] = 1.3302885617493283;
                w[3] = 1.8630639031111309;
                w[4] = 2.4502555580830108;
                w[5] = 3.1227641551351848;
                w[6] = 3.9341526955615240;
                w[7] = 4.9924148721930299;
                w[8] = 6.5722024851307994;
                w[9] = 9.7846958403746243;
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("LAGUERRE_1_SET - Fatal error!");
                Console.WriteLine("  Illegal value of N = " + n + "");
                Console.WriteLine("  Legal values are 1 to 10");
                break;
        }
    }

    private static void legendre_exactness_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEGENDRE_EXACTNESS_TEST tests Gauss-Legendre rules for the Legendre integral.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n;

        Console.WriteLine("");
        Console.WriteLine("LEGENDRE_EXACTNESS_TEST");
        Console.WriteLine("  Gauss-Legendre rules for the Legendre integral.");
        Console.WriteLine("  Density function rho(x) = 1.");
        Console.WriteLine("  Region: -1 <= x <= +1.");
        Console.WriteLine("  Exactness: 2*N-1.");

        for (n = 1; n <= 5; n++)
        {
            double[] x = new double[n];
            double[] w = new double[n];
            legendre_set(n, ref x, ref w);
            int p_max = 2 * n;
            Exactness.legendre_exactness(n, x, w, p_max);
        }
    }

    private static void legendre_set(int n, ref double[] x, ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEGENDRE_SET sets abscissas and weights for Gauss-Legendre quadrature.
        //
        //  Discussion:
        //
        //    The integral:
        //
        //      Integral ( -1 <= X <= 1 ) F(X) dX
        //
        //    Quadrature rule:
        //
        //      Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
        //
        //    The quadrature rule is exact for all polynomials through degree 2*N-1.
        //
        //    The abscissas are the zeroes of the Legendre polynomial P(N)(X).
        //
        //    Mathematica can compute the abscissas and weights of a Gauss-Legendre
        //    rule of order N for the interval [A,B] with P digits of precision
        //    by the commands:
        //
        //    Needs["NumericalDifferentialEquationAnalysis`"]
        //    GaussianQuadratureWeights [ n, a, b, p ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 April 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz, Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    National Bureau of Standards, 1964,
        //    ISBN: 0-486-61272-4,
        //    LC: QA47.A34.
        //
        //    Vladimir Krylov,
        //    Approximate Calculation of Integrals,
        //    Dover, 2006,
        //    ISBN: 0486445798.
        //    LC: QA311.K713.
        //
        //    Arthur Stroud, Don Secrest,
        //    Gaussian Quadrature Formulas,
        //    Prentice Hall, 1966,
        //    LC: QA299.4G3S7.
        //
        //    Stephen Wolfram,
        //    The Mathematica Book,
        //    Fourth Edition,
        //    Cambridge University Press, 1999,
        //    ISBN: 0-521-64314-7,
        //    LC: QA76.95.W65.
        //
        //    Daniel Zwillinger, editor,
        //    CRC Standard Mathematical Tables and Formulae,
        //    30th Edition,
        //    CRC Press, 1996,
        //    ISBN: 0-8493-2479-3,
        //    LC: QA47.M315.
        //
        //  Parameters:
        //
        //    Input, int N, the order.
        //    N must be between 1 and 33 or 63/64/65, 127/128/129, 
        //    255/256/257.
        //
        //    Output, double X[N], the abscissas.
        //
        //    Output, double W[N], the weights.
        //
    {
        switch (n)
        {
            case 1:
                x[0] = 0.000000000000000000000000000000;

                w[0] = 2.000000000000000000000000000000;
                break;
            case 2:
                x[0] = -0.577350269189625764509148780502;
                x[1] = 0.577350269189625764509148780502;

                w[0] = 1.000000000000000000000000000000;
                w[1] = 1.000000000000000000000000000000;
                break;
            case 3:
                x[0] = -0.774596669241483377035853079956;
                x[1] = 0.000000000000000000000000000000;
                x[2] = 0.774596669241483377035853079956;

                w[0] = 0.555555555555555555555555555556;
                w[1] = 0.888888888888888888888888888889;
                w[2] = 0.555555555555555555555555555556;
                break;
            case 4:
                x[0] = -0.861136311594052575223946488893;
                x[1] = -0.339981043584856264802665759103;
                x[2] = 0.339981043584856264802665759103;
                x[3] = 0.861136311594052575223946488893;

                w[0] = 0.347854845137453857373063949222;
                w[1] = 0.652145154862546142626936050778;
                w[2] = 0.652145154862546142626936050778;
                w[3] = 0.347854845137453857373063949222;
                break;
            case 5:
                x[0] = -0.906179845938663992797626878299;
                x[1] = -0.538469310105683091036314420700;
                x[2] = 0.000000000000000000000000000000;
                x[3] = 0.538469310105683091036314420700;
                x[4] = 0.906179845938663992797626878299;

                w[0] = 0.236926885056189087514264040720;
                w[1] = 0.478628670499366468041291514836;
                w[2] = 0.568888888888888888888888888889;
                w[3] = 0.478628670499366468041291514836;
                w[4] = 0.236926885056189087514264040720;
                break;
            case 6:
                x[0] = -0.932469514203152027812301554494;
                x[1] = -0.661209386466264513661399595020;
                x[2] = -0.238619186083196908630501721681;
                x[3] = 0.238619186083196908630501721681;
                x[4] = 0.661209386466264513661399595020;
                x[5] = 0.932469514203152027812301554494;

                w[0] = 0.171324492379170345040296142173;
                w[1] = 0.360761573048138607569833513838;
                w[2] = 0.467913934572691047389870343990;
                w[3] = 0.467913934572691047389870343990;
                w[4] = 0.360761573048138607569833513838;
                w[5] = 0.171324492379170345040296142173;
                break;
            case 7:
                x[0] = -0.949107912342758524526189684048;
                x[1] = -0.741531185599394439863864773281;
                x[2] = -0.405845151377397166906606412077;
                x[3] = 0.000000000000000000000000000000;
                x[4] = 0.405845151377397166906606412077;
                x[5] = 0.741531185599394439863864773281;
                x[6] = 0.949107912342758524526189684048;

                w[0] = 0.129484966168869693270611432679;
                w[1] = 0.279705391489276667901467771424;
                w[2] = 0.381830050505118944950369775489;
                w[3] = 0.417959183673469387755102040816;
                w[4] = 0.381830050505118944950369775489;
                w[5] = 0.279705391489276667901467771424;
                w[6] = 0.129484966168869693270611432679;
                break;
            case 8:
                x[0] = -0.960289856497536231683560868569;
                x[1] = -0.796666477413626739591553936476;
                x[2] = -0.525532409916328985817739049189;
                x[3] = -0.183434642495649804939476142360;
                x[4] = 0.183434642495649804939476142360;
                x[5] = 0.525532409916328985817739049189;
                x[6] = 0.796666477413626739591553936476;
                x[7] = 0.960289856497536231683560868569;

                w[0] = 0.101228536290376259152531354310;
                w[1] = 0.222381034453374470544355994426;
                w[2] = 0.313706645877887287337962201987;
                w[3] = 0.362683783378361982965150449277;
                w[4] = 0.362683783378361982965150449277;
                w[5] = 0.313706645877887287337962201987;
                w[6] = 0.222381034453374470544355994426;
                w[7] = 0.101228536290376259152531354310;
                break;
            case 9:
                x[0] = -0.968160239507626089835576203;
                x[1] = -0.836031107326635794299429788;
                x[2] = -0.613371432700590397308702039;
                x[3] = -0.324253423403808929038538015;
                x[4] = 0.000000000000000000000000000;
                x[5] = 0.324253423403808929038538015;
                x[6] = 0.613371432700590397308702039;
                x[7] = 0.836031107326635794299429788;
                x[8] = 0.968160239507626089835576203;

                w[0] = 0.081274388361574411971892158111;
                w[1] = 0.18064816069485740405847203124;
                w[2] = 0.26061069640293546231874286942;
                w[3] = 0.31234707704000284006863040658;
                w[4] = 0.33023935500125976316452506929;
                w[5] = 0.31234707704000284006863040658;
                w[6] = 0.26061069640293546231874286942;
                w[7] = 0.18064816069485740405847203124;
                w[8] = 0.081274388361574411971892158111;
                break;
            case 10:
                x[0] = -0.973906528517171720077964012;
                x[1] = -0.865063366688984510732096688;
                x[2] = -0.679409568299024406234327365;
                x[3] = -0.433395394129247190799265943;
                x[4] = -0.148874338981631210884826001;
                x[5] = 0.148874338981631210884826001;
                x[6] = 0.433395394129247190799265943;
                x[7] = 0.679409568299024406234327365;
                x[8] = 0.865063366688984510732096688;
                x[9] = 0.973906528517171720077964012;

                w[0] = 0.066671344308688137593568809893;
                w[1] = 0.14945134915058059314577633966;
                w[2] = 0.21908636251598204399553493423;
                w[3] = 0.26926671930999635509122692157;
                w[4] = 0.29552422471475287017389299465;
                w[5] = 0.29552422471475287017389299465;
                w[6] = 0.26926671930999635509122692157;
                w[7] = 0.21908636251598204399553493423;
                w[8] = 0.14945134915058059314577633966;
                w[9] = 0.066671344308688137593568809893;
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("LEGENDRE_SET - Fatal error!");
                Console.WriteLine("  Illegal value of N = " + n + "");
                Console.WriteLine("  Legal values are 1 to 10.");
                break;
        }
    }
}