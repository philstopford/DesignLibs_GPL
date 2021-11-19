using System;
using Burkardt.Chebyshev;
using Burkardt.ChebyshevNS;

namespace TOMS446Test;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TOMS446_TEST.
        //
        //  Discussion:
        //
        //    TOMS446_TEST tests the TOMS446 library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 September 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("TOMS446_TEST");
        Console.WriteLine("  Test the TOMS446 library.");

        cheby_test();
        dfrnt_test();
        echeb_test();
        edcheb_test();
        mltply_test();
        ntgrt_test();

        Console.WriteLine("");
        Console.WriteLine("TOMS446_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void cheby_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHEBY_TEST tests CHEBY, which computes Chebyshev series.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 September 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;
        int nf = 5;
        int npl = 10;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("CHEBY_TEST");
        Console.WriteLine("  CHEBY computes the Chebyshev series for several functions.");

        x = Chebyshev.cheby(nf, npl, functn);

        Console.WriteLine("");
        Console.WriteLine("          Sin(x)          Cos(x)        Sin(2x)         Cos(2x)           X^5");
        Console.WriteLine("");
            
        for (i = 0; i < npl; i++)
        {
            string cout = "";
            for (j = 0; j < nf; j++)
            {
                cout += "  " + x[i + j * npl].ToString().PadLeft(14);
            }

            Console.WriteLine(cout);
        }
    }

    private static void dfrnt_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DFRNT_TEST tests DFRNT, which computes the Chebyshev series of a derivative.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 September 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;
        int nf = 5;
        int npl = 10;
        double[] x;
        double[] x2;
        double[] x3;

        Console.WriteLine("");
        Console.WriteLine("DFRNT_TEST");
        Console.WriteLine("  DFRNT computes the Chebyshev series for the derivative");
        Console.WriteLine("  of several functions.");

        x = Chebyshev.cheby(nf, npl, functn);
        x2 = new double[npl];

        for (j = 0; j < nf; j++)
        {
            for (i = 0; i < npl; i++)
            {
                x2[i] = x[i + j * npl];
            }

            x3 = ChebyshevSeries.dfrnt(x2, npl);
            for (i = 0; i < npl; i++)
            {
                x[i + j * npl] = x3[i];
            }
        }

        Console.WriteLine("");
        Console.WriteLine("  Chebyshev series for d/dx of:");
        Console.WriteLine("");
        Console.WriteLine("        Sin(x)      Cos(x)    Sin(2x)     Cos(2x)       X^5");
        Console.WriteLine("");
            
        for (i = 0; i < npl; i++)
        {
            string cout = "";
            for (j = 0; j < nf; j++)
            {
                cout += "  " + x[i + j * npl].ToString().PadLeft(14);
            }

            Console.WriteLine(cout);
        }
    }

    private static void echeb_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ECHEB_TEST tests ECHEB, which evaluates a Chebyshev series.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 September 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fval;
        double[] fxj;
        int i;
        int j;
        int k;
        int nf = 5;
        int npl = 10;
        int nx;
        double[] x;
        double[] x2;
        double xval;

        nx = 6;

        Console.WriteLine("");
        Console.WriteLine("ECHEB_TEST");
        Console.WriteLine("  ECHEB evaluates a Chebyshev series.");

        x = Chebyshev.cheby(nf, npl, functn);
        x2 = new double[npl];

        for (j = 0; j < nf; j++)
        {
            for (i = 0; i < npl; i++)
            {
                x2[i] = x[i + j * npl];
            }

            Console.WriteLine("");
            switch (j)
            {
                case 0:
                    Console.WriteLine("  Sin(x)");
                    break;
                case 1:
                    Console.WriteLine("  Cos(x)");
                    break;
                case 2:
                    Console.WriteLine("  Sin(2x)");
                    break;
                case 3:
                    Console.WriteLine("  Cos(2x)");
                    break;
                case 4:
                    Console.WriteLine("  x^5");
                    break;
            }

            Console.WriteLine("");

            for (k = 0; k < nx; k++)
            {
                xval = 2.0 * k / (nx - 1) - 1.0;

                fxj = functn(xval);

                fval = ChebyshevSeries.echeb(xval, x2, npl);

                Console.WriteLine("  " + xval.ToString().PadLeft(14)
                                       + "  " + fxj[j].ToString().PadLeft(14)
                                       + "  " + fval.ToString().PadLeft(14) + "");
            }
        }
    }

    private static void edcheb_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EDCHEB_TEST tests EDCHEB, which evaluates the derivative of a Chebyshev series.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 September 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fval;
        double[] fxj;
        int i;
        int j;
        int k;
        int nf = 5;
        int npl = 10;
        int nx;
        double[] x;
        double[] x2;
        double xval;

        nx = 6;

        Console.WriteLine("");
        Console.WriteLine("EDCHEB_TEST");
        Console.WriteLine("  EDCHEB evaluates the derivative of a Chebyshev series.");

        x = Chebyshev.cheby(nf, npl, functn);
        x2 = new double[npl];

        for (j = 0; j < nf; j++)
        {
            for (i = 0; i < npl; i++)
            {
                x2[i] = x[i + j * npl];
            }

            Console.WriteLine("");
            switch (j)
            {
                case 0:
                    Console.WriteLine("  Sin(x)");
                    break;
                case 1:
                    Console.WriteLine("  Cos(x)");
                    break;
                case 2:
                    Console.WriteLine("  Sin(2x)");
                    break;
                case 3:
                    Console.WriteLine("  Cos(2x)");
                    break;
                case 4:
                    Console.WriteLine("  x^5");
                    break;
            }

            Console.WriteLine("");

            for (k = 0; k < nx; k++)
            {
                xval = 2.0 * k / (nx - 1) - 1.0;

                fxj = functn_d(xval);

                fval = ChebyshevSeries.edcheb(xval, x2, npl);

                Console.WriteLine("  " + xval.ToString().PadLeft(14)
                                       + "  " + fxj[j].ToString().PadLeft(14)
                                       + "  " + fval.ToString().PadLeft(14) + "");
            }
        }
    }

    private static void mltply_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MLTPLY_TEST tests MLTPLY, which multiplies two Chebyshev series.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 September 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int nf = 5;
        int npl = 10;
        double[] x;
        double[] x1;
        double[] x2;
        double[] x3;

        Console.WriteLine("");
        Console.WriteLine("MLTPLY_TEST");
        Console.WriteLine("  MLTPLY computes the product of two Chebyshev series.");
        Console.WriteLine("");
        Console.WriteLine("  Multiply series for SIN(X) and COS(X)");
        Console.WriteLine("  and compare with series for 1/2*SIN(2X).");

        x = Chebyshev.cheby(nf, npl, functn);

        x1 = new double[npl];
        x2 = new double[npl];

        for (i = 0; i < npl; i++)
        {
            x1[i] = x[i + 0 * npl];
            x2[i] = x[i + 1 * npl];
            x[i + 2 * npl] = 0.5 * x[i + 2 * npl];
        }

        x3 = ChebyshevSeries.mltply_new(x1, x2, npl);

        Console.WriteLine("");
        Console.WriteLine("          Sin(x)          Cos(x)       1/2*Sin(2x)         RESULT");
        Console.WriteLine("");

        for (i = 0; i < npl; i++)
        {
            Console.WriteLine("  " + x[i + 0 * npl].ToString().PadLeft(14)
                                   + "  " + x[i + 1 * npl].ToString().PadLeft(14)
                                   + "  " + x[i + 2 * npl].ToString().PadLeft(14)
                                   + "  " + x3[i].ToString().PadLeft(14) + "");
        }
    }

    private static void ntgrt_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NTGRT_TEST tests NTGRT, which computes the Chebyshev series of an indefinite integral.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 September 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;
        int nf = 5;
        int npl = 10;
        double[] x;
        double[] x2;
        double[] x3;

        Console.WriteLine("");
        Console.WriteLine("NTGRT_TEST");
        Console.WriteLine("  NTGRT computes the Chebyshev series for the indefinite");
        Console.WriteLine("  integral of several functions.");

        x = Chebyshev.cheby(nf, npl, functn);
        x2 = new double[npl];

        for (j = 0; j < nf; j++)
        {
            for (i = 0; i < npl; i++)
            {
                x2[i] = x[i + j * npl];
            }

            x3 = ChebyshevSeries.ntgrt(x2, npl);
            for (i = 0; i < npl; i++)
            {
                x[i + j * npl] = x3[i];
            }
        }

        Console.WriteLine("");
        Console.WriteLine("  Chebyshev series for indefinite integral of:");
        Console.WriteLine("");
        Console.WriteLine("        Sin(x)      Cos(x)    Sin(2x)     Cos(2x)       X^5");
        Console.WriteLine("");

        for (i = 0; i < npl; i++)
        {
            string cout = "";
            for (j = 0; j < nf; j++)
            {
                cout += "  " + x[i + j * npl].ToString().PadLeft(14);
            }

            Console.WriteLine(cout);
        }
    }

    private static double[] functn(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FUNCTN evaluates several functions at X.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 September 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double FXJ[5], the derivative values.
        //
    {
        double[] fxj;

        fxj = new double[5];

        fxj[0] = Math.Sin(x);
        fxj[1] = Math.Cos(x);
        fxj[2] = Math.Sin(2.0 * x);
        fxj[3] = Math.Cos(2.0 * x);
        fxj[4] = Math.Pow(x, 5);

        return fxj;
    }

    private static double[] functn_d(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FUNCTN_D evaluates the derivatives of several functions at X.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 September 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double FXJ[5], the derivative values.
        //
    {
        double[] fxj;

        fxj = new double[5];

        fxj[0] = Math.Cos(x);
        fxj[1] = -Math.Sin(x);
        fxj[2] = 2.0 * Math.Cos(2.0 * x);
        fxj[3] = -2.0 * Math.Sin(2.0 * x);
        fxj[4] = 5.0 * Math.Pow(x, 4);

        return fxj;
    }
}