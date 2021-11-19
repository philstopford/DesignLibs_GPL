using System;
using System.Globalization;
using Burkardt.AppliedStatistics;

namespace ASA047Test;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for ASA047_TEST.
        //
        //  Discussion:
        //
        //    ASA047_TEST tests the ASA047 library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 February 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("ASA047_TEST:");
        Console.WriteLine("  Test the ASA047 library.");

        test01 ( );
        test02 ( );
        test03 ( );
        test04 ( );

        Console.WriteLine("");
        Console.WriteLine("ASA047_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01 ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 demonstrates the use of NELMIN on ROSENBROCK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 February 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int icount = 0;
        int ifault = 0;
        int numres = 0;

        int n = 2;

        double[] start = new double[n];
        double[] step = new double[n];
        double[] xmin = new double[n];

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Apply NELMIN to ROSENBROCK function.");

        start[0] = -1.2;
        start[1] =  1.0;

        const double reqmin = 1.0E-08;

        step[0] = 1.0;
        step[1] = 1.0;

        const int konvge = 10;
        const int kcount = 500;

        Console.WriteLine("");
        Console.WriteLine("  Starting point X:");
        Console.WriteLine("");
        for (int i = 0; i < n; i++ )
        {
            Console.WriteLine("  " + start[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        double ynewlo = rosenbrock ( start );

        Console.WriteLine("");
        Console.WriteLine("  F(X) = " + ynewlo + "");

        Algorithms.nelmin ( rosenbrock, n, ref start, ref xmin, ref ynewlo, reqmin, step,
            konvge, kcount, ref icount, ref numres, ref ifault );

        Console.WriteLine("");
        Console.WriteLine("  Return code IFAULT = " + ifault + "");
        Console.WriteLine("");
        Console.WriteLine("  Estimate of minimizing value X*:");
        Console.WriteLine("");
        for (int i = 0; i < n; i++ )
        {
            Console.WriteLine("  " + xmin[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  F(X*) = " + ynewlo + "");

        Console.WriteLine("");
        Console.WriteLine("  Number of iterations = " + icount + "");
        Console.WriteLine("  Number of restarts =   " + numres + "");
    }

    private static double rosenbrock ( double[] x )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ROSENBROCK evaluates the Rosenbrock parabolic value function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 February 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    R ONeill,
        //    Algorithm AS 47:
        //    Function Minimization Using a Simplex Procedure,
        //    Applied Statistics,
        //    Volume 20, Number 3, 1971, pages 338-345.
        //
        //  Parameters:
        //
        //    Input, double X[2], the argument.
        //
        //    Output, double ROSENBROCK, the value of the function.
        //
    {
        double fx1 = x[1] - x[0] * x[0];
        double fx2 = 1.0 - x[0];

        double fx = 100.0 * fx1 * fx1 + fx2 * fx2;

        return fx;
    }

    private static void test02 ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 demonstrates the use of NELMIN on POWELL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 February 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int icount = 0;
        int ifault = 0;
        int numres = 0;

        const int n = 4;

        double[] start = new double[n];
        double[] step = new double[n];
        double[] xmin = new double[n];

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  Apply NELMIN to POWELL quartic function.");

        start[0] =   3.0;
        start[1] = - 1.0;
        start[2] =   0.0;
        start[3] =   1.0;

        double reqmin = 1.0E-08;

        step[0] = 1.0;
        step[1] = 1.0;
        step[2] = 1.0;
        step[3] = 1.0;

        int konvge = 10;
        int kcount = 500;

        Console.WriteLine("");
        Console.WriteLine("  Starting point X:");
        Console.WriteLine("");
        for (int i = 0; i < n; i++ )
        {
            Console.WriteLine("  " + start[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        double ynewlo = powell ( start );

        Console.WriteLine("");
        Console.WriteLine("  F(X) = " + ynewlo + "");

        Algorithms.nelmin ( powell, n, ref start, ref xmin, ref ynewlo, reqmin, step,
            konvge, kcount, ref icount, ref numres, ref ifault );

        Console.WriteLine("");
        Console.WriteLine("  Return code IFAULT = " + ifault + "");
        Console.WriteLine("");
        Console.WriteLine("  Estimate of minimizing value X*:");
        Console.WriteLine("");
        for (int i = 0; i < n; i++ )
        {
            Console.WriteLine("  " + xmin[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  F(X*) = " + ynewlo + "");

        Console.WriteLine("");
        Console.WriteLine("  Number of iterations = " + icount + "");
        Console.WriteLine("  Number of restarts =   " + numres + "");

    }

    private static double powell ( double[] x )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POWELL evaluates the Powell quartic function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 February 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    R ONeill,
        //    Algorithm AS 47:
        //    Function Minimization Using a Simplex Procedure,
        //    Applied Statistics,
        //    Volume 20, Number 3, 1971, pages 338-345.
        //
        //  Parameters:
        //
        //    Input, double X[4], the argument.
        //
        //    Output, double POWELL, the value of the function.
        //
    {
        double fx1 = x[0] + 10.0 * x[1];
        double fx2 = x[2] - x[3];
        double fx3 = x[1] - 2.0 * x[2];
        double fx4 = x[0] - x[3];

        double fx = fx1 * fx1
                    +  5.0 * fx2 * fx2
                    +            fx3 * fx3 * fx3 * fx3
                    + 10.0 * fx4 * fx4 * fx4 * fx4;

        return fx;
    }


    private static void test03 ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 demonstrates the use of NELMIN on HELICAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 February 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int icount = 0;
        int ifault = 0;
        int numres = 0;

        const int n = 3;

        double[] start = new double[n];
        double[] step = new double[n];
        double[] xmin = new double[n];

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  Apply NELMIN to the HELICAL function.");

        start[0] = - 1.0;
        start[1] =   0.0;
        start[2] =   0.0;

        double reqmin = 1.0E-08;

        step[0] = 1.0;
        step[1] = 1.0;
        step[2] = 1.0;

        int konvge = 10;
        int kcount = 500;

        Console.WriteLine("");
        Console.WriteLine("  Starting point X:");
        Console.WriteLine("");
        for (int i = 0; i < n; i++ )
        {
            Console.WriteLine("  " + start[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        double ynewlo = helical ( start );

        Console.WriteLine("");
        Console.WriteLine("  F(X) = " + ynewlo + "");

        Algorithms.nelmin ( helical, n, ref start, ref xmin, ref ynewlo, reqmin, step,
            konvge, kcount, ref icount, ref numres, ref ifault );

        Console.WriteLine("");
        Console.WriteLine("  Return code IFAULT = " + ifault + "");
        Console.WriteLine("");
        Console.WriteLine("  Estimate of minimizing value X*:");
        Console.WriteLine("");
        for (int i = 0; i < n; i++ )
        {
            Console.WriteLine("  " + xmin[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  F(X*) = " + ynewlo + "");

        Console.WriteLine("");
        Console.WriteLine("  Number of iterations = " + icount + "");
        Console.WriteLine("  Number of restarts =   " + numres + "");

    }


    private static double helical ( double[] x )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HELICAL evaluates the Fletcher-Powell helical valley function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 February 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    R ONeill,
        //    Algorithm AS 47:
        //    Function Minimization Using a Simplex Procedure,
        //    Applied Statistics,
        //    Volume 20, Number 3, 1971, pages 338-345.
        //
        //  Parameters:
        //
        //    Input, double X[3], the argument.
        //
        //    Output, double HELICAL, the value of the function.
        //
    {
        double theta = x[0] switch
        {
            > 0.0 => Math.Atan2(x[1], x[0]) / 2.0 / Math.PI,
            < 0.0 => 0.5 + Math.Atan2(x[1], x[0]) / 2.0 / Math.PI,
            _ => 0.25
        };

        double fx1 = x[2] - 10.0 * theta;
        double fx2 = Math.Sqrt ( x[0] * x[0] + x[1] * x[1] );
        double fx3 = x[2];

        double fx = 100.0 * fx1 * fx1
                    +         fx2 * fx2
                    +         fx3 * fx3;

        return fx;
    }

    private static void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 demonstrates the use of NELMIN on QUARTIC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 February 2008
//
//  Author:
//
//    John Burkardt
//
    {
        int icount = 0;
        int ifault = 0;
        int numres = 0;

        int n = 10;

        double[] start = new double[n];
        double[] step = new double[n];
        double[] xmin = new double[n];

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  Apply NELMIN to the QUARTIC function.");

        for (int i = 0; i < n; i++ )
        {
            start[i] = 1.0;
        }

        double reqmin = 1.0E-08;

        for (int i = 0; i < n; i++ )
        {
            step[i] = 1.0;
        }

        int konvge = 10;
        int kcount = 500;

        Console.WriteLine("");
        Console.WriteLine("  Starting point X:");
        Console.WriteLine("");
        for (int i = 0; i < n; i++ )
        {
            Console.WriteLine("  " + start[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        double ynewlo = quartic ( start );

        Console.WriteLine("");
        Console.WriteLine("  F(X) = " + ynewlo + "");

        Algorithms.nelmin ( quartic, n, ref start, ref xmin, ref ynewlo, reqmin, step,
            konvge, kcount, ref icount, ref numres, ref ifault );

        Console.WriteLine("");
        Console.WriteLine("  Return code IFAULT = " + ifault + "");
        Console.WriteLine("");
        Console.WriteLine("  Estimate of minimizing value X*:");
        Console.WriteLine("");
        for (int i = 0; i < n; i++ )
        {
            Console.WriteLine("  " + xmin[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  F(X*) = " + ynewlo + "");

        Console.WriteLine("");
        Console.WriteLine("  Number of iterations = " + icount + "");
        Console.WriteLine("  Number of restarts =   " + numres + "");

    }

    private static double quartic ( double[] x )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    QUARTIC evaluates a function defined by a sum of fourth powers.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 February 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    R ONeill,
        //    Algorithm AS 47:
        //    Function Minimization Using a Simplex Procedure,
        //    Applied Statistics,
        //    Volume 20, Number 3, 1971, pages 338-345.
        //
        //  Parameters:
        //
        //    Input, double X[10], the argument.
        //
        //    Output, double QUARTIC, the value of the function.
        //
    {
        int i;

        double fx = 0.0;

        for ( i = 0; i < 10; i++ )
        {
            fx += x[i] * x[i] * x[i] * x[i];
        }

        return fx;
    }
       
}