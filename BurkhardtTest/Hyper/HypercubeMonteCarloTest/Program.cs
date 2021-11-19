using System;
using Burkardt.MonomialNS;
using Burkardt.Types;

namespace HyperCubeMonteCarloTest;

using MonteCarlo = Burkardt.HyperGeometry.Hypercube.MonteCarlo;

internal static class Program
{
    private static void Main()
    {
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for HYPERCUBE_MONTE_CARLO_TEST.
        //
        //  Discussion:
        //
        //    HYPERCUBE_MONTE_CARLO_TEST tests the HYPERCUBE_MONTE_CARLO library.
        //    
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            Console.WriteLine();
            Console.WriteLine("HYPERCUBE_MONTE_CARLO_TEST");
            Console.WriteLine("  Test the HYPERCUBE_MONTE_CARLO library.");

            test01();
            test02();

            Console.WriteLine();
            Console.WriteLine("HYPERCUBE_MONTE_CARLO_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine();
        }
    }


    private static void test01 ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 estimates integrals over the unit hypercube in 3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] e = new int [3];
        int[] e_test = {
            0, 0, 0, 
            1, 0, 0, 
            0, 1, 0, 
            0, 0, 1, 
            2, 0, 0, 
            1, 1, 0, 
            1, 0, 1, 
            0, 2, 0, 
            0, 1, 1, 
            0, 0, 2 };
        int m = 3;

        Console.WriteLine();
        Console.WriteLine("TEST01");
        Console.WriteLine("  Use HYPERCUBE01_SAMPLE to estimate integrals");
        Console.WriteLine("  over the interior of the unit hypercube in 3D.");

        int seed = 123456789;

        Console.WriteLine();
        string cout = "         N";
        cout += "        1";
        cout += "               X";
        cout += "               Y ";
        cout += "              Z";
        cout += "               X^2";
        cout += "              XY";
        cout += "             XZ";
        cout += "              Y^2";
        cout += "             YZ";
        cout += "               Z^2";
        Console.WriteLine(cout);
        Console.WriteLine();

        int n = 1;

        while ( n <= 65536 )
        {
            double[] x = MonteCarlo.hypercube01_sample ( m, n, ref seed );
            cout = "  " + n.ToString().PadLeft(8);
            for (int j = 0; j < 10; j++ )
            {
                for (int i = 0; i < m; i++ )
                {
                    e[i] = e_test[i+j*m];
                }

                double[] value = Monomial.monomial_value ( m, n, e, x );

                double result = MonteCarlo.hypercube01_volume ( m ) * typeMethods.r8vec_sum ( n, value ) / n;
                cout += "  " +result.ToString().PadLeft(14);
            }
            Console.WriteLine(cout);

            n = 2 * n;
        }

        Console.WriteLine();
        cout = "     Exact";
        for (int j = 0; j < 10; j++ )
        {
            for (int i = 0; i < m; i++ )
            {
                e[i] = e_test[i+j*m];
            }
            double exact = MonteCarlo.hypercube01_monomial_integral ( m, e );
            cout += "  " +exact.ToString().PadLeft(14);
        }
        Console.WriteLine(cout);
    }


    private static void test02 ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 estimates integrals over the unit hypercube in 6D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] e = new int[6];
        int[] e_test = {
            0, 0, 0, 0, 0, 0, 
            1, 0, 0, 0, 0, 0, 
            0, 2, 0, 0, 0, 0, 
            0, 2, 2, 0, 0, 0, 
            0, 0, 0, 4, 0, 0, 
            2, 0, 0, 0, 2, 2, 
            0, 0, 0, 0, 0, 6 };
        int m = 6;

        Console.WriteLine();
        Console.WriteLine("TEST02");
        Console.WriteLine("  Use HYPERCUBE01_SAMPLE to estimate integrals");
        Console.WriteLine("  over the interior of the unit hypercube in 6D.");

        int seed = 123456789;

        Console.WriteLine();
        string cout = "         N";
        cout += "        1      ";
        cout += "        U      ";
        cout += "         V^2   ";
        cout += "         V^2W^2";
        cout += "         X^4   ";
        cout += "         Y^2Z^2";
        cout += "         Z^6";
        Console.WriteLine(cout);
        Console.WriteLine();

        int n = 1;

        while ( n <= 65536 )
        {
            double[] x = MonteCarlo.hypercube01_sample ( m, n, ref seed );
            cout = "  " +n.ToString().PadLeft(8);
            for (int j = 0; j < 7; j++ )
            {
                for (int i = 0; i < m; i++ )
                {
                    e[i] = e_test[i+j*m];
                }

                double[] value = Monomial.monomial_value ( m, n, e, x );

                double result = MonteCarlo.hypercube01_volume ( m ) * typeMethods.r8vec_sum ( n, value ) / n;
                cout += "  " + result.ToString().PadLeft(14);
            }
            Console.WriteLine(cout);

            n = 2 * n;
        }

        Console.WriteLine();
        cout = "     Exact";
        for (int j = 0; j < 7; j++ )
        {
            for (int i = 0; i < m; i++ )
            {
                e[i] = e_test[i+j*m];
            }
            double exact = MonteCarlo.hypercube01_monomial_integral ( m, e );
            cout += "  " +exact.ToString().PadLeft(14);
        }
        Console.WriteLine(cout);
    }
        
}