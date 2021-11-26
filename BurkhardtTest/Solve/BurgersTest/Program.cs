using System;
using Burkardt.NavierStokesNS;
using Burkardt.Types;

namespace BurgersTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for BURGERS_SOLUTION_TEST.
        //
        //  Discussion:
        //
        //    BURGERS_SOLUTION_TEST tests the BURGERS_SOLUTION library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 September 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("BURGERS_SOLUTION_TEST");
            
        Console.WriteLine("  Test the BURGERS_SOLUTION library.");

        burgers_viscous_time_exact1_test01();
        burgers_viscous_time_exact1_test02();

        burgers_viscous_time_exact2_test01();
        burgers_viscous_time_exact2_test02();
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("BURGERS_SOLUTION_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void burgers_viscous_time_exact1_test01()

        //****************************************************************************80
        //.
        //  Purpose:
        //
        //    BURGERS_VISCOUS_TIME_EXACT1_TEST01 tests sets up a small test case.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 November 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const string filename = "burgers_solution_test01.txt";
        const int vtn = 11;
        const int vxn = 11;

        double nu = 0.01 / Math.PI;

        Console.WriteLine("");
        Console.WriteLine("BURGERS_VISCOUS_TIME_EXACT1_TEST01");
        Console.WriteLine("  BURGERS_VISCOUS_TIME_EXACT1 evaluates solution #1");
        Console.WriteLine("  to the Burgers equation.");
        Console.WriteLine("");
        Console.WriteLine("  Viscosity NU = " + nu + "");
        Console.WriteLine("  NX = " + vxn + "");
        Console.WriteLine("  NT = " + vtn + "");

        const double xlo = -1.0;
        const double xhi = +1.0;
        double[] vx = typeMethods.r8vec_even_new(vxn, xlo, xhi);
        typeMethods.r8vec_print(vxn, vx, "  X grid points:");

        const double tlo = 0.0;
        const double thi = 3.0 / Math.PI;
        double[] vt = typeMethods.r8vec_even_new(vtn, tlo, thi);
        typeMethods.r8vec_print(vtn, vt, "  T grid points:");

        double[] vu = Burgers.burgers_viscous_time_exact1(nu, vxn, vx, vtn, vt);

        typeMethods.r8mat_print(vxn, vtn, vu, "  U(X,T) at grid points:");

        typeMethods.r8mat_write(filename, vxn, vtn, vu);

        Console.WriteLine("");
        Console.WriteLine("  Data written to file \"" + filename + "\".");

    }

    private static void burgers_viscous_time_exact1_test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BURGERS_VISCOUS_TIME_EXACT1_TEST02 tests sets up a finer test case.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 November 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const string filename = "burgers_solution_test02.txt";
        const int vtn = 41;
        const int vxn = 41;

        const double nu = 0.01 / Math.PI;

        Console.WriteLine("");
        Console.WriteLine("BURGERS_VISCOUS_TIME_EXACT1_TEST02");
        Console.WriteLine("  BURGERS_VISCOUS_TIME_EXACT1 evaluates solution #1");
        Console.WriteLine("  to the Burgers equation.");
        Console.WriteLine("");
        Console.WriteLine("  Viscosity NU = " + nu + "");
        Console.WriteLine("  NX = " + vxn + "");
        Console.WriteLine("  NT = " + vtn + "");

        const double xlo = -1.0;
        const double xhi = +1.0;
        double[] vx = typeMethods.r8vec_even_new(vxn, xlo, xhi);
        typeMethods.r8vec_print(vxn, vx, "  X grid points:");

        const double tlo = 0.0;
        const double thi = 3.0 / Math.PI;
        double[] vt = typeMethods.r8vec_even_new(vtn, tlo, thi);
        typeMethods.r8vec_print(vtn, vt, "  T grid points:");

        double[] vu = Burgers.burgers_viscous_time_exact1(nu, vxn, vx, vtn, vt);

        typeMethods.r8mat_write(filename, vxn, vtn, vu);

        Console.WriteLine("");
        Console.WriteLine("  Data written to file \"" + filename + "\".");

    }

    private static void burgers_viscous_time_exact2_test01()

        //****************************************************************************80
        //.
        //  Purpose:
        //
        //    BURGERS_VISCOUS_TIME_EXACT2_TEST01 tests sets up a small test case.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 September 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const string filename = "burgers_solution_test03.txt";
        const int vtn = 11;
        const int vxn = 11;

        const double nu = 0.5;

        Console.WriteLine("");
        Console.WriteLine("BURGERS_VISCOUS_TIME_EXACT2_TEST01");
        Console.WriteLine("  BURGERS_VISCOUS_TIME_EXACT2 evaluates solution #2");
        Console.WriteLine("  to the Burgers equation.");
        Console.WriteLine("");
        Console.WriteLine("  Viscosity NU = " + nu + "");
        Console.WriteLine("  NX = " + vxn + "");
        Console.WriteLine("  NT = " + vtn + "");

        const double xlo = 0.0;
        const double xhi = 2.0 * Math.PI;
        double[] vx = typeMethods.r8vec_even_new(vxn, xlo, xhi);
        typeMethods.r8vec_print(vxn, vx, "  X grid points:");

        const double tlo = 0.0;
        const double thi = 1.0;
        double[] vt = typeMethods.r8vec_even_new(vtn, tlo, thi);
        typeMethods.r8vec_print(vtn, vt, "  T grid points:");

        double[] vu = Burgers.burgers_viscous_time_exact2(nu, vxn, vx, vtn, vt);

        typeMethods.r8mat_print(vxn, vtn, vu, "  U(X,T) at grid points:");

        typeMethods.r8mat_write(filename, vxn, vtn, vu);

        Console.WriteLine("");
        Console.WriteLine("  Data written to file \"" + filename + "\".");

    }

    private static void burgers_viscous_time_exact2_test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BURGERS_VISCOUS_TIME_EXACT2_TEST02 tests sets up a finer test case.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 September 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const string filename = "burgers_solution_test04.txt";
        const int vtn = 41;
        const int vxn = 41;

        const double nu = 0.5;

        Console.WriteLine("");
        Console.WriteLine("BURGERS_VISCOUS_TIME_EXACT2_TEST02");
        Console.WriteLine("  BURGERS_VISCOUS_TIME_EXACT2 evaluates solution #2");
        Console.WriteLine("  to the Burgers equation.");
        Console.WriteLine("");
        Console.WriteLine("  Viscosity NU = " + nu + "");
        Console.WriteLine("  NX = " + vxn + "");
        Console.WriteLine("  NT = " + vtn + "");

        const double xlo = 0.0;
        const double xhi = 2.0 * Math.PI;
        double[] vx = typeMethods.r8vec_even_new(vxn, xlo, xhi);
        typeMethods.r8vec_print(vxn, vx, "  X grid points:");

        const double tlo = 0.0;
        const double thi = 1.0;
        double[] vt = typeMethods.r8vec_even_new(vtn, tlo, thi);
        typeMethods.r8vec_print(vtn, vt, "  T grid points:");

        double[] vu = Burgers.burgers_viscous_time_exact2(nu, vxn, vx, vtn, vt);

        typeMethods.r8mat_write(filename, vxn, vtn, vu);

        Console.WriteLine("");
        Console.WriteLine("  Data written to file \"" + filename + "\".");

    }
}