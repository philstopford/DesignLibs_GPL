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
        string filename = "burgers_solution_test01.txt";
        double nu;
        double r8_pi = 3.141592653589793;
        double thi;
        double tlo;
        double[] vu;
        double[] vt;
        int vtn = 11;
        double[] vx;
        int vxn = 11;
        double xhi;
        double xlo;

        nu = 0.01 / r8_pi;

        Console.WriteLine("");
        Console.WriteLine("BURGERS_VISCOUS_TIME_EXACT1_TEST01");
        Console.WriteLine("  BURGERS_VISCOUS_TIME_EXACT1 evaluates solution #1");
        Console.WriteLine("  to the Burgers equation.");
        Console.WriteLine("");
        Console.WriteLine("  Viscosity NU = " + nu + "");
        Console.WriteLine("  NX = " + vxn + "");
        Console.WriteLine("  NT = " + vtn + "");

        xlo = -1.0;
        xhi = +1.0;
        vx = typeMethods.r8vec_even_new(vxn, xlo, xhi);
        typeMethods.r8vec_print(vxn, vx, "  X grid points:");

        tlo = 0.0;
        thi = 3.0 / r8_pi;
        vt = typeMethods.r8vec_even_new(vtn, tlo, thi);
        typeMethods.r8vec_print(vtn, vt, "  T grid points:");

        vu = Burgers.burgers_viscous_time_exact1(nu, vxn, vx, vtn, vt);

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
        string filename = "burgers_solution_test02.txt";
        double nu;
        double r8_pi = 3.141592653589793;
        double thi;
        double tlo;
        double[] vu;
        double[] vt;
        int vtn = 41;
        double[] vx;
        int vxn = 41;
        double xhi;
        double xlo;

        nu = 0.01 / r8_pi;

        Console.WriteLine("");
        Console.WriteLine("BURGERS_VISCOUS_TIME_EXACT1_TEST02");
        Console.WriteLine("  BURGERS_VISCOUS_TIME_EXACT1 evaluates solution #1");
        Console.WriteLine("  to the Burgers equation.");
        Console.WriteLine("");
        Console.WriteLine("  Viscosity NU = " + nu + "");
        Console.WriteLine("  NX = " + vxn + "");
        Console.WriteLine("  NT = " + vtn + "");

        xlo = -1.0;
        xhi = +1.0;
        vx = typeMethods.r8vec_even_new(vxn, xlo, xhi);
        typeMethods.r8vec_print(vxn, vx, "  X grid points:");

        tlo = 0.0;
        thi = 3.0 / r8_pi;
        vt = typeMethods.r8vec_even_new(vtn, tlo, thi);
        typeMethods.r8vec_print(vtn, vt, "  T grid points:");

        vu = Burgers.burgers_viscous_time_exact1(nu, vxn, vx, vtn, vt);

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
        string filename = "burgers_solution_test03.txt";
        double nu;
        double r8_pi = 3.141592653589793;
        double thi;
        double tlo;
        double[] vu;
        double[] vt;
        int vtn = 11;
        double[] vx;
        int vxn = 11;
        double xhi;
        double xlo;

        nu = 0.5;

        Console.WriteLine("");
        Console.WriteLine("BURGERS_VISCOUS_TIME_EXACT2_TEST01");
        Console.WriteLine("  BURGERS_VISCOUS_TIME_EXACT2 evaluates solution #2");
        Console.WriteLine("  to the Burgers equation.");
        Console.WriteLine("");
        Console.WriteLine("  Viscosity NU = " + nu + "");
        Console.WriteLine("  NX = " + vxn + "");
        Console.WriteLine("  NT = " + vtn + "");

        xlo = 0.0;
        xhi = 2.0 * r8_pi;
        vx = typeMethods.r8vec_even_new(vxn, xlo, xhi);
        typeMethods.r8vec_print(vxn, vx, "  X grid points:");

        tlo = 0.0;
        thi = 1.0;
        vt = typeMethods.r8vec_even_new(vtn, tlo, thi);
        typeMethods.r8vec_print(vtn, vt, "  T grid points:");

        vu = Burgers.burgers_viscous_time_exact2(nu, vxn, vx, vtn, vt);

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
        string filename = "burgers_solution_test04.txt";
        double nu;
        double r8_pi = 3.141592653589793;
        double thi;
        double tlo;
        double[] vu;
        double[] vt;
        int vtn = 41;
        double[] vx;
        int vxn = 41;
        double xhi;
        double xlo;

        nu = 0.5;

        Console.WriteLine("");
        Console.WriteLine("BURGERS_VISCOUS_TIME_EXACT2_TEST02");
        Console.WriteLine("  BURGERS_VISCOUS_TIME_EXACT2 evaluates solution #2");
        Console.WriteLine("  to the Burgers equation.");
        Console.WriteLine("");
        Console.WriteLine("  Viscosity NU = " + nu + "");
        Console.WriteLine("  NX = " + vxn + "");
        Console.WriteLine("  NT = " + vtn + "");

        xlo = 0.0;
        xhi = 2.0 * r8_pi;
        vx = typeMethods.r8vec_even_new(vxn, xlo, xhi);
        typeMethods.r8vec_print(vxn, vx, "  X grid points:");

        tlo = 0.0;
        thi = 1.0;
        vt = typeMethods.r8vec_even_new(vtn, tlo, thi);
        typeMethods.r8vec_print(vtn, vt, "  T grid points:");

        vu = Burgers.burgers_viscous_time_exact2(nu, vxn, vx, vtn, vt);

        typeMethods.r8mat_write(filename, vxn, vtn, vu);

        Console.WriteLine("");
        Console.WriteLine("  Data written to file \"" + filename + "\".");

    }
}