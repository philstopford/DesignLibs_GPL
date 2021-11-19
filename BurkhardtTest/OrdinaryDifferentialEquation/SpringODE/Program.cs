using System;

namespace SpringODE;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SPRING_ODE.
        //
        //  Discussion:
        //
        //    This is a simple example of how to plot when you don't have a plotter.
        //    This is a particular kind of "ASCII graphics", or "typewriter graphics"
        //    or "lineprinter graphics", and shows you how valuable an illustration 
        //    can be, even when it's as crude as this example.
        //
        //    Hooke's law for a spring observes that the restoring force is
        //    proportional to the displacement: F = - k x
        //
        //    Newton's law relates the force to acceleration: F = m a
        //
        //    Putting these together, we have
        //
        //      m * d^2 x/dt^2 = - k * x
        //
        //    We can add a damping force with coefficient c:
        //
        //      m * d^2 x/dt^2 = - k * x - c * dx/dt
        //
        //    If we write this as a pair of first order equations for (x,v), we have
        //
        //          dx/dt = v
        //      m * dv/dt = - k * x - c * v
        //
        //    and now we can approximate these values for small time steps.
        //
        //    Note that the plotting assumes that the value of X will always be
        //    between -1 and +1.  If the initial condition uses V = 0, and X starts
        //    between -1 and +1, then this will be OK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 May 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    None
        //
    {
        float c;
        float dt;
        int i;
        int j;
        float k;
        float m;
        int n;
        int p = 0;
        float t_final;
        float v;
        float v_old;
        float x;
        float x_old;
        char[] z = new char[21];

        Console.WriteLine("");
        Console.WriteLine("SPRING_ODE");
        Console.WriteLine("  Approximate the solution of a spring equation.");
        Console.WriteLine("  Display the solution with line printer graphics.");
        Console.WriteLine("");
        //
        //  Data
        //
        m = 1.0f;
        k = 1.0f;
        c = 0.3f;
        t_final = 20.0f;
        n = 100;
        dt = t_final / n;
        //
        //  Initial conditions.
        //
        x = 1.0f;
        v = 0.0f;
        //
        //  Compute the approximate solution at equally spaced times.
        //
        for (i = 0; i <= n; i++)
        {
            x_old = x;
            v_old = v;

            //  t = ( float ) ( i ) * t_final / ( float ) ( n );
            x = x_old + dt * v_old;
            v = v_old + dt / m * (-k * x_old - c * v_old);

            p = p switch
            {
                < 0 => 0,
                > 20 => 20,
                //
                //  Approximate the position of X in [-1,+1] to within 1/10.
                //
                _ => (int) (10 * (1.0 + x))
            };

            //
            //  Fill in the next line of the plot, placing 'x' in the p position.
            //
            for (j = 0; j <= 20; j++)
            {
                z[j] = (i % 10) switch
                {
                    0 => '-',
                    _ => ' '
                };
            }

            z[0] = '|';
            z[5] = '.';
            z[10] = '+';
            z[15] = '.';
            z[20] = '|';

            z[p] = 'x';
            string cout = "";
            for (j = 0; j <= 20; j++)
            {
                cout += z[j];
            }

            Console.WriteLine(cout);
        }

        Console.WriteLine("");
        Console.WriteLine("SPRING_ODE:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }
}