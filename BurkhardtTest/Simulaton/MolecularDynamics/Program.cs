using System;
using System.Linq;
using Burkardt.Uniform;

namespace MolecularDynamics;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for MD.
        //
        //  Discussion:
        //
        //    MD implements a simple molecular dynamics simulation.
        //
        //    The velocity Verlet time integration scheme is used. 
        //
        //    The particles interact with a central pair potential.
        //
        //    This program is based on a FORTRAN90 program by Bill Magro.
        //
        //  Usage:
        //
        //    md nd np step_num dt
        //
        //    where
        //
        //    * nd is the spatial dimension (2 or 3);
        //    * np is the number of particles (500, for instance);
        //    * step_num is the number of time steps (500, for instance).
        //    * dt is the time step (0.1 for instance)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 December 2014
        //
        //  Author:
        //
        //    John Burkardt.
        //
    {
        double[] acc;
        DateTime ctime;
        double dt;
        double e0 = 0;
        double[] force;
        double kinetic = 0;
        double mass = 1.0;
        int nd;
        int np;
        double[] pos;
        double potential = 0;
        int step;
        int step_num;
        int step_print;
        int step_print_index;
        int step_print_num;
        double[] vel;

        Console.WriteLine("");
        Console.WriteLine("MD");
        Console.WriteLine("  A molecular dynamics program.");
        //
        //  Get the spatial dimension.
        //
        try
        {
            nd = Convert.ToInt32(args[0]);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("  Enter ND, the spatial dimension (2 or 3).");
            nd = Convert.ToInt32(Console.ReadLine());
        }

        //
        //  Get the number of particles.
        //
        try
        {
            np = Convert.ToInt32(args[1]);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("  Enter NP, the number of particles (500, for instance).");
            np = Convert.ToInt32(Console.ReadLine());
        }

        //
        //  Get the number of time steps.
        //
        try
        {
            step_num = Convert.ToInt32(args[2]);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("  Enter STEP_NUM, the number of time steps (500 or 1000, for instance).");
            step_num = Convert.ToInt32(Console.ReadLine());
        }

        //
        //  Get the time step.
        //
        try
        {
            dt = Convert.ToDouble(args[3]);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("  Enter DT, the time step size (0.1, for instance).");
            dt = Convert.ToDouble(Console.ReadLine());
        }

        //
        //  Report.
        //
        Console.WriteLine("");
        Console.WriteLine("  ND, the spatial dimension, is " + nd + "");
        Console.WriteLine("  NP, the number of particles in the simulation is " + np + "");
        Console.WriteLine("  STEP_NUM, the number of time steps, is " + step_num + "");
        Console.WriteLine("  DT, the size of each time step, is " + dt + "");
        //
        //  Allocate memory.
        //
        acc = new double[nd * np];
        force = new double[nd * np];
        pos = new double[nd * np];
        vel = new double[nd * np];
        //
        //  This is the main time stepping loop:
        //    Compute forces and energies,
        //    Update positions, velocities, accelerations.
        //
        Console.WriteLine("");
        Console.WriteLine("  At each step, we report the potential and kinetic energies.");
        Console.WriteLine("  The sum of these energies should be a constant.");
        Console.WriteLine("  As an accuracy check, we also print the relative error");
        Console.WriteLine("  in the total energy.");
        Console.WriteLine("");
        Console.WriteLine("      Step      Potential       Kinetic        (P+K-E0)/E0");
        Console.WriteLine("                Energy P        Energy K       Relative Energy Error");
        Console.WriteLine("");

        step_print = 0;
        step_print_index = 0;
        step_print_num = 10;

        ctime = DateTime.Now;

        for (step = 0; step <= step_num; step++)
        {
            switch (step)
            {
                case 0:
                    initialize(np, nd, ref pos, ref vel, ref acc);
                    break;
                default:
                    update(np, nd, ref pos, ref vel, force, ref acc, mass, dt);
                    break;
            }

            compute(np, nd, pos, vel, mass, ref force, ref potential, ref kinetic);

            e0 = step switch
            {
                0 => potential + kinetic,
                _ => e0
            };

            if (step == step_print)
            {
                Console.WriteLine("  " + step.ToString().PadLeft(8)
                                       + "  " + potential.ToString().PadLeft(14)
                                       + "  " + kinetic.ToString().PadLeft(14)
                                       + "  " + ((potential + kinetic - e0) / e0).ToString().PadLeft(14) + "");
                step_print_index += 1;
                step_print = step_print_index * step_num / step_print_num;
            }

        }

        //
        //  Report timing.
        //
        Console.WriteLine("");
        Console.WriteLine("  Elapsed cpu time " + (DateTime.Now - ctime).TotalSeconds + " seconds.");

        Console.WriteLine("");
        Console.WriteLine("MD");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void compute(int np, int nd, double[] pos, double[] vel, double mass,
            ref double[] f, ref double pot, ref double kin)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COMPUTE computes the forces and energies.
        //
        //  Discussion:
        //
        //    The computation of forces and energies is fully parallel.
        //
        //    The potential function V(X) is a harmonic well which smoothly
        //    saturates to a maximum value at PI/2:
        //
        //      v(x) = ( Math.Sin ( min ( x, PI2 ) ) )**2
        //
        //    The derivative of the potential is:
        //
        //      dv(x) = 2.0 * Math.Sin ( min ( x, PI2 ) ) * cos ( min ( x, PI2 ) )
        //            = Math.Sin ( 2.0 * min ( x, PI2 ) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 July 2008
        //
        //  Author:
        //
        //    John Burkardt.
        //
        //  Parameters:
        //
        //    Input, int NP, the number of particles.
        //
        //    Input, int ND, the number of spatial dimensions.
        //
        //    Input, double POS[ND*NP], the position of each particle.
        //
        //    Input, double VEL[ND*NP], the velocity of each particle.
        //
        //    Input, double MASS, the mass of each particle.
        //
        //    Output, double F[ND*NP], the forces.
        //
        //    Output, double &POT, the total potential energy.
        //
        //    Output, double &KIN, the total kinetic energy.
        //
    {
        double d;
        double d2;
        int i;
        int j;
        int k;
        double PI2 = Math.PI / 2.0;
        double[] rij = new double[3];

        pot = 0.0;
        kin = 0.0;

        for (k = 0; k < np; k++)
        {
            //
            //  Compute the potential energy and forces.
            //
            for (i = 0; i < nd; i++)
            {
                f[i + k * nd] = 0.0;
            }

            for (j = 0; j < np; j++)
            {
                if (k != j)
                {
                    d = dist(nd, pos.Skip(+k * nd).ToArray(), pos.Skip(+j * nd).ToArray(), ref rij);
                    //
                    //  Attribute half of the potential energy to particle J.
                    //
                    if (d < PI2)
                    {
                        d2 = d;
                    }
                    else
                    {
                        d2 = PI2;
                    }

                    pot += 0.5 * Math.Pow(Math.Sin(d2), 2);

                    for (i = 0; i < nd; i++)
                    {
                        f[i + k * nd] -= rij[i] * Math.Sin(2.0 * d2) / d;
                    }
                }
            }

            //
            //  Compute the kinetic energy.
            //
            for (i = 0; i < nd; i++)
            {
                kin += vel[i + k * nd] * vel[i + k * nd];
            }
        }

        kin = kin * 0.5 * mass;
    }

    private static double dist(int nd, double[] r1, double[] r2, ref double[] dr)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIST computes the displacement (and its norm) between two particles.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 November 2007
        //
        //  Author:
        //
        //    John Burkardt.
        //
        //  Parameters:
        //
        //    Input, int ND, the number of spatial dimensions.
        //
        //    Input, double R1[ND], R2[ND], the positions.
        //
        //    Output, double DR[ND], the displacement vector.
        //
        //    Output, double D, the Euclidean norm of the displacement.
        //
    {
        double d;
        int i;

        d = 0.0;
        for (i = 0; i < nd; i++)
        {
            dr[i] = r1[i] - r2[i];
            d += dr[i] * dr[i];
        }

        d = Math.Sqrt(d);

        return d;
    }

    private static void initialize(int np, int nd, ref double[] pos, ref double[] vel, ref double[] acc)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INITIALIZE initializes the positions, velocities, and accelerations.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    26 December 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NP, the number of particles.
        //
        //    Input, int ND, the number of spatial dimensions.
        //
        //    Output, double POS[ND*NP], the positions.
        //
        //    Output, double VEL[ND*NP], the velocities.
        //
        //    Output, double ACC[ND*NP], the accelerations.
        //
    {
        int i;
        int j;
        int seed;
        //
        //  Set the positions.
        //
        seed = 123456789;
        UniformRNG.r8mat_uniform_ab(nd, np, 0.0, 10.0, ref seed, ref pos);
        //
        //  Set the velocities.
        //
        for (j = 0; j < np; j++)
        {
            for (i = 0; i < nd; i++)
            {
                vel[i + j * nd] = 0.0;
            }
        }

        //
        //  Set the accelerations.
        //
        for (j = 0; j < np; j++)
        {
            for (i = 0; i < nd; i++)
            {
                acc[i + j * nd] = 0.0;
            }
        }
    }

    private static void update(int np, int nd, ref double[] pos, ref double[] vel, double[] f,
            ref double[] acc, double mass, double dt)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UPDATE updates positions, velocities and accelerations.
        //
        //  Discussion:
        //
        //    The time integration is fully parallel.
        //
        //    A velocity Verlet algorithm is used for the updating.
        //
        //    x(t+dt) = x(t) + v(t) * dt + 0.5 * a(t) * dt * dt
        //    v(t+dt) = v(t) + 0.5 * ( a(t) + a(t+dt) ) * dt
        //    a(t+dt) = f(t) / m
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 November 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NP, the number of particles.
        //
        //    Input, int ND, the number of spatial dimensions.
        //
        //    Input/output, double POS[ND*NP], the positions.
        //
        //    Input/output, double VEL[ND*NP], the velocities.
        //
        //    Input, double F[ND*NP], the forces.
        //
        //    Input/output, double ACC[ND*NP], the accelerations.
        //
        //    Input, double MASS, the mass of each particle.
        //
        //    Input, double DT, the time step.
        //
    {
        int i;
        int j;
        double rmass;

        rmass = 1.0 / mass;

        for (j = 0; j < np; j++)
        {
            for (i = 0; i < nd; i++)
            {
                pos[i + j * nd] = pos[i + j * nd] + vel[i + j * nd] * dt + 0.5 * acc[i + j * nd] * dt * dt;
                vel[i + j * nd] += 0.5 * dt * (f[i + j * nd] * rmass + acc[i + j * nd]);
                acc[i + j * nd] = f[i + j * nd] * rmass;
            }
        }

    }
}