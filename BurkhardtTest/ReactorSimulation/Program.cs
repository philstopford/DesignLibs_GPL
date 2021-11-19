using System;
using Burkardt.Uniform;

namespace ReactorSimulation;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for the reactor shielding simulation.
        //
        //  Discussion:
        //
        //    This is a Monte Carlo simulation, using
        //    uniform random numbers, which investigates the
        //    effectiveness of a shield intended to absorb the
        //    neutrons emitted from a nuclear reactor.
        // 
        //    The reactor is modeled as a point source,
        //    located at (0,0,0).
        //   
        //    A particle emitted from the reactor has a random
        //    initial direction, and an energy selected from
        //    [Emin,Emax] with a 1/Sqrt(E) distribution.
        //   
        //    The shield is modeled as a wall of thickness THICK,
        //    extending from 0 to THICK in the X direction, and
        //    extending forever in the Y and Z directions.
        //   
        //    Based on the particle energy, a distance D is computed
        //    which measures how far the particle could travel through
        //    the shield before colliding.
        //   
        //    Based on the particle direction, the position is updated
        //    by D units.
        //   
        //    If the particle is now to the left of the shield, it is
        //    counted as being REFLECTED.
        //   
        //    If the particle is to the right of the shield, it is 
        //    counted as being ABSORBED.
        //   
        //    If the particle is inside the shield, it has COLLIDED.
        //    A particle that collides is either absorbed (end of story)
        //    or SCATTERED with a new random direction and a new (lower)
        //    energy.
        //   
        //    Every particle is followed from origin to its final fate,
        //    which is reflection, transmission, or absorption.
        //    At the end, a summary is printed, giving the number of
        //    particles with each fate, and the average energy of each
        //    group of particles.
        //   
        //    Increasing NTOT, the number of particles used, will improve the
        //    expected reliability of the results.
        //   
        //    Increasing THICK, the thickness of the shield, should 
        //    result in more absorptions and reflections.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 September 2012
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Kahaner, Moler, Nash.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    David Kahaner, Cleve Moler, Steven Nash,
        //    Numerical Methods and Software,
        //    Prentice Hall, 1989,
        //    ISBN: 0-13-627258-4,
        //    LC: TA345.K34.
        //
        //  Local Parameters:
        //
        //    Local, double AZM, the azimuthal angle of the particle's
        //    direction.
        //
        //    Local, double D, the distance that the particle can
        //    travel through the slab, given its current energy.
        //
        //    Local, double E, the energy of the particle.
        //
        //    Local, double EA, energy absorbed by the slab.
        //
        //    Local, double ER, energy reflected by the slab.
        //
        //    Local, double ET, energy transmitted through the slab.
        //
        //    Local, double MU, the cosine of the angle between the
        //    particle's direction and the X axis.
        //
        //    Local, int NA, number of particles absorbed by the slab.
        //
        //    Local, int NPART, the index of the current particle.
        //
        //    Local, int NR, number of particles reflected by the slab.
        //
        //    Local, int NT, number of particles transmitted by the slab.
        //
        //    Local, int NTOT, the total number of particles to be
        //    emitted from the neutron source.
        //
        //    Local, double SA, standard deviation of absorbed energy.
        //
        //    Local, double SR, standard deviation of reflected energy.
        //
        //    Local, double ST, standard deviation of transmitted energy.
        //
        //    Local, double THICK, the thickness of the slab that is
        //    intended to absorb most of the particles.
        //
        //    Local, double X, Y, Z, the current position of the particle.
        //
    {
        double azm = 0;
        double d = 0;
        double e = 0;
        double ea = 0;
        double er = 0;
        double et = 0;
        double mu = 0;
        int na = 0;
        int nr = 0;
        int nt = 0;
        int ntot = 100000;
        int part = 0;
        double sa = 0;
        int seed = 0;
        double sr = 0;
        double st = 0;
        int test = 0;
        int test_num = 5;
        double thick = 2.0;
        double x = 0;
        double y = 0;
        double z = 0;

        Console.WriteLine("");
        Console.WriteLine("REACTOR_SIMULATION");
        Console.WriteLine("  The reactor shielding simulation.");
        Console.WriteLine("");
        Console.WriteLine("  Shield thickness is THICK = " + thick + "");
        Console.WriteLine("  Number of simulated particles is NTOT = " + ntot + "");
        Console.WriteLine("  Number of tests TEST_NUM = " + test_num + "");

        seed = 123456789;

        for (test = 1; test <= test_num; test++)
        {
            Console.WriteLine("");
            Console.WriteLine("  Test # " + test + "");
            Console.WriteLine("  SEED = " + seed + "");
            //
            //  Initialize.
            //
            ea = 0.0;
            er = 0.0;
            et = 0.0;
            na = 0;
            nr = 0;
            nt = 0;
            sa = 0.0;
            sr = 0.0;
            st = 0.0;
            //
            //  Loop over the particles.
            //
            for (part = 1; part <= ntot; part++)
            {
                //
                //  Generate a new particle.
                //
                source(ref seed, ref e, ref mu, ref azm, ref x, ref y, ref z);

                while (true)
                {
                    //
                    //  Compute the distance that the particle can travel through the slab,
                    //  based on its current energy.
                    //
                    d = dist2c(e, ref seed);
                    //
                    //  Update the particle's position by D units.
                    //
                    update(mu, azm, d, ref x, ref y, ref z);
                    //
                    //  The particle was reflected by the shield, and this path is complete.
                    //
                    if (x < 0.0)
                    {
                        nr += 1;
                        er += e;
                        sr += e * e;
                        break;
                    }
                    //
                    //  The particle was transmitted through the shield, and this path is complete.
                    //

                    if (thick < x)
                    {
                        nt += 1;
                        et += e;
                        st += e * e;
                        break;
                    }
                    //
                    //  The particle collided with the shield, and was absorbed.  This path is done.
                    //
                    if (absorb(ref seed))
                    {
                        na += 1;
                        ea += e;
                        sa += e * e;
                        break;
                    }
                    //
                    //  The particle collided with the shield and was scattered.
                    //  Find the scattering angle and energy, and continue along the new path.
                    //
                    scatter(ref seed, ref e, ref mu, ref azm);
                }
            }

            //
            //  Print the results of the simulation.
            //
            output(na, ea, sa, nr, er, sr, nt, et, st, ntot);
        }

        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("REACTOR_SIMULATION:");
        Console.WriteLine("  Normal end of execution.");

        Console.WriteLine("");
    }

    private static bool absorb(ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ABSORB determines if a colliding particle is absorbed.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 September 2012
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Kahaner, Moler, Nash.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    David Kahaner, Cleve Moler, Steven Nash,
        //    Numerical Methods and Software,
        //    Prentice Hall, 1989,
        //    ISBN: 0-13-627258-4,
        //    LC: TA345.K34.
        //
        //  Parameters:
        //
        //    Input/output, int &SEED, a seed for the random
        //    number generator.
        //
        //    Output, logical ABSORB, is TRUE if the particle is absorbed.
        //
        //  Local parameters:
        //
        //    Local, double PA, the probability of absorption.
        //
    {
        double pa = 0.1;
        double u;
        bool value;

        u = UniformRNG.r8_uniform_01(ref seed);

        if (u <= pa)
        {
            value = true;
        }
        else
        {
            value = false;
        }

        return value;
    }

    private static double cross(double e)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CROSS returns the "cross section" of a particle based on its energy.
        //
        //  Discussion:
        //
        //    The particle's cross section is a measure of its likelihood to collide
        //    with the material of the slab.  This quantity typically depends on both
        //    the particle's energy and the kind of medium through which it is traveling.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 September 2012
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Kahaner, Moler, Nash.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    David Kahaner, Cleve Moler, Steven Nash,
        //    Numerical Methods and Software,
        //    Prentice Hall, 1989,
        //    ISBN: 0-13-627258-4,
        //    LC: TA345.K34.
        //
        //  Parameters:
        //
        //    Input, double E, the energy of the particle.
        //
        //    Output, double CROSS, the cross section.
        //
    {
        double s;
        double value = 0;
        double y;

        s = Math.Abs(Math.Sin(100.0 * (Math.Exp(e) - 1.0))
                     + Math.Sin(18.81 * (Math.Exp(e) - 1.0)));

        y = Math.Max(0.02, s);

        value = 10.0 * Math.Exp(-0.1 / y);

        return value;
    }

    private static double dist2c(double e, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIST2C returns the distance to collision.
        //
        //  Discussion:
        //
        //    Assuming the particle has a given energy, and assuming it is currently
        //    somewhere inside the shield, it is possible to determine a typical distance
        //    which the particle can travel before it collides with the material of
        //    the shield.
        //
        //    The computation of the collision distance is made by estimating a
        //    "cross section" (as though having more energy made the particle "bigger"
        //    and hence more likely to collide) and then randomly selecting a distance
        //    that is logarithmically distributed.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 September 2012
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Kahaner, Moler, Nash.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    David Kahaner, Cleve Moler, Steven Nash,
        //    Numerical Methods and Software,
        //    Prentice Hall, 1989,
        //    ISBN: 0-13-627258-4,
        //    LC: TA345.K34.
        //
        //  Parameters:
        //
        //    Input, double E, the energy of the particle.
        //
        //    Input/output, int &SEED, a seed for the random
        //    number generator.
        //
        //    Output, double DIST2C, the distance the particle can travel
        //    through the slab before colliding.
        //
    {
        double u;
        double value = 0;

        u = UniformRNG.r8_uniform_01(ref seed);

        value = -Math.Log(u) / cross(e);

        return value;
    }

    private static double energy(ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ENERGY assigns an energy to an emitted particle.
        //
        //  Discussion:
        //
        //    The energy E is in the range [EMIN,EMAX], with distribution
        //    const/sqrt(energy).
        //
        //    An inverse function approach is used to compute this.
        //
        //    The energies are measured in MeV.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 September 2012
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Kahaner, Moler, Nash.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    David Kahaner, Cleve Moler, Steven Nash,
        //    Numerical Methods and Software,
        //    Prentice Hall, 1989,
        //    ISBN: 0-13-627258-4,
        //    LC: TA345.K34.
        //
        //  Parameters:
        //
        //    Input/output, int &SEED, a seed for the random
        //    number generator.
        //
        //    Output, double ENERGY, a randomly chosen energy that is
        //    distributed as described above.
        //
        //  Local parameters:
        //
        //    Local, double EMIN, EMAX, the minimum and maximum
        //    energies.
        //
    {
        double c;
        double emax = 2.5;
        double emin = 1.0E-03;
        double u;
        double value = 0;

        u = UniformRNG.r8_uniform_01(ref seed);

        c = 1.0 / (2.0 * (Math.Sqrt(emax) - Math.Sqrt(emin)));

        value = u / (2.0 * c) + Math.Sqrt(emin);
        value *= value;

        return value;
    }

    private static void output(int na, double ea, double sa, int nr, double er, double sr,
            int nt, double et, double st, int ntot)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    OUTPUT prints the results of the reactor shielding simulation.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 September 2012
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Kahaner, Moler, Nash.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    David Kahaner, Cleve Moler, Steven Nash,
        //    Numerical Methods and Software,
        //    Prentice Hall, 1989,
        //    ISBN: 0-13-627258-4,
        //    LC: TA345.K34.
        //
        //  Parameters:
        //
        //    Input, int NA, number of particles absorbed by the slab.
        //
        //    Input, double EA, energy absorbed by the slab.
        //
        //    Input, double SA, the sum of the squares of the 
        //    absorbed energies.
        //
        //    Input, int NR, number of particles reflected by the slab.
        //
        //    Input, double ER, energy reflected by the slab.
        //
        //    Input, double SR, the sum of the squares of the 
        //    reflected energies.
        //
        //    Input, int NT, number of particles transmitted by the slab.
        //
        //    Input, double ET, energy transmitted through the slab.
        //
        //    Input, double ST, the sum of the squares of the 
        //    transmitted energies.
        //
        //    Input, int NTOT, the total number of particles.
        //
    {
        double ea_ave;
        double er_ave;
        double et_ave;
        double etot;
        double pa;
        double pr;
        double pt;
        double ptot;

        Console.WriteLine("");
        Console.WriteLine("  The Reactor Shielding Problem:");
        Console.WriteLine("");
        Console.WriteLine("                           Total                   Average");
        Console.WriteLine("                    #      Energy      " +
                          "Percent     Energy         StDev");
        Console.WriteLine("");

        etot = ea + er + et;

        switch (na)
        {
            case > 0:
                ea_ave = ea / na;
                sa = Math.Sqrt(sa / na - ea_ave * ea_ave);
                break;
            default:
                ea_ave = 0.0;
                break;
        }

        pa = na * 100 / (double)ntot;

        Console.WriteLine("Absorbed   "
                          + "  " + na.ToString().PadLeft(8)
                          + "  " + ea.ToString().PadLeft(14)
                          + "  " + pa.ToString().PadLeft(6)
                          + "  " + ea_ave.ToString().PadLeft(14)
                          + "  " + sa.ToString().PadLeft(14) + "");

        switch (nr)
        {
            case > 0:
                er_ave = er / nr;
                sr = Math.Sqrt(sr / nr - er_ave * er_ave);
                break;
            default:
                er_ave = 0.0;
                break;
        }

        pr = nr * 100 / (double)ntot;

        Console.WriteLine("Reflected  "
                          + "  " + nr.ToString().PadLeft(8)
                          + "  " + er.ToString().PadLeft(14)
                          + "  " + pr.ToString().PadLeft(6)
                          + "  " + er_ave.ToString().PadLeft(14)
                          + "  " + sr.ToString().PadLeft(14) + "");

        switch (nt)
        {
            case > 0:
                et_ave = et / nt;
                st = Math.Sqrt(st / nt - et_ave * et_ave);
                break;
            default:
                et_ave = 0.0;
                break;
        }

        pt = nt * 100 / (double)ntot;

        Console.WriteLine("Transmitted  "
                          + "  " + nt.ToString().PadLeft(8)
                          + "  " + et.ToString().PadLeft(14)
                          + "  " + pt.ToString().PadLeft(6)
                          + "  " + et_ave.ToString().PadLeft(14)
                          + "  " + st.ToString().PadLeft(14) + "");

        ptot = 100.0;

        Console.WriteLine("");
        Console.WriteLine("Total      "
                          + "  " + ntot.ToString().PadLeft(8)
                          + "  " + etot.ToString().PadLeft(14)
                          + "  " + ptot.ToString().PadLeft(6) + "");

    }

    private static void scatter(ref int seed, ref double e, ref double mu, ref double azm)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SCATTER returns the new direction and energy of a particle that is scattered.
        //
        //  Discussion:
        //
        //    The scattering direction is chosen uniformly on the sphere.
        //
        //    The energy of the scattered particle is chosen uniformly in
        //    [ 0.3*E, E ].
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 September 2012
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Kahaner, Moler, Nash.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    David Kahaner, Cleve Moler, Steven Nash,
        //    Numerical Methods and Software,
        //    Prentice Hall, 1989,
        //    ISBN: 0-13-627258-4,
        //    LC: TA345.K34.
        //
        //  Parameters:
        //
        //    Input/output, int &SEED, a seed for the random
        //    number generator.
        //
        //    Input/output, double &E.  On input, the particle energy
        //    before collision.  On output, the particle energy after collision
        //    and scattering.
        //
        //    Output, double &MU, the cosine of the angle between the
        //    particle's direction and the X axis.
        //
        //    Output, double &AZM, the azimuthal angle of the particle's
        //    direction.
        //
    {
        double pi = 3.141592653589793;
        double u;

        u = UniformRNG.r8_uniform_01(ref seed);
        mu = -1.0 + 2.0 * u;

        u = UniformRNG.r8_uniform_01(ref seed);
        azm = u * 2.0 * pi;

        u = UniformRNG.r8_uniform_01(ref seed);
        e = (u * 0.7 + 0.3) * e;

    }

    private static void source(ref int seed, ref double e, ref double mu, ref double azm, ref double x,
            ref double y, ref double z)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SOURCE generates a new particle from the neutron source.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 September 2012
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Kahaner, Moler, Nash.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    David Kahaner, Cleve Moler, Steven Nash,
        //    Numerical Methods and Software,
        //    Prentice Hall, 1989,
        //    ISBN: 0-13-627258-4,
        //    LC: TA345.K34.
        //
        //  Parameters:
        //
        //    Input/output, int &SEED, a seed for the random
        //    number generator.
        //
        //    Output, double &E, the initial energy of the particle.
        //
        //    Output, double &MU, the cosine of the angle between the
        //    particle's direction and the X axis.
        //
        //    Output, double &AZM, the azimuthal angle of the particle's
        //    direction.
        //
        //    Output, double &X, &Y, &Z, the initial coordinates of the particle.
        //
    {
        double pi = 3.141592653589793;
        double u;

        u = UniformRNG.r8_uniform_01(ref seed);
        mu = u;

        u = UniformRNG.r8_uniform_01(ref seed);
        azm = u * 2.0 * pi;

        x = 0.0;
        y = 0.0;
        z = 0.0;

        e = energy(ref seed);
    }

    private static void update(double mu, double azm, double d, ref double x, ref double y, ref double z)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UPDATE determines the position of the particle after it has traveled D units.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 September 2012
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Kahaner, Moler, Nash.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    David Kahaner, Cleve Moler, Steven Nash,
        //    Numerical Methods and Software,
        //    Prentice Hall, 1989,
        //    ISBN: 0-13-627258-4,
        //    LC: TA345.K34.
        //
        //  Parameters:
        //
        //    Input, double MU, the cosine of the angle between the
        //    particle's direction and the X axis.
        //
        //    Input, double AZM, the azimuthal angle of the particle's
        //    direction.
        //
        //    Input, double D, the distance the particle traveled.
        //
        //    Input/output, double &X, &Y, &Z.  On input, the previous
        //    coordinates of the particle.  On output, the updated coordinates of the
        //    particle.
        //
    {
        double s;

        s = Math.Sqrt(1.0 - mu * mu);

        x += d * mu;
        y += d * s * Math.Cos(azm);
        z += d * s * Math.Sin(azm);

    }
}