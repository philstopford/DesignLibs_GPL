﻿using System;
using Burkardt.PlotNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Ising2DTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for ISING_2D_SIMULATION.
            //
            //  Usage:
            //
            //    ising_2d_simulation  m  n  iterations  thresh  seed
            //
            //    * M, N, the number of rows and columns.
            //    * ITERATIONS, the number of iterations.
            //    * THRESH, the threshhold.
            //    * SEED, a seed for the random number generator.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    30 June 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int[] c1;
            int i;
            int iterations;
            int m;
            int n;
            string plot_filename = "ising_2d_final.txt";
            string png_filename = "ising_2d_final.png";
            double[] prob = {0.98, 0.85, 0.50, 0.15, 0.02};
            int seed;
            double thresh;
            string title = "Final Configuration";

            Console.WriteLine("");
            Console.WriteLine("ISING_2D_SIMULATION");
            Console.WriteLine("  C++ version");
            Console.WriteLine("  Monte Carlo simulation of a 2D Ising model.");
            //
            //  Get input.
            //
            try
            {
                m = Convert.ToInt32(args[0]);
            }
            catch
            {
                m = 10;
            }

            try
            {
                n = Convert.ToInt32(args[1]);
            }
            catch
            {
                n = 10;
            }

            try
            {
                iterations = Convert.ToInt32(args[2]);
            }
            catch
            {
                iterations = 15;
            }

            try
            {
                thresh = Convert.ToDouble(args[3]);
            }
            catch
            {
                thresh = 0.50;
            }

            try
            {
                seed = Convert.ToInt32(args[4]);
            }
            catch
            {
                seed = 123456789;
            }

            Console.WriteLine("");
            Console.WriteLine("  The number of rows is M = " + m + "");
            Console.WriteLine("  The number of columns is N = " + n + "");
            Console.WriteLine("  The number of iterations taken is ITERATIONS = " + iterations + "");
            Console.WriteLine("  The threshhold THRESH = " + thresh + "");
            Console.WriteLine("  The seed SEED = " + seed + "");
            Console.WriteLine("");
            Console.WriteLine("  The transition probability table, based on the number of");
            Console.WriteLine("  neighbors with the same spin.");
            Console.WriteLine("");
            Console.WriteLine("      1         2         3         4         5");
            Console.WriteLine("");
            string cout = "";
            for (i = 0; i < 5; i++)
            {
                cout += prob[i].ToString().PadLeft(10);
            }

            Console.WriteLine(cout);
            //
            //  Initialize the system.
            //
            c1 = ising_2d_initialize(m, n, thresh, ref seed);
            //
            //  Write the initial state to a gnuplot file.
            //
            Plot.plot_file(m, n, c1, "Initial Configuration", "ising_2d_initial.txt",
                "ising_2d_initial.png");
            //
            //  Do the simulation.
            //
            transition(m, n, iterations, prob, thresh, ref seed, c1);
            //
            //  Write the final state to a gnuplot file.
            //
            Plot.plot_file(m, n, c1, title, plot_filename, png_filename);

            Console.WriteLine("");
            Console.WriteLine("ISING_2D_SIMULATION");
            Console.WriteLine("  Normal end of execution.");

            Console.WriteLine("");
        }

        static void ising_2d_agree(int m, int n, int[] c1, ref int[] c5)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ISING_2D_AGREE returns the number of neighbors agreeing with each cell.
            //
            //  Discussion:
            //
            //    The count includes the cell itself, so it is between 1 and 5.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    22 Noveber 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of cells in each 
            //    spatial dimension.
            //
            //    Input, int C1[M*N], an array of 1's and -1's.
            //
            //    Output, int C5[M*N], the number of neighbors 
            //    that agree.  1, 2, 3, 4, or 5.
            //
        {
            int i;
            int im;
            int ip;
            int j;
            int jm;
            int jp;

            for (j = 0; j < n; j++)
            {
                jp = typeMethods.i4_wrap(j + 1, 0, n - 1);
                jm = typeMethods.i4_wrap(j - 1, 0, n - 1);
                for (i = 0; i < m; i++)
                {
                    ip = typeMethods.i4_wrap(i + 1, 0, m - 1);
                    im = typeMethods.i4_wrap(i - 1, 0, m - 1);
                    c5[i + j * m] = c1[i + j * m] + c1[ip + j * m] + c1[im + j * m] + c1[i + jm * m] + c1[i + jp * m];
                    if (0 < c1[i + j * m])
                    {
                        c5[i + j * m] = (5 + c5[i + j * m]) / 2;
                    }
                    else
                    {
                        c5[i + j * m] = (5 - c5[i + j * m]) / 2;
                    }
                }
            }

            return;
        }

        static int[] ising_2d_initialize(int m, int n, double thresh, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ISING_2D_INITIALIZE initializes the Ising array.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    23 November 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns.
            //
            //    Input, double THRESH, the threshhold.
            //
            //    Input/output, int *SEED, a seed for the random 
            //    number generator.
            //
            //    Output, in ISING_2D_INITIALIZE[M*N], the initial Ising array.
            //
        {
            int[] c1;
            int i;
            int j;
            double[] r;

            r = new double[m * n];

            r = UniformRNG.r8mat_uniform_01(m, n, ref seed);

            c1 = new int[m * n];

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    if (r[i + j * m] <= thresh)
                    {
                        c1[i + j * m] = -1;
                    }
                    else
                    {
                        c1[i + j * m] = +1;
                    }
                }
            }

            return c1;
        }

        static void ising_2d_stats(int step, int m, int n, int[] c1)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ISING_2D_STATS prints information about the current step.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    22 November 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int STEP, the step number.
            //
            //    Input, int M, N, the number of rows and columns.
            //
            //    Input, int C1[M*N], the current state of the system.
            //
        {
            int i;
            int j;
            int pos_count;
            double pos_percent;
            int neg_count;
            double neg_percent;

            if (step == 0)
            {
                Console.WriteLine("");
                Console.WriteLine("  Step     Positives       Negatives");
                Console.WriteLine("             #    %%          #    %%");
                Console.WriteLine("");
            }

            pos_count = 0;
            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    if (0 < c1[i + j * m])
                    {
                        pos_count = pos_count + 1;
                    }
                }
            }

            neg_count = m * n - pos_count;
            pos_percent = (double) (100 * pos_count) / (double) (m * n);
            neg_percent = (double) (100 * neg_count) / (double) (m * n);

            Console.WriteLine("  " + step.ToString().PadLeft(4)
                                   + "  " + pos_count.ToString().PadLeft(6)
                                   + "  " + pos_percent.ToString().PadLeft(6)
                                   + "  " + neg_count.ToString().PadLeft(6)
                                   + "  " + neg_percent.ToString().PadLeft(6) + "");
        }

        static void neighbor_2d_stats ( int step, int m, int n, int[] c1, int[] c5 )

        //****************************************************************************80
        /*
          Purpose:
        
            NEIGHBOR_2D_STATS prints neighbor statistics about the current step.
        
          Licensing:
        
            This code is distributed under the GNU LGPL license.
        
          Modified:
        
            23 November 2011
        
          Author:
        
            John Burkardt
        
          Parameters:
        
            Input, int STEP, the step number.
        
            Input, int M, N, the number of rows and columns.
        
            Input, int C1[M*N], the current state of the system.
        
            Input, int C5[M*N], the number of agreeable neighbors.
        */
        {
            int i;
            int j;
            int[] stats = new int[11];

            if ( step == 0 )
            {
                Console.WriteLine("");
                Console.WriteLine("  Step     Neighborhood Charge:");
                Console.WriteLine("           -5    -4    -3    -2    -1     +1    +2    +3    +4    +5");
                Console.WriteLine("");
            }

            for ( i = - 5; i <= 5; i++ )
            {
                stats[i+5] = 0;
            }

            for (j = 0; j < n; j++ )
            {
                for ( i = 0; i < n; i++ )
                {
                    stats[c5[i+j*m]-1+5] = stats[c5[i+j*m]-1+5] + 1;
                }
            }
            string cout = "  " + step.ToString().PadLeft(4);
            for ( i = - 5; i <= 5; i++ )
            {
                if ( i != 0 )
                {
                    cout += "  " + stats[i+5].ToString().PadLeft(4);
                }
            }
            Console.WriteLine(cout);
        }
        
        static void transition(int m, int n, int iterations, double[] prob,
                double thresh, ref int seed, int[] c1)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRANSITION carries out a Monte Carlo simulation of a 3D Ising model.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    30 June 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns.
            //
            //    Input, int ITERATIONS, the number of iterations.
            //
            //    Input, double PROB[5].  PROB[I-1] represents the probability 
            //    that the spin of a given cell will be reversed, given that it has I 
            //    immediate neighbors (including itself) with spin the same as its own.
            //
            //    Input, double THRESH, the threshhold.
            //
            //    Input/output, int *SEED, a seed for the random number 
            //    generator.
            //
            //    Input/output, int C1[M*N], the current state.
        {
            int[] c5;
            int i;
            int j;
            double[] r;
            int step;

            c5 = new int[m * n];

            r = new double[m * n];

            step = 0;
            ising_2d_stats(step, m, n, c1);

            for (step = 1; step <= iterations; step++)
            {
                //
                //  C5 contains 1 through 5, the number of cells that agree with the center cell.
                //
                ising_2d_agree(m, n, c1, ref c5);

                if (false)
                {
                    neighbor_2d_stats(step, m, n, c1, c5);
                }

                //
                //  Determine the chances of flipping cell (I,J).
                //
                r = UniformRNG.r8mat_uniform_01(m, n, ref seed);

                for (j = 0; j < n; j++)
                {
                    for (i = 0; i < m; i++)
                    {
                        if (r[i + j * m] < prob[c5[i + j * m] - 1])
                        {
                            c1[i + j * m] = -c1[i + j * m];
                        }
                    }
                }

                ising_2d_stats(step, m, n, c1);
            }
        }

    }
}