using System;
using Burkardt.Types;
using Burkardt.Values;

namespace NormalDataset;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for NORMAL_DATASET.
        //
        //  Discussion:
        //
        //    NORMAL_DATASET generates a dataset of multivariate normal random values,
        //    and writes it to a file.
        //
        //  Usage:
        //
        //    normal_dataset m n seed mu a
        //
        //    where
        //
        //    * M, the spatial dimension,
        //    * N, the number of points to generate,
        //    * SEED, the seed, a positive integer.
        //    * MU, the mean vector.
        //    * A, the MxM variance-covariance matrix.
        //
        //    creates an M by N multivariate normal random dataset and writes it 
        //    to the file "normal_M_N.txt"./
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 December 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a;
        string output_filename;
        int i;
        int j;
        int k;
        int m;
        double[] mu;
        string m_ostring;
        int n;
        string n_ostring;
        double[] x;
        int seed;


        Console.WriteLine("");
        Console.WriteLine("NORMAL_DATASET");
        Console.WriteLine("");
        Console.WriteLine("  Generate a multivariate normal random dataset.");
        Console.WriteLine("");
        Console.WriteLine("  The program requests input values from the user:");
        Console.WriteLine("");
        Console.WriteLine("  * M, the spatial dimension,");
        Console.WriteLine("  * N, the number of points to generate,");
        Console.WriteLine("  * SEED, a positive integer.");
        Console.WriteLine("  * MU, the mean vector of length M.");
        Console.WriteLine("  * A, the MxM variance-covariance matrix.");
        Console.WriteLine("");
        Console.WriteLine("  The program generates the data, writes it to the file");
        Console.WriteLine("");
        Console.WriteLine("    normal_M_N.txt");
        Console.WriteLine("");
        Console.WriteLine("  where \"M\" and \"N\" are the numeric values specified");
        Console.WriteLine("  by the user.");
        //
        //  Get the spatial dimension.
        //
        try
        {
            m = Convert.ToInt32(args[0]);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("  Enter the value of M");
            m = Convert.ToInt32(Console.ReadLine());
        }

        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension M = " + m + "");
        //
        //  Get the number of points.
        //
        try
        {
            n = Convert.ToInt32(args[1]);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("  Enter the number of points N");
            n = Convert.ToInt32(Console.ReadLine());
        }

        Console.WriteLine("");
        Console.WriteLine("  Number of points N = " + n + "");
        //
        //  Get the seed.
        //
        try
        {
            seed = Convert.ToInt32(args[2]);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("  Enter the value of SEED");
            seed = Convert.ToInt32(Console.ReadLine());
        }

        Console.WriteLine("");
        Console.WriteLine("  The seed is = " + seed + "");
        //
        //  Get the mean.
        //
        mu = new double[m];

        k = 4;

        try
        {
            for (i = 0; i < m; i++)
            {
                mu[i] = Convert.ToDouble(args[k]);
                k += 1;
            }
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("  Enter MU:");
            for (i = 0; i < m; i++)
            {
                mu[i] = Convert.ToDouble(Console.ReadLine());
            }
        }

        typeMethods.r8vec_print(m, mu, "  The mean vector M:");
        //
        //  Get the variance-covariance matrix.
        //
        a = new double[m * m];

        try
        {
            for (i = 0; i < m; i++)
            {
                for (j = 0; j < m; j++)
                {
                    a[i + j * m] = Convert.ToDouble(args[k]);
                    k += 1;
                }
            }
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("  Enter A:");
            for (i = 0; i < m; i++)
            {
                for (j = 0; j < m; j++)
                {
                    a[i + j * m] = Convert.ToDouble(Console.ReadLine());
                }
            }
        }

        typeMethods.r8mat_print(m, m, a, "  The variance-covariance matrix A:");
        //
        //  Compute the data.
        //
        x = Normal.multinormal_sample(m, n, a, mu, ref seed);
        //
        //  Write it to a file.
        //
        m_ostring = m.ToString();
        n_ostring = n.ToString();

        output_filename = "normal_" + m_ostring + "_"
                          + n_ostring + ".txt";

        typeMethods.r8mat_write(output_filename, m, n, x);

        Console.WriteLine("");
        Console.WriteLine("  The data was written to the file \""
                          + output_filename + "\".");

        Console.WriteLine("");
        Console.WriteLine("NORMAL_DATASET:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }
}