using System;
using Burkardt.Probability;
using Burkardt.RandomNS;

namespace Burkardt.Sampling;

public static class Walker
{
    private static Rand48 rand48 = new();

    public static void walker_build(int n, double[] x, ref double[] y, ref int[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WALKER_BUILD sets up the data for a Walker sampler.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 February 2016
        //
        //  Author:
        //
        //    Original C version by Warren Smith.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Warren Smith,
        //    How to sample from a probability distribution,
        //    April 18, 2002.
        //
        //  Parameters:
        //
        //    Input, int N, indicates the size of X.
        //
        //    Input, double X[N+2], contains in X[1] through X[N] the
        //    probabilities of outcomes 1 through N.  
        //
        //    Output, double Y[N+2], the Walker threshold vector.
        //
        //    Output, int A[N+2], the Walker index vector.
        //
    {
        int[] b;
        int i;
        int j;
        int k;
        //
        //  Initialize A.
        //
        a[0] = 0;
        for (i = 1; i <= n; i++)
        {
            a[i] = i;
        }

        a[n + 1] = n + 1;
        //
        //  Initialize B to the "stay here" value, and set sentinel values at the ends.
        //
        b = new int[n + 2];

        b[0] = 0;
        for (i = 1; i <= n; i++)
        {
            b[i] = i;
        }

        b[n + 1] = n + 1;
        //
        //  Copy Y from X.
        //  Scale the probability vector and set sentinel values at the ends.
        //
        y[0] = 0.0;
        for (i = 1; i <= n; i++)
        {
            y[i] = x[i] * n;
        }

        y[n + 1] = 2.0;

        i = 0;
        j = n + 1;

        for (;;)
        {
            //
            //  Find i so Y[B[i]] needs more.
            //
            do
            {
                i++;
            } while (y[b[i]] < 1.0);

            //
            //  Find j so Y[B[j]] wants less.
            //
            do
            {
                j--;
            } while (1.0 <= y[b[j]]);

            if (j <= i)
            {
                break;
            }

            //
            //  Swap B[i] and B[j].
            //
            k = b[i];
            b[i] = b[j];
            b[j] = k;
        }

        i = j;
        j++;

        while (0 < i)
        {
            //
            //  Find J such that Y[B[j]] needs more.
            //
            while (y[b[j]] <= 1.0)
            {
                j++;
            }

            //
            //  Meanwhile, Y[B[i]] wants less.
            //
            if (n < j)
            {
                break;
            }

            //
            //  B[i] will donate to B[j] to fix up.
            //
            y[b[j]] -= (1.0 - y[b[i]]);
            a[b[i]] = b[j];
            switch (y[b[j]])
            {
                // 
                //  Y[B[j]] now wants less so readjust ordering.
                //
                case < 1.0:
                    k = b[i];
                    b[i] = b[j];
                    b[j] = k;
                    j++;
                    break;
                default:
                    i--;
                    break;
            }
        }
    }

    public static int walker_sampler(int n, double[] y, int[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WALKER_SAMPLER returns a random sample i=1..N with prob X[i].
        //
        //  Discussion:
        //
        //    Implementation of algorithm for sampling from a discrete
        //    probability N-vector X[1], X[2], ..., X[N].  (N>=1.)
        //    Runs on O(1) worst case time per sample,
        //    and uses one integer and one double N-element array for storage.
        //    Preprocessing consumes O(N) time and temporarily uses one 
        //    additional integer array (B[0..N+1]) for bookkeeping. 
        //    X[0] and X[N+1] are also used as sentinels in the Build() algorithm.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 February 2016
        //
        //  Author:
        //
        //    Original C version by Warren Smith.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Warren Smith,
        //    How to sample from a probability distribution,
        //    April 18, 2002.
        //
        //  Parameters:
        //
        //    Input, int N, indicates the size of the probability vector X.
        //
        //    Input, double Y[N+2], the Walker threshold vector.
        //
        //    Input, int A[N+2], the Walker index vector.
        //
        //    Output, int WALKER_SAMPLER, a sample value between 1 and N,
        //    selected according to the probability vector X.
        //
    {
        int i;
        double r;


        // 
        //  Let i = random uniform integer from {1,2,...N};  
        //
        i = 1 + (int) (n * rand48.drand48());

        r = rand48.drand48();

        if (y[i] < r)
        {
            i = a[i];
        }

        return i;
    }

    public static void walker_sampler_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WALKER_SAMPLER_TEST tests WALKER_SAMPLER.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 February 2016
        //
        //  Author:
        //
        //    Original C version by Warren Smith.
        //    C++ version by John Burkardt.
        //
    {
        int[] a;
        int[] count;
        double expval;
        int i;
        int j;
        int n;
        double p;
        int seed;
        int seed2;
        double sum;
        double t;
        double v;
        double[] x;
        double[] y;

        seed = 123456789;
        n = 10;
        p = 2.0;

        Console.WriteLine("");
        Console.WriteLine("WALKER_SAMPLER_TEST:");
        Console.WriteLine("  WALKER_SAMPLER creates Walker sample vectors Y and A");
        Console.WriteLine("  for efficient sampling of a discrete probability vector.");
        Console.WriteLine("  Test the Walker sampler with a Zipf-type probability.");
        //
        //  Initialize the random number generator.
        //
        Console.WriteLine("  Use seed = " + seed + " to initialize srand48():");

        Rand48 rand48 = new();
        rand48.SetSeed((ulong)seed);
        //
        //  "Warm up" the random number generator.
        //
        for (i = 0; i < 100; i++)
        {
            rand48.drand48();
        }

        Console.WriteLine("");
        Console.WriteLine("  After 100 warmup calls, next 3 values of drand48():");

        for (i = 100; i < 103; i++)
        {
            Console.WriteLine("  " + rand48.drand48() + "");
        }

        //
        //  Generate a standard Zipf probability vector for cases 1,...,N,
        //  with parameter P.
        //
        x = Zipf.zipf_probability(n, p);

        Console.WriteLine("");
        Console.WriteLine("  Zipf probabilities");
        Console.WriteLine("  for N = " + n + "");
        Console.WriteLine("  and parameter P = " + p + "");
        Console.WriteLine("");
        Console.WriteLine("     I     X[I]");
        Console.WriteLine("");
        for (i = 1; i <= n; i++)
        {
            Console.WriteLine("  " + i.ToString().PadLeft(4)
                                   + "  " + x[i].ToString().PadLeft(16) + "");
        }

        //
        //  For better testing, randomly scramble the probabilities.
        //
        seed2 = 123456789;
        Permutation.random_permutation(n, ref x, ref seed2);

        Console.WriteLine("");
        Console.WriteLine("  Randomly permuted X:");
        Console.WriteLine("");
        Console.WriteLine("     I     X[I]");
        Console.WriteLine("");
        for (i = 1; i <= n; i++)
        {
            Console.WriteLine("  " + i.ToString().PadLeft(4)
                                   + "  " + x[i].ToString().PadLeft(16) + "");
        }

        //
        //  Build the Walker sampler.
        //
        y = new double[n + 2];
        a = new int[n + 2];

        walker_build(n, x, ref y, ref a);

        Console.WriteLine("");
        Console.WriteLine("  Built the sampler");
        Console.WriteLine("  i Y[i] A[i]:");
        Console.WriteLine("");

        for (i = 1; i <= n; i++)
        {
            Console.WriteLine("  " + i.ToString().PadLeft(3)
                                   + "  " + y[i].ToString().PadLeft(16)
                                   + "  " + a[i].ToString().PadLeft(4) + "");
        }

        //
        //  Prepare to count the frequency of each outcome.
        //
        count = new int[n + 2];
        for (i = 1; i <= n; i++)
        {
            count[i] = 0;
        }

        //
        //  Call the sampler many times.
        //
        for (i = 0; i < 100000; i++)
        {
            j = walker_sampler(n, y, a);
            count[j] += 1;
        }

        //
        //  Compare normalized sample frequencies to the original probabilities in X.
        //
        sum = 0.0;
        Console.WriteLine("");
        Console.WriteLine("  100000 samples:");
        Console.WriteLine("  prob   #samples:");
        Console.WriteLine("");

        for (i = 1; i <= n; i++)
        {
            Console.WriteLine("  " + x[i]
                                   + "  " + count[i] + "");
            expval = x[i] * 100000;
            t = expval - count[i];
            sum += t * t / expval;
        }

        sum /= n;

        Console.WriteLine("");
        Console.WriteLine("  sumvar = " + sum + " (should be about 1)");
        //
        //  Verify the data structure.
        //
        v = walker_verify(n, x, y, a);
        Console.WriteLine("");
        Console.WriteLine("  Verification sum: " + v + "");
        Console.WriteLine("  (Should be close to 0.)");
    }

    public static double walker_verify(int n, double[] x, double[] y,
            int[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WALKER_VERIFY verifies a Walker Sampler structure.
        //
        //  Discussion:
        //
        //    This test applies the sampling algorithms to a Zipfian distribution.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 February 2016
        //
        //  Author:
        //
        //    Original C version by Warren Smith.
        //    C++ version by John Burkardt.
        //
        //  Parameters:
        //
        //    Input, int N, indicates the size of X.
        //
        //    Input, double X[N+2], contains in X[1] through X[N] the
        //    probabilities of outcomes 1 through N.
        //
        //    Input, double Y[N+2], the Walker threshold vector.
        //
        //    Input, int A[N+2], the Walker index vector.
        //
        //    Output, double WALKER_VERIFY, the verification sum, which
        //    should be close to zero.
    {
        int i;
        double v;
        double[] z;

        z = new double[n + 2];
        //
        //  Reverse the scaling.
        //
        for (i = 1; i <= n; i++)
        {
            z[i] = y[i] / n;
        }

        //
        //  Add back the adjustments.
        //
        for (i = 1; i <= n; i++)
        {
            z[a[i]] += (1.0 - y[i]) / n;
        }

        //
        //  Check for discrepancies between Z and X.
        //
        v = 0.0;
        for (i = 1; i <= n; i++)
        {
            v += Math.Abs(z[i] - x[i]);
        }

        return v;
    }
}