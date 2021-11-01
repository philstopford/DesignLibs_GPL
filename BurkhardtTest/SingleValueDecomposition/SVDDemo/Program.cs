using System;
using Burkardt;
using Burkardt.RankingNS;
using Burkardt.SolveNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace SVDDemo
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for SVD_DEMO.
            //
            //  Discussion:
            //
            //    SVD_DEMO demonstrates the singular value decomposition.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 September 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Usage:
            //
            //    svd_demo m n
            //
            //  Command Parameters:
            //
            //    Command parameter, integer M, N, the number of rows and columns
            //    of the matrix.
            //
            //  Local Parameters:
            //
            //    Local, double A[M*N], the matrix whose singular value
            //    decomposition we are investigating.
            //
            //    Local, double S[M*N], the diagonal factor
            //    in the singular value decomposition of A.
            //
            //    Local, int SEED, a seed used to define the random number generator.
            //
            //    Output, double U[M*M], the first orthogonal factor
            //    in the singular value decomposition of A.
            //
            //    Output, double V[N*N], the second orthogonal factor
            //    in the singular value decomposition of A.
            //
        {
            double[] a;
            double[] a_pseudo;
            int i;
            int j;
            int m;
            int n;
            double[] s;
            int seed;
            string string_;
            double[] u;
            double[] v;

            Console.WriteLine("");
            Console.WriteLine("SVD_DEMO:");
            Console.WriteLine("");
            Console.WriteLine("  Demonstrate the singular value decomposition (SVD)");
            Console.WriteLine("");
            Console.WriteLine("  A real MxN matrix A can be factored as:");
            Console.WriteLine("");
            Console.WriteLine("    A = U * S * V'");
            Console.WriteLine("");
            Console.WriteLine("  where");
            Console.WriteLine("");
            Console.WriteLine("    U = MxM orthogonal,");
            Console.WriteLine("    S = MxN zero except for diagonal,");
            Console.WriteLine("    V = NxN orthogonal.");
            Console.WriteLine("");
            Console.WriteLine("  The diagonal of S contains only nonnegative numbers");
            Console.WriteLine("  and these are arranged in descending order.");
            //
            //  If M was not on the command line, get it now.
            //
            try
            {
                string_ = args[0];
            }
            catch (Exception e)
            {
                Console.WriteLine("");
                Console.WriteLine("SVD_DEMO:");
                Console.WriteLine("  Please enter the value of M:");
                Console.WriteLine("  (Number of rows in matrix A).");
                Console.WriteLine("  (We prefer M <= 10!).");

                string_ = Console.ReadLine();
            }

            m = Convert.ToInt32(string_);
            Console.WriteLine("");
            Console.WriteLine("  Matrix row order    M = " + m + "");

            //
            //  If N was not on the command line, get it now.
            //
            try
            {
                string_ = args[1];
            }
            catch (Exception e)
            {
                Console.WriteLine("");
                Console.WriteLine("SVD_DEMO:");
                Console.WriteLine("  Please enter the value of N:");
                Console.WriteLine("  (Number of columns in matrix A).");
                Console.WriteLine("  (We prefer N <= 10!).");

                string_ = Console.ReadLine();
            }

            n = Convert.ToInt32(string_);
            Console.WriteLine("  Matrix column order N = " + n + "");

            //
            //  If SEED was not on the command line, use GET_SEED.
            //
            try
            {
                seed = Convert.ToInt32(args[2]);
                Console.WriteLine("  Random number SEED    = " + seed + "");
                Console.WriteLine("  (Chosen by the user.)");
            }
            catch (Exception e)
            {
                seed = entropyRNG.RNG.nextint();
                Console.WriteLine("  Random number SEED    = " + seed + "");
                Console.WriteLine("  (Chosen by the program.)");
            }

            //
            //  Set aside space for the arrays.
            //
            u = new double[m * m];
            s = new double[m * n];
            v = new double[n * n];
            //
            //  Generate the matrix.
            //
            Console.WriteLine("");
            Console.WriteLine("  We choose a \"random\" matrix A, with integral");
            Console.WriteLine("  values between 0 and 10.");

            a = UniformRNG.r8mat_uniform_01_new(m, n, ref seed);

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    a[i + m * j] = typeMethods.r8_nint(10.0 * a[i + m * j]);
                }
            }

            typeMethods.r8mat_print(m, n, a, "  The matrix A:");
            //
            //  Get the SVD from LINPACK.
            //
            SingleValueDecomposition.r8mat_svd_linpack(m, n, a, ref u, ref s, ref v);
            //
            //  Print the SVD.
            //
            typeMethods.r8mat_print(m, m, u, "  The orthogonal factor U:");

            typeMethods.r8mat_print(m, n, s, "  The diagonal factor S:");

            typeMethods.r8mat_print(n, n, v, "  The orthogonal factor V:");
            //
            //  Check that A = U * S * V'.
            //
            SingleValueDecomposition.svd_product_test(m, n, a, u, s, v);
            //
            //  Compute the norm of the difference between A and the successive
            //  sums of rank one approximants.
            //
            Ranking.rank_one_test(m, n, a, u, s, v);
            //
            //  Actually print the sums of rank one approximants.
            //
            Ranking.rank_one_print_test(m, n, a, u, s, v);
            //
            //  Compute the pseudoinverse.
            //
            a_pseudo = PseudoLinear.pseudo_inverse(m, n, u, s, v);

            typeMethods.r8mat_print(n, m, a_pseudo, "  The pseudoinverse of A:");
            //
            //  Test A*A+ = I+, A+*A = I+
            //
            PseudoLinear.pseudo_product_test(m, n, a, a_pseudo);
            //
            //  Demonstrate the use of the pseudoinverse for linear systems.
            //
            PseudoLinear.pseudo_linear_solve_test(m, n, a, a_pseudo, ref seed);

            Console.WriteLine("");
            Console.WriteLine("SVD_DEMO:");
            Console.WriteLine("  Normal end of execution.");

            Console.WriteLine("");
        }
    }
}