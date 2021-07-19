using System;
using System.Numerics;
using Burkardt.BLAS;
using Burkardt.Linpack;
using Burkardt.Types;
using Burkardt.Uniform;

namespace LinpackZTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for LINPACK_Z_TEST.
            //
            //  Discussion:
            //
            //    LINPACK_Z_TEST tests the LINPACK_Z library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 October 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("LINPACK_Z_TEST");
            Console.WriteLine("  Test the LINPACK_Z library.");

            test01();
            test02();
            test03();
            test04();
            test05();
            test06();
            test07();
            test08();
            test09();

            test10();
            test11();
            test12();
            test13();
            test14();
            test15();
            test16();
            test17();
            test18();
            test19();

            test20();
            test21();
            test22();
            test23();
            test24();
            test25();
            test26();
            zqrdc_test();
            test28();
            test29();

            test30();
            test31();
            test32();
            test33();
            test34();
            test345();
            test35();
            test36();
            test37();
            //
            //  Terminate.
            //
            Console.WriteLine("");
            Console.WriteLine("LINPACK_Z_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void test01()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 tests ZCHDC.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 3;

            Complex[] a = new Complex[N * N];
            Complex[] c = new Complex[N * N];
            int i;
            int info;
            int[] ipvt = new int[N];
            int j;
            int job;
            int k;
            int lda;
            lda = N;
            string cout;

            Console.WriteLine("");
            Console.WriteLine("TEST01");
            Console.WriteLine("  For a complex Hermitian positive definite matrix,");
            Console.WriteLine("  ZCHDC computes the Cholesky decomposition.");
            Console.WriteLine("");
            Console.WriteLine("  The number of equations is N = " + N + "");
            //
            //  Set the values of the matrix A.
            //
            a[0 + 0 * lda] = new Complex(2.5281, 0.0000);
            a[1 + 0 * lda] = new Complex(2.1341, 0.2147);
            a[2 + 0 * lda] = new Complex(2.4187, -0.2932);

            a[0 + 1 * lda] = new Complex(2.1341, -0.2147);
            a[1 + 1 * lda] = new Complex(3.0371, 0.0000);
            a[2 + 1 * lda] = new Complex(2.0905, -1.1505);

            a[0 + 2 * lda] = new Complex(2.4187, 0.2932);
            a[1 + 2 * lda] = new Complex(2.0905, 1.1505);
            a[2 + 2 * lda] = new Complex(2.7638, 0.0000);

            Console.WriteLine("");
            Console.WriteLine("  The matrix:");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                cout = "";
                for (j = 0; j < N; j++)
                {
                    cout += "  " + a[i + j * lda].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }

            //
            //  Decompose the matrix.
            //
            Console.WriteLine("");
            Console.WriteLine("  Decompose the matrix.");

            job = 0;
            for (i = 0; i < N; i++)
            {
                ipvt[i] = 0;
            }

            info = ZCHDC.zchdc(ref a, lda, N, ref ipvt, job);

            if (info != N)
            {
                Console.WriteLine("");
                Console.WriteLine("  ZCHDC returned INFO = " + info + "");
                Console.WriteLine("  The matrix is not Hermitian positive definite.");
                return;
            }

            //
            //  Zero out the lower diagonal.
            //
            for (i = 1; i < N; i++)
            {
                for (j = 0; j < i; j++)
                {
                    a[i + j * lda] = new Complex(0.0, 0.0);
                }
            }

            //
            //  Print the factorization.
            //
            Console.WriteLine("");
            Console.WriteLine("  The Cholesky factor U:");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                cout = "";
                for (j = 0; j < N; j++)
                {
                    cout += "  " + a[i + j * lda].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }

            //
            //  Compute the Cholesky product.
            //
            for (i = 0; i < N; i++)
            {
                for (j = 0; j < N; j++)
                {
                    c[i + j * N] = new Complex(0.0, 0.0);
                    for (k = 0; k < N; k++)
                    {
                        c[i + j * N] = c[i + j * N] + Complex.Conjugate(a[k + i * lda]) * a[k + j * lda];
                    }
                }
            }

            Console.WriteLine("");
            Console.WriteLine("  The product U^H * U:");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                cout = "";
                for (j = 0; j < N; j++)
                {
                    cout += "  " + c[i + j * N].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }

        }

        static void test02()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST02 tests ZCHEX.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 October 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 3;
            int NZ = 1;

            Complex[] a = new Complex[N * N];
            Complex[] b = new Complex[N * N];
            double[] c = new double[N];
            int i;
            int info;
            int[] ipvt = new int[N];
            int j;
            int job;
            int k;
            int l;
            int lda;
            int ldz;
            Complex[] s = new Complex[N];
            Complex[] z = new Complex[N * NZ];
            string cout;

            lda = N;
            ldz = N;

            Console.WriteLine("");
            Console.WriteLine("TEST02");
            Console.WriteLine("  For a complex Hermitian positive definite matrix,");
            Console.WriteLine("  ZCHEX can shift rows and columns in a Cholesky factorization.");
            Console.WriteLine("");
            Console.WriteLine("  The number of equations is N = " + N + "");
            //
            //  Set the values of the matrix A.
            //
            a[0 + 0 * lda] = new Complex(2.5281, 0.0000);
            a[1 + 0 * lda] = new Complex(2.1341, 0.2147);
            a[2 + 0 * lda] = new Complex(2.4187, -0.2932);

            a[0 + 1 * lda] = new Complex(2.1341, -0.2147);
            a[1 + 1 * lda] = new Complex(3.0371, 0.0000);
            a[2 + 1 * lda] = new Complex(2.0905, -1.1505);

            a[0 + 2 * lda] = new Complex(2.4187, 0.2932);
            a[1 + 2 * lda] = new Complex(2.0905, 1.1505);
            a[2 + 2 * lda] = new Complex(2.7638, 0.0000);

            Console.WriteLine("");
            Console.WriteLine("  The matrix A:");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                cout = "";
                for (j = 0; j < N; j++)
                {
                    cout += "  " + a[i + j * lda].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }

            for (i = 0; i < N; i++)
            {
                z[i] = new Complex(i + 1, 0.0);
            }

            Console.WriteLine("");
            Console.WriteLine("  The vector Z:");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                Console.WriteLine("  " + z[i].ToString().PadLeft(20) + "");
            }

            //
            //  Decompose the matrix.
            //
            Console.WriteLine("");
            Console.WriteLine("  Decompose the matrix.");

            job = 0;
            for (i = 0; i < N; i++)
            {
                ipvt[i] = 0;
            }

            info = ZCHDC.zchdc(ref a, lda, N, ref ipvt, job);

            if (info != N)
            {
                Console.WriteLine("");
                Console.WriteLine("  ZCHDC returned INFO = " + info + "");
                Console.WriteLine("  This means the matrix is not positive definite.");
                return;
            }

            //
            //  Zero out the lower diagonal.
            //
            for (i = 1; i < N; i++)
            {
                for (j = 0; j < i; j++)
                {
                    a[i + j * lda] = new Complex(0.0, 0.0);
                }
            }

            //
            //  Print the factorization.
            //
            Console.WriteLine("");
            Console.WriteLine("  The Cholesky factor U:");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                cout = "";
                for (j = 0; j < N; j++)
                {
                    cout += "  " + a[i + j * lda].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }

            //
            //  Right circular shift columns L through K.
            //
            k = 1;
            l = 3;

            Console.WriteLine("");
            Console.WriteLine("  Right circular shift rows and columns K  = " + k
                                                                              + " through L = " + l + "");
            Console.WriteLine("");
            Console.WriteLine("  Logical matrix is now:");
            Console.WriteLine("");
            Console.WriteLine("  33 31 32");
            Console.WriteLine("  13 11 12");
            Console.WriteLine("  23 21 22");

            job = 1;
            ZCHEX.zchex(ref a, lda, N, k, l, ref z, ldz, NZ, ref c, ref s, job);
            //
            //  Left circular shift columns K+1 through L.
            //
            Console.WriteLine("");
            Console.WriteLine("  Left circular shift rows and columns K+1 = " + k + 1
                              + " through L = " + l + "");
            Console.WriteLine("");
            Console.WriteLine("  Logical matrix is now:");
            Console.WriteLine("");
            Console.WriteLine("  33 32 31");
            Console.WriteLine("  23 22 21");
            Console.WriteLine("  13 12 11");

            job = 2;
            ZCHEX.zchex(ref a, lda, N, k + 1, l, ref z, ldz, NZ, ref c, ref s, job);
            //
            //  Print the factorization.
            //
            Console.WriteLine("");
            Console.WriteLine("  The shifted Cholesky factor UU:");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                cout = "";
                for (j = 0; j < N; j++)
                {
                    cout += "  " + a[i + j * lda].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }

            Console.WriteLine("");
            Console.WriteLine("  The shifted vector ZZ:");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                Console.WriteLine("  " + z[i].ToString().PadLeft(20) + "");
            }

            //
            //  Compute the Cholesky product.
            //
            for (i = 0; i < N; i++)
            {
                for (j = 0; j < N; j++)
                {
                    b[i + j * N] = new Complex(0.0, 0.0);
                    for (k = 0; k < N; k++)
                    {
                        b[i + j * N] = b[i + j * N] + Complex.Conjugate(a[k + i * lda]) * a[k + j * lda];
                    }
                }
            }

            Console.WriteLine("");
            Console.WriteLine("  The shifted product AA = UU' * UU:");
            Console.WriteLine("  The rows and columns of the original matrix A reappear,");
            Console.WriteLine("  but in reverse order.");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                cout = "";
                for (j = 0; j < N; j++)
                {
                    cout += "  " + b[i + j * lda].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }
        }

        static void test03()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST03 tests ZCHUD and ZTRSL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int P = 20;
            int NZ = 1;

            Complex[] b = new Complex[P];
            double[] c = new double[P];
            int i;
            int j;
            int job;
            int ldr;
            int ldz;
            Complex[] r = new Complex[P * P];
            double[] rho = new double[NZ];
            Complex[] row;
            Complex[] s = new Complex[P];
            int seed;
            Complex[] x = new Complex[P];
            Complex[] y = new Complex[NZ];
            Complex[] z = new Complex[P * NZ];

            ldr = P;
            ldz = P;

            Console.WriteLine("");
            Console.WriteLine("TEST03");
            Console.WriteLine("  For a complex Hermitian matrix");
            Console.WriteLine("  ZCHUD updates a Cholesky decomposition.");
            Console.WriteLine("  ZTRSL solves a triangular linear system.");
            Console.WriteLine("");
            Console.WriteLine("  In this example, we use ZCHUD to solve a");
            Console.WriteLine("  least squares problem R * b = z.");
            Console.WriteLine("");
            Console.WriteLine("  The number of equations is P = " + P + "");
            //
            //  Initialize.
            //
            for (j = 0; j < P; j++)
            {
                for (i = 0; i < P; i++)
                {
                    r[i + j * ldr] = new Complex(0.0, 0.0);
                }
            }

            for (j = 0; j < NZ; j++)
            {
                for (i = 0; i < P; i++)
                {
                    z[i + j * ldz] = new Complex(0.0, 0.0);
                }
            }

            for (i = 1; i <= P; i++)
            {
                x[i - 1] = new Complex(i, (i % 2));
            }

            //
            //  Use ZCHUD to form R, Z and RHO by adding X and Y a row at a time.
            //  X is a row of the least squares matrix and Y the right hand side.
            //
            seed = 123456789;

            for (i = 1; i <= P; i++)
            {
                row = UniformRNG.c8vec_uniform_01_new(P, ref seed);
                y[0] = BLAS1Z.zdotu(P, row, 1, x, 1);
                rho[0] = 0.0;
                ZCHUD.zchud(ref r, ldr, P, row, ref z, ldz, NZ, ref y, ref rho, ref c, ref s);
            }

            //
            //  Generate the least squares solution, b = inverse ( R ) * Z.
            //
            for (j = 1; j <= NZ; j++)
            {
                for (i = 1; i <= P; i++)
                {
                    b[i - 1] = z[i - 1 + (j - 1) * ldz];
                }

                job = 1;

                ZTRSL.ztrsl(r, ldr, P, ref b, job);

                Console.WriteLine("");
                Console.WriteLine("  Solution vector # " + j + "");
                Console.WriteLine("  (Should be (1,1) (2,0), (3,1) (4,0) ...)");
                Console.WriteLine("");

                for (i = 1; i <= P; i++)
                {
                    if (i <= 5 || P - 5 < i)
                    {
                        Console.WriteLine("  " + i.ToString().PadLeft(8)
                                               + "  " + b[i - 1].ToString().PadLeft(20) + "");
                    }

                    if (i == 5)
                    {
                        Console.WriteLine("  ......  ..............");
                    }
                }
            }
        }

        static void test04()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST04 tests ZGBCO.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 3;
            int ML = 1;
            int MU = 1;

            Complex[] a = new Complex[(2 * ML + MU + 1) * N];
            Complex[] a_save = new Complex[N * N];
            int i;
            int i1;
            int i2;
            int[] ipvt = new int[N];
            int j;
            int k;
            int lda;
            int m;
            double rcond;
            int seed;
            string cout;

            lda = 2 * ML + MU + 1;

            Console.WriteLine("");
            Console.WriteLine("TEST04");
            Console.WriteLine("  For a complex general band storage matrix:");
            Console.WriteLine("  ZGBCO factors the matrix and estimates the");
            Console.WriteLine("  reciprocal condition number.");
            Console.WriteLine("");
            Console.WriteLine("  The matrix order is N = " + N + "");
            Console.WriteLine("  The lower band is ML =  " + ML + "");
            Console.WriteLine("  The upper band is MU =  " + MU + "");
            //
            //  Set the values of the matrix A.
            //
            for (j = 0; j < N; j++)
            {
                for (i = 0; i < N; i++)
                {
                    a_save[i + j * N] = new Complex(0.0, 0.0);
                }
            }

            m = ML + MU + 1;

            seed = 123456789;

            for (j = 1; j <= N; j++)
            {
                i1 = Math.Max(1, j - MU);
                i2 = Math.Min(N, j + ML);
                for (i = i1; i <= i2; i++)
                {
                    k = i - j + m;
                    a[k - 1 + (j - 1) * lda] = UniformRNG.c8_uniform_01(ref seed);
                    a_save[i - 1 + (j - 1) * N] = a[k - 1 + (j - 1) * lda];
                }
            }

            Console.WriteLine("");
            Console.WriteLine("  The matrix:");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                cout = "";
                for (j = 0; j < N; j++)
                {
                    cout += "  " + a_save[i + j * N].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }

            //
            //  Factor the matrix A.
            //
            rcond = ZGBCO.zgbco(ref a, lda, N, ML, MU, ref ipvt);

            Console.WriteLine("");
            Console.WriteLine("  Estimated reciprocal condition RCOND = " + rcond + "");
        }

        static void test05()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST05 tests ZGBFA and ZGBSL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 3;
            int ML = 1;
            int MU = 1;

            Complex[] a = new Complex[(2 * ML + MU + 1) * N];
            Complex[] a_save = new Complex[N * N];
            Complex[] b = new Complex[N];
            int i;
            int i1;
            int i2;
            int info;
            int[] ipvt = new int[N];
            int j;
            int job;
            int k;
            int lda;
            int m;
            int seed;
            Complex[] x;
            string cout;

            lda = 2 * ML + MU + 1;

            Console.WriteLine("");
            Console.WriteLine("TEST05");
            Console.WriteLine("  For a complex general band storage matrix:");
            Console.WriteLine("  ZGBFA factors the matrix;");
            Console.WriteLine("  ZGBSL solves a factored linear system.");
            Console.WriteLine("");
            Console.WriteLine("  The matrix order is N = " + N + "");
            Console.WriteLine("  The lower band is ML =  " + ML + "");
            Console.WriteLine("  The upper band is MU =  " + MU + "");
            //
            //  Set the values of the matrix A.
            //
            for (j = 0; j < N; j++)
            {
                for (i = 0; i < N; i++)
                {
                    a_save[i + j * N] = new Complex(0.0, 0.0);
                }
            }

            m = ML + MU + 1;

            seed = 123456789;

            for (j = 1; j <= N; j++)
            {
                i1 = Math.Max(1, j - MU);
                i2 = Math.Min(N, j + ML);
                for (i = i1; i <= i2; i++)
                {
                    k = i - j + m;
                    a[k - 1 + (j - 1) * lda] = UniformRNG.c8_uniform_01(ref seed);
                    a_save[i - 1 + (j - 1) * N] = a[k - 1 + (j - 1) * lda];
                }
            }

            Console.WriteLine("");
            Console.WriteLine("  The matrix:");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                cout = "";
                for (j = 0; j < N; j++)
                {
                    cout += "  " + a_save[i + j * N].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }

            //
            //  Set the values of the right hand side vector B.
            //
            x = UniformRNG.c8vec_uniform_01_new(N, ref seed);

            for (i = 0; i < N; i++)
            {
                b[i] = 0.0;
                for (j = 0; j < N; j++)
                {
                    b[i] = b[i] + a_save[i + j * N] * x[j];
                }
            }

            Console.WriteLine("");
            Console.WriteLine("  The right hand side:");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                Console.WriteLine("  " + b[i].ToString().PadLeft(20) + "");
            }

            //
            //  Factor the matrix A.
            //
            info = ZGBFA.zgbfa(ref a, lda, N, ML, MU, ref ipvt);

            if (info != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("  ZGBFA returned INFO = " + info + "");
                return;
            }

            //
            //  Solve the system.
            //
            job = 0;
            ZGBSL.zgbsl(a, lda, N, ML, MU, ipvt, ref b, job);

            Console.WriteLine("");
            Console.WriteLine("  Computed                     Exact");
            Console.WriteLine("  Solution                     Solution");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                Console.WriteLine("  " + b[i].ToString().PadLeft(20)
                                       + "  " + x[i].ToString().PadLeft(20) + "");
                ;
            }
        }

        static void test06()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST06 tests ZGBFA and ZGBDI.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 3;
            int ML = 1;
            int MU = 1;

            Complex[] a = new Complex[(2 * ML + MU + 1) * N];
            Complex[] a_save = new Complex[N * N];
            Complex[] det = new Complex[2];
            int i;
            int i1;
            int i2;
            int info;
            int[] ipvt = new int[N];
            int j;
            int k;
            int lda;
            int m;
            int seed;
            string cout;

            lda = 2 * ML + MU + 1;
            Console.WriteLine("");
            Console.WriteLine("TEST06");
            Console.WriteLine("  For a complex general band storage matrix:");
            Console.WriteLine("  ZGBFA factors the matrix.");
            Console.WriteLine("  ZGBDI computes the determinant.");
            Console.WriteLine("");
            Console.WriteLine("  The matrix order is N = " + N + "");
            Console.WriteLine("  The lower band is ML =  " + ML + "");
            Console.WriteLine("  The upper band is MU =  " + MU + "");
            //
            //  Set the values of the matrix A.
            //
            for (j = 0; j < N; j++)
            {
                for (i = 0; i < N; i++)
                {
                    a_save[i + j * N] = new Complex(0.0, 0.0);
                }
            }

            m = ML + MU + 1;

            seed = 123456789;

            for (j = 1; j <= N; j++)
            {
                i1 = Math.Max(1, j - MU);
                i2 = Math.Min(N, j + ML);
                for (i = i1; i <= i2; i++)
                {
                    k = i - j + m;
                    a[k - 1 + (j - 1) * lda] = UniformRNG.c8_uniform_01(ref seed);
                    a_save[i - 1 + (j - 1) * N] = a[k - 1 + (j - 1) * lda];
                }
            }

            Console.WriteLine("");
            Console.WriteLine("  The matrix:");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                cout = "";
                for (j = 0; j < N; j++)
                {
                    cout += "  " + a_save[i + j * N].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }

            //
            //  Factor the matrix A.
            //
            info = ZGBFA.zgbfa(ref a, lda, N, ML, MU, ref ipvt);

            if (info != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("  ZGBFA returned INFO = " + info + "");
                return;
            }

            //
            //  Get the determinant.
            //
            ZGBDI.zgbdi(a, lda, N, ML, MU, ipvt, ref det);

            Console.WriteLine("");
            Console.WriteLine("  Determinant = " + det[0]
                                                 + "  * 10^ " + det[1] + "");
        }

        static void test07()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST07 tests ZGECO.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 3;

            Complex[] a;
            int i;
            int[] ipvt = new int[N];
            int j;
            int lda;
            double rcond;
            int seed;
            string cout;

            lda = N;

            Console.WriteLine("");
            Console.WriteLine("TEST07");
            Console.WriteLine("  For a complex general storage matrix:");
            Console.WriteLine("  ZGECO factors the matrix and estimates the");
            Console.WriteLine("  reciprocal condition number.");
            Console.WriteLine("");
            Console.WriteLine("  The matrix order is N = " + N + "");
            //
            //  Set the values of the matrix A.
            //
            seed = 123456789;

            a = UniformRNG.c8mat_uniform_01_new(N, N, ref seed);

            Console.WriteLine("");
            Console.WriteLine("  The matrix:");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                cout = "";
                for (j = 0; j < N; j++)
                {
                    cout += "  " + a[i + j * N].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }

            //
            //  Factor the matrix A.
            //
            rcond = ZGECO.zgeco(ref a, lda, N, ref ipvt);

            Console.WriteLine("");
            Console.WriteLine("  Estimated reciprocal condition RCOND = " + rcond + "");
        }

        static void test08()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST08 tests ZGEFA and ZGESL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 3;

            Complex[] a;
            Complex[] b = new Complex[N];
            int i;
            int info;
            int[] ipvt = new int[N];
            int j;
            int job;
            int lda;
            int seed;
            Complex[] x;
            string cout;

            lda = N;

            Console.WriteLine("");
            Console.WriteLine("TEST08");
            Console.WriteLine("  For a complex general storage matrix:");
            Console.WriteLine("  ZGEFA factors the matrix.");
            Console.WriteLine("  ZGESL solves a linear system.");
            Console.WriteLine("");
            Console.WriteLine("  The matrix order is N = " + N + "");
            //
            //  Set the values of the matrix A.
            //
            seed = 123456789;

            a = UniformRNG.c8mat_uniform_01_new(N, N, ref seed);

            Console.WriteLine("");
            Console.WriteLine("  The matrix:");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                cout = "";
                for (j = 0; j < N; j++)
                {
                    cout += "  " + a[i + j * N].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }

            //
            //  Set the values of the right hand side vector B.
            //
            x = UniformRNG.c8vec_uniform_01_new(N, ref seed);

            for (i = 0; i < N; i++)
            {
                b[i] = 0.0;
                for (j = 0; j < N; j++)
                {
                    b[i] = b[i] + a[i + j * N] * x[j];
                }
            }

            Console.WriteLine("");
            Console.WriteLine("  The right hand side:");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                Console.WriteLine("  " + b[i].ToString().PadLeft(20) + "");
            }

            //
            //  Factor the matrix A.
            //
            info = ZGEFA.zgefa(ref a, lda, N, ref ipvt);

            if (info != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("  ZGEFA returned an error flag INFO = " + info + "");
                return;
            }

            //
            //  Solve the system.
            //
            job = 0;
            ZGESL.zgesl(a, lda, N, ipvt, ref b, job);

            Console.WriteLine("");
            Console.WriteLine("  Computed                     Exact");
            Console.WriteLine("  Solution                     Solution");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                Console.WriteLine("  " + b[i].ToString().PadLeft(20)
                                       + "  " + x[i].ToString().PadLeft(20) + "");
                ;
            }
        }

        static void test09()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST09 tests ZGEFA and ZGEDI.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 3;

            Complex[] a;
            Complex[] a_save = new Complex[N * N];
            Complex[] c = new Complex[N * N];
            Complex[] det = new Complex[2];
            int i;
            int info;
            int[] ipvt = new int[N];
            int j;
            int job;
            int k;
            int lda;
            int seed;
            string cout;

            lda = N;

            Console.WriteLine("");
            Console.WriteLine("TEST09");
            Console.WriteLine("  For a complex general storage matrix:");
            Console.WriteLine("  ZGEFA factors the matrix.");
            Console.WriteLine("  ZGEDI computes the determinant or inverse.");
            Console.WriteLine("");
            Console.WriteLine("  The matrix order is N = " + N + "");
            //
            //  Set the values of the matrix A.
            //
            seed = 123456789;

            a = UniformRNG.c8mat_uniform_01_new(N, N, ref seed);

            for (i = 0; i < N; i++)
            {
                for (j = 0; j < N; j++)
                {
                    a_save[i + j * N] = a[i + j * N];
                }
            }

            Console.WriteLine("");
            Console.WriteLine("  The matrix:");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                cout = "";
                for (j = 0; j < N; j++)
                {
                    cout += "  " + a[i + j * N].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }

            //
            //  Factor the matrix A.
            //
            info = ZGEFA.zgefa(ref a, lda, N, ref ipvt);

            if (info != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("  ZGEFA returned an error flag INFO = " + info + "");
                return;
            }

            //
            //  Get the determinant.
            //
            job = 10;
            ZGEDI.zgedi(ref a, lda, N, ipvt, ref det, job);

            Console.WriteLine("");
            Console.WriteLine("  Determinant = " + det[0]
                                                 + " * 10^ " + det[1] + "");
            //
            //  Get the inverse.
            //
            job = 01;
            ZGEDI.zgedi(ref a, lda, N, ipvt, ref det, job);

            for (i = 0; i < N; i++)
            {
                for (k = 0; k < N; k++)
                {
                    c[i + k * N] = 0.0;
                    for (j = 0; j < N; j++)
                    {
                        c[i + k * N] = c[i + k * N] + a[i + j * N] * a_save[j + k * N];
                    }
                }
            }

            Console.WriteLine("");
            Console.WriteLine("  The product inv(A) * A is");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                cout = "";
                for (j = 0; j < N; j++)
                {
                    cout += "  " + c[i + j * N].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }
        }

        static void test10()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST10 tests ZGTSL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 10;

            Complex[] b = new Complex[N];
            Complex[] c;
            Complex[] d = new Complex[N];
            Complex[] e;
            int i;
            int seed;
            Complex[] x = new Complex[N];

            Console.WriteLine("");
            Console.WriteLine("TEST10");
            Console.WriteLine("  For a complex tridiagonal matrix:");
            Console.WriteLine("  ZGTSL solves a linear system.");
            Console.WriteLine("");
            Console.WriteLine("  Matrix order N = " + N + "");
            //
            //  Set the matrix.
            //
            seed = 123456789;

            c = UniformRNG.c8vec_uniform_01_new(N, ref seed);
            c[0] = new Complex(0.0, 0.0);

            e = UniformRNG.c8vec_uniform_01_new(N, ref seed);
            e[N - 1] = new Complex(0.0, 0.0);

            for (i = 1; i <= N; i++)
            {
                d[i - 1] = 0.0;
                if (i < N)
                {
                    d[i - 1] = d[i - 1] - new Complex(2.0, 0) * e[i - 1];
                }

                if (1 < i)
                {
                    d[i - 1] = d[i - 1] - new Complex(2.0, 0) * c[i - 1];
                }
            }

            //
            //  Set the desired solution
            //
            for (i = 1; i <= N; i++)
            {
                x[i - 1] = new Complex(i, 10 * i);
            }

            //
            //  Compute the corresponding right hand side.
            //
            b[0] = d[0] * x[0] + e[0] * x[1];
            for (i = 2; i <= N - 1; i++)
            {
                b[i - 1] = c[i - 1] * x[i - 2] + d[i - 1] * x[i - 1] + e[i - 1] * x[i];
            }

            b[N - 1] = c[N - 1] * x[N - 2] + d[N - 1] * x[N - 1];
            //
            //  Solve the tridiagonal system.
            //
            ZGTSL.zgtsl(N, ref c, ref d, ref e, ref b);

            Console.WriteLine("");
            Console.WriteLine("  Computed                     Exact");
            Console.WriteLine("  Solution                     Solution");
            Console.WriteLine("");
            for (i = 0; i < N; i++)
            {
                Console.WriteLine("  " + b[i].ToString().PadLeft(26)
                                       + "  " + x[i].ToString().PadLeft(26) + "");
            }
        }

        static void test11()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST11 tests ZHICO.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 3;

            Complex[] a = new Complex[N * N];
            int i;
            int[] ipvt = new int[N];
            int j;
            int lda;
            double rcond;
            int seed;
            string cout;

            lda = N;

            Console.WriteLine("");
            Console.WriteLine("TEST11");
            Console.WriteLine("  For a complex Hermitian matrix:");
            Console.WriteLine("  ZHICO factors the matrix and estimates");
            Console.WriteLine("  the reciprocal condition number.");
            Console.WriteLine("");
            Console.WriteLine("  The matrix order is N = " + N + "");
            //
            //  Set the values of the matrix A.
            //
            seed = 123456789;

            for (i = 0; i < N; i++)
            {
                a[i + i * lda] = new Complex(UniformRNG.r8_uniform_01(ref seed), 0.0);
                for (j = i + 1; j < N; j++)
                {
                    a[i + j * lda] = UniformRNG.c8_uniform_01(ref seed);
                    a[j + i * lda] = Complex.Conjugate(a[i + j * lda]);
                }
            }

            Console.WriteLine("");
            Console.WriteLine("  The matrix:");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                cout = "";
                for (j = 0; j < N; j++)
                {
                    cout += "  " + a[i + j * N].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }

            //
            //  Factor the matrix A.
            //
            rcond = ZHICO.zhico(ref a, lda, N, ref ipvt);

            Console.WriteLine("");
            Console.WriteLine("  Estimated reciprocal condition RCOND = " + rcond + "");

        }

        static void test12()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST12 tests ZHIFA and ZHISL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 3;

            Complex[] a = new Complex[N * N];
            Complex[] b = new Complex[N];
            int i;
            int info;
            int[] ipvt = new int[N];
            int j;
            int lda;
            int seed;
            Complex[] x;
            string cout;

            lda = N;

            Console.WriteLine("");
            Console.WriteLine("TEST12");
            Console.WriteLine("  For a complex Hermitian matrix:");
            Console.WriteLine("  ZHIFA factors the matrix.");
            Console.WriteLine("  ZHISL solves a linear system.");
            Console.WriteLine("");
            Console.WriteLine("  The matrix order is N = " + N + "");
            //
            //  Set the values of the matrix A.
            //
            seed = 123456789;

            for (i = 0; i < N; i++)
            {
                a[i + i * lda] = new Complex(UniformRNG.r8_uniform_01(ref seed), 0.0);
                for (j = i + 1; j < N; j++)
                {
                    a[i + j * lda] = UniformRNG.c8_uniform_01(ref seed);
                    a[j + i * lda] = Complex.Conjugate(a[i + j * lda]);
                }
            }

            Console.WriteLine("");
            Console.WriteLine("  The matrix:");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                cout = "";
                for (j = 0; j < N; j++)
                {
                    cout += "  " + a[i + j * N].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }

            //
            //  Set the values of the right hand side vector B.
            //
            x = UniformRNG.c8vec_uniform_01_new(N, ref seed);

            for (i = 0; i < N; i++)
            {
                b[i] = 0.0;
                for (j = 0; j < N; j++)
                {
                    b[i] = b[i] + a[i + j * N] * x[j];
                }
            }

            Console.WriteLine("");
            Console.WriteLine("  The right hand side:");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                Console.WriteLine("  " + b[i].ToString().PadLeft(20) + "");
            }

            //
            //  Factor the matrix A.
            //
            info = ZHIFA.zhifa(ref a, lda, N, ref ipvt);

            if (info != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("  ZHIFA returned an error flag INFO = " + info + "");
                return;
            }

            //
            //  Solve the system.
            //
            ZHISL.zhisl(a, lda, N, ipvt, ref b);

            Console.WriteLine("");
            Console.WriteLine("  Computed                     Exact");
            Console.WriteLine("  Solution                     Solution");
            Console.WriteLine("");
            for (i = 0; i < N; i++)
            {
                Console.WriteLine("  " + b[i].ToString().PadLeft(26)
                                       + "  " + x[i].ToString().PadLeft(26) + "");
            }
        }

        static void test13()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST13 tests ZHIFA and ZHIDI.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 3;

            Complex[] a = new Complex[N * N];
            Complex[] a_save = new Complex[N * N];
            Complex[] c = new Complex[N * N];
            double[] det = new double[2];
            int i;
            int[] inert = new int[3];
            int info;
            int[] ipvt = new int[N];
            int j;
            int job;
            int k;
            int lda;
            int seed;
            string cout;

            lda = N;

            Console.WriteLine("");
            Console.WriteLine("TEST13");
            Console.WriteLine("  For a complex hermitian matrix:");
            Console.WriteLine("  ZHIFA factors the matrix.");
            Console.WriteLine("  ZHIDI computes the determinant, inverse,");
            Console.WriteLine("  or inertia.");
            Console.WriteLine("");
            Console.WriteLine("  The matrix order is N = " + N + "");
            //
            //  Set the values of the matrix A.
            //
            seed = 123456789;

            for (i = 0; i < N; i++)
            {
                a[i + i * lda] = new Complex(UniformRNG.r8_uniform_01(ref seed), 0.0);
                for (j = i + 1; j < N; j++)
                {
                    a[i + j * lda] = UniformRNG.c8_uniform_01(ref seed);
                    a[j + i * lda] = Complex.Conjugate(a[i + j * lda]);
                }
            }

            for (i = 0; i < N; i++)
            {
                for (j = 0; j < N; j++)
                {
                    a_save[i + j * lda] = a[i + j * lda];
                }
            }

            Console.WriteLine("");
            Console.WriteLine("  The matrix:");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                cout = "";
                for (j = 0; j < N; j++)
                {
                    cout += "  " + a[i + j * lda].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }

            //
            //  Factor the matrix A.
            //
            info = ZHIFA.zhifa(ref a, lda, N, ref ipvt);

            if (info != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("  ZHIFA returned an error flag INFO = " + info + "");
                return;
            }

            //
            //  Get the determinant.
            //
            job = 10;
            ZHIDI.zhidi(ref a, lda, N, ipvt, ref det, ref inert, job);

            Console.WriteLine("");
            Console.WriteLine("  Determinant = " + det[0]
                                                 + " * 10^ " + det[1] + "");
            //
            //  Get the inertia.
            //
            job = 100;
            ZHIDI.zhidi(ref a, lda, N, ipvt, ref det, ref inert, job);

            Console.WriteLine("");
            Console.WriteLine("  The inertia:");
            Console.WriteLine("");

            for (i = 0; i < 3; i++)
            {
                Console.WriteLine("  " + inert[i] + "");
            }

            //
            //  Get the inverse.
            //
            job = 1;
            ZHIDI.zhidi(ref a, lda, N, ipvt, ref det, ref inert, job);
            //
            //  Only the upper triangle is set, so the user must set up the
            //  lower triangle:
            //
            for (i = 0; i < N; i++)
            {
                for (j = 0; j < i; j++)
                {
                    a[i + j * lda] = Complex.Conjugate(a[j + i * lda]);
                }
            }

            for (i = 0; i < N; i++)
            {
                for (k = 0; k < N; k++)
                {
                    c[i + k * lda] = 0.0;
                    for (j = 0; j < N; j++)
                    {
                        c[i + k * lda] = c[i + k * lda] + a[i + j * lda] * a_save[j + k * lda];
                    }
                }
            }

            Console.WriteLine("");
            Console.WriteLine("  The product inv(A) * A is");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                cout = "";
                for (j = 0; j < N; j++)
                {
                    cout += "  " + c[i + j * lda].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }
        }

        static void test14()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST14 tests ZHPCO.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 3;

            Complex[] a = new Complex[(N * (N + 1)) / 2];
            Complex[] a_save = new Complex[N * N];
            int i;
            int[] ipvt = new int[N];
            int j;
            int k;
            double rcond;
            int seed;
            string cout;

            Console.WriteLine("");
            Console.WriteLine("TEST14");
            Console.WriteLine("  For a complex Hermitian matrix");
            Console.WriteLine("  using packed storage,");
            Console.WriteLine("  ZHPCO factors the matrix and estimates");
            Console.WriteLine("  the reciprocal condition number.");
            Console.WriteLine("");
            Console.WriteLine("  The matrix order is N = " + N + "");
            //
            //  Set the values of the matrix A.
            //
            k = 0;
            seed = 123456789;

            for (j = 0; j < N; j++)
            {
                for (i = 0; i < j; i++)
                {
                    a[k] = UniformRNG.c8_uniform_01(ref seed);
                    a_save[i + j * N] = a[k];
                    a_save[j + i * N] = Complex.Conjugate(a[k]);
                    k = k + 1;
                }

                a[k] = new Complex(UniformRNG.r8_uniform_01(ref seed), 0.0);
                a_save[j + j * N] = a[k];
                k = k + 1;
            }

            Console.WriteLine("");
            Console.WriteLine("  The matrix:");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                cout = "";
                for (j = 0; j < N; j++)
                {
                    cout += "  " + a_save[i + j * N].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }

            //
            //  Factor the matrix A.
            //
            rcond = ZHPCO.zhpco(ref a, N, ref ipvt);

            Console.WriteLine("");
            Console.WriteLine("  Estimated reciprocal condition RCOND = " + rcond + "");

        }

        static void test15()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST15 tests ZHPFA and ZHPSL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 3;

            Complex[] a = new Complex[(N * (N + 1)) / 2];
            Complex[] a_save = new Complex[N * N];
            Complex[] b = new Complex[N];
            int i;
            int info;
            int[] ipvt = new int[N];
            int j;
            int k;
            int seed;
            Complex[] x;
            string cout;

            Console.WriteLine("");
            Console.WriteLine("TEST15");
            Console.WriteLine("  For a complex Hermitian matrix,");
            Console.WriteLine("  using packed storage,");
            Console.WriteLine("  ZHPFA factors the matrix.");
            Console.WriteLine("  ZHPSL solves a linear system.");
            Console.WriteLine("");
            Console.WriteLine("  The matrix order is N = " + N + "");
            //
            //  Set the values of the matrix A.
            //
            k = 0;
            seed = 123456789;

            for (j = 0; j < N; j++)
            {
                for (i = 0; i < j; i++)
                {
                    a[k] = UniformRNG.c8_uniform_01(ref seed);
                    a_save[i + j * N] = a[k];
                    a_save[j + i * N] = Complex.Conjugate(a[k]);
                    k = k + 1;
                }

                a[k] = new Complex(UniformRNG.r8_uniform_01(ref seed), 0.0);
                a_save[j + j * N] = a[k];
                k = k + 1;
            }

            Console.WriteLine("");
            Console.WriteLine("  The matrix:");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                cout = "";
                for (j = 0; j < N; j++)
                {
                    cout += "  " + a_save[i + j * N].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }

            //
            //  Set the values of the right hand side vector B.
            //
            x = UniformRNG.c8vec_uniform_01_new(N, ref seed);

            for (i = 0; i < N; i++)
            {
                b[i] = 0.0;
                for (j = 0; j < N; j++)
                {
                    b[i] = b[i] + a_save[i + j * N] * x[j];
                }
            }

            Console.WriteLine("");
            Console.WriteLine("  The right hand side:");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                Console.WriteLine("  " + b[i].ToString().PadLeft(20) + "");
            }

            //
            //  Factor the matrix A.
            //
            info = ZHPFA.zhpfa(ref a, N, ref ipvt);

            if (info != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("  ZHPFA returned an error flag INFO = " + info + "");
                return;
            }

            //
            //  Solve the system.
            //
            ZHPSL.zhpsl(a, N, ipvt, ref b);

            Console.WriteLine("");
            Console.WriteLine("  Computed                     Exact");
            Console.WriteLine("  Solution                     Solution");
            Console.WriteLine("");
            for (i = 0; i < N; i++)
            {
                Console.WriteLine("  " + b[i].ToString().PadLeft(26)
                                       + "  " + x[i].ToString().PadLeft(26) + "");
            }
        }

        static void test16()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST16 tests ZHPFA and ZHPDI.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 3;

            Complex[] a = new Complex[(N * (N + 1)) / 2];
            Complex[] a_save = new Complex[N * N];
            Complex[] b = new Complex[N * N];
            Complex[] c = new Complex[N * N];
            double[] det = new double[2];
            int i;
            int[] inert = new int[3];
            int info;
            int[] ipvt = new int[N];
            int j;
            int job;
            int k;
            int seed;
            string cout;

            Console.WriteLine("");
            Console.WriteLine("TEST16");
            Console.WriteLine("  For a complex hermitian matrix,");
            Console.WriteLine("  using packed storage,");
            Console.WriteLine("  ZHPFA factors the matrix.");
            Console.WriteLine("  ZHPDI computes the determinant, inverse,");
            Console.WriteLine("  or inertia.");
            Console.WriteLine("");
            Console.WriteLine("  The matrix order is N = " + N + "");
            //
            //  Set the values of the matrix A.
            //
            k = 0;
            seed = 123456789;

            for (j = 0; j < N; j++)
            {
                for (i = 0; i < j; i++)
                {
                    a[k] = UniformRNG.c8_uniform_01(ref seed);
                    a_save[i + j * N] = a[k];
                    a_save[j + i * N] = Complex.Conjugate(a[k]);
                    k = k + 1;
                }

                a[k] = new Complex(UniformRNG.r8_uniform_01(ref seed), 0.0);
                a_save[j + j * N] = a[k];
                k = k + 1;
            }

            Console.WriteLine("");
            Console.WriteLine("  The matrix:");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                cout = "";
                for (j = 0; j < N; j++)
                {
                    cout += "  " + a_save[i + j * N].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }

            //
            //  Factor the matrix A.
            //
            info = ZHPFA.zhpfa(ref a, N, ref ipvt);

            if (info != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("  ZHPFA returned an error flag INFO = " + info + "");
                return;
            }

            //
            //  Get the determinant.
            //
            job = 10;
            ZHPDI.zhpdi(ref a, N, ipvt, ref det, ref inert, job);

            Console.WriteLine("");
            Console.WriteLine("  Determinant = " + det[0]
                                                 + " * 10^ " + det[1] + "");
            //
            //  Get the inertia.
            //
            job = 100;
            ZHPDI.zhpdi(ref a, N, ipvt, ref det, ref inert, job);

            Console.WriteLine("");
            Console.WriteLine("  The inertia:");
            Console.WriteLine("");

            for (i = 0; i < 3; i++)
            {
                Console.WriteLine("  " + inert[i].ToString().PadLeft(8) + "");
            }

            //
            //  Get the inverse.
            //
            job = 1;
            ZHPDI.zhpdi(ref a, N, ipvt, ref det, ref inert, job);
            //
            //  Only the upper triangle is set, so the user must set up the
            //  lower triangle:
            //
            k = 0;
            for (j = 0; j < N; j++)
            {
                for (i = 0; i < j; i++)
                {
                    b[i + j * N] = a[k];
                    b[j + i * N] = Complex.Conjugate(a[k]);
                    k = k + 1;
                }

                b[j + j * N] = a[k];
                k = k + 1;
            }

            for (i = 0; i < N; i++)
            {
                for (k = 0; k < N; k++)
                {
                    c[i + k * N] = 0.0;
                    for (j = 0; j < N; j++)
                    {
                        c[i + k * N] = c[i + k * N] + b[i + j * N] * a_save[j + k * N];
                    }
                }
            }

            Console.WriteLine("");
            Console.WriteLine("  The product inv(A) * A is");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                cout = "";
                for (j = 0; j < N; j++)
                {
                    cout += "  " + c[i + j * N].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }
        }

        static void test17()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST17 tests ZPBCO.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 3;
            int M = 1;

            Complex[] a = new Complex[(M + 1) * N];
            int info = 0;
            int lda;
            double rcond;

            lda = M + 1;

            Console.WriteLine("");
            Console.WriteLine("TEST17");
            Console.WriteLine("  For a complex positive definite hermitian band matrix,");
            Console.WriteLine("  ZPBCO estimates the reciprocal condition number.");
            Console.WriteLine("  The matrix size is N = " + N + "");
            //
            //  Set the value of the superdiagonal and diagonal.
            //
            a[0 + 0 * lda] = new Complex(0.0000, 0.0000);
            a[0 + 1 * lda] = new Complex(2.1341, -0.2147);
            a[0 + 2 * lda] = new Complex(2.0905, 1.1505);

            a[1 + 0 * lda] = new Complex(4.5281, 0.0000);
            a[1 + 1 * lda] = new Complex(5.0371, 0.0000);
            a[1 + 2 * lda] = new Complex(4.7638, 0.0000);
            //
            //  Estimate the condition.
            //
            Console.WriteLine("");
            Console.WriteLine("  Estimate the condition.");

            rcond = ZPBCO.zpbco(a, lda, N, M, ref info);

            if (info != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("  ZPBCO returned INFO = " + info + "");
                Console.WriteLine("  The factorization was not completed.");
                return;
            }

            Console.WriteLine("");
            Console.WriteLine("  Reciprocal condition  = " + rcond + "");

# undef M
        }

        static void test18()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST18 tests ZPBDI.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 3;
            int M = 1;

            Complex[] a = new Complex[(M + 1) * N];
            double[] det = new double[2];
            int info;
            int lda;

            lda = M + 1;

            Console.WriteLine("");
            Console.WriteLine("TEST18");
            Console.WriteLine("  For a complex positive definite hermitian band matrix,");
            Console.WriteLine("  ZPBDI computes the determinant as");
            Console.WriteLine("    det = MANTISSA * 10^EXPONENT");
            Console.WriteLine("");
            //
            //  Set the value of the superdiagonal and diagonal.
            //
            a[0 + 0 * lda] = new Complex(0.0000, 0.0000);
            a[0 + 1 * lda] = new Complex(2.1341, -0.2147);
            a[0 + 2 * lda] = new Complex(2.0905, 1.1505);

            a[1 + 0 * lda] = new Complex(4.5281, 0.0000);
            a[1 + 1 * lda] = new Complex(5.0371, 0.0000);
            a[1 + 2 * lda] = new Complex(4.7638, 0.0000);

            info = ZPBFA.zpbfa(ref a, lda, N, M);

            if (info != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("  Error!  ZPBFA returns INFO = " + info + "");
                return;
            }

            ZPBDI.zpbdi(a, lda, N, M, ref det);

            Console.WriteLine("  Determinant = " + det[0]
                                                 + " * 10^ " + det[1] + "");
        }

        static void test19()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST19 tests ZPBFA and ZPBSL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 3;
            int M = 1;

            Complex[] a = new Complex[(M + 1) * N];
            Complex[] b = new Complex[N];
            int i;
            int info;
            int lda;

            lda = M + 1;

            Console.WriteLine("");
            Console.WriteLine("TEST19");
            Console.WriteLine("  For a complex positive definite hermitian band matrix,");
            Console.WriteLine("  ZPBFA computes the LU factors.");
            Console.WriteLine("  ZPBSL solves a factored linear system.");
            Console.WriteLine("  The matrix size is N = " + N + "");
            //
            //  Set the value of the superdiagonal and diagonal.
            //
            a[0 + 0 * lda] = new Complex(0.0000, 0.0000);
            a[0 + 1 * lda] = new Complex(2.1341, -0.2147);
            a[0 + 2 * lda] = new Complex(2.0905, 1.1505);

            a[1 + 0 * lda] = new Complex(4.5281, 0.0000);
            a[1 + 1 * lda] = new Complex(5.0371, 0.0000);
            a[1 + 2 * lda] = new Complex(4.7638, 0.0000);
            //
            //  Set the right hand side.
            //
            b[0] = new Complex(8.7963, -0.4294);
            b[1] = new Complex(18.4798, 3.6662);
            b[2] = new Complex(18.4724, -2.3010);
            //
            //  Factor the matrix.
            //
            Console.WriteLine("");
            Console.WriteLine("  Factor the matrix.");

            info = ZPBFA.zpbfa(ref a, lda, N, M);

            if (info != 0)
            {
                Console.WriteLine("  Error!  ZPBFA returns INFO = " + info + "");
                return;
            }

            //
            //  Solve the linear system.
            //
            Console.WriteLine("");
            Console.WriteLine("  Solve the linear system.");

            ZPBSL.zpbsl(a, lda, N, M, ref b);
            //
            //  Print the results.
            //
            Console.WriteLine("");
            Console.WriteLine("  The solution:");
            Console.WriteLine("  (Should be roughly (1,2,3)):");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                Console.WriteLine("  " + b[i].ToString().PadLeft(20) + "");
            }
        }

        static void test20()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST20 tests ZPOCO.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 3;

            Complex[] a = new Complex[N * N];
            int info = 0;
            int lda;
            double rcond;

            lda = N;

            Console.WriteLine("");
            Console.WriteLine("TEST20");
            Console.WriteLine("  For a complex Hermitian positive definite matrix,");
            Console.WriteLine("  ZPOCO estimates the reciprocal condition number.");
            Console.WriteLine("  The matrix size is N = " + N + "");
            //
            //  Set the values of the matrix A.
            //
            a[0 + 0 * 3] = new Complex(2.5281, 0.0000);
            a[1 + 0 * 3] = new Complex(2.1341, 0.2147);
            a[2 + 0 * 3] = new Complex(2.4187, -0.2932);

            a[0 + 1 * 3] = new Complex(2.1341, -0.2147);
            a[1 + 1 * 3] = new Complex(3.0371, 0.0000);
            a[2 + 1 * 3] = new Complex(2.0905, -1.1505);

            a[0 + 2 * 3] = new Complex(2.4187, 0.2932);
            a[1 + 2 * 3] = new Complex(2.0905, 1.1505);
            a[2 + 2 * 3] = new Complex(2.7638, 0.0000);
            //
            //  Estimate the condition.
            //
            Console.WriteLine("");
            Console.WriteLine("  Estimate the condition.");

            rcond = ZPOCO.zpoco(ref a, lda, N, ref info);

            Console.WriteLine("");
            Console.WriteLine("  Reciprocal condition  = " + rcond + "");

        }

        static void test21()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST21 tests ZPOFA and ZPODI.
            //
            //  Discussion:
            //
            //    ZPOFA factors a positive definite symmetric matrix,
            //    and ZPODI can compute the determinant or the inverse.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 3;

            Complex[] a = new Complex[N * N];
            double[] det = new double[2];
            int info;
            int j;
            int job;
            int lda;
            string cout;

            lda = N;

            Console.WriteLine("");
            Console.WriteLine("TEST21");
            Console.WriteLine("  For a complex Hermitian positive definite matrix,");
            Console.WriteLine("  ZPOFA computes the LU factors,");
            Console.WriteLine("  ZPODI computes the inverse or determinant.");
            Console.WriteLine("  The matrix size is N = " + N + "");
            //
            //  Set the values of the matrix A.
            //
            a[0 + 0 * 3] = new Complex(2.5281, 0.0000);
            a[1 + 0 * 3] = new Complex(2.1341, 0.2147);
            a[2 + 0 * 3] = new Complex(2.4187, -0.2932);

            a[0 + 1 * 3] = new Complex(2.1341, -0.2147);
            a[1 + 1 * 3] = new Complex(3.0371, 0.0000);
            a[2 + 1 * 3] = new Complex(2.0905, -1.1505);

            a[0 + 2 * 3] = new Complex(2.4187, 0.2932);
            a[1 + 2 * 3] = new Complex(2.0905, 1.1505);
            a[2 + 2 * 3] = new Complex(2.7638, 0.0000);
            //
            //  Factor the matrix.
            //
            Console.WriteLine("");
            Console.WriteLine("  Factor the matrix.");

            info = ZPOFA.zpofa(ref a, lda, N);

            if (info != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("  Error, ZPOFA returns INFO = " + info + "");
                return;
            }

            //
            //  Get the determinant and inverse.
            //
            Console.WriteLine("");
            Console.WriteLine("  Get the determinant and inverse.");

            job = 11;
            ZPODI.zpodi(ref a, lda, N, ref det, job);
            //
            //  Print the results.
            //
            Console.WriteLine("");
            Console.WriteLine("  Determinant  = " + det[0]
                                                  + " * 10^ " + det[1] + "");
            //
            //  ZPODI produces only the 'upper half triangle' of the inverse,
            //  which is actually symmetric.  Thus, the lower half could be
            //  produced by copying from the upper half.  However, the first row
            //  of A, as returned, is exactly the first row of the inverse.
            //
            Console.WriteLine("");
            Console.WriteLine("  First row of inverse:");
            Console.WriteLine("");
            cout = "";
            for (j = 0; j < N; j++)
            {
                cout += "  " + a[0 + j * N].ToString().PadLeft(20);
            }

            Console.WriteLine(cout);

        }

        static void test22()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST22 tests ZPOFA and ZPOSL.
            //
            //  Discussion:
            //
            //    ZPOFA factors a Hermitian positive definite matrix,
            //    and ZPOSL can solve a factored linear system.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 3;

            Complex[] a = new Complex[N * N];
            Complex[] b = new Complex[N];
            int i;
            int info;
            int j;
            int lda;
            Complex[] x = new Complex[N];

            lda = N;

            Console.WriteLine("");
            Console.WriteLine("TEST22");
            Console.WriteLine("  For a complex Hermitian positive definite matrix,");
            Console.WriteLine("  ZPOFA computes the LU factors.");
            Console.WriteLine("  ZPOSL solves a factored linear system.");
            Console.WriteLine("  The matrix size is N = " + N + "");
            //
            //  Set the values of the matrix A.
            //
            a[0 + 0 * 3] = new Complex(2.5281, 0.0000);
            a[1 + 0 * 3] = new Complex(2.1341, 0.2147);
            a[2 + 0 * 3] = new Complex(2.4187, -0.2932);

            a[0 + 1 * 3] = new Complex(2.1341, -0.2147);
            a[1 + 1 * 3] = new Complex(3.0371, 0.0000);
            a[2 + 1 * 3] = new Complex(2.0905, -1.1505);

            a[0 + 2 * 3] = new Complex(2.4187, 0.2932);
            a[1 + 2 * 3] = new Complex(2.0905, 1.1505);
            a[2 + 2 * 3] = new Complex(2.7638, 0.0000);
            //
            //  Set the right hand side.
            //
            for (i = 0; i < N; i++)
            {
                x[i] = new Complex(2 * i + 1, 2 * i + 2);
            }

            for (i = 0; i < N; i++)
            {
                b[i] = 0.0;
                for (j = 0; j < N; j++)
                {
                    b[i] = b[i] + a[i + j * N] * x[j];
                }
            }

            //
            //  Factor the matrix.
            //
            Console.WriteLine("");
            Console.WriteLine("  Factor the matrix.");

            info = ZPOFA.zpofa(ref a, lda, N);

            if (info != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("  Error, ZPOFA returns INFO = " + info + "");
                return;
            }

            //
            //  Solve the linear system.
            //
            Console.WriteLine("");
            Console.WriteLine("  Solve the linear system.");

            ZPOSL.zposl(a, lda, N, ref b);
            //
            //  Print the result.
            //
            Console.WriteLine("");
            Console.WriteLine("  The solution:");
            Console.WriteLine("  (Should be (1+2i),(3+4i),(5+6i):");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                Console.WriteLine("  " + b[i].ToString().PadLeft(20) + "");
            }

        }

        static void test23()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST23 tests ZPPCO.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 3;

            Complex[] a = new Complex[(N * (N + 1)) / 2];
            int info = 0;
            double rcond;

            Console.WriteLine("");
            Console.WriteLine("TEST23");
            Console.WriteLine("  For a complex Hermitian");
            Console.WriteLine("  positive definite packed matrix,");
            Console.WriteLine("  ZPPCO estimates the reciprocal condition number.");
            Console.WriteLine("  The matrix size is N = " + N + "");
            //
            //  Set the values of the matrix A.
            //
            a[0] = new Complex(2.5281, 0.0000);

            a[1] = new Complex(2.1341, -0.2147);
            a[2] = new Complex(3.0371, 0.0000);

            a[3] = new Complex(2.4187, 0.2932);
            a[4] = new Complex(2.0905, 1.1505);
            a[5] = new Complex(2.7638, 0.0000);
            //
            //  Estimate the condition.
            //
            Console.WriteLine("");
            Console.WriteLine("  Estimate the condition number.");

            rcond = ZPPCO.zppco(ref a, N, ref info);

            Console.WriteLine("");
            Console.WriteLine("  Reciprocal condition number = " + rcond + "");

        }

        static void test24()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST24 tests ZPPFA and ZPPDI.
            //
            //  Discussion:
            //
            //    ZPPFA factors a Hermitian positive definite packed matrix,
            //    and ZPPDI can compute the determinant or the inverse.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 3;

            Complex[] a = new Complex[(N * (N + 1)) / 2];
            Complex[] b = new Complex[N * N];
            double[] det = new double[2];
            int i;
            int info;
            int j;
            int job;
            int k;
            string cout;

            Console.WriteLine("");
            Console.WriteLine("TEST24");
            Console.WriteLine("  For a complex Hermitian");
            Console.WriteLine("  positive definite packed matrix,");
            Console.WriteLine("  ZPPFA factors the matrix.");
            Console.WriteLine("  ZPPDI computes the inverse or determinant.");
            Console.WriteLine("  The matrix size is N = " + N + "");
            //
            //  Set the values of the matrix A.
            //
            a[0] = new Complex(2.5281, 0.0000);

            a[1] = new Complex(2.1341, -0.2147);
            a[2] = new Complex(3.0371, 0.0000);

            a[3] = new Complex(2.4187, 0.2932);
            a[4] = new Complex(2.0905, 1.1505);
            a[5] = new Complex(2.7638, 0.0000);
            //
            //  Factor the matrix.
            //
            Console.WriteLine("");
            Console.WriteLine("  Factor the matrix.");

            info = ZPPFA.zppfa(ref a, N);

            if (info != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("  Error, ZPPFA returns INFO = " + info + "");
                return;
            }

            //
            //  Invert the matrix.
            //
            Console.WriteLine("");
            Console.WriteLine("  Get the determinant and inverse.");

            job = 11;
            ZPPDI.zppdi(ref a, N, ref det, job);
            //
            //  Print the results.
            //
            Console.WriteLine("");
            Console.WriteLine("  Determinant  = " + det[0]
                                                  + " * 10^ " + det[1] + "");
            //
            //  ZPPDI produces only the 'upper half triangle' of the inverse,
            //  which is actually symmetric.  Thus, the lower half could be
            //  produced by copying from the upper half.
            //
            k = 0;
            for (j = 0; j < N; j++)
            {
                for (i = 0; i <= j; i++)
                {
                    b[i + j * N] = a[k];
                    b[j + i * N] = Complex.Conjugate(a[k]);
                    k = k + 1;
                }
            }

            Console.WriteLine("");
            Console.WriteLine("  Inverse:");
            Console.WriteLine("");
            for (i = 0; i < N; i++)
            {
                cout = "";
                for (j = 0; j < N; j++)
                {
                    cout += "  " + b[i + j * N].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }

        }

        static void test25()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST25 tests ZPPFA and ZPPSL.
            //
            //  Discussion:
            //
            //    ZPOFA factors a Hermitian positive definite packed matrix,
            //    and ZPOSL can solve a factored linear system.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 3;

            Complex[] a = new Complex[(N * (N + 1)) / 2];
            Complex[] b = new Complex[N];
            int i;
            int info;

            Console.WriteLine("");
            Console.WriteLine("TEST25");
            Console.WriteLine("  For a complex Hermitian");
            Console.WriteLine("  positive definite packed matrix,");
            Console.WriteLine("  ZPPFA factors the matrix.");
            Console.WriteLine("  ZPPSL solves a factored linear system.");
            Console.WriteLine("  The matrix size is N = " + N + "");
            //
            //  Set the values of the matrix A.
            //
            a[0] = new Complex(2.5281, 0.0000);

            a[1] = new Complex(2.1341, -0.2147);
            a[2] = new Complex(3.0371, 0.0000);

            a[3] = new Complex(2.4187, 0.2932);
            a[4] = new Complex(2.0905, 1.1505);
            a[5] = new Complex(2.7638, 0.0000);

            b[0] = new Complex(20.12350, 28.92670);
            b[1] = new Complex(14.36550, 34.92680);
            b[2] = new Complex(27.69760, 26.03750);
            //
            //  Factor the matrix.
            //
            Console.WriteLine("");
            Console.WriteLine("  Factor the matrix.");

            info = ZPPFA.zppfa(ref a, N);

            if (info != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("  Error, ZPPFA returns INFO = " + info + "");
                return;
            }

            //
            //  Solve the linear system.
            //
            Console.WriteLine("");
            Console.WriteLine("  Solve the linear system.");

            ZPPSL.zppsl(a, N, ref b);
            //
            //  Print the result.
            //
            Console.WriteLine("");
            Console.WriteLine("  The solution:");
            Console.WriteLine("  (Should be (1+2i),(3+4i),(5+6i):");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                Console.WriteLine("  " + b[i].ToString().PadLeft(20) + "");
            }

        }

        static void test26()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST26 tests ZPTSL.
            //
            //  Discussion:
            //
            //    ZPTSL factors and solves a Hermitian positive definite
            //    tridiagonal system.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 3;

            Complex[] b = new Complex[N];
            Complex[] d = new Complex[N];
            Complex[] e = new Complex[N];
            int i;

            Console.WriteLine("");
            Console.WriteLine("TEST26");
            Console.WriteLine("  For a complex Hermitian");
            Console.WriteLine("  positive definite tridiagonal matrix,");
            Console.WriteLine("  ZPTSL factors and solves a linear system.");
            Console.WriteLine("  The matrix size is N = " + N + "");
            //
            //  Set the value of the superdiagonal and diagonal.
            //
            e[0] = new Complex(2.1341, -0.2147);
            e[1] = new Complex(2.0905, 1.1505);
            e[2] = new Complex(0.0000, 0.0000);

            d[0] = new Complex(4.5281, 0.0000);
            d[1] = new Complex(5.0371, 0.0000);
            d[2] = new Complex(4.7638, 0.0000);
            //
            //  Set the right hand side.
            //
            b[0] = new Complex(8.7963, -0.4294);
            b[1] = new Complex(18.4798, 3.6662);
            b[2] = new Complex(18.4724, -2.3010);
            //
            //  Factor and solve the system.
            //
            Console.WriteLine("");
            Console.WriteLine("  Factor the matrix and solve the system.");

            ZPTSL.zptsl(N, ref d, ref e, ref b);
            //
            //  Print the result.
            //
            Console.WriteLine("");
            Console.WriteLine("  The solution:");
            Console.WriteLine("  (Should be roughly (1,2,3)):");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                Console.WriteLine("  " + b[i].ToString().PadLeft(26) + "");
            }

        }

        static void zqrdc_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ZQRDC_TEST tests ZQRDC.
            //
            //  Discussion:
            //
            //    ZQRDC computes the QR factorization.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 3;
            int P = 3;

            Complex[] a;
            Complex[] b = new Complex[N * P];
            int i;
            int info;
            int[] ipvt = new int[P];
            int j;
            int job;
            int k;
            int lda;
            Complex[] q = new Complex[N * N];
            Complex[] qraux = new Complex[P];
            Complex[] qty = new Complex[N];
            Complex[] qy = new Complex[N];
            Complex[] r = new Complex[N * P];
            Complex[] rsd = new Complex[N];
            int seed;
            Complex[] xb = new Complex[N];
            Complex[] y = new Complex[N];
            string cout;

            lda = N;

            Console.WriteLine("");
            Console.WriteLine("ZQRDC_TEST");
            Console.WriteLine("  ZQRDC computes the QR decomposition of a rectangular");
            Console.WriteLine("  matrix, but does not return Q and R explicitly.");
            Console.WriteLine("");
            Console.WriteLine("  Show how Q and R can be recovered using ZQRSL.");
            //
            //  Set the values of the matrix A.
            //
            seed = 123456789;

            a = UniformRNG.c8mat_uniform_01_new(N, P, ref seed);

            Console.WriteLine("");
            Console.WriteLine("  The matrix A is");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                cout = "";
                for (j = 0; j < P; j++)
                {
                    cout += "  " + a[i + j * lda].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }

            //
            //  Decompose the matrix.
            //
            Console.WriteLine("");
            Console.WriteLine("  Decompose the matrix.");

            job = 0;
            for (i = 0; i < P; i++)
            {
                ipvt[i] = 0;
            }

            ZQRDC.zqrdc(ref a, lda, N, P, ref qraux, ref ipvt, job);
            //
            //  Print out what ZQRDC has stored in A...
            //
            Console.WriteLine("");
            Console.WriteLine("  The packed matrix A which describes Q and R:");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                cout = "";
                for (j = 0; j < P; j++)
                {
                    cout += "  " + a[i + j * lda].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }

            //
            //  ...and in QRAUX.
            //
            Console.WriteLine("");
            Console.WriteLine("  The QRAUX vector, containing some additional");
            Console.WriteLine("  information defining Q:");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                Console.WriteLine("  " + qraux[i].ToString().PadLeft(20) + "");
            }

            //
            //  Print out the resulting R factor.
            //
            for (i = 0; i < N; i++)
            {
                for (j = 0; j < P; j++)
                {
                    if (j < i)
                    {
                        r[i + j * lda] = new Complex(0.0, 0.0);
                    }
                    else
                    {
                        r[i + j * lda] = a[i + j * lda];
                    }
                }
            }

            Console.WriteLine("");
            Console.WriteLine("  The R factor:");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                cout = "";
                for (j = 0; j < P; j++)
                {
                    cout += "  " + r[i + j * N].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }

            //
            //  Call ZQRSL to extract the information about the Q matrix.
            //  We do this, essentially, by asking ZQRSL to tell us the
            //  value of Q*Y, where Y is a column of the identity matrix.
            //
            job = 10000;

            for (j = 0; j < N; j++)
            {
                //
                //  Set the vector Y.
                //
                for (i = 0; i < N; i++)
                {
                    y[i] = new Complex(0.0, 0.0);
                }

                y[j] = new Complex(1.0, 0.0);
                //
                //  Ask ZQRSL to tell us what Q*Y is.
                //
                info = ZQRSL.zqrsl(a, lda, N, P, qraux, ref y, ref qy, ref qty, ref b, ref rsd, ref xb, job);

                if (info != 0)
                {
                    Console.WriteLine("  Error!  ZQRSL returns INFO = " + info + "");
                    return;
                }

                //
                //  Copy QY into the appropriate column of Q.
                //
                for (i = 0; i < N; i++)
                {
                    q[i + j * N] = qy[i];
                }
            }

            //
            //  Now print out the Q matrix we have extracted.
            //
            Console.WriteLine("");
            Console.WriteLine("  The Q factor:");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                cout = "";
                for (j = 0; j < N; j++)
                {
                    cout += "  " + q[i + j * N].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }

            //
            //  Compute Q*R to verify that it equals A.
            //
            for (i = 0; i < N; i++)
            {
                for (j = 0; j < P; j++)
                {
                    b[i + j * N] = new Complex(0.0, 0.0);
                    for (k = 0; k < N; k++)
                    {
                        b[i + j * N] = b[i + j * N] + q[i + k * N] * r[k + j * N];
                    }
                }

                Console.WriteLine("");
            }

            //
            //  Print the result.
            //
            Console.WriteLine("");
            Console.WriteLine("  The product Q * R:");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                cout = "";
                for (j = 0; j < P; j++)
                {
                    cout += "  " + b[i + j * N].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }
        }

        static void test28()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST28 tests ZSICO.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 3;

            Complex[] a = new Complex[N * N];
            int i;
            int[] ipvt = new int[N];
            int j;
            int lda;
            double rcond;
            int seed;
            string cout;

            lda = N;

            Console.WriteLine("");
            Console.WriteLine("TEST28");
            Console.WriteLine("  For a complex symmetric matrix:");
            Console.WriteLine("  ZSICO factors the matrix and estimates");
            Console.WriteLine("  the reciprocal condition number.");
            Console.WriteLine("");
            Console.WriteLine("  The matrix order is N = " + N + "");
            //
            //  Set the values of the matrix A.
            //
            seed = 123456789;

            for (i = 0; i < N; i++)
            {
                a[i + i * lda] = UniformRNG.c8_uniform_01(ref seed);
                for (j = i + 1; j < N; j++)
                {
                    a[i + j * lda] = UniformRNG.c8_uniform_01(ref seed);
                    a[j + i * lda] = a[i + j * lda];
                }
            }

            Console.WriteLine("");
            Console.WriteLine("  The matrix:");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                cout = "";
                for (j = 0; j < N; j++)
                {
                    cout += "  " + a[i + j * lda].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }

            //
            //  Factor the matrix A.
            //
            rcond = ZSICO.zsico(ref a, lda, N, ref ipvt);

            Console.WriteLine("");
            Console.WriteLine("  Estimated reciprocal condition RCOND = " + rcond + "");

        }

        static void test29()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST29 tests ZSIFA and ZSISL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 3;

            Complex[] a = new Complex[N * N];
            Complex[] b = new Complex[N];
            int i;
            int info;
            int[] ipvt = new int[N];
            int j;
            int lda;
            int seed;
            Complex[] x = new Complex[N];
            string cout;

            lda = N;

            Console.WriteLine("");
            Console.WriteLine("TEST29");
            Console.WriteLine("  For a complex symmetric matrix:");
            Console.WriteLine("  ZSIFA factors the matrix.");
            Console.WriteLine("  ZSISL solves a linear system.");
            Console.WriteLine("");
            Console.WriteLine("  The matrix order is N = " + N + "");
            //
            //  Set the values of the matrix A.
            //
            seed = 123456789;

            for (i = 0; i < N; i++)
            {
                a[i + i * lda] = UniformRNG.c8_uniform_01(ref seed);
                for (j = i + 1; j < N; j++)
                {
                    a[i + j * lda] = UniformRNG.c8_uniform_01(ref seed);
                    a[j + i * lda] = a[i + j * lda];
                }
            }

            Console.WriteLine("");
            Console.WriteLine("  The matrix:");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                cout = "";
                for (j = 0; j < N; j++)
                {
                    cout += "  " + a[i + j * lda].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }

            //
            //  Set the values of the right hand side vector B.
            //
            for (i = 0; i < N; i++)
            {
                x[i] = UniformRNG.c8_uniform_01(ref seed);
            }

            for (i = 0; i < N; i++)
            {
                b[i] = 0.0;
                for (j = 0; j < N; j++)
                {
                    b[i] = b[i] + a[i + j * lda] * x[j];
                }
            }

            Console.WriteLine("");
            Console.WriteLine("  The right hand side:");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                Console.WriteLine("  " + b[i].ToString().PadLeft(20) + "");
            }

            //
            //  Factor the matrix A.
            //
            info = ZSIFA.zsifa(ref a, lda, N, ref ipvt);

            if (info != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("  ZSIFA returned an error flag INFO = " + info + "");
                return;
            }

            //
            //  Solve the system.
            //
            ZSISL.zsisl(a, lda, N, ipvt, ref b);

            Console.WriteLine("");
            Console.WriteLine("  Computed                     Exact");
            Console.WriteLine("  Solution                     Solution");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                Console.WriteLine("  " + b[i].ToString().PadLeft(20)
                                       + "  " + x[i].ToString().PadLeft(20) + "");
                ;
            }

        }

        static void test30()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST30 tests ZSIFA and ZSIDI.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 3;

            Complex[] a = new Complex[N * N];
            Complex[] a_save = new Complex[N * N];
            Complex[] c = new Complex[N * N];
            Complex[] det = new Complex[2];
            int i;
            int info;
            int[] ipvt = new int[N];
            int j;
            int job;
            int k;
            int lda;
            int seed;
            string cout;

            lda = N;

            Console.WriteLine("");
            Console.WriteLine("TEST30");
            Console.WriteLine("  For a complex symmetric matrix:");
            Console.WriteLine("  ZSIFA factors the matrix.");
            Console.WriteLine("  ZSIDI computes the determinant or inverse.");
            Console.WriteLine("");
            Console.WriteLine("  The matrix order is N = " + N + "");
            //
            //  Set the values of the matrix A.
            //
            seed = 123456789;

            for (i = 0; i < N; i++)
            {
                a[i + i * lda] = UniformRNG.c8_uniform_01(ref seed);
                for (j = i + 1; j < N; j++)
                {
                    a[i + j * lda] = UniformRNG.c8_uniform_01(ref seed);
                    a[j + i * lda] = a[i + j * lda];
                }
            }

            for (i = 0; i < N; i++)
            {
                for (j = 0; j < N; j++)
                {
                    a_save[i + j * lda] = a[i + j * lda];
                }
            }

            Console.WriteLine("");
            Console.WriteLine("  The matrix:");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                cout = "";
                for (j = 0; j < N; j++)
                {
                    cout += "  " + a[i + j * lda].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }

            //
            //  Factor the matrix A.
            //
            info = ZSIFA.zsifa(ref a, lda, N, ref ipvt);

            if (info != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("  ZSIFA returned an error flag INFO = " + info + "");
                return;
            }

            //
            //  Get the determinant.
            //
            job = 10;
            ZSIDI.zsidi(ref a, lda, N, ipvt, ref det, job);

            Console.WriteLine("");
            Console.WriteLine("  Determinant = " + det[0]
                                                 + " * 10^ " + det[1] + "");
            //
            //  Get the inverse.
            //
            job = 1;
            ZSIDI.zsidi(ref a, lda, N, ipvt, ref det, job);
            //
            //  Only the upper triangle is set, so the user must set up the
            //  lower triangle:
            //
            for (i = 0; i < N; i++)
            {
                for (j = 0; j < i; j++)
                {
                    a[i + j * lda] = a[j + i * lda];
                }
            }

            for (i = 0; i < N; i++)
            {
                for (k = 0; k < N; k++)
                {
                    c[i + k * N] = 0.0;
                    for (j = 0; j < N; j++)
                    {
                        c[i + k * N] = c[i + k * N] + a[i + j * N] * a_save[j + k * N];
                    }
                }
            }

            Console.WriteLine("");
            Console.WriteLine("  The product inv(A) * A is");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                cout = "";
                for (j = 0; j < N; j++)
                {
                    cout += "  " + c[i + j * N].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }

        }

        static void test31()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST31 tests ZSPCO.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 3;

            Complex[] a = new Complex[(N * (N + 1)) / 2];
            Complex[] a_save = new Complex[N * N];
            int i;
            int[] ipvt = new int[N];
            int j;
            int k;
            double rcond;
            int seed;
            string cout;

            Console.WriteLine("");
            Console.WriteLine("TEST31");
            Console.WriteLine("  For a complex symmetric matrix");
            Console.WriteLine("  in packed storage,");
            Console.WriteLine("  ZSPCO factors the matrix and estimates");
            Console.WriteLine("  the reciprocal condition number.");
            Console.WriteLine("");
            Console.WriteLine("  The matrix order is N = " + N + "");
            //
            //  Set the values of the packed matrix A.
            //
            k = 0;
            seed = 123456789;

            for (j = 0; j < N; j++)
            {
                for (i = 0; i < j; i++)
                {
                    a[k] = UniformRNG.c8_uniform_01(ref seed);
                    k = k + 1;
                }

                a[k] = UniformRNG.c8_uniform_01(ref seed);
                k = k + 1;
            }

            //
            //  Copy the packed matrix into a "normal" matrix.
            //
            k = 0;
            for (j = 0; j < N; j++)
            {
                for (i = 0; i <= j; i++)
                {
                    a_save[i + j * N] = a[k];
                    k = k + 1;
                }
            }

            for (j = 0; j < N; j++)
            {
                for (i = j + 1; i < N; i++)
                {
                    a_save[i + j * N] = a_save[j + i * N];
                }
            }

            Console.WriteLine("");
            Console.WriteLine("  The matrix:");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                cout = "";
                for (j = 0; j < N; j++)
                {
                    cout += "  " + a_save[i + j * N].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }

            //
            //  Factor the matrix A.
            //
            rcond = ZSPCO.zspco(ref a, N, ref ipvt);

            Console.WriteLine("");
            Console.WriteLine("  Estimated reciprocal condition RCOND = " + rcond + "");

        }

        static void test32()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST32 tests ZSPFA and ZSPSL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 3;

            Complex[] a = new Complex[(N * (N + 1)) / 2];
            Complex[] a_save = new Complex[N * N];
            Complex[] b = new Complex[N];
            int i;
            int info;
            int[] ipvt = new int[N];
            int j;
            int k;
            int seed;
            Complex[] x;
            string cout;

            Console.WriteLine("");
            Console.WriteLine("TEST32");
            Console.WriteLine("  For a complex symmetric matrix");
            Console.WriteLine("  in packed storage,");
            Console.WriteLine("  ZSPFA factors the matrix.");
            Console.WriteLine("  ZSPSL solves a linear system.");
            Console.WriteLine("");
            Console.WriteLine("  The matrix order is N = " + N + "");
            //
            //  Set the values of the packed matrix A.
            //
            k = 0;
            seed = 123456789;

            for (j = 0; j < N; j++)
            {
                for (i = 0; i < j; i++)
                {
                    a[k] = UniformRNG.c8_uniform_01(ref seed);
                    k = k + 1;
                }

                a[k] = UniformRNG.c8_uniform_01(ref seed);
                k = k + 1;
            }

            //
            //  Copy the packed matrix into a "normal" matrix.
            //
            k = 0;
            for (j = 0; j < N; j++)
            {
                for (i = 0; i <= j; i++)
                {
                    a_save[i + j * N] = a[k];
                    k = k + 1;
                }
            }

            for (j = 0; j < N; j++)
            {
                for (i = j + 1; i < N; i++)
                {
                    a_save[i + j * N] = a_save[j + i * N];
                }
            }

            Console.WriteLine("");
            Console.WriteLine("  The matrix:");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                cout = "";
                for (j = 0; j < N; j++)
                {
                    cout += "  " + a_save[i + j * N].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }

            //
            //  Set the values of the right hand side vector B.
            //
            x = UniformRNG.c8vec_uniform_01_new(N, ref seed);

            for (i = 0; i < N; i++)
            {
                b[i] = 0.0;
                for (j = 0; j < N; j++)
                {
                    b[i] = b[i] + a_save[i + j * N] * x[j];
                }
            }

            Console.WriteLine("");
            Console.WriteLine("  The right hand side:");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                Console.WriteLine("  " + b[i].ToString().PadLeft(20) + "");
            }

            //
            //  Factor the matrix A.
            //
            info = ZSPFA.zspfa(ref a, N, ref ipvt);

            if (info != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("  ZSPFA returned an error flag INFO = " + info + "");
                return;
            }

            //
            //  Solve the system.
            //
            ZSPSL.zspsl(a, N, ipvt, ref b);

            Console.WriteLine("");
            Console.WriteLine("  Computed                     Exact");
            Console.WriteLine("  Solution                     Solution");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                Console.WriteLine("  " + b[i].ToString().PadLeft(20)
                                       + "  " + x[i].ToString().PadLeft(20) + "");
                ;
            }
        }

        static void test33()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST33 tests ZSPFA and ZSPDI.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 3;

            Complex[] a = new Complex[(N * (N + 1)) / 2];
            Complex[] a_save = new Complex[N * N];
            Complex[] b_save = new Complex[N * N];
            Complex[] c = new Complex[N * N];
            Complex[] det = new Complex[2];
            int i;
            int info;
            int[] ipvt = new int[N];
            int j;
            int job;
            int k;
            int seed;
            string cout;

            Console.WriteLine("");
            Console.WriteLine("TEST33");
            Console.WriteLine("  For a complex symmetric matrix");
            Console.WriteLine("  in packed storage,");
            Console.WriteLine("  ZSPFA factors the matrix.");
            Console.WriteLine("  ZSPDI computes the determinant or inverse.");
            Console.WriteLine("");
            Console.WriteLine("  The matrix order is N = " + N + "");
            //
            //  Set the values of the packed matrix A.
            //
            k = 0;
            seed = 123456789;

            for (j = 0; j < N; j++)
            {
                for (i = 0; i < j; i++)
                {
                    a[k] = UniformRNG.c8_uniform_01(ref seed);
                    k = k + 1;
                }

                a[k] = UniformRNG.c8_uniform_01(ref seed);
                k = k + 1;
            }

            //
            //  Copy the packed matrix into a "normal" matrix.
            //
            k = 0;
            for (j = 0; j < N; j++)
            {
                for (i = 0; i <= j; i++)
                {
                    a_save[i + j * N] = a[k];
                    k = k + 1;
                }
            }

            for (j = 0; j < N; j++)
            {
                for (i = j + 1; i < N; i++)
                {
                    a_save[i + j * N] = a_save[j + i * N];
                }
            }

            Console.WriteLine("");
            Console.WriteLine("  The matrix:");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                cout = "";
                for (j = 0; j < N; j++)
                {
                    cout += "  " + a_save[i + j * N].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }

            //
            //  Factor the matrix A.
            //
            info = ZSPFA.zspfa(ref a, N, ref ipvt);

            if (info != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("  ZSPFA returned an error flag INFO = " + info + "");
                return;
            }

            //
            //  Get the determinant.
            //
            job = 10;
            ZSPDI.zspdi(ref a, N, ipvt, ref det, job);

            Console.WriteLine("");
            Console.WriteLine("  Determinant = " + det[0]
                                                 + " * 10^ " + det[1] + "");
            //
            //  Get the inverse.
            //
            job = 1;
            ZSPDI.zspdi(ref a, N, ipvt, ref det, job);
            //
            //  Copy the packed matrix into a "normal" matrix.
            //
            k = 0;
            for (j = 0; j < N; j++)
            {
                for (i = 0; i <= j; i++)
                {
                    b_save[i + j * N] = a[k];
                    k = k + 1;
                }
            }

            for (j = 0; j < N; j++)
            {
                for (i = j + 1; i < N; i++)
                {
                    b_save[i + j * N] = b_save[j + i * N];
                }
            }

            for (i = 0; i < N; i++)
            {
                for (k = 0; k < N; k++)
                {
                    c[i + k * N] = 0.0;
                    for (j = 0; j < N; j++)
                    {
                        c[i + k * N] = c[i + k * N] + a_save[i + j * N] * b_save[j + k * N];
                    }
                }
            }

            Console.WriteLine("");
            Console.WriteLine("  The product inv(A) * A is");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                cout = "";
                for (j = 0; j < N; j++)
                {
                    cout += "  " + c[i + j * N].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }

        }

        static void test34()

            //**********************************************************************
            //
            //  Purpose:
            //
            //    TEST34 tests ZSVDC.
            //
            //  Discussion:
            //
            //    ZSVDC computes the singular value decomposition:
            //
            //      A = U * S * Complex.Conjugateg-transpose ( V )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int M = 4;
            int N = 3;

            Complex[] a;
            //
            //  E should be dimensioned at least the maximum of M+1 and N.
            //
            Complex[] e = new Complex[M + N];
            int i;
            int info;
            int j;
            int lda;
            int ldu;
            int ldv;
            int job;
            int k;
            //
            //  S should be dimensioned at least the maximum of M+1 and N.
            //
            Complex[] s = new Complex[M + N];
            int seed;
            Complex[] u = new Complex[M * M];
            Complex[] v = new Complex[N * N];
            string cout;

            Console.WriteLine("");
            Console.WriteLine("TEST34");
            Console.WriteLine("  For an MxN matrix A in complex general storage,");
            Console.WriteLine("  ZSVDC computes the singular value decomposition:");
            Console.WriteLine("    A = U * S * V^H");
            Console.WriteLine("");
            Console.WriteLine("  Matrix rows M =    " + M + "");
            Console.WriteLine("  Matrix columns N = " + N + "");
            //
            //  Set A.
            //
            seed = 123456789;
            lda = M;

            a = UniformRNG.c8mat_uniform_01_new(M, N, ref seed);

            Console.WriteLine("");
            Console.WriteLine("  The matrix A:");
            Console.WriteLine("");

            for (i = 0; i < M; i++)
            {
                cout = "";
                for (j = 0; j < N; j++)
                {
                    cout += "  " + a[i + j * lda].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }

            //
            //  Decompose the matrix.
            //
            Console.WriteLine("");
            Console.WriteLine("  Decompose the matrix.");

            job = 11;
            ldu = M;
            ldv = N;

            info = ZSVDC.zsvdc(ref a, lda, M, N, ref s, ref e, ref u, ldu, ref v, ldv, job);

            if (info != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("Warning:");
                Console.WriteLine("  ZSVDC returned nonzero INFO = " + info + "");
                return;
            }

            Console.WriteLine("");
            Console.WriteLine("  Singular values:");
            Console.WriteLine("");

            for (i = 0; i < Math.Min(M, N); i++)
            {
                Console.WriteLine("  " + (i + 1).ToString().PadLeft(8)
                                       + s[i].ToString().PadLeft(20) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("  Left Singular Vector Matrix U:");
            Console.WriteLine("");

            for (i = 0; i < M; i++)
            {
                cout = "";
                for (j = 0; j < M; j++)
                {
                    cout += "  " + u[i + j * ldu].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }

            Console.WriteLine("");
            Console.WriteLine("  Right Singular Vector Matrix V:");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                cout = "";
                for (j = 0; j < N; j++)
                {
                    cout += "  " + v[i + j * ldv].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }

            for (i = 0; i < M; i++)
            {
                for (j = 0; j < N; j++)
                {
                    a[i + j * lda] = new Complex(0.0, 0.0);
                    for (k = 0; k < Math.Min(M, N); k++)
                    {
                        a[i + j * lda] = a[i + j * lda] + u[i + k * ldu] * s[k] * Complex.Conjugate(v[j + k * ldv]);
                    }
                }
            }

            Console.WriteLine("");
            Console.WriteLine("  The product U * S * V^H (should equal A):");
            Console.WriteLine("");

            for (i = 0; i < M; i++)
            {
                cout = "";
                for (j = 0; j < N; j++)
                {
                    cout += "  " + a[i + j * lda].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }
        }

        static void test345()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST345 tests ZSVDC.
            //
            //  Discussion:
            //
            //    ZSVDC computes the singular value decomposition:
            //
            //      A = U * S * Complex.Conjugateg-transpose ( V )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 January 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int M = 4;
            int N = 4;

            Complex[] a;
            //
            //  E must be dimensioned at least the maximum of M+1 and N.
            //
            Complex[] e = new Complex[M + N];
            Complex I;
            int i;
            int info;
            int j;
            int lda;
            int ldu;
            int ldv;
            int job;
            int k;
            //
            //  S must be dimensioned at least the maximum of M+1 and N.
            //
            Complex[] s = new Complex[M + N];
            Complex[] u = new Complex[M * M];
            Complex[] v = new Complex[N * N];

            string cout;

            Console.WriteLine("");
            Console.WriteLine("TEST345");
            Console.WriteLine("  For an MxN matrix A in double complex general storage,");
            Console.WriteLine("  ZSVDC computes the singular value decomposition:");
            Console.WriteLine("    A = U * S * V^H");
            Console.WriteLine("");
            Console.WriteLine("  Matrix rows M =    " + M + "");
            Console.WriteLine("  Matrix columns N = " + N + "");
            //
            //  Set A.
            //
            I = new Complex(0.0, 1.0);

            lda = M;
            a = new Complex[M * N];

            a[0 + 0 * M] = 1.0;
            a[1 + 0 * M] = -I;
            a[2 + 0 * M] = -1.0;
            a[3 + 0 * M] = I;

            a[0 + 1 * M] = 1.0;
            a[1 + 1 * M] = -1.0;
            a[2 + 1 * M] = -1.0;
            a[3 + 1 * M] = 1.0;

            a[0 + 2 * M] = 1.0;
            a[1 + 2 * M] = 1.0;
            a[2 + 2 * M] = 1.0;
            a[3 + 2 * M] = 1.0;

            a[0 + 3 * M] = 1.0;
            a[1 + 3 * M] = I;
            a[2 + 3 * M] = -1.0;
            a[3 + 3 * M] = -I;

            Console.WriteLine("");
            Console.WriteLine("  The matrix A:");
            Console.WriteLine("");

            for (i = 0; i < M; i++)
            {
                cout = "";
                for (j = 0; j < N; j++)
                {
                    cout += "  " + a[i + j * lda].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }

            //
            //  Decompose the matrix.
            //
            Console.WriteLine("");
            Console.WriteLine("  Decompose the matrix.");

            job = 11;
            ldu = M;
            ldv = N;

            info = ZSVDC.zsvdc(ref a, lda, M, N, ref s, ref e, ref u, ldu, ref v, ldv, job);

            if (info != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("Warning:");
                Console.WriteLine("  CSVDC returned nonzero INFO = " + info + "");
                return;
            }

            Console.WriteLine("");
            Console.WriteLine("  Singular values:");
            Console.WriteLine("");

            for (i = 0; i < Math.Min(M, N); i++)
            {
                Console.WriteLine("  " + (i + 1).ToString().PadLeft(8)
                                       + s[i].ToString().PadLeft(20) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("  Left Singular Vector Matrix U:");
            Console.WriteLine("");

            for (i = 0; i < M; i++)
            {
                cout = "";
                for (j = 0; j < M; j++)
                {
                    cout += "  " + u[i + j * ldu].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }

            Console.WriteLine("");
            Console.WriteLine("  Right Singular Vector Matrix V:");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                cout = "";
                for (j = 0; j < N; j++)
                {
                    cout += "  " + v[i + j * ldv].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }

            for (i = 0; i < M; i++)
            {
                for (j = 0; j < N; j++)
                {
                    a[i + j * lda] = new Complex(0.0, 0.0);
                    for (k = 0; k < Math.Min(M, N); k++)
                    {
                        a[i + j * lda] = a[i + j * lda] + u[i + k * ldu] * s[k] * Complex.Conjugate(v[j + k * ldv]);
                    }
                }
            }

            Console.WriteLine("");
            Console.WriteLine("  The product U * S * V^H (should equal A):");
            Console.WriteLine("");

            for (i = 0; i < M; i++)
            {
                cout = "";
                for (j = 0; j < N; j++)
                {
                    cout += "  " + a[i + j * lda].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }
        }

        static void test35()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST35 tests ZTRCO.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 3;

            Complex[] a = new Complex[N * N];
            int i;
            int j;
            int job;
            int lda;
            double rcond;
            int seed;

            lda = N;

            Console.WriteLine("");
            Console.WriteLine("TEST35");
            Console.WriteLine("  For a complex triangular matrix,");
            Console.WriteLine("  ZTRCO estimates the condition.");
            Console.WriteLine("");
            Console.WriteLine("  Matrix order N = " + N + "");
            //
            //  Set the matrix.
            //
            seed = 123456789;

            for (i = 0; i < N; i++)
            {
                for (j = 0; j <= i; j++)
                {
                    a[i + j * lda] = UniformRNG.c8_uniform_01(ref seed);
                }

                for (j = i + 1; j < N; j++)
                {
                    a[i + j * lda] = new Complex(0.0, 0.0);
                }
            }

            //
            //  Get the condition of the lower triangular matrix.
            //
            job = 0;
            rcond = ZTRCO.ztrco(a, lda, N, job);

            Console.WriteLine("");
            Console.WriteLine("  Estimated reciprocal condition RCOND = " + rcond + "");

        }

        static void test36()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST36 tests ZTRDI.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 3;

            Complex[] a = new Complex[N * N];
            Complex[] a_save = new Complex[N * N];
            Complex[] c = new Complex[N * N];
            Complex[] det = new Complex[2];
            int i;
            int j;
            int job;
            int k;
            int lda;
            int seed;
            string cout;

            lda = N;

            Console.WriteLine("");
            Console.WriteLine("TEST36");
            Console.WriteLine("  For a complex triangular matrix,");
            Console.WriteLine("  ZTRDI computes the determinant or inverse.");
            Console.WriteLine("");
            Console.WriteLine("  Matrix order N = " + N + "");
            //
            //  Set the matrix.
            //
            seed = 123456789;

            for (i = 0; i < N; i++)
            {
                for (j = 0; j <= i; j++)
                {
                    a[i + j * lda] = UniformRNG.c8_uniform_01(ref seed);
                }

                for (j = i + 1; j < N; j++)
                {
                    a[i + j * lda] = new Complex(0.0, 0.0);
                }
            }

            for (i = 0; i < N; i++)
            {
                for (j = 0; j < N; j++)
                {
                    a_save[i + j * lda] = a[i + j * lda];
                }
            }

            //
            //  Get the determinant of the lower triangular matrix.
            //
            job = 100;
            ZTRDI.ztrdi(ref a, lda, N, ref det, job);

            Console.WriteLine("");
            Console.WriteLine("  Determinant = " + det[0]
                                                 + " * 10^ " + (det[1].Real) + "");
            //
            //  Get the inverse of the lower triangular matrix.
            //
            job = 10;
            ZTRDI.ztrdi(ref a, lda, N, ref det, job);

            for (i = 0; i < N; i++)
            {
                for (k = 0; k < N; k++)
                {
                    c[i + k * lda] = 0.0;
                    for (j = 0; j < N; j++)
                    {
                        c[i + k * lda] = c[i + k * lda] + a[i + j * lda] * a_save[j + k * lda];
                    }
                }
            }

            Console.WriteLine("");
            Console.WriteLine("  The product inv(A) * A is");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                cout = "";
                for (j = 0; j < N; j++)
                {
                    cout += "  " + c[i + j * lda].ToString().PadLeft(20);
                }

                Console.WriteLine(cout);
            }

        }

        static void test37()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST37 tests ZTRSL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 10;

            Complex[] a = new Complex[N * N];
            Complex[] b = new Complex[N];
            int i;
            int j;
            int job;
            int k;
            int lda;
            int seed;
            Complex[] x = new Complex[N];

            lda = N;

            Console.WriteLine("");
            Console.WriteLine("TEST37");
            Console.WriteLine("  For a complex triangular matrix,");
            Console.WriteLine("  ZTRSL solves a linear system.");
            Console.WriteLine("");
            Console.WriteLine("  Matrix order N = " + N + "");
            //
            //  Set the matrix.
            //
            seed = 123456789;

            k = 0;
            for (i = 0; i < N; i++)
            {
                for (j = 0; j <= i; j++)
                {
                    k = k + 1;
                    a[i + j * lda] = UniformRNG.c8_uniform_01(ref seed);
                }

                for (j = i + 1; j < N; j++)
                {
                    a[i + j * lda] = new Complex(0.0, 0.0);
                }
            }

            //
            //  Set the desired solution
            //
            for (i = 0; i < N; i++)
            {
                x[i] = new Complex(i + 1, 10 * (i + 1));
            }

            //
            //  Compute the corresponding right hand side.
            //
            for (i = 0; i < N; i++)
            {
                b[i] = 0.0;
                for (j = 0; j < N; j++)
                {
                    b[i] = b[i] + a[i + j * N] * x[j];
                }
            }

            //
            //  Solve the lower triangular system.
            //
            job = 0;
            ZTRSL.ztrsl(a, lda, N, ref b, job);

            Console.WriteLine("");
            Console.WriteLine("  Computed                     Exact");
            Console.WriteLine("  Solution                     Solution");
            Console.WriteLine("");
            for (i = 0; i < N; i++)
            {
                Console.WriteLine("  " + b[i].ToString().PadLeft(26)
                                       + "  " + x[i].ToString().PadLeft(26) + "");
            }

        }

    }
}