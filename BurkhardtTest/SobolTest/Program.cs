using System;
using Sobol;
using entropyRNG;

namespace SobolTest
{
    class Program
    {
        static void Main(string[] args)
        {
            Console.WriteLine();
            Console.WriteLine("SOBOL_TEST");
            Console.WriteLine("  Test the SOBOL library.");
 
            or_test ( );
            i4_bit_hi1_test ( );
            i4_bit_lo0_test ( );
            test04 ( );
            test05 ( );
            test055 ( );
            i8_bit_hi1_test ( );
            i8_bit_lo0_test ( );
            test08 ( );
            test09 ( );

            Console.WriteLine();
            Console.WriteLine("SOBOL_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine();
        }
        
        static void or_test ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    OR_TEST tests OR.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            int seed = 123456789;

            Console.WriteLine();
            Console.WriteLine("OR_TEST");
            Console.WriteLine("  The function ^ computes the bitwise exclusive OR");
            Console.WriteLine("  of two integers.");
            Console.WriteLine();
            Console.WriteLine("       I       J     I^J");
            Console.WriteLine();

            for (int test = 1; test <= 10; test++ )
            {
                int i = RNG.nextint( 0, 100, seed );
                int j = RNG.nextint( 0, 100, seed );
                int k = i ^ j;

                string cout = "  ";
                string t = i + "  ";
                cout += t.PadLeft(6); 
                t = j + "  ";
                cout += t.PadLeft(6); 
                t = k + "  ";
                cout += t.PadLeft(6); 
                Console.WriteLine(cout);
            }
        }
        
        
        static void i4_bit_hi1_test ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_BIT_HI1_TEST tests I4_BIT_HI1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            int seed = 123456789;

            Console.WriteLine();;
            Console.WriteLine("I4_BIT_HI1_TEST");
            Console.WriteLine("  I4_BIT_HI1 returns the location of the high bit in an integer.");
            Console.WriteLine();
            Console.WriteLine("       I  I4_BIT_HI1(I)");
            Console.WriteLine();

            for (int test = 1; test <= 10; test++ )
            {
                int i = RNG.nextint( 0, 100, seed );
                int j = SobolSampler.i4_bit_hi1 ( i );

                string cout = "  ";
                string t = i + "  ";
                cout += t.PadLeft(6);
                t = j.ToString();
                cout += t.PadLeft(6);
                Console.WriteLine(cout);
            }
        }
        
        
        static void i4_bit_lo0_test ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_BIT_LO0_TEST tests I4_BIT_LO0.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            int seed = 123456789;

            Console.WriteLine();;
            Console.WriteLine("I4_BIT_LO0_TEST");
            Console.WriteLine("  I4_BIT_LO0 returns the location of the low 0 bit in an integer.");
            Console.WriteLine();
            Console.WriteLine("       I  I4_BIT_LO0(I)");
            Console.WriteLine();

            for (int test = 1; test <= 10; test++ )
            {
                int i = RNG.nextint( 0, 100, seed );
                int j = SobolSampler.i4_bit_lo0 ( i );

                string cout = "  ";
                string t = i + "  ";
                cout += t.PadLeft(6);
                t = j.ToString();
                cout += t.PadLeft(6);
                Console.WriteLine(cout);
            }

            return;
        }

        
        
        static void test04 ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests I4_SOBOL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            int DIM_MAX = 4;
            
            Console.WriteLine();
            Console.WriteLine("TEST04");
            Console.WriteLine("  I4_SOBOL computes the next element of a Sobol sequence.");
            Console.WriteLine();
            Console.WriteLine("  In this test, we call I4_SOBOL repeatedly.");

            for (int dim_num = 2; dim_num <= DIM_MAX; dim_num++ )
            {
                int seed = 0;

                Console.WriteLine();
                Console.WriteLine("  Using dimension DIM_NUM =   " + dim_num);
                Console.WriteLine();
                Console.WriteLine("  Seed  Seed   I4_SOBOL");
                Console.WriteLine("  In    Out");
                Console.WriteLine();

                for (int i = 0; i <= 110; i++ )
                {
                    int seed_in = seed;
                    SobolSampler.SobolVector vec = SobolSampler.i4_sobol ( dim_num, seed);
                    float[] r = vec.quasi;
                    int seed_out = vec.seed;

                    if ( i <= 11 || 95 <= i )
                    {
                        string cout = "";
                        string t = seed_in + "  ";
                        cout += t.PadLeft(6);
                        t = seed_out + "  ";
                        cout += t.PadLeft(6);
                        for (int j = 0; j < dim_num; j++ )
                        {
                            t = r[j] + "  ";
                            cout += t.PadLeft(14);
                        }
                        Console.WriteLine(cout);
                    }
                    else if ( i == 12 )
                    {
                        Console.WriteLine("....................");
                    }
                }
            }
        }
        

        static void test05 ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 tests I4_SOBOL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            int DIM_NUM = 3;
            
            Console.WriteLine();
            Console.WriteLine("TEST05");
            Console.WriteLine("  I4_SOBOL computes the next element of a Sobol sequence.");
            Console.WriteLine();
            Console.WriteLine("  In this test, we demonstrate how the SEED can be");
            Console.WriteLine("  manipulated to skip ahead in the sequence, or");
            Console.WriteLine("  to come back to any part of the sequence.");
            Console.WriteLine();
            Console.WriteLine("  Using dimension DIM_NUM =   " + DIM_NUM);

            int seed = 0;

            Console.WriteLine();
            Console.WriteLine("  Seed  Seed   I4_SOBOL");
            Console.WriteLine("  In    Out");
            Console.WriteLine();

            for (int i = 1; i <= 11; i++ )
            {
                int seed_in = seed;
                SobolSampler.SobolVector vec = SobolSampler.i4_sobol ( DIM_NUM, seed);
                float[] r = vec.quasi;
                int seed_out = vec.seed;
                string cout = "";
                string t = seed_in + "  ";
                cout += t.PadLeft(6);
                t = seed_out + "  ";
                cout += t.PadLeft(6);
                for (int j = 0; j < DIM_NUM; j++ )
                {
                    t = r[j] + "  ";
                    cout += t.PadLeft(14);
                }
                Console.WriteLine(cout);
            }

            Console.WriteLine();
            Console.WriteLine("  Jump ahead by increasing SEED:");

            seed = 100;

            Console.WriteLine();
            Console.WriteLine("  Seed  Seed   I4_SOBOL");
            Console.WriteLine("  In    Out");
            Console.WriteLine();

            for (int i = 1; i <= 5; i++ )
            {
                int seed_in = seed;
                SobolSampler.SobolVector vec = SobolSampler.i4_sobol ( DIM_NUM, seed);
                float[] r = vec.quasi;
                int seed_out = vec.seed;
                string cout = "";
                string t = seed_in + "  ";
                cout += t.PadLeft(6);
                t = seed_out + "  ";
                cout += t.PadLeft(6);
                for (int j = 0; j < DIM_NUM; j++ )
                {
                    t = r[j] + "  ";
                    cout += t.PadLeft(14);
                }
                Console.WriteLine(cout);
            }

            Console.WriteLine();
            Console.WriteLine("  Jump back by decreasing SEED:");

            seed = 3;

            Console.WriteLine();
            Console.WriteLine("  Seed  Seed   I4_SOBOL");
            Console.WriteLine("  In    Out");
            Console.WriteLine();

            for (int i = 1; i <= 11; i++ )
            {
                int seed_in = seed;
                SobolSampler.SobolVector vec = SobolSampler.i4_sobol ( DIM_NUM, seed);
                float[] r = vec.quasi;
                int seed_out = vec.seed;
                string cout = "";
                string t = seed_in + "  ";
                cout += t.PadLeft(6);
                t = seed_out + "  ";
                cout += t.PadLeft(6);
                for (int j = 0; j < DIM_NUM; j++ )
                {
                    t = r[j] + "  ";
                    cout += t.PadLeft(14);
                }
                Console.WriteLine(cout);
            }

            Console.WriteLine();
            Console.WriteLine("  Jump ahead by increasing SEED:");

            seed = 98;

            Console.WriteLine();
            Console.WriteLine("  Seed  Seed   I4_SOBOL");
            Console.WriteLine("  In    Out");
            Console.WriteLine();

            for (int i = 1; i <= 5; i++ )
            {
                int seed_in = seed;
                SobolSampler.SobolVector vec = SobolSampler.i4_sobol ( DIM_NUM, seed);
                float[] r = vec.quasi;
                int seed_out = vec.seed;
                string cout = "";
                string t = seed_in + "  ";
                cout += t.PadLeft(6);
                t = seed_out + "  ";
                cout += t.PadLeft(6);
                for (int j = 0; j < DIM_NUM; j++ )
                {
                    t = r[j] + "  ";
                    cout += t.PadLeft(14);
                }
                Console.WriteLine(cout);
            }
        }
        

        static void test055 ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST055 tests OR on long long ints.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            int seed = 123456789;

            Console.WriteLine();
            Console.WriteLine("TEST055");
            Console.WriteLine("  The C++ function ^ computes the bitwise exclusive OR");
            Console.WriteLine("  of two LONG LONG INT's.");
            Console.WriteLine();
            Console.WriteLine("       I       J     I^J");
            Console.WriteLine();

            for (int  test = 1; test <= 10; test++ )
            {
                long i = ( long ) RNG.nextint( 0, 100, seed );
                long j = ( long ) RNG.nextint( 0, 100, seed );
                long k = i ^ j;

                string cout = "  ";
                string t = i + "  ";
                cout += t.PadLeft(6);
                t = j + "  ";
                cout += t.PadLeft(6);
                t = k.ToString();
                cout += t.PadLeft(6);
                Console.WriteLine(cout);
            }
        }
        
        
        static void i8_bit_hi1_test ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I8_BIT_H1_TEST tests I8_BIT_HI1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            int seed = 123456789;

            Console.WriteLine();
            Console.WriteLine("I8_BIT_HI1_TEST");
            Console.WriteLine("  I8_BIT_HI1 returns the location of the high bit in an integer.");
            Console.WriteLine();
            Console.WriteLine("       I  I8_BIT_HI1(I)");
            Console.WriteLine();

            for (int test = 1; test <= 10; test++ )
            {
                long i = ( long ) RNG.nextint( 0, 100, seed );
                int j = SobolSampler.i8_bit_hi1 ( i );

                string cout = "  ";
                string t = i + "  ";
                cout += t.PadLeft(6);
                t += j.ToString();
                cout += t.PadLeft(6);
                Console.WriteLine(cout);
            }
        }
        
        
        static void i8_bit_lo0_test ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I8_BIT_LO0_TEST tests I8_BIT_LO0.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            int seed = 123456789;

            Console.WriteLine();
            Console.WriteLine("I8_BIT_LO0_TEST");
            Console.WriteLine("  I8_BIT_LO0 returns the location of the low zero bit");
            Console.WriteLine("  in an integer.");
            Console.WriteLine();
            Console.WriteLine("       I  I8_BIT_LO0(I)");
            Console.WriteLine();

            for (int test = 1; test <= 10; test++ )
            {
                long i = ( long ) RNG.nextint( 0, 100, seed );
                int j = SobolSampler.i8_bit_lo0 ( i );

                string cout = "  ";
                string t = i + "  ";
                cout += t.PadLeft(6);
                t = j.ToString();
                cout += t.PadLeft(6);
                Console.WriteLine(cout);
            }
        }

        
        static void test08 ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST08 tests I8_SOBOL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            int DIM_MAX = 4;
            
            Console.WriteLine();
            Console.WriteLine("TEST08");
            Console.WriteLine("  I8_SOBOL computes the next element of a Sobol sequence.");
            Console.WriteLine();
            Console.WriteLine("  In this test, we call I8_SOBOL repeatedly.");

            for (int dim_num = 2; dim_num <= DIM_MAX; dim_num++ )
            {

                long seed = 0;

                Console.WriteLine();
                Console.WriteLine("  Using dimension DIM_NUM =   " + dim_num);
                Console.WriteLine();
                Console.WriteLine("  Seed  Seed   I8_SOBOL");
                Console.WriteLine("  In    Out");
                Console.WriteLine();

                for (int i = 0; i <= 110; i++ )
                {
                    long seed_in = seed;
                    SobolSampler.SobolVectorLarge res = SobolSampler.i8_sobol ( dim_num, seed );
                    double[] r = res.quasi;
                    long seed_out = res.seed;

                    if (i <= 11 || 95 <= i )
                    {
                        string cout = seed_in + "  ";
                        cout = cout.PadLeft(6);
                        string t = seed_out + "  ";
                        cout += t.PadLeft(6);
                        for (int j = 0; j < dim_num; j++ )
                        {
                            t = r[j] + "  ";
                            cout += t.PadLeft(14);
                        }
                        Console.WriteLine(cout);
                    }
                    else if ( i == 12 )
                    {
                        Console.WriteLine("....................");
                    }
                }

            }
        }

        static void test09 ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST09 tests I8_SOBOL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            int DIM_NUM = 3;

            Console.WriteLine();
            Console.WriteLine("TEST09");
            Console.WriteLine("  I8_SOBOL computes the next element of a Sobol sequence.");
            Console.WriteLine();
            Console.WriteLine("  In this test, we demonstrate how the SEED can be");
            Console.WriteLine("  manipulated to skip ahead in the sequence, or");
            Console.WriteLine("  to come back to any part of the sequence.");
            Console.WriteLine();
            Console.WriteLine("  Using dimension DIM_NUM =   " + DIM_NUM);

            long seed = 0;

            Console.WriteLine();
            Console.WriteLine("  Seed  Seed   I8_SOBOL");
            Console.WriteLine("  In    Out");
            Console.WriteLine();

            for (int i = 1; i <= 11; i++ )
            {
                long seed_in = seed;
                SobolSampler.SobolVectorLarge res = SobolSampler.i8_sobol ( DIM_NUM, seed );
                double[] r = res.quasi;
                long seed_out = res.seed;
                string cout = seed_in + "  ";
                cout = cout.PadLeft(6);
                string t = seed_out + "  ";
                cout += t.PadLeft(6);
                for (int j = 0; j < DIM_NUM; j++ )
                {
                    t = r[j] + "  ";
                    cout += t.PadLeft(14);
                }
                Console.WriteLine(cout);
            }

            Console.WriteLine();
            Console.WriteLine("  Jump ahead by increasing SEED:");

            seed = 100;

            Console.WriteLine();
            Console.WriteLine("  Seed  Seed   I8_SOBOL");
            Console.WriteLine("  In    Out");
            Console.WriteLine();

            for (int i = 1; i <= 5; i++ )
            {
                long seed_in = seed;
                SobolSampler.SobolVectorLarge res = SobolSampler.i8_sobol ( DIM_NUM, seed );
                double[] r = res.quasi;
                long seed_out = res.seed;
                string cout = seed_in + "  ";
                cout = cout.PadLeft(6);
                string t = seed_out + "  ";
                cout += t.PadLeft(6);
                for (int j = 0; j < DIM_NUM; j++ )
                {
                    t = r[j] + "  ";
                    cout += t.PadLeft(14);
                }
                Console.WriteLine(cout);
            }

            Console.WriteLine();
            Console.WriteLine("  Jump back by decreasing SEED:");

            seed = 3;

            Console.WriteLine();
            Console.WriteLine("  Seed  Seed   I8_SOBOL");
            Console.WriteLine("  In    Out");
            Console.WriteLine();

            for (int i = 1; i <= 11; i++ )
            {
                long seed_in = seed;
                SobolSampler.SobolVectorLarge res = SobolSampler.i8_sobol ( DIM_NUM, seed );
                double[] r = res.quasi;
                long seed_out = res.seed;
                string cout = seed_in + "  ";
                cout = cout.PadLeft(6);
                string t = seed_out + "  ";
                cout += t.PadLeft(6);
                for (int j = 0; j < DIM_NUM; j++ )
                {
                    t = r[j] + "  ";
                    cout += t.PadLeft(14);
                }
                Console.WriteLine(cout);
            }

            Console.WriteLine();
            Console.WriteLine("  Jump ahead by increasing SEED:");

            seed = 98;

            Console.WriteLine();
            Console.WriteLine("  Seed  Seed   I8_SOBOL");
            Console.WriteLine("  In    Out");
            Console.WriteLine();

            for (int i = 1; i <= 5; i++ )
            {
                long seed_in = seed;
                SobolSampler.SobolVectorLarge res = SobolSampler.i8_sobol ( DIM_NUM, seed );
                double[] r = res.quasi;
                long seed_out = res.seed;
                string cout = seed_in + "  ";
                cout = cout.PadLeft(6);
                string t = seed_out + "  ";
                cout += t.PadLeft(6);
                for (int j = 0; j < DIM_NUM; j++ )
                {
                    t = r[j] + "  ";
                    cout += t.PadLeft(14);
                }
                Console.WriteLine(cout);
            }
        }

    }
}