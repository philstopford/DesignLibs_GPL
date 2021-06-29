using Burkardt.Types;

namespace DREAMTest
{
    public static partial class Problem1c
    {
        class Covariance
        {
            int order;
            double[] array { get; set; }
            double det { get; set; }
            double[] factor { get; set; }
            double[] inv { get; set; }
            double[] mean { get; set; }

            public Covariance(int par_num)
            {
                order = par_num;
                array_set();
                factor_set();
                det_set();
                inv_set();
                mean_set();
            }
            
            public double[] array_get()
            {
                double[] array2 = typeMethods.r8mat_copy_new ( order, order, array );

                return array2;
            }

            public void array_set()
            {
                int i;
                int j;

                array = new double[order*order];

                for ( j = 0; j < order; j++ )
                {
                    for ( i = 0; i < order; i++ )
                    {
                        array[i+j*order] = 0.5;
                    }
                }

                for ( i = 0; i < order; i++ )
                {
                    array[i+i*order] = ( double ) ( i + 1 );
                }
            }

            public double det_get ( )
            {
                double det2 = det;

                return det2;
            }

            public void det_set()
            {
                det = typeMethods.r8mat_podet(order, factor);
            }

            public double[] factor_get ( )
            {
                double[] factor2 = typeMethods.r8mat_copy_new ( order, order, factor );

                return factor2;
            }

            public void factor_set ( )
            {
                factor = typeMethods.r8mat_pofac ( order, array );
            }

            public double[] inv_get ( )
            {
                double[] inv2 = typeMethods.r8mat_copy_new ( order, order, inv );

                return inv2;
            }

            public void inv_set()
            {
                inv = typeMethods.r8mat_poinv(order, factor);
            }

            public double[] mean_get ( )
            {
                double[] mean2 = typeMethods.r8vec_copy_new ( order, mean );

                return mean2;
            }

            public void mean_set()
            {
                int i;

                mean = new double[order];

                for (i = 0; i < order; i++)
                {
                    mean[i] = 0.0;
                }
            }

            public void print(string title)
            {
                typeMethods.r8mat_print ( order, order, array, title );
            }
        }
    }
}