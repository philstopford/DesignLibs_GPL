namespace Burkardt.SolveNS;

public class SuccessiveOverRelaxation
{
    public static double[] sor1 ( int n, double[] a, double[] b, double[] x, double w )

        //****************************************************************************80
    {
        int i;

        double[] x_new = new double[n];
        //
        //  Do the Gauss-Seidel computation.
        //
        for ( i = 0; i < n; i++ )
        {
            x_new[i] = b[i];
            int j;
            for ( j = 0; j < i; j++ )
            {
                x_new[i] -= a[i+j*n] * x_new[j];
            }
            for ( j = i + 1; j < n; j++ )
            {
                x_new[i] -= a[i+j*n] * x[j];
            }
            x_new[i] /= a[i+i*n];
        }
        //
        //  Use W to blend the Gauss-Seidel update with the old solution.
        //
        for ( i = 0; i < n; i++ )
        {
            x_new[i] = ( 1.0 - w ) * x[i] + w * x_new[i];
        }

        return x_new;
    }
}