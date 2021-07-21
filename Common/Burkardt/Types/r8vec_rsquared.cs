using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static double r8vec_rsquared(int n, double[] y_data, double[] y_model )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_RSQUARED computes the R^2 goodness of fit measurement.
        //
        //  Discussion:
        //
        //    We suppose a set of N data values Y_DATA are known, and that a
        //    formula or model has generated a corresponding set of Y_MODEL values.
        //
        //    R^2 measures the extent to which the variation in Y_DATA is captured
        //    by the model data Y_MODEL.  If the model is linear, then a value
        //    of R^2=1 is optimal, and R^2=0 is the worst possible.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 January 2019
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of values.
        //
        //    Input, double Y_DATA[N], Y_MODEL[N], the data and model values.
        //
        //    Output, double R8VEC_RSQUARED, the R^2 value.
        //
        {
            double bot;
            double top;
            int i;
            double value;
            double y_average;

            y_average = 0.0;
            for (i = 0; i < n; i++)
            {
                y_average = y_average + y_data[i];
            }

            y_average = y_average / (double) (n);

            top = 0.0;
            bot = 0.0;
            for (i = 0; i < n; i++)
            {
                top = top + Math.Pow(y_data[i] - y_model[i], 2);
                bot = bot + Math.Pow(y_data[i] - y_average, 2);
            }

            value = 1.0 - top / bot;

            return value;
        }

        public static double r8vec_rsquared_adjusted(int n, double[] y_data, double[] y_model,
        int degree )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_RSQUARED_ADJUSTED computes the adjusted R^2 goodness of fit measurement.
        //
        //  Discussion:
        //
        //    We suppose a set of N data values Y_DATA are known, and that a
        //    formula or model has generated a corresponding set of Y_MODEL values.
        //
        //    R^2 measures the extent to which the variation in Y_DATA is captured
        //    by the model data Y_MODEL.  
        //
        //    The adjusted value of R^2 accounts for the use of a polynomial model
        //    of degree higher than 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 January 2019
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of values.
        //
        //    Input, double Y_DATA[N], Y_MODEL[N], the data and model values.
        //
        //    Input, int DEGREE, the degree of the polynomial model.
        //
        //    Output, double R8VEC_RSQUARED_ADJUSTED, the adjusted R^2 value.
        //
        {
            double bot;
            double top;
            int i;
            double value;
            double y_average;

            y_average = 0.0;
            for (i = 0; i < n; i++)
            {
                y_average = y_average + y_data[i];
            }

            y_average = y_average / (double) (n);

            top = 0.0;
            bot = 0.0;
            for (i = 0; i < n; i++)
            {
                top = top + Math.Pow(y_data[i] - y_model[i], 2);
                bot = bot + Math.Pow(y_data[i] - y_average, 2);
            }

            value = 1.0 - (top / bot)
                * (double) (n - 1) / (double) (n - degree - 1);

            return value;
        }
    }
}