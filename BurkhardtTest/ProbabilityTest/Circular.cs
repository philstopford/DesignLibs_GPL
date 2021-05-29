namespace Burkardt.ProbabilityTest
{
    partial class Program
    {
static void circular_normal_01_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    CIRCULAR_NORMAL_01_SAMPLE_TEST tests CIRCULAR_NORMAL_01_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

int j;
double* mean;
int seed = 123456789;
double* variance;
double x[2 * SAMPLE_NUM];
double* xmax;
double* xmin;
double* y;

Console.WriteLine("");
Console.WriteLine("CIRCULAR_NORMAL_01_SAMPLE_TEST");
Console.WriteLine("  CIRCULAR_NORMAL_01_MEAN computes the Circular Normal 01 mean;");
Console.WriteLine("  CIRCULAR_NORMAL_01_SAMPLE samples the Circular Normal 01 distribution;");
Console.WriteLine("  CIRCULAR_NORMAL_01_VARIANCE computes the Circular Normal 01 variance.");

mean = circular_normal_01_mean();
variance = circular_normal_01_variance();

Console.WriteLine("");
Console.WriteLine("  PDF mean =     "
+ setw(12) + mean[0] + "  "
+ setw(12) + mean[1] + "");
Console.WriteLine("  PDF variance = "
+ setw(12) + variance[0] + "  "
+ setw(12) + variance[1] + "");

delete[] mean;
delete[] variance;

for (j = 0; j < SAMPLE_NUM; j++)
{
y = circular_normal_01_sample(seed);
x[0 + j * 2] = y[0];
x[1 + j * 2] = y[1];
delete[] y;
}

mean = r8row_mean(2, SAMPLE_NUM, x);
variance = r8row_variance(2, SAMPLE_NUM, x);
xmax = r8row_max(2, SAMPLE_NUM, x);
xmin = r8row_min(2, SAMPLE_NUM, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
Console.WriteLine("  Sample mean =     "
+ setw(12) + mean[0] + "  "
+ setw(12) + mean[1] + "");
Console.WriteLine("  Sample variance = "
+ setw(12) + variance[0] + "  "
+ setw(12) + variance[1] + "");
Console.WriteLine("  Sample maximum =  "
+ setw(12) + xmax[0] + "  "
+ setw(12) + xmax[1] + "");
Console.WriteLine("  Sample minimum =  "
+ setw(12) + xmin[0] + "  "
+ setw(12) + xmin[1] + "");

delete[] mean;
delete[] variance;
delete[] xmax;
delete[] xmin;

return;
# undef SAMPLE_NUM
}

static void circular_normal_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    CIRCULAR_NORMAL_SAMPLE_TEST tests CIRCULAR_NORMAL_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 January 2012
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

double a[2];
double b;
int j;
double* mean;
int seed = 123456789;
double* variance;
double x[2 * SAMPLE_NUM];
double* xmax;
double* xmin;
double* y;

a[0] = 1.0;
a[1] = 5.0;
b = 0.75;

Console.WriteLine("");
Console.WriteLine("CIRCULAR_NORMAL_SAMPLE_TEST");
Console.WriteLine("  CIRCULAR_NORMAL_MEAN computes the Circular Normal mean;");
Console.WriteLine("  CIRCULAR_NORMAL_SAMPLE samples the Circular Normal distribution;");
Console.WriteLine("  CIRCULAR_NORMAL_VARIANCE computes the Circular Normal variance.");

mean = circular_normal_mean(a, b);
variance = circular_normal_variance(a, b);

Console.WriteLine("");
Console.WriteLine("  PDF mean =     "
+ setw(12) + mean[0] + "  "
+ setw(12) + mean[1] + "");
Console.WriteLine("  PDF variance = "
+ setw(12) + variance[0] + "  "
+ setw(12) + variance[1] + "");

delete[] mean;
delete[] variance;

for (j = 0; j < SAMPLE_NUM; j++)
{
y = circular_normal_sample(a, b, seed);
x[0 + j * 2] = y[0];
x[1 + j * 2] = y[1];
delete[] y;
}

mean = r8row_mean(2, SAMPLE_NUM, x);
variance = r8row_variance(2, SAMPLE_NUM, x);
xmax = r8row_max(2, SAMPLE_NUM, x);
xmin = r8row_min(2, SAMPLE_NUM, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
Console.WriteLine("  Sample mean =     "
+ setw(12) + mean[0] + "  "
+ setw(12) + mean[1] + "");
Console.WriteLine("  Sample variance = "
+ setw(12) + variance[0] + "  "
+ setw(12) + variance[1] + "");
Console.WriteLine("  Sample maximum =  "
+ setw(12) + xmax[0] + "  "
+ setw(12) + xmax[1] + "");
Console.WriteLine("  Sample minimum =  "
+ setw(12) + xmin[0] + "  "
+ setw(12) + xmin[1] + "");

delete[] mean;
delete[] variance;
delete[] xmax;
delete[] xmin;

return;
# undef SAMPLE_NUM
}
        
    }
}