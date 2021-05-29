namespace Burkardt.ProbabilityTest
{
    partial class Program
    {
static void dirichlet_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    DIRICHLET_SAMPLE_TEST tests DIRICHLET_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 April 2016
//
//  Author:
//
//    John Burkardt
//
{
# define N 3
# define SAMPLE_NUM 1000

double a[N] =  {
0.250, 0.500, 1.250
}
;
int i;
int j;
double* mean;
double* m2;
int seed = 123456789;
double* variance;
double x[N * SAMPLE_NUM];
double* xmax;
double* xmin;
double* y;

Console.WriteLine("");
Console.WriteLine("DIRICHLET_SAMPLE");
Console.WriteLine("  DIRICHLET_MEAN computes the Dirichlet mean;");
Console.WriteLine("  DIRICHLET_SAMPLE samples the Dirichlet distribution;");
Console.WriteLine("  DIRICHLET_VARIANCE computes the Dirichlet variance;");

Console.WriteLine("");
Console.WriteLine("  Number of components N = " + N + "");

r8vec_print(N, a, "  PDF parameter A:");

if (!dirichlet_check(N, a))
{
Console.WriteLine("");
Console.WriteLine("DIRICHLET_SAMPLE - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

mean = dirichlet_mean(N, a);
variance = dirichlet_variance(N, a);
r8vec_print(N, mean, "  PDF mean:");
r8vec_print(N, variance, "  PDF variance:");

delete[] mean;
delete[] variance;

m2 = dirichlet_moment2(N, a);

r8mat_print(N, N, m2, "  Second moment matrix:");

for (j = 0; j < SAMPLE_NUM; j++)
{
y = dirichlet_sample(N, a, seed);
for (i = 0; i < N; i++)
{
x[i + j * N] = y[i];
}

delete[] y;
}

mean = r8row_mean(N, SAMPLE_NUM, x);
variance = r8row_variance(N, SAMPLE_NUM, x);
xmax = r8row_max(N, SAMPLE_NUM, x);
xmin = r8row_min(N, SAMPLE_NUM, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
Console.WriteLine("");
Console.WriteLine("  Component Mean, Variance, Min, Max:");
Console.WriteLine("");

for (i = 0; i < N; i++)
{
Console.WriteLine("  "
+ setw(6) + i + "  "
+ setw(12) + mean[i] + "  "
+ setw(12) + variance[i] + "  "
+ setw(12) + xmax[i] + "  "
+ setw(12) + xmin[i] + "");
}

delete[] mean;
delete[] m2;
delete[] variance;
delete[] xmax;
delete[] xmin;

return;
# undef N
# undef SAMPLE_NUM
}

static void dirichlet_pdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    DIRICHLET_PDF_TEST tests DIRICHLET_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 April 2016
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

double a[N] =  {
0.250, 0.500, 1.250
}
;
double pdf;
double x[N] =  {
0.500, 0.125, 0.375
}
;

Console.WriteLine("");
Console.WriteLine("DIRICHLET_PDF_TEST");
Console.WriteLine("  DIRICHLET_PDF evaluates the Dirichlet PDF;");

Console.WriteLine("");
Console.WriteLine("  Number of components N = " + N + "");

r8vec_print(N, a, "  PDF parameter A:");

if (!dirichlet_check(N, a))
{
Console.WriteLine("");
Console.WriteLine("DIRICHLET_PDF_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

r8vec_print(N, x, "  PDF argument X:");

pdf = dirichlet_pdf(x, N, a);

Console.WriteLine("");
Console.WriteLine("  PDF value = " + pdf + "");

return;
# undef N
}

static void dirichlet_mix_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    DIRICHLET_MIX_SAMPLE_TEST tests DIRICHLET_MIX_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 April 2016
//
//  Author:
//
//    John Burkardt
//
{
# define COMP_NUM 2
# define ELEM_NUM 3
# define SAMPLE_NUM 1000

double a[ELEM_NUM * COMP_NUM] =  {
0.250, 0.500, 1.250,
1.500, 0.500, 2.000
}
;
int comp;
double comp_weight[COMP_NUM] =  {
1.0, 2.0
}
;
int i;
int j;
double* mean;
int seed = 123456789;
double* variance;
double x[ELEM_NUM * SAMPLE_NUM];
double* xmax;
double* xmin;
double* y;

Console.WriteLine("");
Console.WriteLine("DIRICHLET_MIX_SAMPLE_TEST");
Console.WriteLine("  DIRICHLET_MIX_SAMPLE samples the Dirichlet Mix distribution;");
Console.WriteLine("  DIRICHLET_MIX_MEAN computes the Dirichlet Mix mean;");

Console.WriteLine("");
Console.WriteLine("  Number of elements ELEM_NUM =   " + ELEM_NUM + "");
Console.WriteLine("  Number of components COMP_NUM = " + COMP_NUM + "");
r8mat_print(ELEM_NUM, COMP_NUM, a, "  PDF parameters A(ELEM,COMP):");
r8vec_print(COMP_NUM, comp_weight, "  Component weights");

if (!dirichlet_mix_check(COMP_NUM, ELEM_NUM, a, comp_weight))
{
Console.WriteLine("");
Console.WriteLine("DIRICHLET_MIX_SAMPLE_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

mean = dirichlet_mix_mean(COMP_NUM, ELEM_NUM, a, comp_weight);

r8vec_print(ELEM_NUM, mean, "  PDF mean:");

delete[] mean;

for (j = 0; j < SAMPLE_NUM; j++)
{
y = dirichlet_mix_sample(COMP_NUM, ELEM_NUM, a, comp_weight, seed,
&comp);

for (i = 0; i < ELEM_NUM; i++)
{
x[i + j * ELEM_NUM] = y[i];
}

delete[] y;
}

mean = r8row_mean(ELEM_NUM, SAMPLE_NUM, x);
variance = r8row_variance(ELEM_NUM, SAMPLE_NUM, x);
xmax = r8row_max(ELEM_NUM, SAMPLE_NUM, x);
xmin = r8row_min(ELEM_NUM, SAMPLE_NUM, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
Console.WriteLine("");
Console.WriteLine("  Component Mean, Variance, Max, Min:");
Console.WriteLine("");

for (i = 0; i < ELEM_NUM; i++)
{
Console.WriteLine("  "
+ setw(6) + i + "  "
+ setw(12) + mean[i] + "  "
+ setw(12) + variance[i] + "  "
+ setw(12) + xmax[i] + "  "
+ setw(12) + xmin[i] + "");
}

delete[] mean;
delete[] variance;
delete[] xmax;
delete[] xmin;

return;
}

static void dirichlet_mix_pdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    DIRICHLET_MIX_PDF_TEST tests DIRICHLET_MIX_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 April 2016
//
//  Author:
//
//    John Burkardt
//
{
# define COMP_NUM 2
# define ELEM_NUM 3

double a[ELEM_NUM * COMP_NUM] =  {
0.250, 0.500, 1.250,
1.500, 0.500, 2.000
}
;
double comp_weight[COMP_NUM] =  {
1.0, 2.0
}
;
double pdf;
double x[ELEM_NUM] =  {
0.500, 0.125, 0.375
}
;

Console.WriteLine("");
Console.WriteLine("DIRICHLET_MIX_PDF_TEST");
Console.WriteLine("  DIRICHLET_MIX_PDF evaluates the Dirichlet Mix PDF.");

Console.WriteLine(flush;

Console.WriteLine("");
Console.WriteLine("  Number of elements ELEM_NUM =   " + ELEM_NUM + "");
Console.WriteLine("  Number of components COMP_NUM = " + COMP_NUM + "");
r8mat_print(ELEM_NUM, COMP_NUM, a, "  PDF parameters A(ELEM,COMP):");
r8vec_print(COMP_NUM, comp_weight, "  Component weights");

Console.WriteLine(flush;

if (!dirichlet_mix_check(COMP_NUM, ELEM_NUM, a, comp_weight))
{
Console.WriteLine("");
Console.WriteLine("DIRICHLET_MIX_PDF_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

r8vec_print(ELEM_NUM, x, "  PDF argument X:");

Console.WriteLine(flush;

pdf = dirichlet_mix_pdf(x, COMP_NUM, ELEM_NUM, a, comp_weight);

Console.WriteLine("");
Console.WriteLine("  PDF value =           " + pdf + "");

return;
# undef COMP_NUM
# undef ELEM_NUM
}
        
    }
}