using entropyRNG;
using utility;

namespace UnitTests;

public class RNGTests
{
    private static int sampleCount = 25000;

    /// <summary>
    /// High-performance parallel execution optimized for RNG operations
    /// </summary>
    private static void OptimizedParallelFor<T>(T[] array, Func<int, T> generator)
    {
        // Use the new optimized parallel processing utility
        ParallelProcessing.OptimizedParallelFor(array, generator);
    }

    // [SetUp]
    public static void RNGTest()
    {
        CryptoRNGTest();
        MersenneTwisterRNGTest();
        SystemRNGTest();
    }

    // The distinct tests are in case the RNG runs out of samples leading to identical values being returned from calls.
    // The threading also encourages any potential failure.
    // The system should be completely robust against failing.
    [Test]
    public static void CryptoRNGTest()
    {
        double[] values = new double[sampleCount];
        OptimizedParallelFor(values, i => Crypto_RNG.random_gauss3()[0]);

        double[] values2 = new double[sampleCount];
        OptimizedParallelFor(values2, i => Crypto_RNG.random_gauss()[0]);

        double[] values3 = new double[sampleCount];
        OptimizedParallelFor(values3, i => Crypto_RNG.nextdouble());

        int[] ints = new int[sampleCount];
        OptimizedParallelFor(ints, i => Crypto_RNG.nextint());
        
        // Allow for occasional duplicates due to birthday paradox - with 25000 samples from uint.MaxValue space,
        // probability of collision is ~7%. Require at least 99.5% distinct values (at most 125 duplicates).
        int min_unique_samples = (int)(sampleCount * 0.995);


        int[] duplicate_ints = ints
            .GroupBy(x => x)               // group matching items
            .Where(g => g.Skip(1).Any())   // where the group contains more than one item
            .SelectMany(g => g).ToArray();           // re-expand the groups with more than one item
        Assert.That(ints.Distinct().Count(), Is.GreaterThanOrEqualTo(min_unique_samples));

        double[] duplicate_values = values
            .GroupBy(x => x)               // group matching items
            .Where(g => g.Skip(1).Any())   // where the group contains more than one item
            .SelectMany(g => g).ToArray();           // re-expand the groups with more than one item
        Assert.That(values.Distinct().Count(), Is.GreaterThanOrEqualTo(min_unique_samples));

        double[] duplicate_values2 = values2
            .GroupBy(x => x)               // group matching items
            .Where(g => g.Skip(1).Any())   // where the group contains more than one item
            .SelectMany(g => g).ToArray();           // re-expand the groups with more than one item
        Assert.That(values2.Distinct().Count(), Is.GreaterThanOrEqualTo(min_unique_samples));

        double[] duplicate_values3 = values3
            .GroupBy(x => x)               // group matching items
            .Where(g => g.Skip(1).Any())   // where the group contains more than one item
            .SelectMany(g => g).ToArray();           // re-expand the groups with more than one item
        // Allow for occasional duplicates due to birthday paradox - with 25000 samples from uint.MaxValue space,
        // probability of collision is ~7%. Require at least 99.5% distinct values (at most 125 duplicates).
        Assert.That(values3.Distinct().Count(), Is.GreaterThanOrEqualTo((int)(sampleCount * 0.995)));
    }

    [Test]
    public static void MersenneTwisterRNGTest()
    {
        double[] values = new double[sampleCount];
        OptimizedParallelFor(values, i => MersenneTwister_RNG.random_gauss3()[0]);

        double[] values2 = new double[sampleCount];
        OptimizedParallelFor(values2, i => MersenneTwister_RNG.random_gauss()[0]);

        double[] values3 = new double[sampleCount];
        OptimizedParallelFor(values3, i => MersenneTwister_RNG.nextdouble());

        int[] ints = new int[sampleCount];
        OptimizedParallelFor(ints, i => MersenneTwister_RNG.nextint());

        // Allow for occasional duplicates due to birthday paradox - with 25000 samples from uint.MaxValue space,
        // probability of collision is ~7%. Require at least 99.5% distinct values (at most 125 duplicates).
        int min_unique_samples = (int)(sampleCount * 0.995);

        int[] duplicate_ints = ints
            .GroupBy(x => x)               // group matching items
            .Where(g => g.Skip(1).Any())   // where the group contains more than one item
            .SelectMany(g => g).ToArray();           // re-expand the groups with more than one item
        Assert.That(ints.Distinct().Count(), Is.GreaterThanOrEqualTo(min_unique_samples));

        double[] duplicate_values = values
            .GroupBy(x => x)               // group matching items
            .Where(g => g.Skip(1).Any())   // where the group contains more than one item
            .SelectMany(g => g).ToArray();           // re-expand the groups with more than one item
        Assert.That(values.Distinct().Count(), Is.GreaterThanOrEqualTo(min_unique_samples));

        double[] duplicate_values2 = values2
            .GroupBy(x => x)               // group matching items
            .Where(g => g.Skip(1).Any())   // where the group contains more than one item
            .SelectMany(g => g).ToArray();           // re-expand the groups with more than one item
        Assert.That(values2.Distinct().Count(), Is.GreaterThanOrEqualTo(min_unique_samples));

        double[] duplicate_values3 = values3
            .GroupBy(x => x)               // group matching items
            .Where(g => g.Skip(1).Any())   // where the group contains more than one item
            .SelectMany(g => g).ToArray();           // re-expand the groups with more than one item
        Assert.That(values3.Distinct().Count(), Is.GreaterThanOrEqualTo(min_unique_samples));
    }

    [Test]
    public static void SystemRNGTest()
    {
        double[] values = new double[sampleCount];
        OptimizedParallelFor(values, i => RNG.random_gauss3()[0]);

        double[] values2 = new double[sampleCount];
        OptimizedParallelFor(values2, i => RNG.random_gauss()[0]);

        double[] values3 = new double[sampleCount];
        OptimizedParallelFor(values3, i => RNG.nextdouble());

        int[] ints = new int[sampleCount];
        OptimizedParallelFor(ints, i => RNG.nextint());

        // Allow for occasional duplicates due to birthday paradox - with 25000 samples from uint.MaxValue space,
        // probability of collision is ~7%. Require at least 99.5% distinct values (at most 125 duplicates).
        int min_unique_samples = (int)(sampleCount * 0.995);

        int[] duplicate_ints = ints
            .GroupBy(x => x)               // group matching items
            .Where(g => g.Skip(1).Any())   // where the group contains more than one item
            .SelectMany(g => g).ToArray();           // re-expand the groups with more than one item
        Assert.That(ints.Distinct().Count(), Is.GreaterThanOrEqualTo(min_unique_samples));

        double[] duplicate_values = values
            .GroupBy(x => x)               // group matching items
            .Where(g => g.Skip(1).Any())   // where the group contains more than one item
            .SelectMany(g => g).ToArray();           // re-expand the groups with more than one item
        Assert.That(values.Distinct().Count(), Is.GreaterThanOrEqualTo(min_unique_samples));

        double[] duplicate_values2 = values2
            .GroupBy(x => x)               // group matching items
            .Where(g => g.Skip(1).Any())   // where the group contains more than one item
            .SelectMany(g => g).ToArray();           // re-expand the groups with more than one item
        Assert.That(values2.Distinct().Count(), Is.GreaterThanOrEqualTo(min_unique_samples));

        double[] duplicate_values3 = values3
            .GroupBy(x => x)               // group matching items
            .Where(g => g.Skip(1).Any())   // where the group contains more than one item
            .SelectMany(g => g).ToArray();           // re-expand the groups with more than one item
        Assert.That(values3.Distinct().Count(), Is.GreaterThanOrEqualTo(min_unique_samples));
    }
}
