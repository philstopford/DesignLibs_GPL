using entropyRNG;

namespace UnitTests;

public class RNGTests
{
    private static int sampleCount = 25000;
    [SetUp]
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
        Parallel.For(0, sampleCount, i => { values[i] = Crypto_RNG.random_gauss3()[0]; });
        Assert.AreEqual(sampleCount, values.Distinct().Count());
        
        double[] values2 = new double[sampleCount];
        Parallel.For(0, sampleCount, i => { values2[i] = Crypto_RNG.random_gauss()[0]; });
        Assert.AreEqual(sampleCount, values2.Distinct().Count());

        double[] values3 = new double[sampleCount];
        Parallel.For(0, sampleCount, i => { values3[i] = Crypto_RNG.nextdouble(); });
        Assert.AreEqual(sampleCount, values3.Distinct().Count());

        int[] ints = new int[sampleCount];
        Parallel.For(0, sampleCount, i => { ints[i] = Crypto_RNG.nextint(); });
        Assert.AreEqual(sampleCount, ints.Distinct().Count());
    }
    
    [Test]
    public static void MersenneTwisterRNGTest()
    {
        double[] values = new double[sampleCount];
        Parallel.For(0, sampleCount, i => { values[i] = MersenneTwister_RNG.random_gauss3()[0]; });
        Assert.AreEqual(sampleCount, values.Distinct().Count());

        double[] values2 = new double[sampleCount];
        Parallel.For(0, sampleCount, i => { values2[i] = MersenneTwister_RNG.random_gauss()[0]; });
        Assert.AreEqual(sampleCount, values2.Distinct().Count());

        double[] values3 = new double[sampleCount];
        Parallel.For(0, sampleCount, i => { values3[i] = MersenneTwister_RNG.nextdouble(); });
        Assert.AreEqual(sampleCount, values3.Distinct().Count());

        int[] ints = new int[sampleCount];
        Parallel.For(0, sampleCount, i => { ints[i] = MersenneTwister_RNG.nextint(); });
        Assert.AreEqual(sampleCount, ints.Distinct().Count());
    }
    
    [Test]
    public static void SystemRNGTest()
    {
        double[] values = new double[sampleCount];
        Parallel.For(0, sampleCount, i => { values[i] = RNG.random_gauss3()[0]; });
        Assert.AreEqual(sampleCount, values.Distinct().Count());
        
        double[] values2 = new double[sampleCount];
        Parallel.For(0, sampleCount, i => { values2[i] = RNG.random_gauss()[0]; });
        Assert.AreEqual(sampleCount, values2.Distinct().Count());

        double[] values3 = new double[sampleCount];
        Parallel.For(0, sampleCount, i => { values3[i] = RNG.nextdouble(); });
        Assert.AreEqual(sampleCount, values3.Distinct().Count());

        int[] ints = new int[sampleCount];
        Parallel.For(0, sampleCount, i => { ints[i] = RNG.nextint(); });
        Assert.AreEqual(sampleCount, ints.Distinct().Count());
    }
}