using Tral.Randomness;
using Tral.Randomness.Algorithms;
using NUnit.Framework;

namespace UnitTests;

/// <summary>
/// Tests for the Tral.Randomness library functionality, covering random number generation,
/// algorithm implementations, and statistical properties.
/// </summary>
[TestFixture]
public class TralRandomnessTests
{
    [Test]
    public void RandomGenerator_Constructor_ShouldCreateGenerator()
    {
        var generator = new RandomGenerator(false); // Don't randomize for predictable testing

        Assert.That(generator, Is.Not.Null);
        Assert.That(generator.AlgorithmName, Is.Not.Null.And.Not.Empty);
    }

    [Test]
    public void RandomGenerator_ConstructorWithRandomize_ShouldCreateRandomizedGenerator()
    {
        var generator = new RandomGenerator(true);

        Assert.That(generator, Is.Not.Null);
        Assert.That(generator.AlgorithmName, Is.Not.Null.And.Not.Empty);
    }

    [Test]
    public void RandomGenerator_Next64_ShouldGenerateValues()
    {
        var generator = new RandomGenerator(false);

        ulong value1 = generator.Next64();
        ulong value2 = generator.Next64();

        // Values should be different (extremely unlikely to be the same)
        Assert.That(value1, Is.Not.EqualTo(value2));
    }

    [Test]
    public void RandomGenerator_NextInt_ShouldGenerateValues()
    {
        var generator = new RandomGenerator(false);

        int value1 = generator.NextInt(100);
        int value2 = generator.NextInt(100);

        Assert.That(value1, Is.GreaterThanOrEqualTo(0));
        Assert.That(value2, Is.GreaterThanOrEqualTo(0));
        Assert.That(value1, Is.LessThan(100));
        Assert.That(value2, Is.LessThan(100));
    }

    [Test]
    public void RandomGenerator_NextInt_WithMaxValue_ShouldRespectBounds()
    {
        var generator = new RandomGenerator(false);
        int maxValue = 100;

        for (int i = 0; i < 50; i++)
        {
            int value = generator.NextInt(maxValue);
            Assert.That(value, Is.GreaterThanOrEqualTo(0).And.LessThan(maxValue));
        }
    }

    [Test]
    public void RandomGenerator_NextInt_WithMinMaxValues_ShouldRespectBounds()
    {
        var generator = new RandomGenerator(false);
        int minValue = 10;
        int maxValue = 20;

        for (int i = 0; i < 50; i++)
        {
            int value = generator.NextInt(minValue, maxValue);
            Assert.That(value, Is.GreaterThanOrEqualTo(minValue).And.LessThan(maxValue));
        }
    }

    [Test]
    public void RandomGenerator_NextDouble_ShouldGenerateValuesInRange()
    {
        var generator = new RandomGenerator(false);

        for (int i = 0; i < 50; i++)
        {
            double value = generator.NextDouble();
            Assert.That(value, Is.GreaterThanOrEqualTo(0.0).And.LessThan(1.0));
        }
    }

    [Test]
    public void RandomGenerator_NextFlip_ShouldGenerateBothValues()
    {
        var generator = new RandomGenerator(false);
        bool foundTrue = false;
        bool foundFalse = false;

        // Generate many values to ensure we get both true and false
        for (int i = 0; i < 100; i++)
        {
            bool value = generator.NextFlip();
            if (value) foundTrue = true;
            else foundFalse = true;

            if (foundTrue && foundFalse) break;
        }

        Assert.That(foundTrue, Is.True);
        Assert.That(foundFalse, Is.True);
    }

    [Test]
    public void RandomGenerator_Global_ShouldProvideThreadLocalInstance()
    {
        var global = RandomGenerator.Global;

        Assert.That(global, Is.Not.Null);
        Assert.That(global, Is.SameAs(RandomGenerator.Global)); // Should be same instance for same thread
    }

    [Test]
    public void RandomGenerator_Randomize_ShouldChangeState()
    {
        var generator = new RandomGenerator(false);

        // Get initial state  
        ulong value1 = generator.Next64();

        // Randomize
        generator.Randomize();

        // Should produce different sequence after randomize
        ulong value2 = generator.Next64();
        Assert.That(value1, Is.Not.EqualTo(value2));
    }

    [Test]
    public void Xoshiro256pp_Algorithm_ShouldHaveCorrectProperties()
    {
        var algorithm = new Xoshiro256pp();

        Assert.That(algorithm.AlgorithmName, Is.EqualTo("xoshiro256++ 1.0"));
        Assert.That(algorithm.MaxNext, Is.EqualTo(ulong.MaxValue));
    }

    [Test]
    public void SplitMix64_Algorithm_ShouldHaveCorrectProperties()
    {
        var algorithm = new SplitMix64();

        Assert.That(algorithm.AlgorithmName, Is.EqualTo("SplitMix64"));
        Assert.That(algorithm.MaxNext, Is.EqualTo(ulong.MaxValue));
    }

    [Test]
    public void RandomGenerator_NextBytes_ShouldFillArray()
    {
        var generator = new RandomGenerator(false);
        byte[] bytes = new byte[10];

        generator.NextBytes(bytes);

        // Check that at least some bytes are non-zero (extremely unlikely to be all zeros)
        bool hasNonZero = false;
        foreach (byte b in bytes)
        {
            if (b != 0)
            {
                hasNonZero = true;
                break;
            }
        }
        Assert.That(hasNonZero, Is.True);
    }

    [Test]
    public void RandomGenerator_NextString_ShouldGenerateStringOfCorrectLength()
    {
        var generator = new RandomGenerator(false);

        string result = generator.NextString(10, "ABCDEF");

        Assert.That(result.Length, Is.EqualTo(10));
        Assert.That(result, Does.Match("^[ABCDEF]*$"));
    }

    [Test]
    public void RandomGenerator_StatisticalDistribution_ShouldBeReasonableForNextDouble()
    {
        var generator = new RandomGenerator(false);
        int sampleCount = 10000;
        double sum = 0;

        for (int i = 0; i < sampleCount; i++)
        {
            sum += generator.NextDouble();
        }

        double average = sum / sampleCount;

        // For uniform distribution [0,1), average should be around 0.5
        Assert.That(average, Is.GreaterThan(0.4).And.LessThan(0.6));
    }
}