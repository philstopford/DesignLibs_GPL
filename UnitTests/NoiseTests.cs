using Noise;

namespace UnitTests;

public class NoiseTests
{
    // [SetUp]
    public static void NoiseSetup()
    {
        OpenSimplexNoiseTest();
        SimplexNoiseTest();
        PerlinNoiseTest();
    }

    [Test]
    public static void OpenSimplexNoiseTest()
    {
        OpenSimplexNoise noise = new(1234);
        double value_x0_y0 = noise.Evaluate(0, 0);
        double value_x0_y0_z0 = noise.Evaluate(0, 0, 0);
        double value_x0_y1 = noise.Evaluate(0, 1);
        double value_x0_y1_z0 = noise.Evaluate(0, 1, 0);
        double value_x0_y0_z1 = noise.Evaluate(0, 0, 1);
        double value_x1_y0 = noise.Evaluate(1, 0);
        double value_x1_y0_z0 = noise.Evaluate(1, 0, 0);
        double value_x1_y1 = noise.Evaluate(1, 1);
        double value_x1_y1_z0 = noise.Evaluate(1, 1, 0);
        double value_x1_y0_z1 = noise.Evaluate(1, 0, 1);
        double value_x1_y1_z1 = noise.Evaluate(1, 1, 1);
        
        Assert.AreEqual(0, value_x0_y0);
        Assert.AreEqual(2.1510994893307009E-67, value_x0_y0_z0);
        Assert.AreEqual(-0.21026991281233706, value_x0_y1);
        Assert.AreEqual(-0.031883015701785991, value_x0_y1_z0);
        Assert.AreEqual(0.036078149346757737, value_x0_y0_z1);
        Assert.AreEqual(0.27409970004637962, value_x1_y0);
        Assert.AreEqual(-0.23584641815494023, value_x1_y0_z0);
        Assert.AreEqual(0.080930776649610278, value_x1_y1);
        Assert.AreEqual(-0.3311358823764432, value_x1_y1_z0);
        Assert.AreEqual(0.62775180790283258, value_x1_y0_z1);
        Assert.AreEqual(0.025889967637540496, value_x1_y1_z1);
    }
    
    [Test]
    public static void SimplexNoiseTest()
    {
        SimplexNoise noise = new SimplexNoise(1234);
        double value_x0_y0 = noise.GetNoise(0, 0);
        double value_x0_y1 = noise.GetNoise(0, 1);
        double value_x1_y0 = noise.GetNoise(1, 0);
        double value_x1_y1 = noise.GetNoise(1, 1);
        
        Assert.AreEqual(0, value_x0_y0);
        Assert.AreEqual(-0.30292468861274469, value_x0_y1);
        Assert.AreEqual(0.29083935882358936, value_x1_y0);
        Assert.AreEqual(0.19624229464828005, value_x1_y1);
    }
    
    [Test]
    public static void PerlinNoiseTest()
    {
        PerlinNoise noise = new(1234);
        double value_x0_y0_z0 = noise.Noise(0, 0, 0);
        double value_x0_y1_z0 = noise.Noise(0, 0.1, 0);
        double value_x0_y0_z1 = noise.Noise(0, 0, 0.1);
        double value_x1_y0_z0 = noise.Noise(0.1, 0, 0);
        double value_x1_y1_z0 = noise.Noise(0.1, 0.1, 0);
        double value_x1_y0_z1 = noise.Noise(0.1, 0, 0.1);
        double value_x1_y1_z1 = noise.Noise(0.1, 0.1, 0.1);
        Assert.AreEqual(0, value_x0_y0_z0);
        Assert.AreEqual(-0.026731927582850909, value_x0_y1_z0);
        Assert.AreEqual(-0.065324357781058356, value_x0_y0_z1);
        Assert.AreEqual(0.066910458280449669, value_x1_y0_z0);
        Assert.AreEqual(0.04097485526117766, value_x1_y1_z0);
        Assert.AreEqual(-0.0010061032904174726, value_x1_y0_z1);
        Assert.AreEqual(-0.022851565377773237, value_x1_y1_z1);
        
        double value_x0_y0_z0_1 = noise.Noise(0, 0, 0);
        double value_x0_y1_z0_1 = noise.Noise(0, 1.1, 0);
        double value_x0_y0_z1_1 = noise.Noise(0, 0, 1.1);
        double value_x1_y0_z0_1 = noise.Noise(1.1, 0, 0);
        double value_x1_y1_z0_1 = noise.Noise(1.1, 1.1, 0);
        double value_x1_y0_z1_1 = noise.Noise(1.1, 0, 1.1);
        double value_x1_y1_z1_1 = noise.Noise(1.1, 1.1, 1.1);
        Assert.AreEqual(0, value_x0_y0_z0_1);
        Assert.AreEqual(0.015682294387549912, value_x0_y1_z0_1);
        Assert.AreEqual(0.097720348513927566, value_x0_y0_z1_1);
        Assert.AreEqual(0.080157433296947633, value_x1_y0_z0_1);
        Assert.AreEqual(-0.13671986449997819, value_x1_y1_z0_1);
        Assert.AreEqual(0.11427810320822852, value_x1_y0_z1_1);
        Assert.AreEqual(-0.12862142076516839, value_x1_y1_z1_1);
    }
}