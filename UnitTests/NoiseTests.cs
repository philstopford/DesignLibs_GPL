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

        Assert.That(value_x0_y0, Is.EqualTo(0));
        Assert.That(value_x0_y0_z0, Is.EqualTo(2.1510994893307009E-67));
        Assert.That(value_x0_y1, Is.EqualTo(-0.21026991281233706));
        Assert.That(value_x0_y1_z0, Is.EqualTo(-0.031883015701785991));
        Assert.That(value_x0_y0_z1, Is.EqualTo(0.036078149346757737));
        Assert.That(value_x1_y0, Is.EqualTo(0.27409970004637962));
        Assert.That(value_x1_y0_z0, Is.EqualTo(-0.23584641815494023));
        Assert.That(value_x1_y1, Is.EqualTo(0.080930776649610278));
        Assert.That(value_x1_y1_z0, Is.EqualTo(-0.3311358823764432));
        Assert.That(value_x1_y0_z1, Is.EqualTo(0.62775180790283258));
        Assert.That(value_x1_y1_z1, Is.EqualTo(0.025889967637540496));
    }
    
    [Test]
    public static void SimplexNoiseTest()
    {
        SimplexNoise noise = new SimplexNoise(1234);
        double value_x0_y0 = noise.GetNoise(0, 0);
        double value_x0_y1 = noise.GetNoise(0, 1);
        double value_x1_y0 = noise.GetNoise(1, 0);
        double value_x1_y1 = noise.GetNoise(1, 1);

        Assert.That(value_x0_y0, Is.EqualTo(0));
        Assert.That(value_x0_y1, Is.EqualTo(-0.30292468861274469));
        Assert.That(value_x1_y0, Is.EqualTo(0.29083935882358936));
        Assert.That(value_x1_y1, Is.EqualTo(0.19624229464828005));
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
        Assert.That(value_x0_y0_z0, Is.EqualTo(0));
        Assert.That(value_x0_y1_z0, Is.EqualTo(-0.026731927582850909));
        Assert.That(value_x0_y0_z1, Is.EqualTo(-0.065324357781058356));
        Assert.That(value_x1_y0_z0, Is.EqualTo(0.066910458280449669));
        Assert.That(value_x1_y1_z0, Is.EqualTo(0.04097485526117766));
        Assert.That(value_x1_y0_z1, Is.EqualTo(-0.0010061032904174726));
        Assert.That(value_x1_y1_z1, Is.EqualTo(-0.022851565377773237));
        
        double value_x0_y0_z0_1 = noise.Noise(0, 0, 0);
        double value_x0_y1_z0_1 = noise.Noise(0, 1.1, 0);
        double value_x0_y0_z1_1 = noise.Noise(0, 0, 1.1);
        double value_x1_y0_z0_1 = noise.Noise(1.1, 0, 0);
        double value_x1_y1_z0_1 = noise.Noise(1.1, 1.1, 0);
        double value_x1_y0_z1_1 = noise.Noise(1.1, 0, 1.1);
        double value_x1_y1_z1_1 = noise.Noise(1.1, 1.1, 1.1);
        Assert.That(value_x0_y0_z0_1, Is.EqualTo(0));
        Assert.That(value_x0_y1_z0_1, Is.EqualTo(0.015682294387549912));
        Assert.That(value_x0_y0_z1_1, Is.EqualTo(0.097720348513927566));
        Assert.That(value_x1_y0_z0_1, Is.EqualTo(0.080157433296947633));
        Assert.That(value_x1_y1_z0_1, Is.EqualTo(-0.13671986449997819));
        Assert.That(value_x1_y0_z1_1, Is.EqualTo(0.11427810320822852));
        Assert.That(value_x1_y1_z1_1, Is.EqualTo(-0.12862142076516839));
    }
}