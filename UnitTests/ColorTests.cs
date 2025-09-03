using color;

namespace UnitTests;

/// <summary>
/// Tests for the Color library, which provides color representations and common variables
/// used by Eto.Forms-based viewports and end-user tools.
/// </summary>
public class ColorTests
{
    [Test]
    public static void MyColor_IntConstructor_CreatesValidColor()
    {
        // Test integer constructor
        MyColor color = new(255, 128, 64);

        Assert.That(color.R, Is.EqualTo(255));
        Assert.That(color.G, Is.EqualTo(128));
        Assert.That(color.B, Is.EqualTo(64));
        Assert.That(color.A, Is.EqualTo(255)); // Default alpha

        // Test float values are correctly calculated
        Assert.That(color.Rf, Is.EqualTo(1.0f).Within(0.001f));
        Assert.That(color.Gf, Is.EqualTo(128 / 255.0f).Within(0.001f));
        Assert.That(color.Bf, Is.EqualTo(64 / 255.0f).Within(0.001f));
        Assert.That(color.Af, Is.EqualTo(1.0f).Within(0.001f));
    }

    [Test]
    public static void MyColor_FloatConstructor_CreatesValidColor()
    {
        // Test float constructor
        MyColor color = new(1.0f, 0.5f, 0.25f);

        Assert.That(color.Rf, Is.EqualTo(1.0f).Within(0.001f));
        Assert.That(color.Gf, Is.EqualTo(0.5f).Within(0.001f));
        Assert.That(color.Bf, Is.EqualTo(0.25f).Within(0.001f));
        Assert.That(color.Af, Is.EqualTo(1.0f).Within(0.001f));

        // Test integer values are correctly calculated
        Assert.That(color.R, Is.EqualTo(255));
        Assert.That(color.G, Is.EqualTo(127)); // 0.5 * 255 = 127.5, truncated to 127
        Assert.That(color.B, Is.EqualTo(63));  // 0.25 * 255 = 63.75, truncated to 63
        Assert.That(color.A, Is.EqualTo(255));
    }

    [Test]
    public static void MyColor_CopyConstructor_CreatesIdenticalColor()
    {
        // Test copy constructor
        MyColor original = new(200, 150, 100);
        MyColor copy = new(original);

        Assert.That(copy.R, Is.EqualTo(original.R));
        Assert.That(copy.G, Is.EqualTo(original.G));
        Assert.That(copy.B, Is.EqualTo(original.B));
        Assert.That(copy.A, Is.EqualTo(original.A));
        Assert.That(copy.Rf, Is.EqualTo(original.Rf).Within(0.001f));
        Assert.That(copy.Gf, Is.EqualTo(original.Gf).Within(0.001f));
        Assert.That(copy.Bf, Is.EqualTo(original.Bf).Within(0.001f));
        Assert.That(copy.Af, Is.EqualTo(original.Af).Within(0.001f));
    }

    [Test]
    public static void MyColor_ToArgb_ReturnsCorrectValue()
    {
        // Test ARGB conversion
        MyColor color = new(255, 128, 64);
        int argb = color.toArgb();

        // ARGB should be 0xFF FF 80 40 = -8355776 (as signed int)
        // Breaking it down: A=255, R=255, G=128, B=64
        string expected = "FF" + "FF" + "80" + "40"; // ARGB format
        int expectedValue = Convert.ToInt32(expected, 16);

        Assert.That(argb, Is.EqualTo(expectedValue));
    }

    [Test]
    public static void MyColor_ToHtml_ReturnsCorrectValue()
    {
        // Test HTML color conversion
        MyColor color = new(255, 128, 64);
        string html = color.ToHtml();

        Assert.That(html, Is.EqualTo("#FF8040"));

        // Test with black color
        MyColor black = new(0, 0, 0);
        Assert.That(black.ToHtml(), Is.EqualTo("#000000"));

        // Test with white color
        MyColor white = new(255, 255, 255);
        Assert.That(white.ToHtml(), Is.EqualTo("#FFFFFF"));
    }

    [Test]
    public static void MyColor_PredefinedColors_HaveCorrectValues()
    {
        // Test some key predefined colors
        Assert.That(MyColor.Black.R, Is.EqualTo(0));
        Assert.That(MyColor.Black.G, Is.EqualTo(0));
        Assert.That(MyColor.Black.B, Is.EqualTo(0));

        Assert.That(MyColor.White.R, Is.EqualTo(255));
        Assert.That(MyColor.White.G, Is.EqualTo(255));
        Assert.That(MyColor.White.B, Is.EqualTo(255));

        Assert.That(MyColor.Red.R, Is.EqualTo(255));
        Assert.That(MyColor.Red.G, Is.EqualTo(0));
        Assert.That(MyColor.Red.B, Is.EqualTo(0));

        Assert.That(MyColor.Green.R, Is.EqualTo(0));
        Assert.That(MyColor.Green.G, Is.EqualTo(128));
        Assert.That(MyColor.Green.B, Is.EqualTo(0));

        Assert.That(MyColor.Blue.R, Is.EqualTo(0));
        Assert.That(MyColor.Blue.G, Is.EqualTo(0));
        Assert.That(MyColor.Blue.B, Is.EqualTo(255));

        // Test custom color
        Assert.That(MyColor.Algae.R, Is.EqualTo(78));
        Assert.That(MyColor.Algae.G, Is.EqualTo(161));
        Assert.That(MyColor.Algae.B, Is.EqualTo(71));
    }

    [Test]
    public static void Colors_Constructor_InitializesDefaults()
    {
        // Test Colors class initialization
        Colors colors = new();

        // Verify some default colors are set
        Assert.That(colors.default_selected_Color, Is.EqualTo(MyColor.Black));
        Assert.That(colors.default_deselected_Color, Is.EqualTo(MyColor.Algae));
        Assert.That(colors.default_enabled_Color, Is.EqualTo(MyColor.Black));
        Assert.That(colors.default_subshape1_Color, Is.EqualTo(MyColor.Red));
        Assert.That(colors.default_subshape2_Color, Is.EqualTo(MyColor.Blue));
        Assert.That(colors.default_subshape3_Color, Is.EqualTo(MyColor.Green));

        // Verify layer colors
        Assert.That(colors.default_layer1_Color, Is.EqualTo(MyColor.Blue));
        Assert.That(colors.default_layer2_Color, Is.EqualTo(MyColor.Orange));
        Assert.That(colors.default_background_Color, Is.EqualTo(MyColor.White));
    }

    [Test]
    public static void Colors_Reset_RestoresDefaults()
    {
        Colors colors = new();

        // Modify some colors
        colors.selected_Color = MyColor.Red;
        colors.layer1_Color = MyColor.Yellow;

        // Reset and verify they're back to defaults
        colors.reset();

        Assert.That(colors.selected_Color.R, Is.EqualTo(MyColor.Black.R));
        Assert.That(colors.selected_Color.G, Is.EqualTo(MyColor.Black.G));
        Assert.That(colors.selected_Color.B, Is.EqualTo(MyColor.Black.B));

        Assert.That(colors.layer1_Color.R, Is.EqualTo(MyColor.Blue.R));
        Assert.That(colors.layer1_Color.G, Is.EqualTo(MyColor.Blue.G));
        Assert.That(colors.layer1_Color.B, Is.EqualTo(MyColor.Blue.B));
    }

    [Test]
    public static void Colors_RebuildLists_CreatesCorrectCollections()
    {
        Colors colors = new();

        // Verify preview colors list has correct number of items
        Assert.That(colors.simPreviewColors.Count, Is.EqualTo(16));
        Assert.That(colors.simPreviewColors[0], Is.EqualTo(colors.layer1_Color));
        Assert.That(colors.simPreviewColors[1], Is.EqualTo(colors.layer2_Color));
        Assert.That(colors.simPreviewColors[15], Is.EqualTo(colors.layer16_Color));

        // Verify output colors list
        Assert.That(colors.simOutputColors.Count, Is.EqualTo(4));
        Assert.That(colors.simOutputColors[0], Is.EqualTo(colors.result_Color));
        Assert.That(colors.simOutputColors[3], Is.EqualTo(colors.result4_Color));

        // Verify result colors array
        Assert.That(colors.resultColors.Length, Is.EqualTo(4));
        Assert.That(colors.resultColors[0], Is.EqualTo(colors.result_Color));
        Assert.That(colors.resultColors[3], Is.EqualTo(colors.result4_Color));
    }

    [Test]
    public static void MyColor_EdgeCases_HandleCorrectly()
    {
        // Test boundary values
        MyColor minColor = new(0, 0, 0);
        Assert.That(minColor.Rf, Is.EqualTo(0.0f));
        Assert.That(minColor.Gf, Is.EqualTo(0.0f));
        Assert.That(minColor.Bf, Is.EqualTo(0.0f));

        MyColor maxColor = new(255, 255, 255);
        Assert.That(maxColor.Rf, Is.EqualTo(1.0f).Within(0.001f));
        Assert.That(maxColor.Gf, Is.EqualTo(1.0f).Within(0.001f));
        Assert.That(maxColor.Bf, Is.EqualTo(1.0f).Within(0.001f));

        // Test float boundary values
        MyColor minFloatColor = new(0.0f, 0.0f, 0.0f);
        Assert.That(minFloatColor.R, Is.EqualTo(0));
        Assert.That(minFloatColor.G, Is.EqualTo(0));
        Assert.That(minFloatColor.B, Is.EqualTo(0));

        MyColor maxFloatColor = new(1.0f, 1.0f, 1.0f);
        Assert.That(maxFloatColor.R, Is.EqualTo(255));
        Assert.That(maxFloatColor.G, Is.EqualTo(255));
        Assert.That(maxFloatColor.B, Is.EqualTo(255));
    }
}