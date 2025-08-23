using System;

namespace color;

/// <summary>
/// Custom color implementation to avoid platform-specific quirks of System.Drawing.
/// Provides both integer (0-255) and floating-point (0.0-1.0) color component access.
/// Based on web color standards: https://en.wikipedia.org/wiki/Web_colors
/// </summary>
public class MyColor
{
    // Reimplementation to avoid platform-quirks of System.Drawing.
    // https://en.wikipedia.org/wiki/Web_colors
    public static readonly MyColor AliceBlue = new(240, 248, 255);
    public static readonly MyColor AntiqueWhite = new(250, 235, 215);
    public static readonly MyColor Aqua = new(0, 255, 255);
    public static readonly MyColor Aquamarine = new(127, 255, 212);
    public static readonly MyColor Azure = new(240, 255, 255);
    public static readonly MyColor Beige = new(245, 245, 220);
    public static readonly MyColor Bisque = new(255, 228, 196);
    public static readonly MyColor Black = new(0, 0, 0);
    public static readonly MyColor BlanchedAlmond = new(255, 255, 205);
    public static readonly MyColor Blue = new(0, 0, 255);
    public static readonly MyColor BlueViolet = new(138, 43, 226);
    public static readonly MyColor Brown = new(165, 42, 42);
    public static readonly MyColor BurlyWood = new(222, 184, 135);
    public static readonly MyColor CadetBlue = new(95, 158, 160);
    public static readonly MyColor Chartreuse = new(127, 255, 0);
    public static readonly MyColor Chocolate = new(210, 105, 30);
    public static readonly MyColor Coral = new(210, 105, 30);
    public static readonly MyColor CornflowerBlue = new(100, 149, 237);
    public static readonly MyColor Cornsilk = new(100, 149, 237);
    public static readonly MyColor Crimson = new(220, 20, 60);
    public static readonly MyColor DarkGoldenrod = new(184, 134, 11);
    public static readonly MyColor DarkKhaki = new(189, 183, 107);
    public static readonly MyColor DarkOrange = new(255, 140, 0);
    public static readonly MyColor DarkRed = new(139, 0, 0);
    public static readonly MyColor DarkSalmon = new(233, 150, 122);
    public static readonly MyColor DarkGray = new(169, 169, 169);
    public static readonly MyColor DeepPink = new(255, 20, 147);
    public static readonly MyColor FireBrick = new(178, 34, 34);
    public static readonly MyColor Gold = new(255, 215, 0);
    public static readonly MyColor Goldenrod = new(218, 165, 32);
    public static readonly MyColor Green = new(0, 128, 0);
    public static readonly MyColor HotPink = new(255, 0, 255);
    public static readonly MyColor IndianRed = new(205, 92, 92);
    public static readonly MyColor Indigo = new(75, 0, 130);
    public static readonly MyColor Khaki = new(240, 230, 140);
    public static readonly MyColor LemonChiffon = new(255, 250, 205);
    public static readonly MyColor LightCoral = new(240, 128, 128);
    public static readonly MyColor LightGoldenrodYellow = new(250, 250, 210);
    public static readonly MyColor LightGray = new(211, 211, 211);
    public static readonly MyColor LightPink = new(255, 182, 193);
    public static readonly MyColor LightSalmon = new(255, 160, 122);
    public static readonly MyColor LightYellow = new(255, 255, 224);
    public static readonly MyColor MediumVioletRed = new(199, 21, 133);
    public static readonly MyColor Moccasin = new(255, 228, 181);
    public static readonly MyColor NavajoWhite = new(255, 228, 196);
    public static readonly MyColor Orange = new(255, 165, 0);
    public static readonly MyColor OrangeRed = new(255, 69, 0);
    public static readonly MyColor PaleGoldenrod = new(238, 232, 170);
    public static readonly MyColor PaleVioletRed = new(219, 112, 147);
    public static readonly MyColor PapayaWhip = new(255, 239, 213);
    public static readonly MyColor PeachPuff = new(255, 218, 185);
    public static readonly MyColor Pink = new(255, 192, 203);
    public static readonly MyColor Purple = new(128, 0, 128);
    public static readonly MyColor Red = new(255, 0, 0);
    public static readonly MyColor RosyBrown = new(188, 143, 143);
    public static readonly MyColor Salmon = new(250, 128, 114);
    public static readonly MyColor SandyBrown = new(244, 164, 96);
    public static readonly MyColor SlateGray = new(112, 128, 144);
    public static readonly MyColor Tan = new(210, 180, 140);
    public static readonly MyColor Tomato = new(255, 90, 71);
    public static readonly MyColor Wheat = new(245, 222, 179);
    public static readonly MyColor White = new(255, 255, 255);
    public static readonly MyColor Yellow = new(255, 255, 0);

    public static readonly MyColor Algae = new(78, 161, 71);

    /// <summary>
    /// Gets or sets the red component as a floating-point value (0.0 to 1.0).
    /// </summary>
    public float Rf { get; set; }
    
    /// <summary>
    /// Gets or sets the green component as a floating-point value (0.0 to 1.0).
    /// </summary>
    public float Gf { get; set; }
    
    /// <summary>
    /// Gets or sets the blue component as a floating-point value (0.0 to 1.0).
    /// </summary>
    public float Bf { get; set; }
    
    /// <summary>
    /// Gets or sets the alpha (transparency) component as a floating-point value (0.0 to 1.0).
    /// </summary>
    public float Af { get; set; }

    /// <summary>
    /// Gets or sets the red component as an integer value (0 to 255).
    /// </summary>
    public int R { get; set; }
    
    /// <summary>
    /// Gets or sets the green component as an integer value (0 to 255).
    /// </summary>
    public int G { get; set; }
    
    /// <summary>
    /// Gets or sets the blue component as an integer value (0 to 255).
    /// </summary>
    public int B { get; set; }
    
    /// <summary>
    /// Gets or sets the alpha (transparency) component as an integer value (0 to 255).
    /// </summary>
    public int A { get; set; }

    /// <summary>
    /// Initializes the color components from integer RGB values.
    /// Sets alpha to full opacity (255/1.0).
    /// </summary>
    /// <param name="R">Red component (0-255)</param>
    /// <param name="G">Green component (0-255)</param>
    /// <param name="B">Blue component (0-255)</param>
    private void init(int R, int G, int B)
    {
        this.R = R;
        Rf = R / 255.0f;

        this.G = G;
        Gf = G / 255.0f;

        this.B = B;
        Bf = B / 255.0f;

        A = 255;
        Af = 1.0f;
    }

    /// <summary>
    /// Initializes a new MyColor instance with integer RGB values.
    /// Alpha is set to full opacity (255).
    /// </summary>
    /// <param name="R">Red component (0-255)</param>
    /// <param name="G">Green component (0-255)</param>
    /// <param name="B">Blue component (0-255)</param>
    public MyColor(int R, int G, int B)
    {
        init(R, G, B);
    }

    /// <summary>
    /// Initializes a new MyColor instance with floating-point RGB values.
    /// Alpha is set to full opacity (1.0).
    /// </summary>
    /// <param name="R">Red component (0.0-1.0)</param>
    /// <param name="G">Green component (0.0-1.0)</param>
    /// <param name="B">Blue component (0.0-1.0)</param>
    public MyColor(float R, float G, float B)
    {
        init(R, G, B);
    }

    /// <summary>
    /// Initializes the color components from floating-point RGB values.
    /// Sets alpha to full opacity (255/1.0).
    /// </summary>
    /// <param name="R">Red component (0.0-1.0)</param>
    /// <param name="G">Green component (0.0-1.0)</param>
    /// <param name="B">Blue component (0.0-1.0)</param>
    private void init(float R, float G, float B)
    {
        Rf = R;
        this.R = (int)(R * 255);

        Gf = G;
        this.G = (int)(G * 255);

        Bf = B;
        this.B = (int)(B * 255);

        A = 255;
        Af = 1.0f;
    }

    /// <summary>
    /// Initializes a new MyColor instance by copying values from another MyColor.
    /// </summary>
    /// <param name="sourceColor">The source color to copy from</param>
    public MyColor(MyColor sourceColor)
    {
        init(sourceColor.R, sourceColor.G, sourceColor.B);
    }

    /// <summary>
    /// Converts the color to an ARGB integer value.
    /// The format is 0xAARRGGBB where AA=alpha, RR=red, GG=green, BB=blue.
    /// </summary>
    /// <returns>32-bit ARGB integer representation</returns>
    public int toArgb()
    {
        string aHex = A.ToString("X2");  // gives you hex
        string rHex = R.ToString("X2");  // gives you hex
        string gHex = G.ToString("X2");  // gives you hex
        string bHex = B.ToString("X2");  // gives you hex
        int myNewInt = Convert.ToInt32(aHex + rHex + gHex + bHex, 16);  // back to int again.

        return myNewInt;
    }

    /// <summary>
    /// Converts the color to HTML/CSS hex color format.
    /// </summary>
    /// <returns>HTML color string in #RRGGBB format</returns>
    public string ToHtml()
    {
        return "#" + R.ToString("X2") + G.ToString("X2") + B.ToString("X2");
    }
}