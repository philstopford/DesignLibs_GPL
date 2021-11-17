using System;

namespace color;

public class MyColor
{
    // Reimplementation to avoid platform-quirks of System.Drawing.
    // https://en.wikipedia.org/wiki/Web_colors
    public static MyColor AliceBlue = new(240, 248, 255);
    public static MyColor AntiqueWhite = new(250, 235, 215);
    public static MyColor Aqua = new(0, 255, 255);
    public static MyColor Aquamarine = new(127, 255, 212);
    public static MyColor Azure = new(240, 255, 255);
    public static MyColor Beige = new(245, 245, 220);
    public static MyColor Bisque = new(255, 228, 196);
    public static MyColor Black = new(0, 0, 0);
    public static MyColor BlanchedAlmond = new(255, 255, 205);
    public static MyColor Blue = new(0, 0, 255);
    public static MyColor BlueViolet = new(138, 43, 226);
    public static MyColor Brown = new(165, 42, 42);
    public static MyColor BurlyWood = new(222, 184, 135);
    public static MyColor CadetBlue = new(95, 158, 160);
    public static MyColor Chartreuse = new(127, 255, 0);
    public static MyColor Chocolate = new(210, 105, 30);
    public static MyColor Coral = new(210, 105, 30);
    public static MyColor CornflowerBlue = new(100, 149, 237);
    public static MyColor Cornsilk = new(100, 149, 237);
    public static MyColor Crimson = new(220, 20, 60);
    public static MyColor DarkGoldenrod = new(184, 134, 11);
    public static MyColor DarkKhaki = new(189, 183, 107);
    public static MyColor DarkOrange = new(255, 140, 0);
    public static MyColor DarkRed = new(139, 0, 0);
    public static MyColor DarkSalmon = new(233, 150, 122);
    public static MyColor DarkGray = new(169, 169, 169);
    public static MyColor DeepPink = new(255, 20, 147);
    public static MyColor FireBrick = new(178, 34, 34);
    public static MyColor Gold = new(255, 215, 0);
    public static MyColor Goldenrod = new(218, 165, 32);
    public static MyColor Green = new(0, 128, 0);
    public static MyColor HotPink = new(255, 0, 255);
    public static MyColor IndianRed = new(205, 92, 92);
    public static MyColor Indigo = new(75, 0, 130);
    public static MyColor Khaki = new(240, 230, 140);
    public static MyColor LemonChiffon = new(255, 250, 205);
    public static MyColor LightCoral = new(240, 128, 128);
    public static MyColor LightGoldenrodYellow = new(250, 250, 210);
    public static MyColor LightGray = new(211, 211, 211);
    public static MyColor LightPink = new(255, 182, 193);
    public static MyColor LightSalmon = new(255, 160, 122);
    public static MyColor LightYellow = new(255, 255, 224);
    public static MyColor MediumVioletRed = new(199, 21, 133);
    public static MyColor Moccasin = new(255, 228, 181);
    public static MyColor NavajoWhite = new(255, 228, 196);
    public static MyColor Orange = new(255, 165, 0);
    public static MyColor OrangeRed = new(255, 69, 0);
    public static MyColor PaleGoldenrod = new(238, 232, 170);
    public static MyColor PaleVioletRed = new(219, 112, 147);
    public static MyColor PapayaWhip = new(255, 239, 213);
    public static MyColor PeachPuff = new(255, 218, 185);
    public static MyColor Pink = new(255, 192, 203);
    public static MyColor Purple = new(128, 0, 128);
    public static MyColor Red = new(255, 0, 0);
    public static MyColor RosyBrown = new(188, 143, 143);
    public static MyColor Salmon = new(250, 128, 114);
    public static MyColor SandyBrown = new(244, 164, 96);
    public static MyColor SlateGray = new(112, 128, 144);
    public static MyColor Tan = new(210, 180, 140);
    public static MyColor Tomato = new(255, 90, 71);
    public static MyColor Wheat = new(245, 222, 179);
    public static MyColor White = new(255, 255, 255);
    public static MyColor Yellow = new(255, 255, 0);

    public static MyColor Algae = new(78, 161, 71);

    public float Rf { get; set; }
    public float Gf { get; set; }
    public float Bf { get; set; }
    public float Af { get; set; }

    public int R { get; set; }
    public int G { get; set; }
    public int B { get; set; }
    public int A { get; set; }

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

    public MyColor(int R, int G, int B)
    {
        init(R, G, B);
    }

    public MyColor(float R, float G, float B)
    {
        init(R, G, B);
    }

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

    public MyColor(MyColor sourceColor)
    {
        init(sourceColor.R, sourceColor.G, sourceColor.B);
    }

    public int toArgb()
    {
        string aHex = A.ToString("X2");  // gives you hex
        string rHex = R.ToString("X2");  // gives you hex
        string gHex = G.ToString("X2");  // gives you hex
        string bHex = B.ToString("X2");  // gives you hex
        int myNewInt = Convert.ToInt32(aHex + rHex + gHex + bHex, 16);  // back to int again.

        return myNewInt;
    }

    public string ToHtml()
    {
        return "#" + R.ToString("X2") + G.ToString("X2") + B.ToString("X2");
    }
}