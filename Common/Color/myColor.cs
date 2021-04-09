using System;

namespace color
{
    public class MyColor
    {
        // Reimplementation to avoid platform-quirks of System.Drawing.
        // https://en.wikipedia.org/wiki/Web_colors
        public static MyColor AliceBlue = new MyColor(240, 248, 255);
        public static MyColor AntiqueWhite = new MyColor(250, 235, 215);
        public static MyColor Aqua = new MyColor(0, 255, 255);
        public static MyColor Aquamarine = new MyColor(127, 255, 212);
        public static MyColor Azure = new MyColor(240, 255, 255);
        public static MyColor Beige = new MyColor(245, 245, 220);
        public static MyColor Bisque = new MyColor(255, 228, 196);
        public static MyColor Black = new MyColor(0, 0, 0);
        public static MyColor BlanchedAlmond = new MyColor(255, 255, 205);
        public static MyColor Blue = new MyColor(0, 0, 255);
        public static MyColor BlueViolet = new MyColor(138, 43, 226);
        public static MyColor Brown = new MyColor(165, 42, 42);
        public static MyColor BurlyWood = new MyColor(222, 184, 135);
        public static MyColor CadetBlue = new MyColor(95, 158, 160);
        public static MyColor Chartreuse = new MyColor(127, 255, 0);
        public static MyColor Chocolate = new MyColor(210, 105, 30);
        public static MyColor Coral = new MyColor(210, 105, 30);
        public static MyColor CornflowerBlue = new MyColor(100, 149, 237);
        public static MyColor Cornsilk = new MyColor(100, 149, 237);
        public static MyColor Crimson = new MyColor(220, 20, 60);
        public static MyColor DarkGoldenrod = new MyColor(184, 134, 11);
        public static MyColor DarkKhaki = new MyColor(189, 183, 107);
        public static MyColor DarkOrange = new MyColor(255, 140, 0);
        public static MyColor DarkRed = new MyColor(139, 0, 0);
        public static MyColor DarkSalmon = new MyColor(233, 150, 122);
        public static MyColor DarkGray = new MyColor(169, 169, 169);
        public static MyColor DeepPink = new MyColor(255, 20, 147);
        public static MyColor FireBrick = new MyColor(178, 34, 34);
        public static MyColor Gold = new MyColor(255, 215, 0);
        public static MyColor Goldenrod = new MyColor(218, 165, 32);
        public static MyColor Green = new MyColor(0, 128, 0);
        public static MyColor HotPink = new MyColor(255, 0, 255);
        public static MyColor IndianRed = new MyColor(205, 92, 92);
        public static MyColor Indigo = new MyColor(75, 0, 130);
        public static MyColor Khaki = new MyColor(240, 230, 140);
        public static MyColor LemonChiffon = new MyColor(255, 250, 205);
        public static MyColor LightCoral = new MyColor(240, 128, 128);
        public static MyColor LightGoldenrodYellow = new MyColor(250, 250, 210);
        public static MyColor LightGray = new MyColor(211, 211, 211);
        public static MyColor LightPink = new MyColor(255, 182, 193);
        public static MyColor LightSalmon = new MyColor(255, 160, 122);
        public static MyColor LightYellow = new MyColor(255, 255, 224);
        public static MyColor MediumVioletRed = new MyColor(199, 21, 133);
        public static MyColor Moccasin = new MyColor(255, 228, 181);
        public static MyColor NavajoWhite = new MyColor(255, 228, 196);
        public static MyColor Orange = new MyColor(255, 165, 0);
        public static MyColor OrangeRed = new MyColor(255, 69, 0);
        public static MyColor PaleGoldenrod = new MyColor(238, 232, 170);
        public static MyColor PaleVioletRed = new MyColor(219, 112, 147);
        public static MyColor PapayaWhip = new MyColor(255, 239, 213);
        public static MyColor PeachPuff = new MyColor(255, 218, 185);
        public static MyColor Pink = new MyColor(255, 192, 203);
        public static MyColor Purple = new MyColor(128, 0, 128);
        public static MyColor Red = new MyColor(255, 0, 0);
        public static MyColor RosyBrown = new MyColor(188, 143, 143);
        public static MyColor Salmon = new MyColor(250, 128, 114);
        public static MyColor SandyBrown = new MyColor(244, 164, 96);
        public static MyColor SlateGray = new MyColor(112, 128, 144);
        public static MyColor Tan = new MyColor(210, 180, 140);
        public static MyColor Tomato = new MyColor(255, 90, 71);
        public static MyColor Wheat = new MyColor(245, 222, 179);
        public static MyColor White = new MyColor(255, 255, 255);
        public static MyColor Yellow = new MyColor(255, 255, 0);

        public static MyColor Algae = new MyColor(78, 161, 71);

        public float Rf { get; set; }
        public float Gf { get; set; }
        public float Bf { get; set; }
        public float Af { get; set; }

        public Int32 R { get; set; }
        public Int32 G { get; set; }
        public Int32 B { get; set; }
        public Int32 A { get; set; }

        void init(Int32 R, Int32 G, Int32 B)
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

        public MyColor(Int32 R, Int32 G, Int32 B)
        {
            init(R, G, B);
        }

        public MyColor(float R, float G, float B)
        {
            init(R, G, B);
        }

        void init(float R, float G, float B)
        {
            Rf = R;
            this.R = (Int32)(R * 255);

            Gf = G;
            this.G = (Int32)(G * 255);

            Bf = B;
            this.B = (Int32)(B * 255);

            A = 255;
            Af = 1.0f;
        }

        public MyColor(MyColor sourceColor)
        {
            init(sourceColor.R, sourceColor.G, sourceColor.B);
        }

        public Int32 toArgb()
        {
            string aHex = A.ToString("X2");  // gives you hex
            string rHex = R.ToString("X2");  // gives you hex
            string gHex = G.ToString("X2");  // gives you hex
            string bHex = B.ToString("X2");  // gives you hex
            Int32 myNewInt = Convert.ToInt32(aHex + rHex + gHex + bHex, 16);  // back to int again.

            return myNewInt;
        }

        public string ToHtml()
        {
            return "#" + R.ToString("X2") + G.ToString("X2") + B.ToString("X2");
        }
    }
}
