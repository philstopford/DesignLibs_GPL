using System.Collections.Generic;

namespace color
{
    public class Colors
    {
        public MyColor subshape1_Color { get; set; }
        public MyColor subshape2_Color { get; set; }
        public MyColor subshape3_Color { get; set; }
        public MyColor enabled_Color { get; set; }

        public MyColor selected_Color { get; set; }

        public MyColor deselected_Color { get; set; }

        public MyColor layer1_Color { get; set; }
        public MyColor layer2_Color { get; set; }
        public MyColor layer3_Color { get; set; }
        public MyColor layer4_Color { get; set; }
        public MyColor layer5_Color { get; set; }
        public MyColor layer6_Color { get; set; }
        public MyColor layer7_Color { get; set; }
        public MyColor layer8_Color { get; set; }
        public MyColor layer9_Color { get; set; }
        public MyColor layer10_Color { get; set; }
        public MyColor layer11_Color { get; set; }
        public MyColor layer12_Color { get; set; }
        public MyColor layer13_Color { get; set; }
        public MyColor layer14_Color { get; set; }
        public MyColor layer15_Color { get; set; }
        public MyColor layer16_Color { get; set; }
        public MyColor result_Color { get; set; }
        public MyColor result2_Color { get; set; }
        public MyColor result3_Color { get; set; }
        public MyColor result4_Color { get; set; }

        public MyColor major_Color { get; set; }
        public MyColor minor_Color { get; set; }
        public MyColor background_Color { get; set; }
        public MyColor axis_Color { get; set; }

        public MyColor implantResist_Color { get; set; }
        public MyColor implantMin_Color { get; set; }
        public MyColor implantMean_Color { get; set; }
        public MyColor implantMax_Color { get; set; }

        public MyColor extents_Color { get; set; }

        public MyColor default_subshape1_Color { get; set; }
        public MyColor default_subshape2_Color { get; set; }
        public MyColor default_subshape3_Color { get; set; }
        public MyColor default_enabled_Color { get; set; }

        public MyColor default_selected_Color { get; set; }

        public MyColor default_deselected_Color { get; set; }

        public MyColor default_layer1_Color { get; set; }
        public MyColor default_layer2_Color { get; set; }
        public MyColor default_layer3_Color { get; set; }
        public MyColor default_layer4_Color { get; set; }
        public MyColor default_layer5_Color { get; set; }
        public MyColor default_layer6_Color { get; set; }
        public MyColor default_layer7_Color { get; set; }
        public MyColor default_layer8_Color { get; set; }
        public MyColor default_layer9_Color { get; set; }
        public MyColor default_layer10_Color { get; set; }
        public MyColor default_layer11_Color { get; set; }
        public MyColor default_layer12_Color { get; set; }
        public MyColor default_layer13_Color { get; set; }
        public MyColor default_layer14_Color { get; set; }
        public MyColor default_layer15_Color { get; set; }
        public MyColor default_layer16_Color { get; set; }
        public MyColor default_result_Color { get; set; }
        public MyColor default_result2_Color { get; set; }
        public MyColor default_result3_Color { get; set; }
        public MyColor default_result4_Color { get; set; }

        public MyColor default_major_Color { get; set; }
        public MyColor default_minor_Color { get; set; }
        public MyColor default_axis_Color { get; set; }
        public MyColor default_background_Color { get; set; }

        public List<MyColor> simOutputColors { get; set; }
        public List<MyColor> simPreviewColors { get; set; }
        public MyColor[] resultColors { get; set; }

        public MyColor default_implantResist_Color { get; set; }
        public MyColor default_implantMin_Color { get; set; }
        public MyColor default_implantMean_Color { get; set; }
        public MyColor default_implantMax_Color { get; set; }

        public MyColor default_extents_Color { get; set; }

        public Colors()
        {
            pColors();
        }

        void pColors()
        {
            default_selected_Color = MyColor.Black;

            default_deselected_Color = MyColor.Algae;

            default_enabled_Color = MyColor.Black;
            default_subshape1_Color = MyColor.Red;
            default_subshape2_Color = MyColor.Blue;
            default_subshape3_Color = MyColor.Green;

            default_layer1_Color = MyColor.Blue;
            default_layer2_Color = MyColor.Orange;
            default_layer3_Color = MyColor.Indigo;
            default_layer4_Color = MyColor.Cornsilk;

            default_layer5_Color = MyColor.Chocolate;
            default_layer6_Color = MyColor.Moccasin;
            default_layer7_Color = MyColor.CornflowerBlue;
            default_layer8_Color = MyColor.MediumVioletRed;

            default_layer9_Color = MyColor.Green;
            default_layer10_Color = MyColor.Purple;
            default_layer11_Color = MyColor.CadetBlue;
            default_layer12_Color = MyColor.BlueViolet;

            default_layer13_Color = MyColor.DarkGray;
            default_layer14_Color = MyColor.RosyBrown;
            default_layer15_Color = MyColor.DarkKhaki;
            default_layer16_Color = MyColor.SandyBrown;

            default_result_Color = MyColor.Red;
            default_result2_Color = MyColor.DarkRed;
            default_result3_Color = MyColor.Salmon;
            default_result4_Color = MyColor.DarkSalmon;

            default_major_Color = MyColor.SlateGray;
            default_minor_Color = MyColor.LightGray;
            default_axis_Color = MyColor.Black;
            default_background_Color = MyColor.White;

            default_implantResist_Color = MyColor.Blue;
            default_implantMax_Color = MyColor.Red;
            default_implantMean_Color = MyColor.Green;
            default_implantMin_Color = MyColor.DarkRed;

            default_extents_Color = MyColor.HotPink;

            reset();
        }

        public void rebuildLists()
        {
            pRebuildLists();
        }

        void pRebuildLists()
        {
            simPreviewColors = new List<MyColor>();
            simPreviewColors.Add(layer1_Color);
            simPreviewColors.Add(layer2_Color);
            simPreviewColors.Add(layer3_Color);
            simPreviewColors.Add(layer4_Color);
            simPreviewColors.Add(layer5_Color);
            simPreviewColors.Add(layer6_Color);
            simPreviewColors.Add(layer7_Color);
            simPreviewColors.Add(layer8_Color);
            simPreviewColors.Add(layer9_Color);
            simPreviewColors.Add(layer10_Color);
            simPreviewColors.Add(layer11_Color);
            simPreviewColors.Add(layer12_Color);
            simPreviewColors.Add(layer13_Color);
            simPreviewColors.Add(layer14_Color);
            simPreviewColors.Add(layer15_Color);
            simPreviewColors.Add(layer16_Color);

            simOutputColors = new List<MyColor>();
            simOutputColors.Add(result_Color);
            simOutputColors.Add(result2_Color);
            simOutputColors.Add(result3_Color);
            simOutputColors.Add(result4_Color);

            resultColors = new MyColor[4];
            resultColors[0] = result_Color;
            resultColors[1] = result2_Color;
            resultColors[2] = result3_Color;
            resultColors[3] = result4_Color;
        }

        public void reset()
        {
            pReset();
        }

        void pReset()
        {
            selected_Color = new MyColor(default_selected_Color);
            deselected_Color = new MyColor(default_deselected_Color);
            enabled_Color = new MyColor(default_enabled_Color);
            subshape1_Color = new MyColor(default_subshape1_Color);
            subshape2_Color = new MyColor(default_subshape2_Color);
            subshape3_Color = new MyColor(default_subshape3_Color);

            layer1_Color = new MyColor(default_layer1_Color);
            layer2_Color = new MyColor(default_layer2_Color);
            layer3_Color = new MyColor(default_layer3_Color);
            layer4_Color = new MyColor(default_layer4_Color);
            layer5_Color = new MyColor(default_layer5_Color);
            layer6_Color = new MyColor(default_layer6_Color);
            layer7_Color = new MyColor(default_layer7_Color);
            layer8_Color = new MyColor(default_layer8_Color);

            layer9_Color = new MyColor(default_layer9_Color);
            layer10_Color = new MyColor(default_layer10_Color);
            layer11_Color = new MyColor(default_layer11_Color);
            layer12_Color = new MyColor(default_layer12_Color);
            layer13_Color = new MyColor(default_layer13_Color);
            layer14_Color = new MyColor(default_layer14_Color);
            layer15_Color = new MyColor(default_layer15_Color);
            layer16_Color = new MyColor(default_layer16_Color);

            result_Color = new MyColor(default_result_Color);
            result2_Color = new MyColor(default_result2_Color);
            result3_Color = new MyColor(default_result3_Color);
            result4_Color = new MyColor(default_result4_Color);

            axis_Color = new MyColor(default_axis_Color);
            major_Color = new MyColor(default_major_Color);
            minor_Color = new MyColor(default_minor_Color);
            background_Color = new MyColor(default_background_Color);

            implantResist_Color = new MyColor(default_implantResist_Color);
            implantMin_Color = new MyColor(default_implantMin_Color);
            implantMean_Color = new MyColor(default_implantMean_Color);
            implantMax_Color = new MyColor(default_implantMax_Color);

            extents_Color = new MyColor(default_extents_Color);

            rebuildLists();
        }
    }
}
