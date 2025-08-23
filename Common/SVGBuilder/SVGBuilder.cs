using System.Globalization;
using Clipper2Lib;
using color;

namespace OLD__SVGBuilder;

/// <summary>
/// A simple class that builds an SVG file with any number of polygons of the specified formats.
/// Provides functionality to create scalable vector graphics from polygon data.
/// </summary>
public class OLD__SVGBuilder
{
    /// <summary>
    /// Defines the styling information for rendering polygons in SVG format.
    /// </summary>
    public class StyleInfo
    {
        /// <summary>
        /// Gets or sets the polygon fill type.
        /// </summary>
        public int pft { get; init; }
        
        /// <summary>
        /// Gets or sets the brush color for filling polygons.
        /// </summary>
        public MyColor brushClr { get; init; }
        
        /// <summary>
        /// Gets or sets the pen color for polygon outlines.
        /// </summary>
        public MyColor penClr { get; init; }
        
        /// <summary>
        /// Gets or sets the pen width for polygon outlines.
        /// </summary>
        public double penWidth { get; init; }
        
        /// <summary>
        /// Gets or sets the dash array pattern for polygon outlines.
        /// </summary>
        public int[] dashArray { get; init; }
        
        /// <summary>
        /// Gets or sets a value indicating whether to show coordinates.
        /// </summary>
        public bool showCoords { get; init; }
        
        /// <summary>
        /// Creates a deep copy of this StyleInfo instance.
        /// </summary>
        /// <returns>A new StyleInfo instance with the same property values.</returns>
        public StyleInfo Clone()
        {
            StyleInfo si = new()
            {
                pft = pft,
                brushClr = brushClr,
                dashArray = dashArray,
                penClr = penClr,
                penWidth = penWidth,
                showCoords = showCoords
            };
            return si;
        }
        
        /// <summary>
        /// Initializes a new instance of the StyleInfo class with default values.
        /// </summary>
        public StyleInfo()
        {
            pft = 0; // PolyFillType.pftNonZero;
            brushClr = MyColor.AntiqueWhite;
            dashArray = [];
            penClr = MyColor.Black;
            penWidth = 0.8;
            showCoords = false;
        }
    }

    /// <summary>
    /// Container for polygon data and associated styling information.
    /// </summary>
    public class PolyInfo
    {
        /// <summary>
        /// Gets or sets the collection of polygon paths.
        /// </summary>
        public PathsD polygons { get; init; }
        
        /// <summary>
        /// Gets or sets the styling information for the polygons.
        /// </summary>
        public StyleInfo si { get; init; }

        /// <summary>
        /// Initializes a new instance of the PolyInfo class.
        /// </summary>
        public PolyInfo()
        {
            polygons = [];
            si = new StyleInfo();
        }
    }

    /// <summary>
    /// Represents a bounding rectangle with top, bottom, left, and right coordinates.
    /// </summary>
    public class BoundingRect
    {
        /// <summary>
        /// Gets or sets the bottom coordinate of the bounding rectangle.
        /// </summary>
        public double bottom { get; set; }
        
        /// <summary>
        /// Gets or sets the top coordinate of the bounding rectangle.
        /// </summary>
        public double top { get; set; }
        
        /// <summary>
        /// Gets or sets the left coordinate of the bounding rectangle.
        /// </summary>
        public double left { get; set; }
        
        /// <summary>
        /// Gets or sets the right coordinate of the bounding rectangle.
        /// </summary>
        public double right { get; set; }
    }

    /// <summary>
    /// Gets or sets the current style configuration for newly added polygons.
    /// </summary>
    public StyleInfo style;
    
    /// <summary>
    /// Internal list containing all polygon information and their associated styles.
    /// </summary>
    private List<PolyInfo> PolyInfoList;

    /// <summary>
    /// SVG file header template with placeholders for width, height, and viewBox dimensions.
    /// </summary>
    private const string svg_header = "<?xml version=\"1.0\" standalone=\"no\"?>\n" +
                                      "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.0//EN\"\n" +
                                      "\"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd\">\n\n" +
                                      "<svg width=\"{0}px\" height=\"{1}px\" viewBox=\"0 0 {2} {3}\" " +
                                      "version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\">\n\n";

    /// <summary>
    /// SVG path formatting template for styling polygon elements.
    /// </summary>
    private const string svg_path_format = "\"\n style=\"fill:{0};" +
                                           " fill-opacity:0; fill-rule:{2}; stroke:{3};" +
                                           " stroke-opacity:{4:f2}; stroke-width:{5:f2};\"/>\n\n";

    /// <summary>
    /// Initializes a new instance of the OLD__SVGBuilder class with default styling.
    /// </summary>
    public OLD__SVGBuilder()
    {
        PolyInfoList = [];
        style = new StyleInfo();
    }
    
    /// <summary>
    /// Adds a single polygon path to the SVG builder using the current style settings.
    /// </summary>
    /// <param name="pointArray">The polygon path to add.</param>
    public void AddPolygons(PathD pointArray)
    {
        PathsD tempPolygonsList = [];
        PathD tempPolygon = new (pointArray);
        tempPolygonsList.Add(new PathD(tempPolygon));
        AddPolygons(tempPolygonsList);
    }

    /// <summary>
    /// Adds a collection of polygon paths to the SVG builder using the current style settings.
    /// </summary>
    /// <param name="poly">The collection of polygon paths to add. If empty, no action is taken.</param>
    public void AddPolygons(PathsD poly)
    {
        if (poly.Count == 0)
        {
            return;
        }

        PolyInfo pi = new() {polygons = poly, si = style.Clone()};
        PolyInfoList.Add(pi);
    }

    /// <summary>
    /// Saves the SVG content to a file with the specified filename, scale, and margin.
    /// </summary>
    /// <param name="filename">The filename to save the SVG content to.</param>
    /// <param name="scale">The scale factor to apply to coordinates (default: 10.0).</param>
    /// <param name="margin">The margin to add around the content in pixels (default: 10).</param>
    /// <returns>True if the file was saved successfully; otherwise, false.</returns>
    public bool SaveToFile(string filename, double scale = 10.0, int margin = 10)
    {
        scale = Math.Max(1, scale);
        margin = Math.Max(0, margin);
        
        //calculate the bounding rect ...
        int i = 0, j = 0;
        while (i < PolyInfoList.Count)
        {
            j = 0;
            while (j < PolyInfoList[i].polygons.Count &&
                   PolyInfoList[i].polygons[j].Count == 0)
            {
                j++;
            }

            if (j < PolyInfoList[i].polygons.Count)
            {
                break;
            }

            i++;
        }
        if (i == PolyInfoList.Count)
        {
            return false;
        }

        BoundingRect rec = new() {left = PolyInfoList[i].polygons[j][0].x};
        rec.right = rec.left;
        rec.top = PolyInfoList[0].polygons[j][0].y;
        rec.bottom = rec.top;

        for (; i < PolyInfoList.Count; i++)
        {
            foreach (PointD pt in PolyInfoList[i].polygons.SelectMany(pg => pg))
            {
                if (pt.x < rec.left)
                {
                    rec.left = pt.x;
                }
                else if (pt.x > rec.right)
                {
                    rec.right = pt.x;
                }

                if (pt.y < rec.top)
                {
                    rec.top = pt.y;
                }
                else if (pt.y > rec.bottom)
                {
                    rec.bottom = pt.y;
                }
            }
        }

        rec.left *= scale;
        rec.top *= scale;
        rec.right *= scale;
        rec.bottom *= scale;
        double offsetX = -rec.left + margin;
        double offsetY = -rec.top + margin;

        using StreamWriter writer = new(filename);
        writer.Write(svg_header,
            rec.right - rec.left + margin * 2,
            rec.bottom - rec.top + margin * 2,
            rec.right - rec.left + margin * 2,
            rec.bottom - rec.top + margin * 2);

        foreach (PolyInfo pi in PolyInfoList)
        {
            writer.Write(" <path d=\"");
            foreach (PathD p in pi.polygons.Where(p => p.Count >= 3))
            {
                writer.Write(string.Format(NumberFormatInfo.InvariantInfo, " M {0:f2} {1:f2}",
                    p[0].x * scale + offsetX,
                    p[0].y * scale + offsetY));
                for (int k = 1; k < p.Count; k++)
                {
                    writer.Write(string.Format(NumberFormatInfo.InvariantInfo, " L {0:f2} {1:f2}",
                        p[k].x * scale + offsetX,
                        p[k].y * scale + offsetY));
                }
                writer.Write(" z");
            }

            writer.Write(string.Format(NumberFormatInfo.InvariantInfo, svg_path_format,
                pi.si.brushClr.ToHtml(),
                (float)pi.si.brushClr.A / 255,
                pi.si.pft == 0,
                pi.si.penClr.ToHtml(),
                (float)pi.si.penClr.A / 255,
                pi.si.penWidth));

            switch (pi.si.showCoords)
            {
                case true:
                {
                    writer.Write("<g font-family=\"Verdana\" font-size=\"11\" fill=\"black\">\n\n");
                    foreach (PathD p in pi.polygons)
                    {
                        foreach (PointD pt in p)
                        {
                            double x = pt.x;
                            double y = pt.y;
                            writer.Write($"<text x=\"{x * scale + offsetX}\" y=\"{y * scale + offsetY}\">{x},{y}</text>\n");

                        }
                        writer.Write("\n");
                    }
                    writer.Write("</g>\n");
                    break;
                }
            }
        }
        writer.Write("</svg>\n");
        return true;
    }
}