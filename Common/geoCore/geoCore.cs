using gds;
using geoLib;
using geoWrangler;
using oasis;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using Clipper2Lib;

namespace geoCoreLib;

/// <summary>
/// Core geometric design library that provides file format handling and geometric operations.
/// GeoCore is a hybrid design inspired by multiple libraries including layouteditor and KLayout.
/// It handles GDSII and OASIS format parsing and provides integer-based geometry representation
/// to avoid precision issues.
/// </summary>
public class GeoCore
{
    /*
     * GeoCore is a hybrid design. It takes inspiration from a number of other libraries.
     * The inheritance/hierarchical model is inspired by layouteditor (now long since closed source)
     * The parsing code for GDSII and OASIS comes from multiple inspirations (layouteditor, KLayout, gdstk).
     * By design, this tool is all integer based internally. Mapping to the layout grid is deferred to avoid precision issues.
     */
    
    /// <summary>
    /// The main drawing field containing all geometric data and layout information.
    /// </summary>
    private GCDrawingfield drawingField;
    
    /// <summary>
    /// Default tolerance value used for geometric operations and comparisons.
    /// Set to 0.001 to provide appropriate precision for layout operations.
    /// </summary>
    public static double tolerance = 0.001;
    
    /// <summary>
    /// Collection of error messages generated during file parsing or geometric operations.
    /// Used to provide feedback about issues encountered during processing.
    /// </summary>
    public List<string> error_msgs;
    
    /// <summary>
    /// Internal flag indicating whether the GeoCore instance contains valid data.
    /// </summary>
    private bool valid;

    /// <summary>
    /// Sets the validity state of the GeoCore instance.
    /// </summary>
    /// <param name="val">True if the instance contains valid data, false otherwise</param>
    public void setValid(bool val)
    {
        pSetValid(val);
    }

    /// <summary>
    /// Internal implementation for setting the validity state.
    /// </summary>
    /// <param name="val">The validity state to set</param>
    private void pSetValid(bool val)
    {
        valid = val;
    }

    /// <summary>
    /// Gets the current validity state of the GeoCore instance.
    /// </summary>
    /// <returns>True if the instance contains valid data, false otherwise</returns>
    public bool isValid()
    {
        return pIsValid();
    }

    /// <summary>
    /// Internal implementation for getting the validity state.
    /// </summary>
    /// <returns>The current validity state</returns>
    private bool pIsValid()
    {
        return valid;
    }

    /// <summary>
    /// The filename of the currently loaded file.
    /// </summary>
    private string filename;
    
    /// <summary>
    /// Internal identifier for the file format type.
    /// </summary>
    private int fileFormat;
    
    /// <summary>
    /// Enumeration of supported file format types.
    /// </summary>
    public enum fileType 
    { 
        /// <summary>GDSII Stream format</summary>
        gds, 
        /// <summary>OASIS format</summary>
        oasis 
    }

    /// <summary>
    /// Dictionary mapping layer numbers to their string names for improved readability.
    /// </summary>
    private Dictionary<string, string> layerNames;

    /// <summary>
    /// Internal list of structure names in the loaded file.
    /// </summary>
    private List<string> pStructureList;

    /// <summary>
    /// Gets the list of structure names from the loaded file.
    /// </summary>
    /// <returns>List of structure names</returns>
    public List<string> getStructureList()
    {
        return pGetStructureList();
    }

    /// <summary>
    /// Internal implementation for getting the structure list.
    /// </summary>
    /// <returns>List of structure names</returns>
    private List<string> pGetStructureList()
    {
        return pStructureList;
    }

    /// <summary>
    /// Internal list of layer/datatype combinations for the active structure.
    /// </summary>
    private List<string> pActiveStructure_LDList;

    /// <summary>
    /// Gets the layer/datatype list for the currently active structure.
    /// </summary>
    /// <returns>List of layer/datatype combinations</returns>
    public List<string> getActiveStructureLDList()
    {
        return pGetActiveStructureLDList();
    }

    /// <summary>
    /// Internal implementation for getting the active structure's layer/datatype list.
    /// </summary>
    /// <returns>List of layer/datatype combinations</returns>
    private List<string> pGetActiveStructureLDList()
    {
        return pActiveStructure_LDList;
    }

    /// <summary>
    /// Observable collection of structure names for UI binding support.
    /// </summary>
    public ObservableCollection<string> structureList_ { get; set; }
    
    /// <summary>
    /// Observable collection of layer/datatype combinations for the active structure for UI binding.
    /// </summary>
    public ObservableCollection<string> activeStructure_LayerDataTypeList_ { get; set; }

    /// <summary>
    /// Index of the currently active structure.
    /// </summary>
    public int activeStructure { get; set; }
    
    /// <summary>
    /// Index of the currently active layer/datatype combination.
    /// </summary>
    public int activeLD { get; set; }
    
    /// <summary>
    /// List of layer/datatype combinations for each structure in the file.
    /// </summary>
    private List<List<string>> structure_LayerDataTypeList;
    
    /// <summary>
    /// Gets the complete list of layer/datatype combinations for all structures.
    /// </summary>
    /// <returns>Nested list where each inner list contains layer/datatype combinations for a structure</returns>
    public List<List<string>> getStructureLayerDataTypeList()
    {
        return pGetStructureLayerDataTypeList();
    }

    /// <summary>
    /// Internal implementation for getting the structure layer/datatype list.
    /// </summary>
    /// <returns>Nested list of layer/datatype combinations</returns>
    private List<List<string>> pGetStructureLayerDataTypeList()
    {
        return structure_LayerDataTypeList;
    }

    /// <summary>
    /// Internal class representing a layout structure containing geometric elements.
    /// </summary>
    private class Structure
    {
        /// <summary>
        /// List of layer/datatype combinations present in this structure.
        /// </summary>
        public List<string> layerDataTypes;

        /// <summary>
        /// List of geometric elements (polygons, paths, cell references) in this structure.
        /// </summary>
        public List<Element> elements;

        private List<string> bakedGeo_LD;
        private List<BakedGeo> bakedGeo;

        public class BakedGeo
        {
            public string LD { get; set; }
            public PathsD fgeo { get; set; }
            public List<bool> isText { get; set; }

            public BakedGeo()
            {
                fgeo = [];
                isText = [];
            }

            public BakedGeo(PathsD source, List<bool> text, string ld)
            {
                fgeo = new PathsD(source);
                isText = text;
                LD = ld;
            }
        }

        public void addBakedGeo(PathsD source, List<bool> text, string ld)
        {
            for (int i = 0; i < source.Count; i++)
            {
                addBakedGeo(source[i], text[i], ld);
            }
        }

        public void addBakedGeo(PathD source, bool text, string ld)
        {
            int index = bakedGeo_LD.IndexOf(ld);
            switch (index)
            {
                case -1:
                    bakedGeo.Add(new BakedGeo([new PathD(source)], [text], ld));
                    bakedGeo_LD.Add(ld);
                    break;
                default:
                    bakedGeo[index].fgeo.Add(new PathD(source));
                    bakedGeo[index].isText.Add(text);
                    break;
            }
        }

        public PathsD getBakedGeo(string ld)
        {
            int index = bakedGeo_LD.IndexOf(ld);
            return index switch
            {
                -1 => [],
                _ => new PathsD(bakedGeo[index].fgeo)
            };
        }

        public List<GeoLibArray> getArrayData(string ld)
        {
            List<GeoLibArray> ret = [];
            foreach (Element t in elements)
            {
                //if (elements[i].LD == ld)
                {
                    ret.Add(t.arrayData);
                }
            }

            return ret;
        }

        public List<bool> isText(string ld)
        {
            List<bool> ret = [];
            for (int i = 0; i < bakedGeo_LD.Count; i++)
            {
                if (bakedGeo_LD[i] == ld)
                {
                    ret.AddRange(bakedGeo[i].isText);
                }
            }

            return ret;
        }

        public class Element
        {
            public string LD; // might be useful to track originating LD.

            public bool isText { get; set; }
            public string isCellRefArray { get; set; }

            public PathD geometry { get; set; }

            public GeoLibArray arrayData { get; set; }

            public string name { get; set; }

            public Element()
            {
                init();
            }

            private void init()
            {
                geometry = [];
                isText = false;
                name = "";
                isCellRefArray = "";
                arrayData = null;
            }

            public Element(PathD sourceGeo, bool text)
            {
                init(sourceGeo, text);
            }

            private void init(PathD sourceGeo, bool text)
            {
                geometry = new PathD(sourceGeo);
                isText = text;
            }
        }

        public void addPoly(PathD poly, string ldString)
        {
            pAddPoly(poly, ldString);
        }

        private void pAddPoly(PathD poly, string ldString)
        {
            elements.Add(new Element(poly, false));
            elements[^1].LD = ldString;
            switch (layerDataTypes.IndexOf(ldString))
            {
                case -1:
                    layerDataTypes.Add(ldString);
                    break;
            }
        }

        public void addText(string name, PathD poly, string ldString)
        {
            pAddText(name, poly, ldString);
        }

        private void pAddText(string text, PathD poly, string ldString)
        {
            elements.Add(new Element(poly, true));
            elements[^1].name = text;
            elements[^1].LD = ldString;
            switch (layerDataTypes.IndexOf(ldString))
            {
                case -1:
                    layerDataTypes.Add(ldString);
                    break;
            }
        }

        public Structure()
        {
            init();
        }

        private void init()
        {
            elements = [];
            layerDataTypes = [];
            bakedGeo = [];
            bakedGeo_LD = [];
        }

        public void addElement()
        {
            pAddElement();
        }

        private void pAddElement()
        {
            elements.Add(new Element());
        }
    }

    private List<Structure> structures;

    public void updateCollections()
    {
        pUpdateCollections();
    }

    private void pUpdateCollections()
    {
        structureList_.Clear();
        foreach (string t in pStructureList)
        {
            structureList_.Add(t);
        }
        activeStructure_LayerDataTypeList_.Clear();
        foreach (string t in pActiveStructure_LDList)
        {
            activeStructure_LayerDataTypeList_.Add(t);
        }
    }

    public void updateGeoCore(string filename_, fileType type)
    {
        pUpdateGeoCore(filename_, type);
    }

    private void pUpdateGeoCore(string filename_, fileType type)
    {
        reset();
        filename = filename_;

        switch (filename)
        {
            case "":
                return;
        }

        if (!File.Exists(filename))
        {
            throw new Exception("File not accessible: " + filename);
        }

        /*
        structureList.Clear();
        structure_LayerDataTypeList[0].Clear();
        structure_Element_Poly_PointList.Clear();
        */
        switch (type)
        {
            case fileType.gds:
                gdsReader gdsFileData = new(filename);
                valid = gdsFileData.load(ref drawingField);
                fileFormat = (int)fileType.gds;
                error_msgs = gdsFileData.error_msgs;
                switch (valid)
                {
                    case true:
                        processGeometry(ref drawingField, gdsFileData.layerNames);
                        layerNames = gdsFileData.layerNames;
                        activeStructure = drawingField.active_cell;
                        break;
                }
                break;
            case fileType.oasis:
                oasReader oasFileData = new(filename);
                valid = oasFileData.load(ref drawingField);
                fileFormat = (int)fileType.oasis;
                error_msgs = oasFileData.error_msgs;
                switch (valid)
                {
                    case true:
                        processGeometry(ref drawingField, oasFileData.layerNames);
                        layerNames = oasFileData.layerNames;
                        activeStructure = drawingField.active_cell;
                        break;
                }
                break;
            default:
                valid = false;
                break;
        }
        genLDList();
        updateCollections();
    }

    public GCDrawingfield getDrawing()
    {
        return pGetDrawing();
    }

    private GCDrawingfield pGetDrawing()
    {
        return valid switch
        {
            true => drawingField,
            _ => null
        };
    }

    public void setDrawing(GCDrawingfield source)
    {
        drawingField = source.copy();
    }

    public Dictionary<string, string> getLayerNames()
    {
        return pGetLayerNames();
    }

    private Dictionary<string, string> pGetLayerNames()
    {
        return layerNames;
    }

    public void addLayerName(string key, string value)
    {
        pAddLayerName(key, value);
    }

    private void pAddLayerName(string key, string value)
    {
        layerNames.TryAdd(key, value);
    }

    public GeoCore()
    {
        pGeoCore();
    }

    private void pGeoCore()
    {
        drawingField = new GCDrawingfield("");
        pStructureList = [""];
        pActiveStructure_LDList = [""];
        structureList_ = [];
        activeStructure_LayerDataTypeList_ = [];
        error_msgs = [];
        updateCollections();
        reset();
        genLDList();
    }

    public void reset()
    {
        pReset();
    }

    private void pReset()
    {
        drawingField.reset();

        valid = false;
        error_msgs.Clear();
        activeStructure = 0;
        filename = "";

        layerNames = new Dictionary<string, string>();

        pStructureList.Clear();
        pStructureList.Add("");
        pActiveStructure_LDList.Clear();
        pActiveStructure_LDList.Add("");

        structure_LayerDataTypeList = [new List<string>()];
        structure_LayerDataTypeList[0].Add("");

        structures = [new Structure()];
        //structures[0].addElement();

        pUpdateCollections();
    }

    private void genLDList()
    {
        pActiveStructure_LDList.Clear();
        for (int ld = 0; ld < structure_LayerDataTypeList[activeStructure].Count; ld++)
        {
            pActiveStructure_LDList.Add(structure_LayerDataTypeList[activeStructure][ld]);
        }
    }

    public void readValues(GeoCore sourceGeoCore)
    {
        pReadValues(sourceGeoCore);
    }

    private void pReadValues(GeoCore sourceGeoCore)
    {
        reset();
        filename = sourceGeoCore.filename;
        int source = 0;

        pStructureList[source] = sourceGeoCore.pStructureList[source];
        structure_LayerDataTypeList = sourceGeoCore.structure_LayerDataTypeList.ToList();
        source++;

        while (source < sourceGeoCore.pStructureList.Count)
        {
            pStructureList.Add(sourceGeoCore.pStructureList[source]);
            source++;
        }

        drawingField = new GCDrawingfield(sourceGeoCore.getDrawing());

        structures = sourceGeoCore.structures.ToList();
        valid = sourceGeoCore.valid;
        activeStructure = sourceGeoCore.activeStructure;
        layerNames = sourceGeoCore.layerNames;
        genLDList();
    }

    private void processGeometry(ref GCDrawingfield drawing_, Dictionary<string, string> layerNames_)
    {
        // Should have been reset before this call. Remove the defaults.
        pStructureList.Clear();
        structure_LayerDataTypeList[0].Clear();
        structures.Clear();
        
        for (int cell = 0; cell < drawing_.cellList.Count; cell++)
        {
            if (drawing_.cellList[cell].elementList == null)
            {
                continue;
            }

            // Check whether our cell is known already.
            string cellName = drawing_.cellList[cell].cellName;

            int cellIndex = -1;

            try
            {
                cellIndex = pStructureList.IndexOf(cellName);
            }
            catch (Exception)
            {
            }

            switch (cellIndex)
            {
                case -1:
                    pStructureList.Add(cellName); // new cell.
                    cellIndex = pStructureList.Count - 1; // index.

                    // Cell marker is null so we need to skip in this case.
                    structure_LayerDataTypeList.Add([]);

                    structures.Add(new Structure());
                    break;
            }

            // Check whether our layer / datatype combination is known.
            // We have to be careful here - Oasis has the ability to have strings, GDS is numeric.

            List<string> hashList = [];
            for (int element = 0; element < drawing_.cellList[cell].elementList.Count; element++)
            {
                int layer = drawing_.cellList[cell].elementList[element].layer_nr;
                int datatype = drawing_.cellList[cell].elementList[element].datatype_nr;

                // See if our layer/datatype combination is known to us already.

                string searchString = "L" + layer + "D" + datatype;

                // Query our dictionary.
                try
                {
                    if (layerNames_.TryGetValue(searchString, out string resultString))
                    {
                        searchString = resultString;
                    }
                }
                catch (Exception)
                {

                }

                int ldIndex = -1;

                try
                {
                    ldIndex = structure_LayerDataTypeList[cellIndex].IndexOf(searchString);
                }
                catch (Exception)
                {
                }

                if (ldIndex == -1 && searchString != "L-1D-1")
                {
                    structure_LayerDataTypeList[cellIndex].Add(searchString);
                    //structures[cellIndex].addElement();
                    // ldIndex = structure_LayerDataTypeList[cellIndex].Count - 1;
                }

                getGeometry(ref drawing_, cell, element, hashList, cellIndex, searchString);
            }
        }

        // Final clean-up
        structure_LayerDataTypeList.RemoveAt(structure_LayerDataTypeList.Count - 1);
        // Need to remove any layer/datatypes that don't have any entries.

        for (int structure = 0; structure < pStructureList.Count; structure++)
        {
            List<int> indicesToRemove = [];
            for (int ld = structure_LayerDataTypeList[structure].Count - 1; ld > -1; ld--)
            {
                switch (structure_LayerDataTypeList[structure][ld])
                {
                    case "L-1D-1":
                        indicesToRemove.Add(ld);
                        break;
                }
            }

            for (int index = 0; index < indicesToRemove.Count; index++)
            {
                structure_LayerDataTypeList[structure].RemoveAt(index);
                structures[structure].elements.RemoveAt(index);
            }
        }
    }

    private void getGeometry(ref GCDrawingfield drawing_, int cell, int element, List<string> hashList, int cellIndex, string searchString)
    {
        getGeometry(drawing_.cellList[cell], element, hashList, cellIndex, searchString);
    }

    private GeoData getGeometry_simple(GCCell gcCell, int element, List<string> hashList, int cellIndex)
    {
        PathsD ret = [];
        List<string> lds = [];
        List<GCPolygon> lp = gcCell.elementList[element].convertToPolygons(drawingField.getDrawingScale());
        List<bool> isText = [];
        List<string> names = [];
        foreach (GCPolygon p in lp)
        {
            isText.Add(p.isText());
            names.Add(p.getName());

            string ldString = "L" + p.layer_nr + "D" + p.datatype_nr;

            // We should remove identical polygons here in case of doubled-up input geometry.
            string crP_Hash = utility.Utils.GetMD5Hash(p.pointarray);

            if (hashList.IndexOf(crP_Hash) == -1)
            {
                hashList.Add(crP_Hash);
                PathD t = new(p.pointarray.Select(t1 => new PointD(t1.X, t1.Y)));
                if (gcCell.elementList[element].isText())
                {
                    string text = gcCell.elementList[element].getName();
                    structures[cellIndex].addText(text, t, ldString);
                }
                else
                {
                    structures[cellIndex].addPoly(t, ldString);
                }

                ret.Add(t);
                lds.Add(ldString);
            }
        }

        return new GeoData(ret, isText, lds, names);
    }

    private GeoData getGeometry_complex(GCCell gcCell, int element, List<string> hashList, int cellIndex)
    {
        GeoData ret = new();
        // Need to de-reference these cases.
        double angle = 0.0f;
        Point64 point = new();
        double mag = 1.0f;
        GCCell tmpCel = gcCell;
        double xSpace = 0;
        double ySpace = 0;
        int xCount = 1;
        int yCount = 1;

        try
        {
            if (gcCell.elementList[element].isCellref())
            {
                GCCellref refCell = (GCCellref)gcCell.elementList[element];
                point = refCell.getPos();
                mag = refCell.trans.mag;
                angle = refCell.trans.angle;
                tmpCel = refCell.cell_ref;
            }
            else
            {
                GCCellRefArray refCell = (GCCellRefArray)gcCell.elementList[element];
                point = refCell.getPos();
                mag = refCell.trans.mag;
                angle = refCell.trans.angle;
                tmpCel = refCell.cell_ref;
                xSpace = refCell.repetition.rowVector.X + refCell.repetition.colVector.X;
                ySpace = refCell.repetition.colVector.Y + refCell.repetition.rowVector.Y;
                xCount = refCell.repetition.columns;
                yCount = refCell.repetition.rows;
            }
        }
        catch (Exception)
        {
        }

        if (tmpCel != null) // guard against broken cellref
        {
            for (int referenceElement = 0; referenceElement < tmpCel.elementList.Count; referenceElement++)
            {
                // ReSharper disable once CheckForReferenceEqualityInstead.1
                // ReSharper disable once CheckForReferenceEqualityInstead.3
                if (!tmpCel.elementList[referenceElement].isCellref() && !tmpCel.elementList[referenceElement].GetType().Equals(typeof(GCCellRefArray)))
                {
                    ret = getGeometry_2(gcCell, element, hashList, tmpCel, cellIndex, referenceElement, point, xCount, yCount, xSpace, ySpace, angle, mag);
                }
                else
                {
                    ret = getGeometry_complex(tmpCel, referenceElement, hashList, cellIndex);
                }
            }
        }

        // Need to construct the array of our geometry here!
        GeoData postArray = new();
        for (int i = 0; i < ret.geo.Count; i++)
        {
            postArray.geo = GeoWrangler.makeArray(ret.geo[i], xCount, xSpace, yCount, ySpace);
            postArray.ld = [];
            postArray.isText = [];
            for (int j = 0; j < postArray.geo.Count; j++)
            {
                postArray.ld.Add(ret.ld[i]);
                postArray.isText.Add(ret.isText[i]);
            }
        }

        return postArray;

        // Might need to track nested  array configurations here to handle recursive settings.
    }


    private GeoData getGeometry_2(GCCell gcCell, int element, List<string> hashList, GCCell tmpCel, int cellIndex, int referenceElement, Point64 point, int xCount, int yCount, double xSpace, double ySpace, double angle, double mag)
    {
        PathsD ret = [];
        int crLayer = tmpCel.elementList[referenceElement].layer_nr;
        int crDatatype = tmpCel.elementList[referenceElement].datatype_nr;

        List<string> lds = [];

        List<bool> text = [];

        List<string> names = [];

        // See if our layer/datatype combination is known to us already.

        string crSearchString = "L" + crLayer + "D" + crDatatype;
            
        int crLDIndex = -1;

        try
        {
            crLDIndex = structure_LayerDataTypeList[cellIndex].IndexOf(crSearchString);
        }
        catch (Exception)
        {
        }

        switch (crLDIndex)
        {
            case -1:
            {
                structure_LayerDataTypeList[cellIndex].Add(crSearchString);
                structures[cellIndex].addElement();
                int adIndex = structures[cellIndex].elements.Count - 1;

                // ReSharper disable once CheckForReferenceEqualityInstead.1
                if (gcCell.elementList[element].GetType().Equals(typeof(GCCellRefArray)))
                {
                    structures[cellIndex].elements[adIndex].isCellRefArray = gcCell.elementList[element].depend().cellName;
                }

                GeoLibArray tmpArray = new()
                {
                    count = new Point64(xCount, yCount),
                    point = new Point64(point.X, point.Y),
                    pitch = new Point64(xSpace, ySpace)
                };
                structures[cellIndex].elements[adIndex].arrayData = tmpArray;
                break;
            }
        }

        try
        {
            List<GCPolygon> lCRP = tmpCel.elementList[referenceElement].convertToPolygons(drawingField.getDrawingScale());
            foreach (GCPolygon crP in lCRP)
            {
                crP.move(point);
                crP.rotate(angle, point);
                crP.scale(point, mag);

                string ldString = "L" + crP.layer_nr + "D" + crP.datatype_nr;

                // We should remove identical polygons here in case of doubled-up input geometry.
                string crP_Hash = utility.Utils.GetMD5Hash(crP.pointarray);

                switch (hashList.IndexOf(crP_Hash))
                {
                    case -1:
                    {
                        hashList.Add(crP_Hash);

                        const int x = 0;
                        const int y = 0;
                        PathD t = new (crP.pointarray.Select(t1 => new PointD((t1.X + x * xSpace), (t1.Y + y * ySpace))));
                        structures[cellIndex].addPoly(t, ldString);
                        ret.Add(t);
                        text.Add(crP.isText());
                        lds.Add(ldString);
                        names.Add(crP.getName());
                        break;
                    }
                }
            }
        }
        catch (Exception)
        {
            // Exception is not a big deal.
        }

        return new GeoData(ret, text, lds, names);
    }


    public class GeoData
    {
        public PathsD geo;
        public List<string> ld;
        public List<bool> isText;
        public List<string> name;

        public GeoData()
        {
            geo = [];
            ld = [];
            isText = [];
            name = [];
        }

        public GeoData(PathsD poly, List<bool> text, List<string> LDs, List<string> names)
        {
            geo = new PathsD(poly);
            ld = LDs.ToList();
            isText = text.ToList();
            name = names.ToList();
        }
    }

    private void getGeometry(GCCell gcCell, int element, List<string> hashList, int cellIndex, string searchString)
    {
        GeoData ret;
        // Now we have to process our geometry into the right place.
        if (!gcCell.elementList[element].isCellref() && !gcCell.elementList[element].isCellrefArray())
        {
            ret = getGeometry_simple(gcCell, element, hashList, cellIndex);
        }
        else
        {
            ret = getGeometry_complex(gcCell, element, hashList, cellIndex);
        }
            
        for (int i = 0; i < ret.geo.Count; i++)
        {
            structures[cellIndex].addBakedGeo(new PathD(ret.geo[i]), ret.isText[i], ret.ld[i]);
        }
    }

    public List<GeoLibArray> getArrayParameters()
    {
        return pGetArrayParameters();
    }

    private List<GeoLibArray> pGetArrayParameters()
    {
        return structures[activeStructure].getArrayData(pActiveStructure_LDList[activeLD]);
    }

    public List<bool> isText()
    {
        return pIsText();
    }

    private List<bool> pIsText()
    {
        //string text = pActiveStructure_LDList[activeLD];
        //return structures[activeStructure].isText(text);

        List<GCPolygon> tmp = convertToPolygons(true);

        return tmp.Select(t => t.isText()).ToList();

    }

    public List<string> names()
    {
        return pNames();
    }

    private List<string> pNames()
    {
        //string text = pActiveStructure_LDList[activeLD];
        //return structures[activeStructure].isText(text);

        List<GCPolygon> tmp = convertToPolygons(true);

        return tmp.Select(t => t.getName()).ToList();

    }

    public PathsD points(bool flatten)
    {
        return pPoints(flatten);
    }

    private PathsD pPoints(bool flatten)
    {
        PathsD points = [];

        switch (flatten)
        {
            case true:
                List<GCPolygon> tmp = convertToPolygons(true);
                points = new PathsD(tmp.Select(t => GeoWrangler.PathDFromPath64(t.pointarray)));
                break;
            default:
                // Do we ever get here?

                Path64 array_count = [];
                PathD array_pitch = [];

                points.Add(new PathD(structures[activeStructure].elements[activeLD].geometry));
                if (structures[activeStructure].elements[activeLD].arrayData != null)
                {
                    array_count.Add(new Point64(structures[activeStructure].elements[activeLD].arrayData.count.X,structures[activeStructure].elements[activeLD].arrayData.count.Y));
                    array_pitch.Add(new PointD(structures[activeStructure].elements[activeLD].arrayData.pitch.X, structures[activeStructure].elements[activeLD].arrayData.pitch.Y));
                }
                else
                {
                    array_count.Add(new Point64(1, 1));
                    array_pitch.Add(new PointD(0, 0));
                }

                break;
        }

        double resizeFactor = fileFormat switch
        {
            (int)fileType.gds => drawingField.userunits / 1E-3,
            (int)fileType.oasis => 1000.0 / drawingField.databaseunits,
            _ => 1.0
        };

        points = GeoWrangler.resize(points, resizeFactor);
        // points = GeoWrangler.removeDuplicates(points);
        return points;
    }

    public string getName(int index)
    {
        return pGetName(index);
    }

    private string pGetName(int index)
    {
        return structures[activeStructure].elements[index].name;
    }

    public void updateGeometry(int structure, int layerdatatype)
    {
        pUpdateGeometry(structure, layerdatatype);
    }

    private void pUpdateGeometry(int structure, int layerdatatype)
    {
        switch (valid)
        {
            case true:
                drawingField.active_cell = structure;
                activeStructure = structure;
                activeLD = layerdatatype;
                genLDList();
                break;
        }
    }

    public List<GCPolygon> convertToPolygons(int layer, int datatype)
    {
        return pConvertToPolygons(layer: layer, datatype: datatype);
    }

    public List<GCPolygon> convertToPolygons(bool activeLDOnly = false)
    {
        return pConvertToPolygons(activeLDOnly);
    }

    private List<GCPolygon> pConvertToPolygons(bool activeLDOnly = false)
    {
        int layer = -1;
        int datatype = -1;
        switch (activeLDOnly)
        {
            case true:
            {
                // Recover our layer and datatype from the string representation.
                string ld = layerNames.FirstOrDefault(x => x.Value == pActiveStructure_LDList[activeLD]).Key;
                string[] temp = ld.Split(['L'])[1].Split(['D']);
                layer = Convert.ToInt32(temp[0]);
                datatype = Convert.ToInt32(temp[1]);
                break;
            }
        }
        return pConvertToPolygons(layer: layer, datatype: datatype);
    }

    private List<GCPolygon> pConvertToPolygons(int layer = -1, int datatype = -1)
    {
        return drawingField.cellList[activeStructure].convertToPolygons(drawingField.getDrawingScale(), layer: layer, datatype: datatype);
    }

    public bool nestedCellRef(int cellIndex)
    {
        bool ret = false;

        for (int i = 0; i < drawingField.cellList[activeStructure].elementList.Count; i++)
        {
            ret = nestedCellRef(cellIndex, i);

            if (ret)
            {
                break;
            }
        }

        return ret;
    }

    public bool nestedCellRef(int cellIndex, int elementIndex)
    {
        return nestedCellRef(drawingField.cellList[cellIndex], elementIndex);
    }

    private static bool nestedCellRef(GCCell cell, int elementIndex)
    {
        bool ret = false;
        if (!cell.elementList[elementIndex].isCellref() && !cell.elementList[elementIndex].isCellrefArray())
        {
            return ret;
        }

        GCCell rCell = cell.elementList[elementIndex].getCellref();
        if (rCell.elementList.Any(t => t.isCellref() || t.isCellrefArray()))
        {
            ret = true;
        }

        return ret;
    }

}