using geoLib;
using geoWrangler;

namespace shapeEngine;

[Serializable]
public class ShapeSettings
{
    public static List<string> getAvailableTipsLocations()
    {
        return pGetAvailableTipsLocations();
    }

    private static List<string> pGetAvailableTipsLocations()
    {
        return availableTipsLocations;
    }

    public enum tipLocations { none, L, R, LR, T, B, TB, TL, TR, TLR, BL, BR, BLR, TBL, TBR, all }
    
    public static List<string> getAvailableSubShapePositions()
    {
        return pGetAvailableSubShapePositions();
    }

    private static List<string> pGetAvailableSubShapePositions()
    {
        return availableSubShapePositions;
    }

    public enum subShapeLocations { TL, TR, BL, BR, TS, RS, BS, LS, C }
    
    public static List<string> getPolyFillTypes()
    {
        return pGetPolyFillTypes();
    }

    private static List<string> pGetPolyFillTypes()
    {
        return polyFillTypes;
    }

    private static readonly List<string> availableSubShapePositions = new List<string>
    { "Corner: Top Left", "Corner: Top Right", "Corner: Bottom Left", "Corner: Bottom Right",
        "Middle: Top Side", "Middle: Right Side", "Middle: Bottom Side", "Middle: Left Side",
        "Center"};

    private static readonly List<string> availableTipsLocations = new List<string>
    { "None", "Left", "Right", "Left & Right",
        "Top", "Bottom", "Top & Bottom",
        "Top & Left", "Top & Right", "Top & Left & Right",
        "Bottom & Left", "Bottom & Right", "Bottom & Left & Right",
        "Top & Bottom & Left", "Top & Bottom & Right",
        "All"};

    private static readonly List<string> polyFillTypes = new List<string> { "Even/Odd", "Non-zero", "Positive", "Negative" };

    public enum PolyFill { pftEvenOdd, pftNonZero, pftPositive, pftNegative }

    public enum typeShapes { none, rectangle, L, T, X, U, S, GEOCORE, BOOLEAN }

    private static readonly int default_subShapeTipLocIndex = 0;
    private static readonly int default_subShape2TipLocIndex = 0;

    private static readonly int default_subShape3TipLocIndex = 0;
    private static readonly int default_enabled = 0;
    private static readonly int default_geoCoreShapeEngine = 0;

    private static readonly int default_shapeIndex = (int)typeShapes.none;
    private static readonly int default_subShapeRefIndex = 0;
    private static readonly int default_posInSubShapeIndex = (int)subShapeLocations.BL;

    private static readonly int default_edgeSlide = 0;
    private static readonly int default_proximitySideRays = 2;
    private static readonly int default_proximitySideRaysFallOff = 0;

    private static readonly decimal default_LWR = 0;
    private static readonly int default_LWRNoiseType = (int)NoiseC.noiseIndex.perlin;
    private static readonly int default_LWRNoisePreview = 0;
    private static readonly decimal default_LWRNoiseFreq = 0.2m;

    private int enabled = default_enabled;
    private int shapeIndex = default_shapeIndex;
    private int gCSEngine = default_geoCoreShapeEngine;
    private int subShapeTipLocIndex = default_subShapeTipLocIndex;
    private int subShape2TipLocIndex = default_subShape2TipLocIndex;
    private int subShape3TipLocIndex = default_subShape3TipLocIndex;
    private int subShapeRefIndex = default_subShapeRefIndex;
    private int posInSubShapeIndex = default_posInSubShapeIndex;
    private int edgeSlide = default_edgeSlide;
    private int proximitySideRays = default_proximitySideRays;
    private int LWRNoiseType = default_LWRNoiseType;
    private int LWR2NoiseType = default_LWRNoiseType;
    private int proxSideRaysFallOff = default_proximitySideRaysFallOff;

    private static string default_layerName = "";
    
    public enum properties_i
    {
        enabled,
        gCSEngine,
        shapeIndex,
        subShapeTipLocIndex, subShape2TipLocIndex, subShape3TipLocIndex,
        subShapeRefIndex,posInSubShapeIndex,
        edgeSlide,
        proxRays,proxSideRaysFallOff,
        lwrType, lwr2Type
    }

    public int getInt(properties_i p)
    {
        return pGetInt(p);
    }

    int pGetInt(properties_i p)
    {
        int ret = 0;
        switch (p)
        {
            case properties_i.enabled:
                ret = enabled;
                break;
            case properties_i.gCSEngine:
                ret = gCSEngine;
                break;
            case properties_i.shapeIndex:
                ret = shapeIndex;
                break;
            case properties_i.subShapeTipLocIndex:
                ret = subShapeTipLocIndex;
                break;
            case properties_i.subShape2TipLocIndex:
                ret = subShape2TipLocIndex;
                break;
            case properties_i.subShape3TipLocIndex:
                ret = subShape3TipLocIndex;
                break;
            case properties_i.subShapeRefIndex:
                ret = subShapeRefIndex;
                break;
            case properties_i.posInSubShapeIndex:
                ret = posInSubShapeIndex;
                break;
            case properties_i.edgeSlide:
                ret = edgeSlide;
                break;
            case properties_i.proxRays:
                ret = proximitySideRays;
                break;
            case properties_i.proxSideRaysFallOff:
                ret = proxSideRaysFallOff;
                break;
            case properties_i.lwrType:
                ret = LWRNoiseType;
                break;
            case properties_i.lwr2Type:
                ret = LWR2NoiseType;
                break;
        }

        return ret;
    }

    public void setInt(properties_i p, int val)
    {
        pSetInt(p, val);
    }

    void pSetInt(properties_i p, int val)
    {
        switch (p)
        {
            case properties_i.enabled:
                enabled = val;
                break;
            case properties_i.gCSEngine:
                gCSEngine = val;
                break;
            case properties_i.shapeIndex:
                shapeIndex = val;
                break;
            case properties_i.subShapeTipLocIndex:
                subShapeTipLocIndex = val;
                break;
            case properties_i.subShape2TipLocIndex:
                subShape2TipLocIndex = val;
                break;
            case properties_i.subShape3TipLocIndex:
                subShape3TipLocIndex = val;
                break;
            case properties_i.subShapeRefIndex:
                subShapeRefIndex = val;
                break;
            case properties_i.posInSubShapeIndex:
                posInSubShapeIndex = val;
                break;
            case properties_i.edgeSlide:
                edgeSlide = val;
                break;
            case properties_i.proxRays:
                proximitySideRays = val;
                break;
            case properties_i.proxSideRaysFallOff:
                proxSideRaysFallOff = val;
                break;
            case properties_i.lwrType:
                LWRNoiseType = val;
                break;
            case properties_i.lwr2Type:
                LWR2NoiseType = val;
                break;
        }
    }

    public void defaultInt(properties_i p)
    {
        pDefaultInt(p);
    }

    private void pDefaultInt(properties_i p)
    {
        switch (p)
        {
            case properties_i.enabled:
                enabled = default_enabled;
                break;
            case properties_i.gCSEngine:
                gCSEngine = default_geoCoreShapeEngine;
                break;
            case properties_i.posInSubShapeIndex:
                posInSubShapeIndex = default_posInSubShapeIndex;
                break;
            case properties_i.subShapeTipLocIndex:
                subShapeTipLocIndex = default_subShapeTipLocIndex;
                break;
            case properties_i.subShape2TipLocIndex:
                subShape2TipLocIndex = default_subShape2TipLocIndex;
                break;
            case properties_i.subShape3TipLocIndex:
                subShape3TipLocIndex = default_subShape3TipLocIndex;
                break;
            case properties_i.shapeIndex:
                shapeIndex = default_shapeIndex;
                break;
            case properties_i.subShapeRefIndex:
                subShapeRefIndex = default_subShapeRefIndex;
                break;
            case properties_i.edgeSlide:
                edgeSlide = default_edgeSlide;
                break;
            case properties_i.proxRays:
                proximitySideRays = default_proximitySideRays;
                break;
            case properties_i.proxSideRaysFallOff:
                proxSideRaysFallOff = default_proximitySideRaysFallOff;
                break;
            case properties_i.lwrType:
                LWRNoiseType = default_LWRNoiseType;
                break;
            case properties_i.lwr2Type:
                LWR2NoiseType = default_LWRNoiseType;
                break;
        }
    }

    public static int getDefaultInt(properties_i p)
    {
        return pGetDefaultInt(p);
    }

    private static int pGetDefaultInt(properties_i p)
    {
        int val = 0;
        switch (p)
        {
            case properties_i.enabled:
                val = default_enabled;
                break;
            case properties_i.gCSEngine:
                val = default_geoCoreShapeEngine;
                break;
            case properties_i.posInSubShapeIndex:
                val = default_posInSubShapeIndex;
                break;
            case properties_i.subShapeTipLocIndex:
                val = default_subShapeTipLocIndex;
                break;
            case properties_i.subShape2TipLocIndex:
                val = default_subShape2TipLocIndex;
                break;
            case properties_i.subShape3TipLocIndex:
                val = default_subShape3TipLocIndex;
                break;
            case properties_i.shapeIndex:
                val = default_shapeIndex;
                break;
            case properties_i.subShapeRefIndex:
                val = default_subShapeRefIndex;
                break;
            case properties_i.edgeSlide:
                val = default_edgeSlide;
                break;
            case properties_i.proxRays:
                val = default_proximitySideRays;
                break;
            case properties_i.proxSideRaysFallOff:
                val = default_proximitySideRaysFallOff;
                break;
            case properties_i.lwrType:
                val = default_LWRNoiseType;
                break;
            case properties_i.lwr2Type:
                val = default_LWRNoiseType;
                break;
        }

        return val;
    }

    private static readonly decimal default_subShapeHorLength = 0;
    private static readonly decimal default_subShapeHorOffset = 0;
    private static readonly decimal default_subShapeVerLength = 0;
    private static readonly decimal default_subShapeVerOffset = 0;
    private static readonly decimal default_subShape2HorLength = 0;
    private static readonly decimal default_subShape2HorOffset = 0;
    private static readonly decimal default_subShape2VerLength = 0;
    private static readonly decimal default_subShape2VerOffset = 0;
    private static readonly decimal default_subShape3HorLength = 0;
    private static readonly decimal default_subShape3HorOffset = 0;
    private static readonly decimal default_subShape3VerLength = 0;
    private static readonly decimal default_subShape3VerOffset = 0;
    private static readonly decimal default_sideBias = 0;
    private static readonly decimal default_horTipBias = 0;
    private static readonly decimal default_verTipBias = 0;
    private static readonly decimal default_innerCRR = 0;
    private static readonly decimal default_outerCRR = 0;
    private static readonly decimal default_rotation = 0;
    private static readonly decimal default_proximityBias = 0;
    private static readonly decimal default_proximityIsoDistance = 0;
    private static readonly decimal default_proximitySideRaysFallOffMultiplier = 1.0m;
    private static readonly decimal default_edgeSlideTension = 0.35m;
    private static readonly decimal default_horTipBiasNVar = 0;
    private static readonly decimal default_horTipBiasPVar = 0;
    private static readonly decimal default_verTipBiasNVar = 0;
    private static readonly decimal default_verTipBiasPVar = 0;

    private static readonly decimal default_globalHorOffset = 0;
    private static readonly decimal default_globalVerOffset = 0;

    private static readonly decimal default_rayExtension = 1.03m;

    private decimal subShapeHorLength = default_subShapeHorLength;
    private decimal subShapeHorOffset = default_subShapeHorOffset;
    private decimal subShapeVerLength = default_subShapeVerLength;
    private decimal subShapeVerOffset = default_subShapeVerOffset;
    private decimal subShape2HorLength = default_subShape2HorLength;
    private decimal subShape2HorOffset = default_subShape2HorOffset;
    private decimal subShape2VerLength = default_subShape2VerLength;
    private decimal subShape2VerOffset = default_subShape2VerOffset;
    private decimal subShape3HorLength = default_subShape3HorLength;
    private decimal subShape3HorOffset = default_subShape3HorOffset;
    private decimal subShape3VerLength = default_subShape3VerLength;
    private decimal subShape3VerOffset = default_subShape3HorOffset;
    private decimal innerCRR = default_innerCRR;
    private decimal outerCRR = default_outerCRR;
    private decimal globalHorOffset = default_globalHorOffset;
    private decimal globalVerOffset = default_globalVerOffset;
    private decimal sideBias = default_sideBias;
    private decimal horTipBias = default_horTipBias;
    private decimal horTipBiasPVar = default_horTipBias;
    private decimal horTipBiasNVar;
    private decimal verTipBias;
    private decimal verTipBiasPVar;
    private decimal verTipBiasNVar;
    private decimal proximityBias = default_proximityBias;
    private decimal proximityIsoDistance = default_proximityIsoDistance;
    private decimal rotation = default_rotation;
    private decimal edgeSlideTension= default_edgeSlideTension;
    private decimal LWR;
    private decimal LWR2;
    private decimal LWRNoiseFreq = default_LWRNoiseFreq;
    private decimal LWR2NoiseFreq = default_LWRNoiseFreq;
    private decimal proxSideRaysMultiplier = default_proximitySideRaysFallOffMultiplier;
    private decimal rayExtension = default_rayExtension;
    private decimal gcRayExtension;

    public enum properties_decimal
    {
        s0HorLength, s1HorLength, s2HorLength,
        s0VerLength, s1VerLength, s2VerLength,
        s0HorOffset, s1HorOffset, s2HorOffset,
        s0VerOffset, s1VerOffset, s2VerOffset,
        gHorOffset, gVerOffset,
        iCR, oCR,
        sBias, hTBias, hTNVar, hTPVar, vTBias, vTNVar, vTPVar,
        lwr, lwrFreq, lwr2, lwr2Freq,
        eTension, rot,
        pBias, pBiasDist, proxSideRaysMultiplier,
        rayExtension, keyhole_factor
    }

    public decimal getDecimal(properties_decimal p)
    {
        return pGetDecimal(p);
    }

    decimal pGetDecimal(properties_decimal p)
    {
        decimal ret = 0;
        switch (p)
        {
            case properties_decimal.s0HorLength:
                ret = subShapeHorLength;
                break;
            case properties_decimal.s1HorLength:
                ret = subShape2HorLength;
                break;
            case properties_decimal.s2HorLength:
                ret = subShape3HorLength;
                break;
            case properties_decimal.s0VerLength:
                ret = subShapeVerLength;
                break;
            case properties_decimal.s1VerLength:
                ret = subShape2VerLength;
                break;
            case properties_decimal.s2VerLength:
                ret = subShape3VerLength;
                break;
            case properties_decimal.s0HorOffset:
                ret = subShapeHorOffset;
                break;
            case properties_decimal.s1HorOffset:
                ret = subShape2HorOffset;
                break;
            case properties_decimal.s2HorOffset:
                ret = subShape3HorOffset;
                break;
            case properties_decimal.s0VerOffset:
                ret = subShapeVerOffset;
                break;
            case properties_decimal.s1VerOffset:
                ret = subShape2VerOffset;
                break;
            case properties_decimal.s2VerOffset:
                ret = subShape3VerOffset;
                break;
            case properties_decimal.gHorOffset:
                ret = globalHorOffset;
                break;
            case properties_decimal.gVerOffset:
                ret = globalVerOffset;
                break;
            case properties_decimal.iCR:
                ret = innerCRR;
                break;
            case properties_decimal.oCR:
                ret = outerCRR;
                break;
            case properties_decimal.sBias:
                ret = sideBias;
                break;
            case properties_decimal.hTBias:
                ret = horTipBias;
                break;
            case properties_decimal.hTNVar:
                ret = horTipBiasNVar;
                break;
            case properties_decimal.hTPVar:
                ret = horTipBiasPVar;
                break;
            case properties_decimal.vTBias:
                ret = verTipBias;
                break;
            case properties_decimal.vTNVar:
                ret = verTipBiasNVar;
                break;
            case properties_decimal.vTPVar:
                ret = verTipBiasPVar;
                break;
            case properties_decimal.lwr:
                ret = LWR;
                break;
            case properties_decimal.lwr2:
                ret = LWR2;
                break;
            case properties_decimal.lwrFreq:
                ret = LWRNoiseFreq;
                break;
            case properties_decimal.lwr2Freq:
                ret = LWR2NoiseFreq;
                break;
            case properties_decimal.eTension:
                ret = edgeSlideTension;
                break;
            case properties_decimal.rot:
                ret = rotation;
                break;
            case properties_decimal.pBias:
                ret = proximityBias;
                break;
            case properties_decimal.pBiasDist:
                ret = proximityIsoDistance;
                break;
            case properties_decimal.proxSideRaysMultiplier:
                ret = proxSideRaysMultiplier;
                break;
            case properties_decimal.keyhole_factor:
                ret = gcRayExtension;
                break;
            case properties_decimal.rayExtension:
                ret = rayExtension;
                break;
        }

        return ret;
    }

    public void setDecimal(properties_decimal p, decimal val)
    {
        pSetDecimal(p, val);
    }

    void pSetDecimal(properties_decimal p, decimal val)
    {
        switch (p)
        {
            case properties_decimal.s0HorLength:
                subShapeHorLength = val;
                break;
            case properties_decimal.s1HorLength:
                subShape2HorLength = val;
                break;
            case properties_decimal.s2HorLength:
                subShape3HorLength = val;
                break;
            case properties_decimal.s0VerLength:
                subShapeVerLength = val;
                break;
            case properties_decimal.s1VerLength:
                subShape2VerLength = val;
                break;
            case properties_decimal.s2VerLength:
                subShape3VerLength = val;
                break;
            case properties_decimal.s0HorOffset:
                subShapeHorOffset = val;
                break;
            case properties_decimal.s1HorOffset:
                subShape2HorOffset = val;
                break;
            case properties_decimal.s2HorOffset:
                subShape3HorOffset = val;
                break;
            case properties_decimal.s0VerOffset:
                subShapeVerOffset = val;
                break;
            case properties_decimal.s1VerOffset:
                subShape2VerOffset = val;
                break;
            case properties_decimal.s2VerOffset:
                subShape3VerOffset = val;
                break;
            case properties_decimal.gHorOffset:
                globalHorOffset = val;
                break;
            case properties_decimal.gVerOffset:
                globalVerOffset = val;
                break;
            case properties_decimal.iCR:
                innerCRR = val;
                break;
            case properties_decimal.oCR:
                outerCRR = val;
                break;
            case properties_decimal.sBias:
                sideBias = val;
                break;
            case properties_decimal.hTBias:
                horTipBias = val;
                break;
            case properties_decimal.hTNVar:
                horTipBiasNVar = val;
                break;
            case properties_decimal.hTPVar:
                horTipBiasPVar = val;
                break;
            case properties_decimal.vTBias:
                verTipBias = val;
                break;
            case properties_decimal.vTNVar:
                verTipBiasNVar = val;
                break;
            case properties_decimal.vTPVar:
                verTipBiasPVar = val;
                break;
            case properties_decimal.lwr:
                LWR = val;
                break;
            case properties_decimal.lwr2:
                LWR2 = val;
                break;
            case properties_decimal.lwrFreq:
                LWRNoiseFreq = val;
                break;
            case properties_decimal.lwr2Freq:
                LWR2NoiseFreq = val;
                break;
            case properties_decimal.eTension:
                edgeSlideTension = val;
                break;
            case properties_decimal.rot:
                rotation = val;
                break;
            case properties_decimal.pBias:
                proximityBias = val;
                break;
            case properties_decimal.pBiasDist:
                proximityIsoDistance = val;
                break;
            case properties_decimal.proxSideRaysMultiplier:
                proxSideRaysMultiplier = val;
                break;
            case properties_decimal.rayExtension:
                rayExtension = val;
                break;
            case properties_decimal.keyhole_factor:
                gcRayExtension = val;
                break;
        }
    }

    public void defaultDecimal(properties_decimal p)
    {
        pDefaultDecimal(p);
    }

    void pDefaultDecimal(properties_decimal p)
    {
        switch (p)
        {
            case properties_decimal.s0HorLength:
                subShapeHorLength = default_subShapeHorLength;
                break;
            case properties_decimal.s1HorLength:
                subShape2HorLength = default_subShape2HorLength;
                break;
            case properties_decimal.s2HorLength:
                subShape3HorLength = default_subShape3HorLength;
                break;
            case properties_decimal.s0VerLength:
                subShapeVerLength = default_subShapeVerLength;
                break;
            case properties_decimal.s1VerLength:
                subShape2VerLength = default_subShape2VerLength;
                break;
            case properties_decimal.s2VerLength:
                subShape3VerLength = default_subShape3VerLength;
                break;
            case properties_decimal.s0HorOffset:
                subShapeHorOffset = default_subShapeHorOffset;
                break;
            case properties_decimal.s1HorOffset:
                subShape2HorOffset = default_subShape2HorOffset;
                break;
            case properties_decimal.s2HorOffset:
                subShape3HorOffset = default_subShape3HorOffset;
                break;
            case properties_decimal.s0VerOffset:
                subShapeVerOffset = default_subShapeVerOffset;
                break;
            case properties_decimal.s1VerOffset:
                subShape2VerOffset = default_subShape2VerOffset;
                break;
            case properties_decimal.s2VerOffset:
                subShape3VerOffset = default_subShape3VerOffset;
                break;
            case properties_decimal.gHorOffset:
                globalHorOffset = default_globalHorOffset;
                break;
            case properties_decimal.gVerOffset:
                globalVerOffset = default_globalVerOffset;
                break;
            case properties_decimal.iCR:
                innerCRR = default_innerCRR;
                break;
            case properties_decimal.oCR:
                outerCRR = default_outerCRR;
                break;
            case properties_decimal.sBias:
                sideBias = default_sideBias;
                break;
            case properties_decimal.hTBias:
                horTipBias = default_horTipBias;
                break;
            case properties_decimal.hTNVar:
                horTipBiasNVar = default_horTipBiasNVar;
                break;
            case properties_decimal.hTPVar:
                horTipBiasPVar = default_horTipBiasPVar;
                break;
            case properties_decimal.vTBias:
                verTipBias = default_verTipBias;
                break;
            case properties_decimal.vTNVar:
                verTipBiasNVar = default_verTipBiasNVar;
                break;
            case properties_decimal.vTPVar:
                verTipBiasPVar = default_verTipBiasPVar;
                break;
            case properties_decimal.lwr:
                LWR = default_LWR;
                break;
            case properties_decimal.lwr2:
                LWR2 = default_LWR;
                break;
            case properties_decimal.lwrFreq:
                LWRNoiseFreq = default_LWRNoiseFreq;
                break;
            case properties_decimal.lwr2Freq:
                LWR2NoiseFreq = default_LWRNoiseFreq;
                break;
            case properties_decimal.eTension:
                edgeSlideTension = default_edgeSlideTension;
                break;
            case properties_decimal.rot:
                rotation = default_rotation;
                break;
            case properties_decimal.pBias:
                proximityBias = default_proximityBias;
                break;
            case properties_decimal.pBiasDist:
                proximityIsoDistance = default_proximityIsoDistance;
                break;
            case properties_decimal.proxSideRaysMultiplier:
                proxSideRaysMultiplier = default_proximitySideRaysFallOffMultiplier;
                break;
            case properties_decimal.rayExtension:
                rayExtension = default_rayExtension;
                break;
            case properties_decimal.keyhole_factor:
                gcRayExtension = default_rayExtension;
                break;
        }
    }

    public static decimal getDefaultDecimal(properties_decimal p)
    {
        return pGetDefaultDecimal(p);
    }

    static decimal pGetDefaultDecimal(properties_decimal p)
    {
        decimal ret = 0;
        switch (p)
        {
            case properties_decimal.s0HorLength:
                ret = default_subShapeHorLength;
                break;
            case properties_decimal.s1HorLength:
                ret = default_subShape2HorLength;
                break;
            case properties_decimal.s2HorLength:
                ret = default_subShape3HorLength;
                break;
            case properties_decimal.s0VerLength:
                ret = default_subShapeVerLength;
                break;
            case properties_decimal.s1VerLength:
                ret = default_subShape2VerLength;
                break;
            case properties_decimal.s2VerLength:
                ret = default_subShape3VerLength;
                break;
            case properties_decimal.s0HorOffset:
                ret = default_subShapeHorOffset;
                break;
            case properties_decimal.s1HorOffset:
                ret = default_subShape2HorOffset;
                break;
            case properties_decimal.s2HorOffset:
                ret = default_subShape3HorOffset;
                break;
            case properties_decimal.s0VerOffset:
                ret = default_subShapeVerOffset;
                break;
            case properties_decimal.s1VerOffset:
                ret = default_subShape2VerOffset;
                break;
            case properties_decimal.s2VerOffset:
                ret = default_subShape3VerOffset;
                break;
            case properties_decimal.gHorOffset:
                ret = default_globalHorOffset;
                break;
            case properties_decimal.gVerOffset:
                ret = default_globalVerOffset;
                break;
            case properties_decimal.iCR:
                ret = default_innerCRR;
                break;
            case properties_decimal.oCR:
                ret = default_outerCRR;
                break;
            case properties_decimal.sBias:
                ret = default_sideBias;
                break;
            case properties_decimal.hTBias:
                ret = default_horTipBias;
                break;
            case properties_decimal.hTNVar:
                ret = default_horTipBiasNVar;
                break;
            case properties_decimal.hTPVar:
                ret = default_horTipBiasPVar;
                break;
            case properties_decimal.vTBias:
                ret = default_verTipBias;
                break;
            case properties_decimal.vTNVar:
                ret = default_verTipBiasNVar;
                break;
            case properties_decimal.vTPVar:
                ret = default_verTipBiasPVar;
                break;
            case properties_decimal.lwr:
                ret = default_LWR;
                break;
            case properties_decimal.lwr2:
                ret = default_LWR;
                break;
            case properties_decimal.lwrFreq:
                ret = default_LWRNoiseFreq;
                break;
            case properties_decimal.lwr2Freq:
                ret = default_LWRNoiseFreq;
                break;
            case properties_decimal.eTension:
                ret = default_edgeSlideTension;
                break;
            case properties_decimal.rot:
                ret = default_rotation;
                break;
            case properties_decimal.pBias:
                ret = default_proximityBias;
                break;
            case properties_decimal.pBiasDist:
                ret = default_proximityIsoDistance;
                break;
            case properties_decimal.proxSideRaysMultiplier:
                ret = default_proximitySideRaysFallOffMultiplier;
                break;
            case properties_decimal.rayExtension:
                ret = default_rayExtension;
                break;
            case properties_decimal.keyhole_factor:
                ret = default_rayExtension;
                break;
        }

        return ret;
    }

    private string layerName;

    public enum properties_s
    {
        s_name
    }

    public void setString(properties_s p, string val)
    {
        pSetString(p, val);
    }

    void pSetString(properties_s p, string val)
    {
        switch (p)
        {
            case properties_s.s_name:
                layerName = val;
                break;
        }
        
    }

    public string getString(properties_s p)
    {
        return pGetString(p);
    }

    string pGetString(properties_s p)
    {
        string ret = "";
        switch (p)
        {
            case properties_s.s_name:
                ret = layerName;
                break;
        }

        return ret;
    }

    public void defaultString(properties_s p)
    {
        pDefaultString(p);
    }

    void pDefaultString(properties_s p)
    {
        switch (p)
        {
            case properties_s.s_name:
                layerName = default_layerName;
                break;
        }
    }
    
    public static string getDefaultString(properties_s p)
    {
        return pGetDefaultString(p);
    }

    static string pGetDefaultString(properties_s p)
    {
        string ret = "";
        switch (p)
        {
            case properties_s.s_name:
                ret = default_layerName;
                break;
        }

        return ret;
    }
    
}

public class ShapeLibrary
{
    private static readonly List<string> availableShapes= new () { "(None)", "Rectangle/Square", "L-shape", "T-shape", "X-shape", "U-shape", "S-shape", "GDS/Oasis", "Boolean" };

    public static List<string> getAvailableShapes()
    {
        return pGetAvailableShapes();
    }

    private static List<string> pGetAvailableShapes()
    {
        return availableShapes;
    }

    public enum shapeNames { none, rect, Lshape, Tshape, Xshape, Ushape, Sshape, GEOCORE, BOOLEAN }

    private int shapeIndex;
    public bool shapeValid { get; private set; }
    public bool geoCoreShapeOrthogonal { get; private set; }
    public MyVertex[] Vertex { get; private set; }
    public MyRound[] round1 { get; private set; }
    public bool[] tips { get; private set; }
    private ShapeSettings layerSettings;

    public ShapeLibrary(ShapeSettings mcLayerSettings)
    {
        pShapeLibrary(mcLayerSettings);
    }

    private void pShapeLibrary(ShapeSettings mcLayerSettings)
    {
        shapeValid = false;
        layerSettings = mcLayerSettings;
    }

    public ShapeLibrary(int shapeIndex_, ShapeSettings mcLayerSettings)
    {
        pShapeLibrary(shapeIndex_, mcLayerSettings);
    }

    private void pShapeLibrary(int shapeIndex_, ShapeSettings mcLayerSettings)
    {
        shapeIndex = shapeIndex_;
        shapeValid = false;
        layerSettings = mcLayerSettings;
        pSetShape(shapeIndex);
    }

    public void setShape(int shapeIndex_, GeoLibPointF[] sourcePoly = null)
    {
        pSetShape(shapeIndex_, sourcePoly);
    }

    private void pSetShape(int shapeIndex_, GeoLibPointF[] sourcePoly = null)
    {
        try
        {
            shapeIndex = shapeIndex_;
            switch (shapeIndex)
            {
                case (int)shapeNames.rect:
                    rectangle();
                    break;
                case (int)shapeNames.Lshape:
                    Lshape();
                    break;
                case (int)shapeNames.Tshape:
                    Tshape();
                    break;
                case (int)shapeNames.Xshape:
                    crossShape();
                    break;
                case (int)shapeNames.Ushape:
                    Ushape();
                    break;
                case (int)shapeNames.Sshape:
                    Sshape();
                    break;
                case (int)shapeNames.GEOCORE:
                    if (layerSettings.getInt(ShapeSettings.properties_i.gCSEngine) == 1)
                    {
                        customShape(sourcePoly);
                    }
                    break;
            }
        }
        catch (Exception)
        {

        }
    }

    private void configureArrays()
    {
        double vertexCount = 0;
        switch (shapeIndex)
        {
            case (int)shapeNames.rect: // rectangle
                vertexCount = 9;
                break;
            case (int)shapeNames.Lshape: // L
                vertexCount = 13;
                break;
            case (int)shapeNames.Tshape: // T
                vertexCount = 17;
                break;
            case (int)shapeNames.Xshape: // Cross
                vertexCount = 25;
                break;
            case (int)shapeNames.Ushape: // U
                vertexCount = 17;
                break;
            case (int)shapeNames.Sshape: // S
                vertexCount = 25;
                break;
        }
        int arrayLength = (int)vertexCount;
        Vertex = new MyVertex[arrayLength];
        tips = new bool[arrayLength];
#if !SHAPELIBSINGLETHREADED
        Parallel.For(0, arrayLength, i => 
#else
            for (Int32 i = 0; i < arrayLength; i++)
#endif
            {
                tips[i] = false;
            }
#if !SHAPELIBSINGLETHREADED
        );
#endif
        arrayLength = (int)Math.Floor(vertexCount / 2) + 1;
        round1 = new MyRound[arrayLength];
#if !SHAPELIBSINGLETHREADED
        Parallel.For(0, arrayLength, i => 
#else
            for (Int32 i = 0; i < arrayLength; i++)
#endif
            {
                round1[i] = new MyRound();
            }
#if !SHAPELIBSINGLETHREADED
        );
#endif
    }

    private void rectangle()
    {
        configureArrays();

        // Sort out the tips by setting 'true' to center vertex that defines tip in each case.
        switch (layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex))
        {
            case (int)ShapeSettings.tipLocations.none: // None
                break;
            case (int)ShapeSettings.tipLocations.L: // Left
                tips[1] = true;
                tips[8] = true;
                break;
            case (int)ShapeSettings.tipLocations.R: // Right
                tips[5] = true;
                break;
            case (int)ShapeSettings.tipLocations.LR: // Left and right
                tips[1] = true;
                tips[8] = true;
                tips[5] = true;
                break;
            case (int)ShapeSettings.tipLocations.T: // Top
                tips[3] = true;
                break;
            case (int)ShapeSettings.tipLocations.B: // Bottom
                tips[7] = true;
                break;
            case (int)ShapeSettings.tipLocations.TB: // Top and Bottom
                tips[3] = true;
                tips[7] = true;
                break;
            case (int)ShapeSettings.tipLocations.TL: // Top and Left
                tips[1] = true;
                tips[3] = true;
                tips[8] = true;
                break;
            case (int)ShapeSettings.tipLocations.TR: // Top and Right
                tips[3] = true;
                tips[5] = true;
                break;
            case (int)ShapeSettings.tipLocations.TLR: // Top and Left and Right
                tips[1] = true;
                tips[3] = true;
                tips[5] = true;
                tips[8] = true;
                break;
            case (int)ShapeSettings.tipLocations.BL: // Bottom and Left
                tips[1] = true;
                tips[7] = true;
                tips[8] = true;
                break;
            case (int)ShapeSettings.tipLocations.BR: // Bottom and Right
                tips[5] = true;
                tips[7] = true;
                break;
            case (int)ShapeSettings.tipLocations.BLR: // Bottom and Left and Right
                tips[1] = true;
                tips[5] = true;
                tips[7] = true;
                tips[8] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBL: // Top and Bottom and Left
                tips[3] = true;
                tips[7] = true;
                tips[1] = true;
                tips[8] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBR: // Top and Bottom and Right
                tips[3] = true;
                tips[7] = true;
                tips[5] = true;
                break;
            default: // All
                tips[3] = true;
                tips[7] = true;
                tips[1] = true;
                tips[8] = true;
                tips[5] = true;
                break;
        }

        // NOTE:
        // Subshape offsets are applied later to simplify ellipse generation
        // Horizontal and vertical global offsets are applied in callsite.
        // Repositioning with respect to subshape reference is also applied in callsite.

        double tmpX = 0.0;
        double tmpY = 0.0;
        Vertex[0] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0VerLength)) / 2;
        Vertex[1] = new MyVertex(tmpX, tmpY, typeDirection.left1, true, false, typeVertex.center);

        tmpY = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0VerLength));
        Vertex[2] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0HorLength)) / 2;
        Vertex[3] = new MyVertex(tmpX, tmpY, typeDirection.up1, false, false, typeVertex.center);

        tmpX = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0HorLength));
        Vertex[4] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY -= Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0VerLength)) / 2;
        Vertex[5] = new MyVertex(tmpX, tmpY, typeDirection.right1, true, false, typeVertex.center);

        tmpY = 0.0;
        Vertex[6] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0HorLength)) / 2;
        tmpY = 0.0;
        Vertex[7] = new MyVertex(tmpX, tmpY, typeDirection.down1, false, false, typeVertex.center);

        processEdgesForRounding();

        shapeValid = true;
    }

    private void Tshape()
    {
        configureArrays();
        // Sort out the tips by setting 'true' to center vertex that defines tip in each case.
        switch (layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex))
        {
            case (int)ShapeSettings.tipLocations.none: // None
                break;
            case (int)ShapeSettings.tipLocations.L: // Left
                tips[1] = true;
                break;
            case (int)ShapeSettings.tipLocations.R: // Right
                tips[5] = true;
                tips[13] = true;
                break;
            case (int)ShapeSettings.tipLocations.LR: // Left and right
                tips[1] = true;
                tips[5] = true;
                tips[13] = true;
                break;
            case (int)ShapeSettings.tipLocations.T: // Top
                tips[3] = true;
                break;
            case (int)ShapeSettings.tipLocations.B: // Bottom
                tips[15] = true;
                break;
            case (int)ShapeSettings.tipLocations.TB: // Top and Bottom
                tips[3] = true;
                tips[15] = true;
                break;
            case (int)ShapeSettings.tipLocations.TL: // Top and Left
                tips[3] = true;
                tips[1] = true;
                break;
            case (int)ShapeSettings.tipLocations.TR: // Top and Right
                tips[3] = true;
                tips[5] = true;
                tips[13] = true;
                break;
            case (int)ShapeSettings.tipLocations.TLR: // Top and Left and Right
                tips[1] = true;
                tips[3] = true;
                tips[5] = true;
                tips[13] = true;
                break;
            case (int)ShapeSettings.tipLocations.BL: // Bottom and Left
                tips[15] = true;
                tips[1] = true;
                break;
            case (int)ShapeSettings.tipLocations.BR: // Bottom and Right
                tips[15] = true;
                tips[5] = true;
                tips[13] = true;
                break;
            case (int)ShapeSettings.tipLocations.BLR: // Bottom and Left and Right
                tips[1] = true;
                tips[15] = true;
                tips[5] = true;
                tips[13] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBL: // Top and Bottom and Left
                tips[1] = true;
                tips[3] = true;
                tips[15] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBR: // Top and Bottom and Right
                tips[3] = true;
                tips[5] = true;
                tips[13] = true;
                tips[15] = true;
                break;
            default: // All
                tips[1] = true;
                tips[3] = true;
                tips[5] = true;
                tips[13] = true;
                tips[15] = true;
                break;
        }
        switch (layerSettings.getInt(ShapeSettings.properties_i.subShape2TipLocIndex))
        {
            case (int)ShapeSettings.tipLocations.none: // None
                break;
            case (int)ShapeSettings.tipLocations.L: // Left
                break;
            case (int)ShapeSettings.tipLocations.R: // Right
                tips[9] = true;
                break;
            case (int)ShapeSettings.tipLocations.LR: // Left and right
                tips[9] = true;
                break;
            case (int)ShapeSettings.tipLocations.T: // Top
                tips[7] = true;
                break;
            case (int)ShapeSettings.tipLocations.B: // Bottom
                tips[11] = true;
                break;
            case (int)ShapeSettings.tipLocations.TB: // Top and Bottom
                tips[7] = true;
                tips[11] = true;
                break;
            case (int)ShapeSettings.tipLocations.TL: // Top and Left
                tips[7] = true;
                break;
            case (int)ShapeSettings.tipLocations.TR: // Top and Right
                tips[7] = true;
                tips[9] = true;
                break;
            case (int)ShapeSettings.tipLocations.TLR: // Top and Left and Right
                tips[7] = true;
                tips[9] = true;
                break;
            case (int)ShapeSettings.tipLocations.BL: // Bottom and Left
                tips[11] = true;
                break;
            case (int)ShapeSettings.tipLocations.BR: // Bottom and Right
                tips[11] = true;
                tips[9] = true;
                break;
            case (int)ShapeSettings.tipLocations.BLR: // Bottom and Left and Right
                tips[11] = true;
                tips[9] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBL: // Top and Bottom and Left
                tips[7] = true;
                tips[11] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBR: // Top and Bottom and Right
                tips[7] = true;
                tips[9] = true;
                tips[11] = true;
                break;
            default: // All
                tips[7] = true;
                tips[9] = true;
                tips[11] = true;
                break;
        }

        // NOTE:
        // Subshape offsets are applied later to simplify ellipse generation
        // Horizontal and vertical global offsets are applied in callsite.

        double tmpX = 0.0;
        double tmpY = 0.0;
        Vertex[0] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0VerLength)) / 2;
        Vertex[1] = new MyVertex(tmpX, tmpY, typeDirection.left1, true, false, typeVertex.center);

        tmpY = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0VerLength));
        Vertex[2] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0HorLength)) / 2;
        Vertex[3] = new MyVertex(tmpX, tmpY, typeDirection.up1, false, false, typeVertex.center);

        tmpX += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0HorLength)) / 2;
        Vertex[4] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY = (Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0VerLength)) - (Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1VerLength)) +
            (Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1VerOffset)) - Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0VerOffset))))) / 2;
        Vertex[5] = new MyVertex(tmpX, tmpY, typeDirection.right1, true, false, typeVertex.center);

        tmpY = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1VerLength)) + Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1VerOffset));
        Vertex[6] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1HorLength) / 2);
        Vertex[7] = new MyVertex(tmpX, tmpY, typeDirection.up1, false, false, typeVertex.center);

        tmpX += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1HorLength) / 2);
        Vertex[8] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY -= Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1VerLength) / 2);
        Vertex[9] = new MyVertex(tmpX, tmpY, typeDirection.right1, true, false, typeVertex.center);

        tmpY -= Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1VerLength) / 2);
        Vertex[10] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX -= Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1HorLength) / 2);
        Vertex[11] = new MyVertex(tmpX, tmpY, typeDirection.down1, false, false, typeVertex.center);

        tmpX -= Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1HorLength) / 2);
        Vertex[12] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1VerOffset)) / 2;
        Vertex[13] = new MyVertex(tmpX, tmpY, typeDirection.right1, true, false, typeVertex.center);

        tmpY = 0;
        Vertex[14] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX -= Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0HorLength) / 2);
        Vertex[15] = new MyVertex(tmpX, tmpY, typeDirection.down1, false, false, typeVertex.center);

        processEdgesForRounding();

        shapeValid = true;
    }

    private void Lshape()
    {
        configureArrays();
        // Sort out the tips by setting 'true' to center vertex that defines tip in each case.
        switch (layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex))
        {
            case (int)ShapeSettings.tipLocations.none: // None
                break;
            case (int)ShapeSettings.tipLocations.L: // Left
                tips[1] = true;
                break;
            case (int)ShapeSettings.tipLocations.R: // Right
                tips[5] = true;
                break;
            case (int)ShapeSettings.tipLocations.LR: // Left and right
                tips[1] = true;
                tips[5] = true;
                break;
            case (int)ShapeSettings.tipLocations.T: // Top
                tips[3] = true;
                break;
            case (int)ShapeSettings.tipLocations.B: // Bottom
                tips[11] = true;
                break;
            case (int)ShapeSettings.tipLocations.TB: // Top and Bottom
                tips[11] = true;
                tips[3] = true;
                break;
            case (int)ShapeSettings.tipLocations.TL: // Top and Left
                tips[1] = true;
                tips[3] = true;
                break;
            case (int)ShapeSettings.tipLocations.TR: // Top and Right
                tips[3] = true;
                tips[5] = true;
                break;
            case (int)ShapeSettings.tipLocations.TLR: // Top and Left and Right
                tips[1] = true;
                tips[3] = true;
                tips[5] = true;
                break;
            case (int)ShapeSettings.tipLocations.BL: // Bottom and Left
                tips[11] = true;
                tips[1] = true;
                break;
            case (int)ShapeSettings.tipLocations.BR: // Bottom and Right
                tips[11] = true;
                tips[5] = true;
                break;
            case (int)ShapeSettings.tipLocations.BLR: // Bottom and Left and Right
                tips[1] = true;
                tips[5] = true;
                tips[11] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBL: // Top and Bottom and Left
                tips[1] = true;
                tips[11] = true;
                tips[3] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBR: // Top and Bottom and Right
                tips[11] = true;
                tips[3] = true;
                tips[5] = true;
                break;
            default: // All
                tips[1] = true;
                tips[5] = true;
                tips[11] = true;
                tips[3] = true;
                break;
        }
        switch (layerSettings.getInt(ShapeSettings.properties_i.subShape2TipLocIndex))
        {
            case (int)ShapeSettings.tipLocations.none: // None
                break;
            case (int)ShapeSettings.tipLocations.L: // Left
                break;
            case (int)ShapeSettings.tipLocations.R: // Right
                tips[9] = true;
                break;
            case (int)ShapeSettings.tipLocations.LR: // Left and right
                tips[9] = true;
                break;
            case (int)ShapeSettings.tipLocations.T: // Top
                tips[7] = true;
                break;
            case (int)ShapeSettings.tipLocations.B: // Bottom
                tips[11] = true;
                break;
            case (int)ShapeSettings.tipLocations.TB: // Top and Bottom
                tips[7] = true;
                tips[11] = true;
                break;
            case (int)ShapeSettings.tipLocations.TL: // Top and Left
                tips[7] = true;
                break;
            case (int)ShapeSettings.tipLocations.TR: // Top and Right
                tips[7] = true;
                tips[9] = true;
                break;
            case (int)ShapeSettings.tipLocations.TLR: // Top and Left and Right
                tips[7] = true;
                tips[9] = true;
                break;
            case (int)ShapeSettings.tipLocations.BL: // Bottom and Left
                tips[11] = true;
                break;
            case (int)ShapeSettings.tipLocations.BR: // Bottom and Right
                tips[11] = true;
                tips[9] = true;
                break;
            case (int)ShapeSettings.tipLocations.BLR: // Bottom and Left and Right
                tips[11] = true;
                tips[9] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBL: // Top and Bottom and Left
                tips[7] = true;
                tips[11] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBR: // Top and Bottom and Right
                tips[7] = true;
                tips[11] = true;
                tips[9] = true;
                break;
            default: // All
                tips[9] = true;
                tips[7] = true;
                tips[11] = true;
                break;
        }

        // NOTE:
        // Subshape offsets are applied later to simplify ellipse generation
        // Horizontal and vertical global offsets are applied in callsite.
        // Repositioning with respect to subshape reference is also applied in callsite.
        double tmpX = 0.0;
        double tmpY = 0.0;
        Vertex[0] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0VerLength)) / 2;
        Vertex[1] = new MyVertex(tmpX, tmpY, typeDirection.left1, true, false, typeVertex.center);

        tmpY += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0VerLength)) / 2;
        Vertex[2] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0HorLength)) / 2;
        Vertex[3] = new MyVertex(tmpX, tmpY, typeDirection.up1, false, false, typeVertex.center);

        tmpX += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0HorLength)) / 2;
        Vertex[4] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY -= (Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0VerLength)) - Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1VerLength))) / 2;
        Vertex[5] = new MyVertex(tmpX, tmpY, typeDirection.right1, true, false, typeVertex.center);

        tmpY -= (Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0VerLength)) - Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1VerLength))) / 2;
        Vertex[6] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1HorLength)) / 2;
        Vertex[7] = new MyVertex(tmpX, tmpY, typeDirection.up1, false, false, typeVertex.center);

        tmpX += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1HorLength)) / 2;
        Vertex[8] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY -= Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1VerLength)) / 2;
        Vertex[9] = new MyVertex(tmpX, tmpY, typeDirection.right1, true, false, typeVertex.center);

        tmpY -= Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1VerLength)) / 2;
        Vertex[10] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX -= (Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0HorLength)) + Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1HorLength))) / 2;
        Vertex[11] = new MyVertex(tmpX, tmpY, typeDirection.down1, false, false, typeVertex.center);

        processEdgesForRounding();

        shapeValid = true;
    }

    private void Ushape()
    {
        configureArrays();
        // Sort out the tips by setting 'true' to center vertex that defines tip in each case.
        switch (layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex))
        {
            case (int)ShapeSettings.tipLocations.none: // None
                break;
            case (int)ShapeSettings.tipLocations.L: // Left
                tips[1] = true;
                break;
            case (int)ShapeSettings.tipLocations.R: // Right
                tips[13] = true;
                break;
            case (int)ShapeSettings.tipLocations.LR: // Left and right
                tips[1] = true;
                tips[13] = true;
                break;
            case (int)ShapeSettings.tipLocations.T: // Top
                tips[3] = true;
                tips[11] = true;
                break;
            case (int)ShapeSettings.tipLocations.B: // Bottom
                tips[15] = true;
                break;
            case (int)ShapeSettings.tipLocations.TB: // Top and Bottom
                tips[3] = true;
                tips[11] = true;
                tips[15] = true;
                break;
            case (int)ShapeSettings.tipLocations.TL: // Top and Left
                tips[1] = true;
                tips[3] = true;
                tips[11] = true;
                break;
            case (int)ShapeSettings.tipLocations.TR: // Top and Right
                tips[3] = true;
                tips[11] = true;
                tips[13] = true;
                break;
            case (int)ShapeSettings.tipLocations.TLR: // Top and Left and Right
                tips[1] = true;
                tips[3] = true;
                tips[11] = true;
                tips[13] = true;
                break;
            case (int)ShapeSettings.tipLocations.BL: // Bottom and Left
                tips[1] = true;
                tips[15] = true;
                break;
            case (int)ShapeSettings.tipLocations.BR: // Bottom and Right
                tips[13] = true;
                tips[15] = true;
                break;
            case (int)ShapeSettings.tipLocations.BLR: // Bottom and Left and Right
                tips[1] = true;
                tips[13] = true;
                tips[15] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBL: // Top and Bottom and Left
                tips[1] = true;
                tips[3] = true;
                tips[11] = true;
                tips[15] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBR: // Top and Bottom and Right
                tips[3] = true;
                tips[11] = true;
                tips[13] = true;
                tips[15] = true;
                break;
            default: // All
                tips[1] = true;
                tips[3] = true;
                tips[11] = true;
                tips[13] = true;
                tips[15] = true;
                break;
        }
        switch (layerSettings.getInt(ShapeSettings.properties_i.subShape2TipLocIndex))
        {
            case (int)ShapeSettings.tipLocations.none: // None
                break;
            case (int)ShapeSettings.tipLocations.L: // Left
                tips[5] = true;
                break;
            case (int)ShapeSettings.tipLocations.R: // Right
                tips[9] = true;
                break;
            case (int)ShapeSettings.tipLocations.LR: // Left and right
                tips[5] = true;
                tips[9] = true;
                break;
            case (int)ShapeSettings.tipLocations.T: // Top
                break;
            case (int)ShapeSettings.tipLocations.B: // Bottom
                tips[7] = true;
                break;
            case (int)ShapeSettings.tipLocations.TB: // Top and Bottom
                tips[7] = true;
                break;
            case (int)ShapeSettings.tipLocations.TL: // Top and Left
                tips[5] = true;
                break;
            case (int)ShapeSettings.tipLocations.TR: // Top and Right
                tips[9] = true;
                break;
            case (int)ShapeSettings.tipLocations.TLR: // Top and Left and Right
                tips[5] = true;
                tips[9] = true;
                break;
            case (int)ShapeSettings.tipLocations.BL: // Bottom and Left
                tips[5] = true;
                tips[7] = true;
                break;
            case (int)ShapeSettings.tipLocations.BR: // Bottom and Right
                tips[7] = true;
                tips[9] = true;
                break;
            case (int)ShapeSettings.tipLocations.BLR: // Bottom and Left and Right
                tips[5] = true;
                tips[7] = true;
                tips[9] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBL: // Top and Bottom and Left
                tips[7] = true;
                tips[5] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBR: // Top and Bottom and Right
                tips[7] = true;
                tips[9] = true;
                break;
            default: // All
                tips[5] = true;
                tips[7] = true;
                tips[9] = true;
                break;
        }

        // NOTE:
        // Subshape offsets are applied later to simplify ellipse generation
        // Horizontal and vertical global offsets are applied in callsite.
        // Repositioning with respect to subshape reference is also applied in callsite.
        double tmpX = 0.0;
        double tmpY = 0.0;
        Vertex[0] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0VerLength)) / 2;
        Vertex[1] = new MyVertex(tmpX, tmpY, typeDirection.left1, true, false, typeVertex.center);

        tmpY = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0VerLength));
        Vertex[2] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1HorOffset)) / 2;
        Vertex[3] = new MyVertex(tmpX, tmpY, typeDirection.up1, false, false, typeVertex.center);

        tmpX += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1HorOffset)) / 2;
        Vertex[4] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY -= Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1VerLength)) / 2;
        Vertex[5] = new MyVertex(tmpX, tmpY, typeDirection.right1, true, false, typeVertex.center);

        tmpY -= Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1VerLength)) / 2;
        Vertex[6] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1HorLength)) / 2;
        Vertex[7] = new MyVertex(tmpX, tmpY, typeDirection.up1, false, false, typeVertex.center);

        tmpX += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1HorLength)) / 2;
        Vertex[8] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1VerLength)) / 2;
        Vertex[9] = new MyVertex(tmpX, tmpY, typeDirection.left1, true, false, typeVertex.center);

        tmpY = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0VerLength));
        Vertex[10] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX += (Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0HorLength)) - (Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1HorOffset)) + Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1HorLength)))) / 2;
        Vertex[11] = new MyVertex(tmpX, tmpY, typeDirection.up1, false, false, typeVertex.center);

        tmpX = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0HorLength));
        Vertex[12] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY /= 2;
        Vertex[13] = new MyVertex(tmpX, tmpY, typeDirection.right1, true, false, typeVertex.center);

        tmpY = 0;
        Vertex[14] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX /= 2;
        Vertex[15] = new MyVertex(tmpX, tmpY, typeDirection.down1, false, false, typeVertex.center);

        processEdgesForRounding();

        shapeValid = true;
    }

    private void crossShape()
    {
        configureArrays();
        // Sort out the tips.
        switch (layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex))
        {
            case (int)ShapeSettings.tipLocations.none: // None
                break;
            case (int)ShapeSettings.tipLocations.L: // Left
                tips[1] = true;
                tips[9] = true;
                break;
            case (int)ShapeSettings.tipLocations.R: // Right
                tips[13] = true;
                tips[21] = true;
                break;
            case (int)ShapeSettings.tipLocations.LR: // Left and right
                tips[1] = true;
                tips[9] = true;
                tips[13] = true;
                tips[21] = true;
                break;
            case (int)ShapeSettings.tipLocations.T: // Top
                tips[11] = true;
                break;
            case (int)ShapeSettings.tipLocations.B: // Bottom
                tips[23] = true;
                break;
            case (int)ShapeSettings.tipLocations.TB: // Top and Bottom
                tips[11] = true;
                tips[23] = true;
                break;
            case (int)ShapeSettings.tipLocations.TL: // Top and Left
                tips[1] = true;
                tips[9] = true;
                tips[11] = true;
                break;
            case (int)ShapeSettings.tipLocations.TR: // Top and Right
                tips[11] = true;
                tips[13] = true;
                tips[21] = true;
                break;
            case (int)ShapeSettings.tipLocations.TLR: // Top and Left and Right
                tips[1] = true;
                tips[9] = true;
                tips[11] = true;
                tips[13] = true;
                tips[21] = true;
                break;
            case (int)ShapeSettings.tipLocations.BL: // Bottom and Left
                tips[1] = true;
                tips[9] = true;
                tips[23] = true;
                break;
            case (int)ShapeSettings.tipLocations.BR: // Bottom and Right
                tips[13] = true;
                tips[21] = true;
                tips[23] = true;
                break;
            case (int)ShapeSettings.tipLocations.BLR: // Bottom and Left and Right
                tips[1] = true;
                tips[9] = true;
                tips[13] = true;
                tips[21] = true;
                tips[23] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBL: // Top and Bottom and Left
                tips[11] = true;
                tips[1] = true;
                tips[9] = true;
                tips[23] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBR: // Top and Bottom and Right
                tips[13] = true;
                tips[11] = true;
                tips[23] = true;
                tips[21] = true;
                break;
            default: // All
                tips[1] = true;
                tips[9] = true;
                tips[13] = true;
                tips[21] = true;
                tips[11] = true;
                tips[23] = true;
                break;
        }
        switch (layerSettings.getInt(ShapeSettings.properties_i.subShape2TipLocIndex))
        {
            case (int)ShapeSettings.tipLocations.none: // None
                break;
            case (int)ShapeSettings.tipLocations.L: // Left
                tips[5] = true;
                break;
            case (int)ShapeSettings.tipLocations.R: // Right
                tips[17] = true;
                break;
            case (int)ShapeSettings.tipLocations.LR: // Left and right
                tips[5] = true;
                tips[17] = true;
                break;
            case (int)ShapeSettings.tipLocations.T: // Top
                tips[7] = true;
                tips[15] = true;
                break;
            case (int)ShapeSettings.tipLocations.B: // Bottom
                tips[3] = true;
                tips[19] = true;
                break;
            case (int)ShapeSettings.tipLocations.TB: // Top and Bottom
                tips[3] = true;
                tips[7] = true;
                tips[15] = true;
                tips[19] = true;
                break;
            case (int)ShapeSettings.tipLocations.TL: // Top and Left
                tips[5] = true;
                tips[7] = true;
                tips[15] = true;
                break;
            case (int)ShapeSettings.tipLocations.TR: // Top and Right
                tips[7] = true;
                tips[15] = true;
                tips[17] = true;
                break;
            case (int)ShapeSettings.tipLocations.TLR: // Top and Left and Right
                tips[5] = true;
                tips[7] = true;
                tips[15] = true;
                tips[17] = true;
                break;
            case (int)ShapeSettings.tipLocations.BL: // Bottom and Left
                tips[3] = true;
                tips[5] = true;
                tips[19] = true;
                break;
            case (int)ShapeSettings.tipLocations.BR: // Bottom and Right
                tips[3] = true;
                tips[17] = true;
                tips[19] = true;
                break;
            case (int)ShapeSettings.tipLocations.BLR: // Bottom and Left and Right
                tips[3] = true;
                tips[5] = true;
                tips[17] = true;
                tips[19] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBL: // Top and Bottom and Left
                tips[5] = true;
                tips[3] = true;
                tips[7] = true;
                tips[15] = true;
                tips[19] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBR: // Top and Bottom and Right
                tips[3] = true;
                tips[7] = true;
                tips[15] = true;
                tips[17] = true;
                tips[19] = true;
                break;
            default: // All
                tips[5] = true;
                tips[17] = true;
                tips[3] = true;
                tips[7] = true;
                tips[15] = true;
                tips[19] = true;
                break;
        }

        // NOTE:
        // Subshape offsets are applied later to simplify ellipse generation
        // Horizontal and vertical global offsets are applied in callsite.

        double tmpX = 0.0;
        double tmpY = 0.0;
        Vertex[0] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1VerOffset)) / 2;
        Vertex[1] = new MyVertex(tmpX, tmpY, typeDirection.left1, true, false, typeVertex.center);

        tmpY = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1VerOffset));
        Vertex[2] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1HorOffset));
        tmpX /= 2;
        Vertex[3] = new MyVertex(tmpX, tmpY, typeDirection.down1, false, false, typeVertex.center);

        tmpX *= 2;
        Vertex[4] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1VerLength)) / 2;
        Vertex[5] = new MyVertex(tmpX, tmpY, typeDirection.left1, true, false, typeVertex.center);

        tmpY += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1VerLength)) / 2;
        Vertex[6] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX /= 2;
        Vertex[7] = new MyVertex(tmpX, tmpY, typeDirection.up1, false, false, typeVertex.center);

        tmpX = 0.0;
        Vertex[8] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0VerLength)) - (Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1VerLength)) + Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1VerOffset))) / 2;
        Vertex[9] = new MyVertex(tmpX, tmpY, typeDirection.left1, true, false, typeVertex.center);

        tmpY = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0VerLength));
        Vertex[10] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0HorLength)) / 2;
        Vertex[11] = new MyVertex(tmpX, tmpY, typeDirection.up1, false, false, typeVertex.center);

        tmpX = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0HorLength));
        Vertex[12] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0VerLength)) - (Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1VerLength)) + Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1VerOffset))) / 2;
        Vertex[13] = new MyVertex(tmpX, tmpY, typeDirection.right1, true, false, typeVertex.center);

        tmpY = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1VerOffset)) + Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1VerLength));
        Vertex[14] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        // Need midpoint of edge
        tmpX = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0HorLength));
        tmpX += (Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1HorLength)) - Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0HorLength))) / 2 / 2; // midpoint
        tmpX -= Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1HorOffset));
        Vertex[15] = new MyVertex(tmpX, tmpY, typeDirection.up1, false, false, typeVertex.center);

        tmpX = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1HorLength));
        tmpX += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1HorOffset));
        Vertex[16] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY -= Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1VerLength)) / 2;
        Vertex[17] = new MyVertex(tmpX, tmpY, typeDirection.right1, true, false, typeVertex.center);

        tmpY -= Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1VerLength)) / 2;
        Vertex[18] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        // Need midpoint of edge
        tmpX = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0HorLength));
        tmpX += (Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1HorLength)) - Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0HorLength))) / 2 / 2; // midpoint
        tmpX -= Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1HorOffset));
        Vertex[19] = new MyVertex(tmpX, tmpY, typeDirection.down1, false, false, typeVertex.center);

        tmpX = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0HorLength));
        Vertex[20] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY -= Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1VerOffset)) / 2;
        Vertex[21] = new MyVertex(tmpX, tmpY, typeDirection.right1, true, false, typeVertex.center);

        tmpY = 0.0;
        Vertex[22] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0HorLength)) / 2;
        Vertex[23] = new MyVertex(tmpX, tmpY, typeDirection.down1, true, false, typeVertex.center);

        processEdgesForRounding();

        shapeValid = true;
    }

    private void Sshape()
    {
        configureArrays();
        // Sort out the tips by setting 'true' to center vertex that defines tip in each case.
        switch (layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex))
        {
            case (int)ShapeSettings.tipLocations.none: // None
                break;
            case (int)ShapeSettings.tipLocations.L: // Left
                tips[1] = true;
                tips[9] = true;
                break;
            case (int)ShapeSettings.tipLocations.R: // Right
                tips[13] = true;
                tips[21] = true;
                break;
            case (int)ShapeSettings.tipLocations.LR: // Left and right
                tips[1] = true;
                tips[9] = true;
                tips[13] = true;
                tips[21] = true;
                break;
            case (int)ShapeSettings.tipLocations.T: // Top
                tips[11] = true;
                break;
            case (int)ShapeSettings.tipLocations.B: // Bottom
                tips[23] = true;
                break;
            case (int)ShapeSettings.tipLocations.TB: // Top and Bottom
                tips[11] = true;
                tips[23] = true;
                break;
            case (int)ShapeSettings.tipLocations.TL: // Top and Left
                tips[1] = true;
                tips[9] = true;
                tips[11] = true;
                break;
            case (int)ShapeSettings.tipLocations.TR: // Top and Right
                tips[11] = true;
                tips[13] = true;
                tips[21] = true;
                break;
            case (int)ShapeSettings.tipLocations.TLR: // Top and Left and Right
                tips[1] = true;
                tips[9] = true;
                tips[11] = true;
                tips[13] = true;
                tips[21] = true;
                break;
            case (int)ShapeSettings.tipLocations.BL: // Bottom and Left
                tips[1] = true;
                tips[9] = true;
                tips[23] = true;
                break;
            case (int)ShapeSettings.tipLocations.BR: // Bottom and Right
                tips[13] = true;
                tips[21] = true;
                tips[23] = true;
                break;
            case (int)ShapeSettings.tipLocations.BLR: // Bottom and Left and Right
                tips[1] = true;
                tips[9] = true;
                tips[13] = true;
                tips[21] = true;
                tips[23] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBL: // Top and Bottom and Left
                tips[1] = true;
                tips[9] = true;
                tips[11] = true;
                tips[23] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBR: // Top and Bottom and Right
                tips[11] = true;
                tips[13] = true;
                tips[21] = true;
                tips[23] = true;
                break;
            default: // All
                tips[1] = true;
                tips[9] = true;
                tips[11] = true;
                tips[13] = true;
                tips[21] = true;
                tips[23] = true;
                break;
        }
        switch (layerSettings.getInt(ShapeSettings.properties_i.subShape2TipLocIndex)) // Bottom notch
        {
            case (int)ShapeSettings.tipLocations.none: // None
                break;
            case (int)ShapeSettings.tipLocations.L: // Left
                tips[5] = true;
                break;
            case (int)ShapeSettings.tipLocations.R: // Right
                break;
            case (int)ShapeSettings.tipLocations.LR: // Left and right
                tips[5] = true;
                break;
            case (int)ShapeSettings.tipLocations.T: // Top
                tips[3] = true;
                break;
            case (int)ShapeSettings.tipLocations.B: // Bottom
                tips[7] = true;
                break;
            case (int)ShapeSettings.tipLocations.TB: // Top and Bottom
                tips[3] = true;
                tips[7] = true;
                break;
            case (int)ShapeSettings.tipLocations.TL: // Top and Left
                tips[3] = true;
                tips[5] = true;
                break;
            case (int)ShapeSettings.tipLocations.TR: // Top and Right
                tips[3] = true;
                break;
            case (int)ShapeSettings.tipLocations.TLR: // Top and Left and Right
                tips[3] = true;
                tips[5] = true;
                break;
            case (int)ShapeSettings.tipLocations.BL: // Bottom and Left
                tips[7] = true;
                break;
            case (int)ShapeSettings.tipLocations.BR: // Bottom and Right
                tips[5] = true;
                tips[7] = true;
                break;
            case (int)ShapeSettings.tipLocations.BLR: // Bottom and Left and Right
                tips[5] = true;
                tips[7] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBL: // Top and Bottom and Left
                tips[3] = true;
                tips[7] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBR: // Top and Bottom and Right
                tips[3] = true;
                tips[5] = true;
                tips[7] = true;
                break;
            default: // All
                tips[3] = true;
                tips[5] = true;
                tips[7] = true;
                break;
        }
        switch (layerSettings.getInt(ShapeSettings.properties_i.subShape3TipLocIndex)) // Top notch
        {
            case (int)ShapeSettings.tipLocations.none: // None
                break;
            case (int)ShapeSettings.tipLocations.L: // Left
                break;
            case (int)ShapeSettings.tipLocations.R: // Right
                tips[17] = true;
                break;
            case (int)ShapeSettings.tipLocations.LR: // Left and right
                tips[17] = true;
                break;
            case (int)ShapeSettings.tipLocations.T: // Top
                tips[19] = true;
                break;
            case (int)ShapeSettings.tipLocations.B: // Bottom
                tips[15] = true;
                break;
            case (int)ShapeSettings.tipLocations.TB: // Top and Bottom
                tips[15] = true;
                tips[19] = true;
                break;
            case (int)ShapeSettings.tipLocations.TL: // Top and Left
                tips[19] = true;
                break;
            case (int)ShapeSettings.tipLocations.TR: // Top and Right
                tips[17] = true;
                tips[19] = true;
                break;
            case (int)ShapeSettings.tipLocations.TLR: // Top and Left and Right
                tips[17] = true;
                tips[19] = true;
                break;
            case (int)ShapeSettings.tipLocations.BL: // Bottom and Left
                tips[15] = true;
                break;
            case (int)ShapeSettings.tipLocations.BR: // Bottom and Right
                tips[15] = true;
                tips[17] = true;
                break;
            case (int)ShapeSettings.tipLocations.BLR: // Bottom and Left and Right
                tips[15] = true;
                tips[17] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBL: // Top and Bottom and Left
                tips[17] = true;
                tips[19] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBR: // Top and Bottom and Right
                tips[15] = true;
                tips[17] = true;
                tips[19] = true;
                break;
            default: // All
                tips[15] = true;
                tips[17] = true;
                tips[19] = true;
                break;
        }

        // NOTE:
        // Subshape offsets are applied later to simplify ellipse generation
        // Horizontal and vertical global offsets are applied in callsite.
        // Repositioning with respect to subshape reference is also applied in callsite.
        double tmpX = 0.0;
        double tmpY = 0.0;
        Vertex[0] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1VerOffset)) / 2;
        Vertex[1] = new MyVertex(tmpX, tmpY, typeDirection.left1, true, false, typeVertex.center);

        tmpY = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1VerOffset));
        Vertex[2] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1HorLength)) / 2;
        Vertex[3] = new MyVertex(tmpX, tmpY, typeDirection.up1, false, false, typeVertex.center);

        tmpX = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1HorLength));
        Vertex[4] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1VerLength)) / 2;
        Vertex[5] = new MyVertex(tmpX, tmpY, typeDirection.left1, true, false, typeVertex.center);

        tmpY = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1VerOffset)) + Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1VerLength));
        Vertex[6] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s1HorLength)) / 2;
        Vertex[7] = new MyVertex(tmpX, tmpY, typeDirection.down1, false, false, typeVertex.center);

        tmpX = 0;
        Vertex[8] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY = (Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0VerLength)) - tmpY) / 2;
        Vertex[9] = new MyVertex(tmpX, tmpY, typeDirection.left1, true, false, typeVertex.center);

        tmpY = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0VerLength));
        Vertex[10] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0HorLength)) / 2;
        Vertex[11] = new MyVertex(tmpX, tmpY, typeDirection.up1, false, false, typeVertex.center);
        // Center so no rounding definition

        tmpX = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0HorLength));
        Vertex[12] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY -= Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s2VerOffset)) / 2;
        Vertex[13] = new MyVertex(tmpX, tmpY, typeDirection.right1, true, false, typeVertex.center);
        // Center so no rounding definition

        tmpY = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0VerLength)) - Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s2VerOffset));
        Vertex[14] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX -= Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s2HorLength) / 2);
        Vertex[15] = new MyVertex(tmpX, tmpY, typeDirection.down1, false, false, typeVertex.center);

        tmpX = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s2HorOffset));
        Vertex[16] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY -= Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s2VerLength) / 2);
        Vertex[17] = new MyVertex(tmpX, tmpY, typeDirection.right1, false, false, typeVertex.center);

        tmpY = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0VerLength)) - Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s2VerLength)) - Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s2VerOffset));
        Vertex[18] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s2HorLength) / 2);
        Vertex[19] = new MyVertex(tmpX, tmpY, typeDirection.up1, false, false, typeVertex.center);

        tmpX = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0HorLength));
        Vertex[20] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY /= 2;
        Vertex[21] = new MyVertex(tmpX, tmpY, typeDirection.right1, false, false, typeVertex.center);

        tmpY = 0;
        Vertex[22] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.s0HorLength)) / 2;
        Vertex[23] = new MyVertex(tmpX, tmpY, typeDirection.down1, false, false, typeVertex.center);

        processEdgesForRounding();

        shapeValid = true;
    }

    private void processEdgesForRounding()
    {
        int horEdge = Vertex.Length - 2; // deal with padding.
        int verEdge = 1;
        for (int r = 0; r < round1.Length - 1; r++)
        {
            try
            {
                round1[r].index = r * 2;
                round1[r].horFace = horEdge;
                round1[r].verFace = verEdge;

                // Figure out our corner type. First is a special case.
                if (r == 0)
                {
                    round1[r].direction = typeRound.exter;
                    round1[r].MaxRadius = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.oCR));
                }
                else
                {
                    if (
                        Vertex[round1[r].verFace].direction == typeDirection.right1 && Vertex[round1[r].horFace].direction == typeDirection.up1 && horEdge < verEdge ||
                        Vertex[round1[r].verFace].direction == typeDirection.left1 && Vertex[round1[r].horFace].direction == typeDirection.up1 && horEdge > verEdge ||
                        Vertex[round1[r].verFace].direction == typeDirection.right1 && Vertex[round1[r].horFace].direction == typeDirection.down1 && horEdge > verEdge ||
                        Vertex[round1[r].verFace].direction == typeDirection.left1 && Vertex[round1[r].horFace].direction == typeDirection.down1 && horEdge < verEdge
                    )
                    {
                        round1[r].direction = typeRound.exter;
                        round1[r].MaxRadius = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.oCR));
                    }
                    else
                    {
                        round1[r].direction = typeRound.inner;
                        round1[r].MaxRadius = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.iCR));
                    }
                }

                // Small fudge for the 0 case
                if (r == 0)
                {
                    horEdge = -1;
                }

                // Change our edge configuration for the next loop. We need to handle overflow references as well
                if (r % 2 == 0)
                {
                    horEdge += 4;
                    horEdge %= Vertex.Length;
                }
                else
                {
                    verEdge += 4;
                    verEdge %= Vertex.Length;
                }
            }
            catch (Exception)
            {
            }
        }

        // First and last are the same.
        round1[^1] = round1[0];
    }

    // Intended to take geometry from an external source and map it into our shape engine.
    private void customShape(GeoLibPointF[] sourcePoly)
    {
        if (sourcePoly == null)
        {
            shapeValid = false;
            return;
        }

        // Note that we assume the point order matches our general primitives; might need upstream review to ensure this is being
        // fed correctly.
        // Upstream should trim array to ensure end point is different from start point, but we'll force the issue here for robustness.
        sourcePoly = GeoWrangler.stripTerminators(sourcePoly, true);
        sourcePoly = GeoWrangler.stripColinear(sourcePoly);
        //  Strip the terminator again to meet the requirements below.
        sourcePoly = GeoWrangler.stripTerminators(sourcePoly, false);
        sourcePoly = GeoWrangler.clockwise(sourcePoly);

        // We need to look at our incoming shape to see whether it's orthogonal and suitable for contouring.
        geoCoreShapeOrthogonal = GeoWrangler.orthogonal(sourcePoly, angularTolerance:0.003);

        if (!geoCoreShapeOrthogonal)
        {
            customShape_nonOrthogonal(sourcePoly);
        }
        else
        {
            customShape_orthogonal(sourcePoly);
        }

        shapeValid = true;
    }

    private void customShape_nonOrthogonal(GeoLibPointF[] sourcePoly)
    {
        int sCount = sourcePoly.Length;
        Vertex = new MyVertex[sCount + 1]; // add one to close.
        // Assign shape vertices to Vertex and move on. EntropyShape will know what to do.
#if !SHAPELIBSINGLETHREADED
        Parallel.For(0, sCount, pt => 
#else
            for (int pt = 0; pt < sCount; pt++)
#endif
            {
                Vertex[pt] = new MyVertex(sourcePoly[pt].X, sourcePoly[pt].Y, typeDirection.tilt1, false, false, typeVertex.corner);
            }
#if !SHAPELIBSINGLETHREADED
        );
#endif
        // Close the shape.
        Vertex[^1] = new MyVertex(Vertex[0]);
    }

    private void customShape_orthogonal(GeoLibPointF[] sourcePoly)
    {
        int sCount = sourcePoly.Length;
        int vertexCount = sCount * 2 + 1; // assumes no point in midpoint of edges, and 1 to close.
        Vertex = new MyVertex[vertexCount];
        tips = new bool[vertexCount];
        int vertexCounter = 0; // set up our vertex counter.

#if !SHAPELIBSINGLETHREADED
        Parallel.For(0, vertexCount, i =>
#else
            for (Int32 i = 0; i < vertexCount; i++)
#endif
            {
                tips[i] = false;
            }
#if !SHAPELIBSINGLETHREADED
        );
#endif

        int roundCount = sourcePoly.Length + 1;
        round1 = new MyRound[roundCount];
#if !SHAPELIBSINGLETHREADED
        Parallel.For(0, roundCount, i =>
#else
            for (Int32 i = 0; i < roundCount; i++)
#endif
            {
                round1[i] = new MyRound();
            }
#if !SHAPELIBSINGLETHREADED
        );
#endif
        // Set up first rounding entry
        round1[0].direction = typeRound.exter;
        round1[0].MaxRadius = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.oCR));
        round1[0].verFace = 1;
        round1[0].horFace = vertexCount - 2;
        round1[^1] = round1[0]; // close the loop

        // Set up first vertex.
        Vertex[0] = new MyVertex(sourcePoly[0].X, sourcePoly[0].Y, typeDirection.tilt1, false, false, typeVertex.corner);
        vertexCounter++;
        // Set up first midpoint.
        Vertex[1] = new MyVertex((sourcePoly[0].X + sourcePoly[1].X) / 2.0f, (sourcePoly[0].Y + sourcePoly[1].Y) / 2.0f, typeDirection.left1, true, false, typeVertex.center);
        if (layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) == (int)ShapeSettings.tipLocations.L ||
            layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) == (int)ShapeSettings.tipLocations.LR ||
            layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) == (int)ShapeSettings.tipLocations.BL ||
            layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) == (int)ShapeSettings.tipLocations.TL ||
            layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) == (int)ShapeSettings.tipLocations.TBL ||
            layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) == (int)ShapeSettings.tipLocations.BLR ||
            layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) == (int)ShapeSettings.tipLocations.TLR ||
            layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) == (int)ShapeSettings.tipLocations.all)
        {
            tips[vertexCounter] = true;
        }
        vertexCounter++;

        // Also set our end points
        Vertex[vertexCount - 2] = new MyVertex((sourcePoly[0].X + sourcePoly[^1].X) / 2.0f,
            (sourcePoly[0].Y + sourcePoly[^1].Y) / 2.0f, typeDirection.down1, false, false, typeVertex.center);

        // Figure out our rounding characteristics.

        // First edge is always vertical, left facing.
        bool left = true;
        bool up = false;

        for (int pt = 1; pt < roundCount - 1; pt++)
        {
            // Link to our vertical/horizontal edges
            round1[pt].index = vertexCounter;
            if (pt % 2 == 1)
            {
                round1[pt].verFace = vertexCounter - 1;
                round1[pt].horFace = vertexCounter + 1;
            }
            else
            {
                round1[pt].verFace = vertexCounter + 1;
                round1[pt].horFace = vertexCounter - 1;
            }

            // Register our corner point into the vertex array.
            Vertex[vertexCounter] = new MyVertex(sourcePoly[pt].X, sourcePoly[pt].Y, typeDirection.tilt1, false, false, typeVertex.corner);
            vertexCounter++;

            // Now we have to wrangle the midpoint.

            int next = (pt + 1) % sourcePoly.Length; // wrap to polygon length

            // Find the normal for the edge to the next point.

            double dx = sourcePoly[next].X - sourcePoly[pt].X;
            double dy = sourcePoly[next].Y - sourcePoly[pt].Y;

            // Set up our midpoint for convenience.
            GeoLibPointF midPt = new(sourcePoly[pt].X + dx / 2.0f, sourcePoly[pt].Y + dy / 2.0f);

            // The normal, to match convention in the distance calculation is assessed from this point to the next point.

            // Get average angle for this vertex based on angles from line segments.
            // http://stackoverflow.com/questions/1243614/how-do-i-calculate-the-normal-vector-of-a-line-segmen
            GeoLibPointF normalPt = new(-dy, dx);

            // Vertical edge has a normal with an X value non-zero and Y value ~0.
            // treating a 0.01 difference as being ~0
            bool vertical = Math.Abs(normalPt.X) > 0.01;

            // Assess the normal to establish direction
            if (vertical)
            {
                // left facing vertical edge has normal with negative X value.
                left = normalPt.X < 0;
            }
            else
            {
                // down facing horizontal edge has normal with negative Y value.
                up = !(normalPt.Y < 0);
            }

            if (!vertical)
            {
                if (up)
                {
                    Vertex[vertexCounter] = new MyVertex(midPt.X, midPt.Y, typeDirection.up1, vertical, false, typeVertex.center);
                    if (layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) == (int)ShapeSettings.tipLocations.T ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) == (int)ShapeSettings.tipLocations.TB ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) == (int)ShapeSettings.tipLocations.TL ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) == (int)ShapeSettings.tipLocations.TR ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) == (int)ShapeSettings.tipLocations.TBL ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) == (int)ShapeSettings.tipLocations.TBR ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) == (int)ShapeSettings.tipLocations.TLR ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) == (int)ShapeSettings.tipLocations.all)
                    {
                        tips[vertexCounter] = true;
                    }
                }
                else
                {
                    Vertex[vertexCounter] = new MyVertex(midPt.X, midPt.Y, typeDirection.down1, vertical, false, typeVertex.center);
                    if (layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) == (int)ShapeSettings.tipLocations.B ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) == (int)ShapeSettings.tipLocations.TB ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) == (int)ShapeSettings.tipLocations.BL ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) == (int)ShapeSettings.tipLocations.BR ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) == (int)ShapeSettings.tipLocations.TBL ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) == (int)ShapeSettings.tipLocations.TBR ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) == (int)ShapeSettings.tipLocations.BLR ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) == (int)ShapeSettings.tipLocations.all)
                    {
                        tips[vertexCounter] = true;
                    }
                }
            }
            else
            {
                if (left)
                {
                    Vertex[vertexCounter] = new MyVertex(midPt.X, midPt.Y, typeDirection.left1, vertical, false, typeVertex.center);
                    if (layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) == (int)ShapeSettings.tipLocations.L ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) == (int)ShapeSettings.tipLocations.LR ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) == (int)ShapeSettings.tipLocations.BL ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) == (int)ShapeSettings.tipLocations.TL ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) == (int)ShapeSettings.tipLocations.TBL ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) == (int)ShapeSettings.tipLocations.BLR ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) == (int)ShapeSettings.tipLocations.TLR ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) == (int)ShapeSettings.tipLocations.all)
                    {
                        tips[vertexCounter] = true;
                    }
                }
                else
                {
                    Vertex[vertexCounter] = new MyVertex(midPt.X, midPt.Y, typeDirection.right1, vertical, false, typeVertex.center);
                    if (layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) == (int)ShapeSettings.tipLocations.R ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) == (int)ShapeSettings.tipLocations.LR ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) == (int)ShapeSettings.tipLocations.BR ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) == (int)ShapeSettings.tipLocations.TR ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) == (int)ShapeSettings.tipLocations.TBR ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) == (int)ShapeSettings.tipLocations.BLR ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) == (int)ShapeSettings.tipLocations.TLR ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) == (int)ShapeSettings.tipLocations.all)
                    {
                        tips[vertexCounter] = true;
                    }
                }
            }
            vertexCounter++;
        }

        // Reprocess our corners for inner/outer rounding based on horFace/verFace directions
#if !SHAPELIBSINGLETHREADED
        Parallel.For(0, roundCount, pt => 
#else
            for (int pt = 0; pt < roundCount; pt++)
#endif
            {
                // Only certain changes in direction correspond to an outer vertex, for a clockwise ordered series of points.
                bool outerVertex = pt == 0 || pt == round1.Length - 1 ||
                                   round1[pt].verFace < round1[pt].horFace && Vertex[round1[pt].verFace].direction == typeDirection.left1 && Vertex[round1[pt].horFace].direction == typeDirection.up1 ||
                                   round1[pt].verFace > round1[pt].horFace && Vertex[round1[pt].horFace].direction == typeDirection.up1 && Vertex[round1[pt].verFace].direction == typeDirection.right1 ||
                                   round1[pt].verFace < round1[pt].horFace && Vertex[round1[pt].verFace].direction == typeDirection.right1 && Vertex[round1[pt].horFace].direction == typeDirection.down1 ||
                                   round1[pt].verFace > round1[pt].horFace && Vertex[round1[pt].horFace].direction == typeDirection.down1 && Vertex[round1[pt].verFace].direction == typeDirection.left1;

                if (outerVertex)
                {
                    round1[pt].direction = typeRound.exter;
                    round1[pt].MaxRadius = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.oCR));
                }
                else
                {
                    round1[pt].direction = typeRound.inner;
                    round1[pt].MaxRadius = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.iCR));
                }

                Vertex[round1[pt].index].inner = !outerVertex;
            }
#if !SHAPELIBSINGLETHREADED
        );
#endif
    }
}