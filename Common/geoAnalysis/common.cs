namespace geoAnalysis;

/// <summary>
/// Constants used throughout the geoAnalysis library for geometric calculations.
/// Provides consistent tolerance values and calculation modes.
/// </summary>
public static class ConstantsGA
{
    /// <summary>
    /// Default geometric tolerance value used for floating-point comparisons and calculations.
    /// Set to 0.001 to balance precision with numerical stability.
    /// </summary>
    public const double tolerance = 0.001;
}

/// <summary>
/// Defines supported calculation modes for geometric analysis operations.
/// These modes determine the type of geometric analysis to be performed.
/// </summary>
public static class Supported
{
    /// <summary>
    /// Enumeration of supported geometric calculation modes.
    /// </summary>
    public enum calcModes 
    { 
        /// <summary>Area calculation between polygon sets</summary>
        area, 
        /// <summary>Enclosure, spacing, and overlap analysis</summary>
        enclosure_spacing_overlap, 
        /// <summary>Chord length analysis</summary>
        chord, 
        /// <summary>Angle analysis between intersecting edges</summary>
        angle 
    }
}