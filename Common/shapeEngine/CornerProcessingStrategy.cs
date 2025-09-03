using System;
using System.Diagnostics;
using System.Reflection;
using Clipper2Lib;
using geoLib;
using utility;

namespace shapeEngine
{
    /// <summary>
    /// Strategy selector for choosing optimal corner processing approach based on input complexity
    /// and performance requirements. Provides hybrid approach combining legacy and new algorithms.
    /// </summary>
    public static class CornerProcessingStrategy
    {
        /// <summary>
        /// Configuration options for corner processing performance vs quality trade-offs
        /// </summary>
        public class ProcessingConfig
        {
            /// <summary>
            /// If true, always use legacy approach regardless of complexity analysis
            /// </summary>
            public bool ForceLegacy { get; set; } = false;
            
            /// <summary>
            /// If true, always use new approach regardless of complexity analysis  
            /// </summary>
            public bool ForceNew { get; set; } = false;
            
            /// <summary>
            /// Enable parallel processing for applicable algorithms
            /// </summary>
            public bool EnableParallelProcessing { get; set; } = true;
            
            /// <summary>
            /// Maximum number of corners before switching to legacy approach for simple cases
            /// </summary>
            public int ComplexityThreshold { get; set; } = 8;
            
            /// <summary>
            /// Minimum resolution before considering high-quality approach
            /// </summary>
            public double HighQualityResolutionThreshold { get; set; } = 1.0;
            
            /// <summary>
            /// Enable performance profiling and logging
            /// </summary>
            public bool EnableProfiling { get; set; } = false;
        }

        /// <summary>
        /// Geometric complexity analysis result
        /// </summary>
        public class ComplexityAnalysis
        {
            public int CornerCount { get; set; }
            public double AverageCornerRadius { get; set; }
            public double MinEdgeLength { get; set; }
            public double MaxEdgeLength { get; set; }
            public bool HasSharpCorners { get; set; }
            public bool IsSimpleGeometry { get; set; }
            public string RecommendedApproach { get; set; } = "";
            public string Reasoning { get; set; } = "";
        }

        /// <summary>
        /// Analyzes geometric complexity to determine optimal processing approach
        /// </summary>
        public static ComplexityAnalysis AnalyzeComplexity(ShapeLibrary shapeLib, double resolution, 
            double cornerRadius, ProcessingConfig config)
        {
            var analysis = new ComplexityAnalysis();
            
            // Get vertex count as a proxy for complexity
            analysis.CornerCount = shapeLib.Vertex?.Length ?? 0;
            analysis.AverageCornerRadius = cornerRadius;
            
            // Calculate edge statistics if vertices are available
            if (shapeLib.Vertex != null && shapeLib.Vertex.Length > 2)
            {
                double minLength = double.MaxValue;
                double maxLength = 0;
                
                for (int i = 0; i < shapeLib.Vertex.Length - 1; i++)
                {
                    double dx = shapeLib.Vertex[i + 1].X - shapeLib.Vertex[i].X;
                    double dy = shapeLib.Vertex[i + 1].Y - shapeLib.Vertex[i].Y;
                    double length = Math.Sqrt(dx * dx + dy * dy);
                    
                    minLength = Math.Min(minLength, length);
                    maxLength = Math.Max(maxLength, length);
                }
                
                analysis.MinEdgeLength = minLength;
                analysis.MaxEdgeLength = maxLength;
                analysis.HasSharpCorners = cornerRadius < minLength * 0.1;
            }

            // Determine if geometry is simple
            analysis.IsSimpleGeometry = analysis.CornerCount <= config.ComplexityThreshold &&
                                       !analysis.HasSharpCorners &&
                                       resolution >= config.HighQualityResolutionThreshold;

            // Make recommendation
            if (config.ForceLegacy)
            {
                analysis.RecommendedApproach = "Legacy";
                analysis.Reasoning = "Forced legacy mode via configuration";
            }
            else if (config.ForceNew)
            {
                analysis.RecommendedApproach = "New";
                analysis.Reasoning = "Forced new mode via configuration";
            }
            else if (analysis.IsSimpleGeometry && analysis.CornerCount <= 4)
            {
                analysis.RecommendedApproach = "Legacy";
                analysis.Reasoning = "Simple rectangular geometry - legacy approach is efficient";
            }
            else if (analysis.HasSharpCorners || resolution < config.HighQualityResolutionThreshold)
            {
                analysis.RecommendedApproach = "New";
                analysis.Reasoning = "High quality requirements or sharp corners - new approach provides better results";
            }
            else if (analysis.CornerCount > config.ComplexityThreshold)
            {
                analysis.RecommendedApproach = "NewParallel";
                analysis.Reasoning = "Complex geometry benefits from parallel processing";
            }
            else
            {
                analysis.RecommendedApproach = "New";
                analysis.Reasoning = "Balanced case - new approach provides good quality/performance trade-off";
            }

            return analysis;
        }

        /// <summary>
        /// Process corners using hybrid strategy selection based on complexity analysis
        /// </summary>
        public static PathD ProcessCornersHybrid(ShapeLibrary shapeLib, bool previewMode, bool cornerCheck, 
            int cornerSegments, int optimizeCorners, double resolution, 
            bool iCPA = false, bool oCPA = false, double iCV = 0, double iCVariation_scalar = 0, 
            double oCV = 0, double oCVariation_scalar = 0, double shortEdgeLength = 0, 
            double maxShortEdgeLength = 0, ProcessingConfig? config = null)
        {
            config ??= new ProcessingConfig();
            
            var stopwatch = config.EnableProfiling ? Stopwatch.StartNew() : null;
            
            // Get corner radius for complexity analysis
            double iCR = Convert.ToDouble(shapeLib.LayerSettings.getDecimal(ShapeSettings.properties_decimal.iCR));
            double oCR = Convert.ToDouble(shapeLib.LayerSettings.getDecimal(ShapeSettings.properties_decimal.oCR));
            double avgRadius = (iCR + oCR) / 2.0;
            
            var complexity = AnalyzeComplexity(shapeLib, resolution, avgRadius, config);
            
            PathD result;
            
            switch (complexity.RecommendedApproach)
            {
                case "Legacy":
                    result = ProcessCornersLegacy(shapeLib, previewMode, cornerCheck, cornerSegments, 
                        optimizeCorners, resolution, iCPA, oCPA, iCV, iCVariation_scalar, oCV, oCVariation_scalar);
                    break;
                    
                case "NewParallel":
                    result = ProcessCornersNewWithParallel(shapeLib, previewMode, cornerCheck, cornerSegments,
                        optimizeCorners, resolution, iCPA, oCPA, iCV, iCVariation_scalar, oCV, oCVariation_scalar,
                        shortEdgeLength, maxShortEdgeLength, config);
                    break;
                    
                case "New":
                default:
                    result = shapeLib.processCorners(previewMode, cornerCheck, cornerSegments, optimizeCorners,
                        resolution, iCPA, oCPA, iCV, iCVariation_scalar, oCV, oCVariation_scalar,
                        shortEdgeLength, maxShortEdgeLength);
                    break;
            }
            
            if (config.EnableProfiling && stopwatch != null)
            {
                stopwatch.Stop();
                Console.WriteLine($"Corner processing: {complexity.RecommendedApproach} approach took {stopwatch.ElapsedMilliseconds}ms");
                Console.WriteLine($"Reasoning: {complexity.Reasoning}");
                Console.WriteLine($"Corners: {complexity.CornerCount}, Output points: {result?.Count ?? 0}");
            }
            
            return result;
        }

        /// <summary>
        /// Process corners using legacy approach via reflection
        /// </summary>
        private static PathD ProcessCornersLegacy(ShapeLibrary shapeLib, bool previewMode, bool cornerCheck,
            int cornerSegments, int optimizeCorners, double resolution, bool iCPA, bool oCPA, 
            double iCV, double iCVariation_scalar, double oCV, double oCVariation_scalar)
        {
            var method = typeof(ShapeLibrary).GetMethod("legacy_processCorners_actual",
                BindingFlags.NonPublic | BindingFlags.Instance);
                
            if (method == null)
            {
                throw new InvalidOperationException("Legacy corner processing method not found");
            }
            
            var result = method.Invoke(shapeLib, new object[]
            {
                previewMode, cornerCheck, cornerSegments, optimizeCorners, resolution,
                iCPA, oCPA, iCV, iCVariation_scalar, oCV, oCVariation_scalar
            });
            
            return (PathD)result;
        }

        /// <summary>
        /// Process corners using new approach with parallel processing optimizations
        /// </summary>
        private static PathD ProcessCornersNewWithParallel(ShapeLibrary shapeLib, bool previewMode, bool cornerCheck,
            int cornerSegments, int optimizeCorners, double resolution, bool iCPA, bool oCPA,
            double iCV, double iCVariation_scalar, double oCV, double oCVariation_scalar,
            double shortEdgeLength, double maxShortEdgeLength, ProcessingConfig config)
        {
            if (!config.EnableParallelProcessing)
            {
                return shapeLib.processCorners(previewMode, cornerCheck, cornerSegments, optimizeCorners,
                    resolution, iCPA, oCPA, iCV, iCVariation_scalar, oCV, oCVariation_scalar,
                    shortEdgeLength, maxShortEdgeLength);
            }

            // For now, delegate to standard approach
            // Future enhancement: implement parallel corner processing in contourGen
            return shapeLib.processCorners(previewMode, cornerCheck, cornerSegments, optimizeCorners,
                resolution, iCPA, oCPA, iCV, iCVariation_scalar, oCV, oCVariation_scalar,
                shortEdgeLength, maxShortEdgeLength);
        }

        /// <summary>
        /// Creates a configuration optimized for performance over quality
        /// </summary>
        public static ProcessingConfig CreatePerformanceConfig()
        {
            return new ProcessingConfig
            {
                ForceLegacy = false,
                ForceNew = false,
                EnableParallelProcessing = true,
                ComplexityThreshold = 6, // Lower threshold favors legacy for simple cases
                HighQualityResolutionThreshold = 2.0, // Higher threshold favors legacy
                EnableProfiling = false
            };
        }

        /// <summary>
        /// Creates a configuration optimized for quality over performance
        /// </summary>
        public static ProcessingConfig CreateQualityConfig()
        {
            return new ProcessingConfig
            {
                ForceLegacy = false,
                ForceNew = false,
                EnableParallelProcessing = true,
                ComplexityThreshold = 12, // Higher threshold favors new approach
                HighQualityResolutionThreshold = 0.5, // Lower threshold favors new approach
                EnableProfiling = false
            };
        }
        
        /// <summary>
        /// Creates a configuration that always uses the legacy approach
        /// </summary>
        public static ProcessingConfig CreateLegacyConfig()
        {
            return new ProcessingConfig
            {
                ForceLegacy = true,
                EnableProfiling = false
            };
        }
    }
}