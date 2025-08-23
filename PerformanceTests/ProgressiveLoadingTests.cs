using System;
using System.Threading.Tasks;
using Eto.Drawing;
using VeldridEto;
using Xunit;

namespace ProgressiveLoadingTests
{
    public class ProgressivePolygonLoadingTests
    {
        [Fact]
        public async Task BatchedPolygonLoader_ProcessesPolygonsInBatches()
        {
            // Arrange
            var settings = new OVPSettings();
            settings.setProgressiveLoading(true);
            settings.setBatchSize(10);
            
            // Add test polygons
            for (int i = 0; i < 25; i++)
            {
                var poly = CreateTestPolygon(i * 10, i * 10);
                settings.addPolygon(poly, Color.Red, 1.0f, false, i);
            }
            
            var loader = new BatchedPolygonLoader(settings);
            var processedBatches = 0;
            var totalPolygonsProcessed = 0;
            
            // Act
            var result = await loader.LoadPolygonsProgressivelyAsync(async batch =>
            {
                processedBatches++;
                totalPolygonsProcessed += batch.Polygons.Count;
                await Task.Delay(1); // Simulate processing time
            });
            
            // Assert
            Assert.True(result.Success);
            Assert.Equal(3, processedBatches); // 25 polygons / 10 batch size = 3 batches
            Assert.Equal(25, totalPolygonsProcessed);
            Assert.Equal(25, result.ProcessedCount);
        }
        
        [Fact]
        public void ProgressivePolygonRenderer_CreatesVertexData()
        {
            // Arrange
            var settings = new OVPSettings();
            var renderer = new ProgressivePolygonRenderer(settings);
            
            var batch = new PolygonBatch
            {
                Polygons = new()
                {
                    new ovp_Poly(CreateTestPolygon(0, 0), Color.Blue, 1.0f)
                },
                PolygonType = PolygonType.Foreground,
                BatchIndex = 0,
                StartIndex = 0,
                EndIndex = 0,
                IsLastBatch = true
            };
            
            // Act
            var batchData = renderer.ProcessBatch(batch);
            
            // Assert
            Assert.NotNull(batchData.PolyVertices);
            Assert.True(batchData.PolyVertices.Length > 0);
            Assert.NotNull(batchData.PolyIndices);
            Assert.True(batchData.PolyIndices.Length > 0);
        }
        
        [Fact]
        public void OVPSettings_ProgressiveLoadingConfiguration()
        {
            // Arrange
            var settings = new OVPSettings();
            
            // Act & Assert - Default values
            Assert.False(settings.progressiveLoadingEnabled());
            Assert.Equal(200, settings.getBatchSize());
            Assert.Equal(16, settings.getMaxBatchProcessingTimeMs());
            
            // Act & Assert - Set values
            settings.setProgressiveLoading(true);
            settings.setBatchSize(100);
            settings.setMaxBatchProcessingTimeMs(8);
            
            Assert.True(settings.progressiveLoadingEnabled());
            Assert.Equal(100, settings.getBatchSize());
            Assert.Equal(8, settings.getMaxBatchProcessingTimeMs());
        }
        
        [Fact]
        public void OVPSettings_BackwardCompatibility()
        {
            // Arrange & Act
            var settings = new OVPSettings();
            
            // Assert - Progressive loading disabled by default
            Assert.False(settings.progressiveLoadingEnabled());
            
            // Should be able to add polygons normally
            var poly = CreateTestPolygon(0, 0);
            settings.addPolygon(poly, Color.Green, 1.0f, false, 1);
            
            Assert.Equal(1, settings.polyList?.Count);
        }
        
        private static PointF[] CreateTestPolygon(float x, float y)
        {
            return new PointF[]
            {
                new(x, y),
                new(x + 10, y),
                new(x + 10, y + 10),
                new(x, y + 10),
                new(x, y) // Close polygon
            };
        }
    }
}