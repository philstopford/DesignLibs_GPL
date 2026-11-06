# Layout regression fixtures

`doetest.gds` and `doetest.oas` are the original DOE layouts from Variance.

The JSON fixtures capture geometry and settings at the shared-library boundary
from `ilb_4_test.variance`, `ilb_16.variance`, `ilb_16b.variance`, and `ilb_17.variance`. They avoid a
dependency on Variance project loading, UI, or a sibling checkout.

- ILB4: the original 50 nm bar produces an open notch; a 30 nm bar encloses a hole.
- ILB16/16b/17: suffixes 0 and 1 select modern and legacy rounding respectively.
- Source is layer 5's Boolean geometry before rounding; Clip is layer 6's geometry.
- ILB17 shares ILB16 geometry with rotation changed from -55 to -60 degrees,
  as in the original project XML.
- Expected areas include the finite-width slits chosen by the global planner;
  physical-hole counts and connectivity are checked independently.

The tests exercise ShapeLibrary rounding and rotation, overlap processing,
proximity bias, and the layer 7 Boolean directly, retaining area, topology,
spike, input-preservation, and seam-relocation assertions.

`stair.GDS` is the supplied 350 nm staircase layout (1 nm database units).
Its step edges range from 50 to 100 nm. GEOCORE short-edge regressions compare
modern contouring with the direct contour generator, including the hole path,
and verify that 0.05/2 nm thresholds leave the larger steps intact.
