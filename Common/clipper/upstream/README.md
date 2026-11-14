# Upstream alignment and maintenance

The maintained source is [PolyClipping/Clipper2](https://github.com/PolyClipping/Clipper2).

- **C++ reference:** `884ba821abc002d895ce72dcc21fe139c207229f` (`884ba82`, “Add Y-stripe point-in-polygon accelerator”).
- **Final C# baseline:** `8b130712b3999ec4d7d591687c88d89a721e20ce`, before C# removal in `27a64a6`.
- **Native comparison configuration:** C++17, `USINGZ`, default intersection precision, `-ffp-contract=off`.
- **Managed configuration:** .NET 10, existing `Clipper2Lib` namespace, `USINGZ` in all configurations.

Read [the library guide](../GUIDE.md) for concepts and examples, and [the developer guide](../DEVELOPER_GUIDE.md) for source interoperation, state invariants, debugging, and application integration.

## Compatibility policy

`Clipper64`, `ClipperD`, and `ClipperOffset` default to `Compatibility = ClipperCompatibility.Legacy`. Existing callers and static `Clipper` convenience methods retain established intersection truncation and offset square-join calculations. Offset output follows the first polygon group's input orientation by default in this mode.

Set `Compatibility = ClipperCompatibility.Upstream` on an instance to select upstream intersection rounding, canonical segment evaluation, the native right-edge join threshold, and square-join intersection construction. Offsets then default to positive outer winding. `PreserveInputOrientation` can explicitly override that winding choice in either mode, and `ReverseSolution` still reverses the selected winding.

The compatibility mode is local to an instance. A finishing clipper created by an offsetter inherits its mode. There is no global mutable switch and no change to the existing public method signatures, namespace, Z fields, callback types, or enum numeric values. The previous `InternalClipper.GetSegmentIntersectPt` name and `MultiplyUInt64Result` type remain available. Path-list string representations retain their previous newline contract.

Legacy mode is not a frozen old engine: fixes to orientation predicates, invalid contour handling, duplicate stripping, state cleanup, and other algorithms apply in both modes. Existing regression results pass unchanged, but untested historical behavior is not guaranteed. Migrate whole application pipelines using representative design files before selecting upstream mode in production.

Zero/insignificant polygon offsets now undergo cleanup; insignificant open strokes disappear. Empty paths are ignored; singleton paths survive duplicate stripping and can grow under an offset. These are intentional fixes in both modes.

## Port audit

The final C# baseline supplies the later C# engine fixes, exact-sign uses, offset and rectangle-clipping updates, and reusable object pools. Repository-specific helpers, namespace, Z build support, and public compatibility names are retained. The discontinued triangulation API was not added: current upstream removed that module too.

| C++ change after C# removal | Managed treatment |
|---|---|
| `f97d65d` cyclic split ownership | Mark emptied split records before recursion to prevent cycles |
| `9b460cb` skip sorted intersection passes | Detect already sorted active edges; reuse top X values |
| `7813a71`, `b496b45` local minima scanline merge | Merge the sorted pending minima with the scanline list |
| `32b0a72` cleanup on exception | Managed execution cleanup in `finally`; retain the existing boolean failure contract |
| `6c7d083` offset per-path state | Local end type and delta, fresh normals, shared arc-step setup for callbacks |
| `ca30beb` empty offset paths | Remove empty paths after duplicate stripping |
| `f770e28` square-join Z | Preserve source Z on both square-join vertices |
| `ec5f2d4` shared-vertex edge order | Compare directions and restore AEL ordering at shared vertices |
| `f774daf` rounded/stable intersections | Upstream mode uses ties-to-even and canonical segment direction; stable event ordering |
| `872bfa0` floating contraction | Ordinary managed arithmetic; native fixture generators explicitly disable contraction |
| `90f96cc` integer product/rounding intrinsics | Use .NET `Int128` product comparison and `Math.Round`; no native-intrinsic dependency |
| `4a3bf64` scanbeam clamp | Project onto the flatter edge's scanbeam portion and clamp Y after rounding |
| `0812644` arena allocation | Keep managed vertex/output/join pools and active-edge reuse; clear index caches on record reuse |
| `64a9ecc` offset orientation | Exact 128-bit area ranking, mixed-group nesting, stray-hole rejection, degenerate growth, positive winding in upstream mode |
| `884ba82` stripe PIP acceleration | Snapshot stripe index with fallback, cached owner queries, configurable threshold, and reuse cleanup |

Additional differences discovered by exact native comparison were aligned: local minima sort by descending Y then ascending X with stable ties, invalid/tiny closed-path checks, the split-triangle area boundary, and the asymmetric right-edge join distance tolerance (selected by compatibility mode).

C++ export/DLL signatures, compiler flags, platform-specific intrinsics, and build-system changes have no direct managed API equivalent. Native object arenas are represented by managed pools rather than translated byte-for-byte. The optional `CLIPPER2_HI_PRECISION` C++ variant and floating expansion predicates are **not implemented** by this C# port; comparisons target the default native build. No claim is made of identical performance, allocation patterns, or output for all possible inputs.

## Validation

The original 469 DesignLibs tests were retained without changing their expected geometry. The suite now contains 1,180 tests:

- 195 polygon and 16 open-path native fixtures, compared by complete coordinates and winding. Every case also checks PolyTree output with forced indexed, forced brute-force, and default nesting thresholds.
- 480 offset fixtures covering eight input arrangements, four joins, five end types, and three signed/zero deltas. Both flat and tree output are compared.
- 20 focused managed regression cases for compatibility winding, callbacks, pooled state, empty/degenerate paths, exact area, and index rebuilding.

Closed ring comparison ignores starting vertex and list order but does not reverse winding or remove vertices. Open paths retain their vertex order. The fixtures include complete input and output geometry rather than only counts or area tolerances.

Validation in this workspace:

| Check | Result |
|---|---|
| Complete DesignLibs Release suite | 1,180 passed, 0 failed |
| Guide C# examples | All 9 compile with 0 warnings/errors |
| Quilt layout/decomposition tests | All 10 pass in an isolated project referencing the unchanged application and test source |
| Quilt complete test project | Follow-up test-harness repair: all 44 pass. Exclude historical duplicate/standalone sources and update stale test API references and expectations to exercise production implementations. |
| Variance complete unit suite, updated library | 193 pass, 25 fail |
| Variance complete unit suite, original library | 193 pass, the identical 25 tests fail |

The Variance baseline was built from the original tracked Clipper sources and substituted into an isolated copy of the same application test output, with its checked-in test data available. The comparison found no new failing tests. Existing failures include null references, assertion problems, and timing benchmarks. This does not certify untested interactive application workflows. Application source files were not changed for these checks.

## Reproducing native reference fixtures

Obtain a checkout of upstream and check out the pinned C++ commit. Then, from the DesignLibs repository root:

```sh
sh Common/clipper/upstream/refresh-reference.sh /path/to/Clipper2
```

The script verifies the commit, builds two small reference generators in a temporary directory, and writes:

- `UnitTests/TestData/Clipper2/polygons.json`
- `UnitTests/TestData/Clipper2/lines.json`
- `UnitTests/TestData/Clipper2/offsets.json`

It requires a C++17 compiler and Python 3. It does not fetch from the network. The generators use the upstream polygon/line test loader and native clipping/offsetting implementations. The upstream Boost license is included as [LICENSE](LICENSE).

Run the checks:

```sh
dotnet test UnitTests/UnitTests.csproj -m:1 --filter FullyQualifiedName~ClipperUpstreamTests
dotnet test UnitTests/UnitTests.csproj -m:1 --configuration Release
```

For a future upstream update, inspect every change in `clipper.core.h`, `clipper.engine.*`, `clipper.offset.*`, `clipper.rectclip.*`, `clipper.minkowski.h`, and `clipper.h`, along with the native tests. Port relevant behavior, preserve deliberate compatibility decisions, then update both the pin and this audit. Regenerate fixtures only after reviewing the upstream changes; do not replace expected outputs merely to make failing tests pass. Re-run the application regression suites before promoting a new mode or changing defaults.
