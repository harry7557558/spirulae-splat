# SfM unused-feature compaction

`--compact-unused-features` removes, in memory only, feature rows that no
stored match record references. It is a representation change: it does not rank
features, remove pairs or correspondences, build tracks, or change mapper,
camera, bundle-adjustment, assembly, or merge options. Input feature and match
files remain unchanged. `auto` and `map` both take the switch; it defaults off.

On `auto` the pass runs after `matches.bin` has been written and after the
descriptors are released, so the file on disk keeps indexing the feature files
and the compaction sees the same arrays `map` would have loaded.

## Mapping invariants

Mapping uses feature counts to size per-feature arrays and correspondence-graph
rows. Registration, triangulation, visibility, audit, split, merge, and bundle
adjustment consume supported features, stored correspondences, or reconstructed
observations; they do not normalize a decision by the total number of extracted
features. Removing a row with no stored correspondence therefore changes the
representation, not a mapping decision.

An image with no referenced feature remains in the image table and camera
group. It has empty feature arrays, contributes no correspondence-graph edge,
and cannot be selected for registration.

## Transformation

1. Read `matches.bin` and allocate one `std::vector<uint32_t>` old-to-new map
   per image, with one entry per original feature and `UINT32_MAX` as the
   unreferenced sentinel.
2. Validate both image IDs and both feature indices of every stored match,
   without filtering on pair configuration.
3. Assign retained indices in increasing original-row order using a `uint64_t`
   counter checked before conversion to the stored `uint32_t` endpoint type.
4. Compact each feature set's keypoint and optional RGB rows -- during the
   parallel load in `map`, in place after matching in `auto`. The primitive also
   preserves descriptors when they are present.
5. Validate image feature counts, image names/order, camera metadata, pair
   identity/order/configuration, and pair/correspondence counts before remapping
   both endpoints in place.
6. Update only the in-memory image feature counts and endpoint indices.
7. Destroy the plan before camera setup and `Mapper` construction.

## What the written model inherits

`Image::points2D` holds one entry per feature, so a compacted run writes only
the matched keypoints into `images.bin`. Poses, cameras, 3D points and track
image IDs are unaffected, and the training-side reader skips the POINTS2D
block entirely (`src/data/parsers/ColmapParser.cpp`). What changes is that a
`point2D_idx` in such a model indexes the compacted feature array, not the
feature file's rows.

Two consequences:

- `--resume` and `--audit` read those indices back and apply them to the
  feature arrays of the *current* run, so the flag has to match the run that
  wrote the model. `Mapper::adopt` compares each image's keypoint count against
  this run's features, drops the observations of any image that disagrees, and
  warns; poses still adopt.
- `merge` refuses to align two models whose shared images disagree on keypoint
  count (`sharedImages` in `src/sfm/map/Merge.h`), so a compacted and an
  uncompacted model of one capture will not merge.

An external COLMAP consumer that expects the full keypoint list per image, or
that cross-references feature-file rows against the model, sees the compacted
indexing.

## Feature files must match the database

`compactFeatureSet` throws when a feature file's keypoint count disagrees with
the `num_features` recorded in `matches.bin`. Without the flag,
`CorrespondenceGraph::build` silently skips out-of-range endpoints instead, so
a stale feature directory that used to map badly now fails outright.

## Temporary storage and lifetime

The only temporary proportional to the original feature count is
`old_to_new`: one `uint32_t` entry per extracted feature. `compact_counts` and
the image snapshots hold one small entry per image; the pair snapshots hold two
`uint32_t` image IDs, one `int32_t` configuration, and one `size_t` match count
per stored pair. The compact feature sets coexist with their individual loaded
source set only while that image is compacted. The complete plan is reset before
camera setup and before Mapper allocates its own per-feature structures.

No separate peak-memory claim is inferred from this lifetime description; the
measured end-to-end validation results belong in the pull-request evidence.
