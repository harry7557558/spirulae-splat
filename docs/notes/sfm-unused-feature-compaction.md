# SfM unused-feature compaction

`spirula sfm map --compact-unused-features` removes, in memory only, feature
rows that no stored match record references. It is a representation change:
it does not rank features, remove pairs or correspondences, build tracks, or
change mapper, camera, bundle-adjustment, assembly, or merge options. Input
feature and match files remain unchanged. The map-only switch defaults off.

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
4. Read each feature file without descriptors in the map CLI and compact its
   keypoint and optional RGB rows. The reusable primitive also preserves
   descriptors when they are present.
5. Validate image feature counts, image names/order, camera metadata, pair
   identity/order/configuration, and pair/correspondence counts before remapping
   both endpoints in place.
6. Update only the in-memory image feature counts and endpoint indices.
7. Destroy the plan before camera setup and `Mapper` construction.

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
