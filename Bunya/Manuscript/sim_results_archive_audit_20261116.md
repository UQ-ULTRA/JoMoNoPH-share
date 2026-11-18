# Sim-Results archive audit

> Historical inventory, not current folder-layout or running instructions.
> The archive was subsequently reorganised into the audit package; old paths
> and temporary files below may no longer exist. See the
> [current audit README](../../reproducibility-BiomJ/README.md).

Checked 2026-11-16. Source: `Bunya/ENZAMET/Sim-Results/`.

## Conclusion

The archive contains the three simulation collections needed for the revised
manuscript: the main proposed sAFT-JM simulation, three additional alpha-prior
sensitivity simulations, and the JMbayes2 comparator simulations. Older
collections, failed placeholders and superseded batches must not be pooled
with these results.

No source data, production scripts or archived results were changed. This is an
inventory and identity check, not a new model fit or a new full diagnostics audit.

## Verified collections

| Collection | Archive location | Configurations | Batches | Replicates per batch | Datasets/fits represented | Selected uncompressed size |
|---|---|---:|---:|---:|---:|---:|
| Main sAFT-JM, default alpha SD log(2)/1.96 | 20260526 account slices plus rescue folder | 40 | 10000 | 4 | 40000 | 19.621 GiB |
| Alpha sensitivity, SD log(1.25)/1.96 | 20260529, maxF125 slices plus rescue | 10 | 2500 | 4 | 10000 | 4.907 GiB |
| Alpha sensitivity, SD log(1.5)/1.96 | 20260529, maxF150 slices plus rescue | 10 | 2500 | 4 | 10000 | 4.907 GiB |
| Alpha sensitivity, SD log(3)/1.96 | 20260529, maxF300 slices plus rescue | 10 | 2500 | 4 | 10000 | 4.907 GiB |
| JMbayes2 comparators | 20260602 main output ZIP plus rescue output ZIP | 40 | 4000 | 10 | 40000 datasets | 2.443 GiB |

The main/comparator grid crosses four distributions (LL120_230, W090_380,
W130_380, W210_380), five coefficient scenarios, and two censoring mechanisms
(controlled 50% and administrative). Main sAFT-JM settings include n = 1100,
five BP weights, and maximum follow-up 120 months. Additional alpha-sensitivity
fits cover the ten log-logistic configurations only. The default alpha prior
setting is already represented by the log-logistic subset of the main run.

The comparator records contain Cox-JM, Weibull-JM, log-logistic-JM and LMM
results; these do not contain the proposed sAFT-JM fits.

## Verification performed

- Main sAFT-JM: selected 10000 unique canonical batch filenames after applying
  19 rescue replacements. All 10000 files have identical MD5 checksums to the
  existing consolidated `output_20260525_ENZAMET-altDM-gamma5-1000rep-zeroprior-initinfo` files.
- Alpha maxF125: 2500 unique batches, applying six rescue replacements; all
  2500 are byte-identical by MD5 to the matching consolidated output folder.
- Alpha maxF150: 2500 unique batches, applying one rescue replacement; all
  2500 are byte-identical by MD5 to the matching consolidated output folder.
- Alpha maxF300: 2500 unique batches, applying six rescue replacements; all
  2500 are byte-identical by MD5 to the matching consolidated output folder.
- Comparator main ZIP: 3999 batches. Its rescue ZIP contains 23 batches,
  including the missing filename and 22 replacements. The union is 4000 unique
  batches; all 4000 match the consolidated comparator files by CRC-32 and
  uncompressed file size.
- Sampled 102 sAFT-JM batch contents, covering every configuration in each set
  and every selected rescue batch. All sampled batches contain four non-NULL
  results. Metadata confirm response centring, zero beta prior means, and raw
  alpha prior SDs 0.3536465, 0.1138488, 0.2068700 and 0.5605165 for the main,
  maxF125, maxF150 and maxF300 sets, respectively.
- The 19 files in the archive's main summary folder match the local main
  summary folder byte-for-byte. The parameter summary reports 1000 replicates
  for each of 40 configurations.
- The 22 files in the archive's comparator summary folder match the local
  comparator summary folder byte-for-byte.

File identity was checked for every selected batch. Replicate content was
sampled, rather than rereading every sAFT-JM object or recomputing summaries.

## Replacement and exclusion rules

1. Main results: combine AJM, DM, MZ and PM under `20260526`, then replace
   matching filenames with those in
   `rescue 20 timeout replicates 20260527/output_20260525_ENZAMET-altDM-gamma5-1000rep-zeroprior-initinfo`.
   The 19 rescue batches contain the previously identified 20 timeout replicates.
2. Alpha sensitivity: combine the four account slices separately for each
   maxF setting, then replace matching filenames from the corresponding
   `20260529/rescue run` output folder.
3. Exclude the 1870 `-AMARTIN-7400.rds` conflict files in the alpha AJM slices:
   623 for maxF125, 623 for maxF150, and 624 for maxF300. Every one was read and
   confirmed to contain four NULL entries; each file is 50 bytes. They are not
   additional valid replicates and must not enter aggregation.
4. Comparator: unpack the main output ZIP, then apply the rescue output ZIP
   with replacement. The missing original filename is
   `JMbayes2_All_k5_n1100_cens50_LL120_230_scen5_FU120_batch079.rds`.
5. Exclude `20260509`, `20260601 abandon`, and `version 1 - original submission`
   when preparing the current revised-manuscript results. Representative older
   files show different batching/metadata; original-submission filenames also
   include FU72 and cens0 settings.

## Recommended source for packaging

Use the already consolidated folders immediately under `Bunya/ENZAMET/` as
the simplest source for the reproducibility package. They have already applied
the account merges and rescue replacements, and the checks above establish
that the relevant archive selections match them. Keep `Sim-Results` unchanged
as the underlying archive.

These are original saved batch result objects, containing summaries and
diagnostics. Do not describe them as full posterior-draw archives without
checking their contents. File coverage also does not imply that every model fit
was successful: the existing comparator summary records seven Weibull-JM
failures, as noted in the previous inventory. Preserve their status information.

Detailed per-file comparison manifests and sampled metadata were written to
`/private/tmp/sim-results-audit-20261116/` during this audit. No reproducibility
package has been created yet.
