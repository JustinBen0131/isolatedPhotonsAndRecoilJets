# Photon-trigger sampled luminosity

[photon_luminosity.py](photon_luminosity.py) uses one function for pp and AuAu:

```text
L(run) = MB_live / sigma_MB
         × photon_scaled / photon_live
         × pileup × mb_coverage
L(total) = sum over selected runs of L(run)
```

The MB live count measures the exposure before MB prescaling. Multiplying by
the photon scaled/live ratio retains the fraction sampled by the photon
trigger. No additional photon prescale factor should be applied to this
luminosity. Trigger efficiency/turn-on and photon acceptance are separate.

This is different from the narrow-vertex database MB reference made from
scaled MB counts. Do not compare the two totals as if they used the same
vertex acceptance and counter convention.

## Inputs and units

Supply one CSV row per selected run, with these mandatory columns:

```csv
run,mb_live,photon_live,photon_scaled,pileup,mb_coverage
1,1000000,10000,2500,1,1
```

The row above is synthetic. Use independently measured counters for physics.
Live counts precede prescaling but include DAQ live-time gating; scaled counts
follow prescaling. Check the trigger identity for each run and use matching
acquisition intervals for all counters. The script rejects missing numerical
inputs and duplicate runs; it does not infer the run list or missing coverage.

```bash
python3 scripts/data_prep/recoiljets/photon_luminosity.py pp.csv \
  --sigma-barn 0.0252 --unit 'pb^-1' > pp_luminosity.json
python3 scripts/data_prep/recoiljets/photon_luminosity.py auau.csv \
  --sigma-barn 6.324 --unit 'nb^-1' > auau_luminosity.json
```

The cross section is explicitly supplied in barns. The function returns
inverse barns; reporting divides by 10^12 for pb^-1 or 10^9 for nb^-1.
Output contains every run's inputs/result and the input file SHA-256.

## Reproducing the pp PPG12 convention

For the checked PPG12 sample, MB is bit 10 and photon is bit 30. From each
upstream luminosity counter summary, use:

| CSV input | Counter quantity |
| --- | --- |
| mb_live | sumgoodlive[10] |
| photon_live | sumgoodlive[30] |
| photon_scaled | sumgoodscaled[30] |
| pileup | nmbdc / nocorN |
| mb_coverage | tottrigcounts[5][10] / sumgoodscaled[10] |

Here `tottrigcounts[5][10]` is the observed all-vertex MB count. Because
`avgPS[bit] = sumgoodlive[bit] / sumgoodscaled[bit]`, the shared expression is
algebraically equivalent to the PPG12 expression:

```text
observed_MB_allz × avgPS[10] × (nmbdc/nocorN) / avgPS[30] / sigma_MB
```

`mb_coverage` is a measured count-convention adjustment, not a factor fitted
to the published luminosity. Setting it to one instead selects the direct
good-segment scaler convention. The stored pileup integrals are reused;
this cross-check does not independently recalibrate the pileup model.

The 1,428 positive-reference runs give **64.3717886308 pb^-1**, compared with
the reference sum **64.3717887 pb^-1**. Every run agrees within 10^-7 pb^-1
(maximum difference 5.66 × 10^-8 pb^-1). The reference table selects the
comparison roster and supplies comparison values, never correction factors.
With `mb_coverage=1`, the result is **64.5491481542 pb^-1**, or **0.275524%**
higher. This difference is the MB counting convention, not a pp-to-AuAu
correction to transfer.

## AuAu configuration and scope

The checked 2,448-run calculation uses wide-vertex MB bit 14 and photon bit
22 counters, with `sigma_MB=6.324 b`, `pileup=1`, and `mb_coverage=1`.
It gives **7.8825203013 nb^-1** from 2,061 positive-exposure runs and matches
the previous per-run calculation within 10^-12 nb^-1 for every run.
Trigger bit numbers are menu-dependent, not universal identifiers.

This verifies shared arithmetic. The AuAu cross-section/vertex calibration,
unity correction assumptions, and correspondence to processed segments must
be established for final yield normalization. Run-level exposure alone does
not establish the exposure of a partially processed analysis sample.

## Tests

```bash
python3 -m unittest discover -s scripts/data_prep/recoiljets/tests -p test_photon_luminosity.py
```

Tests cover the PPG12 algebra, units, zero exposure, invalid inputs, and
duplicate-run rejection. Full measured-input comparisons require the original
counter summaries; those data are not distributed with this script.
