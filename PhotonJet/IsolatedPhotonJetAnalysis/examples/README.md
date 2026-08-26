# Examples

`fixture_walkthrough.sh` runs the complete implemented downstream vertical
slice on the two-event p+p fixture: tree projection and validation, event-
leading histogram and skim, ABCD occupancy, a provenance-compiled PNG, and a
generated ROOT annotation include.

```bash
./examples/fixture_walkthrough.sh outputs/example_pp
```

The fixture is for engineering regression only. Its plot deliberately says
`Engineering fixture`; it is not a physics result.

The script contains analysis arguments but no duplicated plot labels. Every
visible cut is compiled from the executable selection receipt. The same
pattern applies to real dataset manifests once their upstream products and
exposure receipts are available.
