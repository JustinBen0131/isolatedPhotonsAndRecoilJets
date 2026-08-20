# Photon-jet balance diagnostics

This directory contains read-only diagnostics for the reconstructed and
response-level photon-jet balance observable. The helpers consume explicit
JSON/CSV/ROOT inputs, emit plots plus provenance sidecars, and do not submit
jobs or mutate source products.

The scripts are validation tools, not a declaration that a plotted candidate
is unfolded or publication-ready. Interpret each output using the status and
science-boundary fields written by its generator.
