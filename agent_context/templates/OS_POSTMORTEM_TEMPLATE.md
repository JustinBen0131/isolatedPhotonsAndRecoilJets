# OS Postmortem Template

Use this for agent/OS failures, near misses, or repeated friction. Keep it
short enough to write immediately.

## Summary

- Date/time:
- Workstream:
- Severity: low / medium / high / critical
- What happened:
- User-visible impact:

## Evidence

- Exact command, file, email, job ID, deck, plot, or log:
- First known bad state:
- Detection point:

## Classification

Choose all that apply:

- missing context
- stale state
- missing duplicate guard
- missing scope gate
- weak policy
- missing executable check
- connector/session drift
- human-agent handoff gap
- plot/slide provenance gap
- SDCC/Condor operational issue

## Root Cause

What made the wrong action possible?

## Recovery

What fixed or contained it?

## Prevention

At least one required:

- policy patch:
- script/doctor check:
- guard rule:
- artifact/provenance rule:
- task/register correction:
- Linear follow-up:
- accepted risk with reason:

Prevented by next time:

- exact policy/script/guard/dashboard mechanism:

## Follow-Up Owner

- Owner:
- Due/stale-after:
- Link to Linear/register:
