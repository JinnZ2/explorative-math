# Run schema and audit protocol

A dataset consists of `manifest.json`, `records.jsonl`, `observers.jsonl`, `divergence.jsonl`, and `notes.jsonl`. JSONL means one complete JSON object per newline-terminated line. UTF-8 JSON is canonicalized with sorted keys; nonfinite numbers and duplicate object keys are rejected. Initial observer, divergence, and note files are empty.

## Base runs

Every base run has `id`, `params`, `outputs`, `human_interest`, `ai_interest`, `compressibility`, `divergence_note`, and `status`. Optional identifying fields are `run_id`, `family`, and `kind`. An ID identifies one parameterized run, not a unique output value. Repeated outputs and separately enumerated equivalent rules remain distinct.

`params` records the rule or experiment, finite bounds, dimension, and seed as applicable. `outputs` retains the complete generated run, normally as an ordered list. Every trajectory contains step zero and all requested updates. Geometry pairs retain their vectors; shell samples retain their points; direction packings retain every vector and every unordered pair angle. Projection frames retain full shadows, composition details, degeneracies, and all catalog comparisons.

All four observer slots are JSON `null` in `records.jsonl` forever. Changes are explicit reports in a separate log, not destructive edits to the experiment. The `show` command overlays the latest report for each slot in memory.

Interest slots accept only booleans or `null`. `compressibility` accepts a nonnegative finite number or `null`. A numeric value requires a method description specifying the measurement convention; ratios may exceed one, and values from different methods are not inherently comparable. `divergence_note` is a string or `null`, filled only by an explicit report.

## Status

`OBSERVED` means an executed result or the observation that a named observer reported something. It does not make an observer's claim objectively true, and floating-point measurements are not exact mathematical proofs.

`DERIVED` means a result calculated from recorded inputs, an aggregate over retained outputs, a mathematical formula with a stated derivation, or a disagreement mechanically derived from explicit interest reports. Geometry run records are observed executions; their summary and comparison outputs are labeled `DERIVED` individually.

`PROPOSED` means an untested claim or suggested next experiment. The `note` command requires a deliberate status and does not promote proposals into observations. Notes may reference known runs, observer reports, divergence events, and earlier notes.

## Observer events

An annotation includes `id`, `entry_id`, `slot`, `value`, `observer`, `reason`, `method`, `status: OBSERVED`, and a UTC `timestamp`. The observer name and reason are required. Existing values can be replaced or cleared to `null`, but every previous report remains in `observers.jsonl`. A later report for a slot is the effective value; the archive does not pretend to maintain a consensus across many human or AI observers.

The supplied observer name is attribution, not authenticated identity. No model is called, and writing code with AI does not fill `ai_interest`.

## Divergence events

The derivation rule is `explicit_nonnull_boolean_disagreement:v1`:

```text
both human_interest and ai_interest are booleans
AND human_interest != ai_interest
```

Only explicit observer events can activate the rule. A derived event references the run, exact current and previous values, causal annotation, latest interest annotation IDs, timestamp, and derivation rule. `opened` establishes disagreement. `resolved` records agreement or a slot being cleared to unknown. `changed` is supported for a change in an existing disagreement. Identical effective values do not produce duplicate divergence events, although the repeated observer report itself is kept.

Resolution entries remain in the historical log so an old disagreement is not erased. This log is independent of `divergence_note`; the latter remains blank unless someone supplies text. It is also independent of a shape's numerical catalog matches and misses.

## Manifest and durability

The manifest includes schema/generator versions, source, complete configuration, expected and recorded run/output counts, `state`, `completed`, and failure information. Counts include every top-level output entry, including derived summaries, not merely samples. Each entire run is serialized, appended, flushed, and fsynced before proceeding. Memory use includes the largest run and its full outputs; there is no silent truncation or interest-based retention policy.

Dataset creation refuses every existing output path. Cooperating commands use `.write.lock` created exclusively. There is no automatic lock stealing. After an abrupt process kill, inspect the process and archive before manually removing a stale lock. A caught exception attempts to record an observed execution error and leaves `completed: false`. A hard kill may leave the `started` manifest and lose the unfinished in-memory run; such a dataset is not declared complete.

Observer and divergence appends are individually durable but not one cross-file transaction. A crash between them can leave a missing derived divergence event. Validation detects inconsistent replay, and subsequent annotations refuse to proceed rather than quietly rewriting history. Archive editing outside the CLI can also invalidate the logs.

Validation for supported generator sources deterministically replays records and checks coverage. It detects malformed/truncated JSON, undeclared fields, unexpected base observer values, duplicate IDs, invalid status labels, mismatched counts, bad references, and divergence inconsistencies. It establishes implementation consistency, not scientific validity or optimality.
