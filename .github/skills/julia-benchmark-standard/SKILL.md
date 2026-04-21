---
name: julia-benchmark-standard
description: 'Standardize Julia BenchmarkTools setup, warmup policy, and before/after performance reporting for repeatable runtime comparisons in this package.'
argument-hint: 'Target function/module, benchmark inputs, and acceptance threshold (runtime, allocations, regressions)'
user-invocable: true
disable-model-invocation: false
---

# Julia Benchmark Standard

## What This Skill Produces
- A reproducible BenchmarkTools benchmark setup.
- A fixed warmup policy used consistently across runs.
- A before/after report in [BENCHMARK_REPORT_TEMPLATE.md](./assets/BENCHMARK_REPORT_TEMPLATE.md).
- A concise pass/fail decision against defined performance goals.

## When To Use
- Before and after performance refactors.
- When comparing algorithm or data-structure alternatives.
- When validating runtime regressions in CI or release preparation.

## Inputs
- Target function(s) and realistic input generator.
- Primary KPI: runtime (mandatory).
- Secondary KPIs: allocations, memory, compile-time impact (optional).
- Constraints: acceptable regression threshold and numerical equivalence tolerance.

## Standard Setup
1. Environment hygiene
- Use a stable Julia session and avoid unrelated background load.
- Keep package versions fixed for baseline and candidate runs.
- If possible, pin CPU governor/power profile and note platform details.

2. BenchmarkTools configuration
- Use `BenchmarkTools` with explicit interpolation (`$`) for external values.
- Benchmark function kernels, not setup logic; move setup into `setup=` or pre-step.
- Use identical inputs for baseline and candidate unless sensitivity testing is explicit.

3. Warmup policy (mandatory)
- Run at least one non-measured warmup execution per benchmark case.
- Run one BenchmarkTools trial as calibration warmup before recorded trial.
- Only compare recorded trials gathered after warmup.

4. Measurement policy
- Capture baseline report first.
- Apply code change.
- Capture after report under same conditions.
- Report median runtime and allocation metrics as primary comparison.

5. Correctness guard
- Validate output equivalence before interpreting speedups.
- For floating-point outputs, define tolerance and verify acceptance bounds.

6. Reporting
- Fill [BENCHMARK_REPORT_TEMPLATE.md](./assets/BENCHMARK_REPORT_TEMPLATE.md).
- Include absolute numbers and percentage deltas.
- Include clear decision: improved, neutral, or regressed.

## Decision Rules
- Runtime improves and correctness holds: accept.
- Runtime neutral but complexity/readability significantly better: case-by-case accept.
- Runtime regresses beyond threshold without strong justification: reject or iterate.
- Any correctness failure: reject until fixed.

## Quality Gates
- Same benchmark harness used before and after.
- Warmup policy followed and documented.
- Correctness checked prior to final interpretation.
- Before/after report completed with deltas and decision.

## Example Prompts
- Benchmark src/IF/build.jl constructors with standard warmup policy and produce before/after report.
- Compare interpolation kernels in src/MTk/Interpolations.jl and gate on runtime median delta.
- Run standardized benchmark for LG polynomial evaluation and include allocation trend.
