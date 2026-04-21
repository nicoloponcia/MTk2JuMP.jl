---
name: julia-performance-refactor
description: 'Refactor Julia code for type stability, runtime-first optimization, API/interface redesign, data structure refactoring, and aggressive code debloating while preserving mathematical robustness.'
argument-hint: 'Target module/function, runtime goal, and constraints (mathematical robustness, migration impact, readability)'
user-invocable: true
disable-model-invocation: false
---

# Julia Performance Refactor

## What This Skill Produces
- A prioritized refactor plan for target code paths.
- Concrete code changes that improve type stability and runtime performance.
- API/interface, naming, and data structure changes where they reduce runtime cost.
- Aggressively debloated implementations that preserve mathematical behavior while reducing complexity.
- Validation artifacts: correctness checks, timing/allocation comparisons, and rollback-safe notes.
- A complete downstream migration summary in [INTERFACE_MIGRATION_SUMMARY.txt](./INTERFACE_MIGRATION_SUMMARY.txt).

## When To Use
- You need code refactoring plus measurable speedups.
- `@code_warntype` reveals unstable return or local variable types.
- Profiling shows allocation-heavy hot loops.
- Data models are over-abstracted or container choices are slowing execution.
- The codebase has duplicated logic, dead branches, or unnecessary wrappers.

## Inputs
- Target scope: file/module/function or call chain.
- Performance objective: runtime throughput/latency (primary), plus supporting metrics.
- Constraints: mathematical robustness, readability, dependency policy.

## Workflow
1. Baseline and Scope
- Identify the hot path and define a measurable target.
- Capture baseline with representative inputs using `BenchmarkTools` (`@btime`, `@benchmark`).
- Record type stability evidence (`@code_warntype`) for suspected call sites.

2. Diagnose Root Causes
- Classify issues by type:
- Type instability (abstract containers, `Any`, union growth, non-concrete fields).
- Allocation pressure (temporary arrays, slicing, closure capture, boxing).
- Dispatch overhead (runtime polymorphism in tight loops).
- Structural bloat (duplicate transforms, unnecessary indirection, stale paths).

3. Select Refactor Strategy (Decision Points)
- If type instability dominates: introduce concrete field/container types, function barriers, and stable return contracts.
- If allocations dominate: preallocate buffers, use in-place operations (`!`), avoid temporary views/copies unless needed.
- If dispatch dominates: move dynamic behavior outside loops; specialize kernels with parametric types.
- If data layout dominates: replace generic nested structures with cache-friendly, concrete, domain-shaped structs.
- If complexity dominates: remove dead code and merge duplicated logic; when comments indicate "future enhancement" or "other applications", preserve the path behind a clear extension hook or isolate it into a compatibility module instead of silent deletion.

4. Implement Incrementally
- Apply one coherent change set at a time.
- API, function signatures, naming, and data structures may change when they provide meaningful runtime gains.
- Add concise comments only where non-obvious invariants are introduced.

5. Verify Correctness
- Run existing tests.
- Add focused tests for changed behavior, edge cases, and invariants.
- Confirm numerical equivalence/tolerance for floating-point code.
- Validate mathematical robustness with invariant checks and representative domain cases.

6. Verify Performance
- Re-run benchmarks under same conditions.
- Compare runtime, memory, and allocations against baseline.
- Confirm no regression in neighboring call paths.

7. Finalize
- Summarize what changed and why each change helps.
- Note residual risks and follow-up opportunities.
- Provide before/after metrics with representative inputs.
- Produce and update [INTERFACE_MIGRATION_SUMMARY.txt](./INTERFACE_MIGRATION_SUMMARY.txt) with old/new APIs, renamed symbols, data layout changes, and concrete call-site migration examples.

## Quality Gates (Completion Checks)
- Functional parity: all relevant tests pass.
- Type stability: target functions show stable inferred types where expected.
- Performance gain: measurable improvement against baseline on agreed metrics.
- Simplicity gain: fewer redundant abstractions/branches without hidden coupling.
- Mathematical robustness: invariants and domain-specific numerical behavior remain valid.
- Migration safety: downstream projects have a complete, actionable migration guide.

## Default Checklist
- Baseline captured and documented.
- Hot path narrowed to concrete functions.
- Root-cause category identified.
- Refactor strategy selected with explicit rationale.
- Correctness tests passed.
- Benchmarks re-run and compared.
- Interface migration summary updated for downstream callers.
- Summary + residual risks documented.

## Example Prompts
- Refactor `src/MTk/Interpolations.jl` for type stability and allocation reduction in the interpolation hot path.
- Optimize data structures in `src/IF/LGPoly.jl` to reduce dynamic dispatch and improve cache locality.
- Debloat `src/IF/build.jl` by removing redundant branches and simplifying construction logic while preserving behavior.
- Improve performance in `src/MTk.jl` with function barriers and in-place kernels, then benchmark before/after.
