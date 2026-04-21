---
name: julia-global-correctness-integration
description: 'Ensure global correctness and integration across all package modules by validating interfaces, invariants, cross-module behavior, and end-to-end workflows.'
argument-hint: 'Target change set or release range and required confidence level (smoke, standard, strict)'
user-invocable: true
disable-model-invocation: false
---

# Julia Global Correctness Integration

## What This Skill Produces
- A package-wide integration validation plan.
- Cross-module interface consistency checks.
- Invariant and numerical robustness validation matrix.
- A consolidated integration report using [INTEGRATION_VALIDATION_TEMPLATE.md](./assets/INTEGRATION_VALIDATION_TEMPLATE.md).

## When To Use
- Before release or major merge.
- After API/data-structure refactors.
- After changes touching multiple modules in src.
- When failures appear only in composed end-to-end flows.

## Inputs
- Target range: commit(s) or release candidate.
- Confidence level: smoke, standard, strict.
- Critical workflows that must remain correct.

## Validation Scope
- Module entrypoints and exports.
- Cross-module call contracts.
- Data flow consistency (type and shape expectations).
- End-to-end package workflows.
- Numerical invariants and tolerance-sensitive outputs.

## Workflow
1. Build Integration Map
- Enumerate modules under src and their primary interfaces.
- Record upstream/downstream dependencies between modules.
- Identify high-risk integration seams changed by the diff.

2. Contract Consistency Checks
- Verify function names, signatures, and return assumptions match across module boundaries.
- Verify struct field expectations remain coherent for producer/consumer modules.
- Verify changed defaults do not silently alter downstream behavior.

3. Correctness and Invariant Checks
- Run unit and integration tests (`test/runtests.jl`).
- Add focused integration tests for changed seams when missing.
- Validate mathematical invariants and tolerance-bound equivalence on representative cases.

4. End-to-End Scenarios
- Execute at least one end-to-end path per critical workflow.
- Compare outputs to baseline/known-good references.
- Record any acceptable drift with rationale.

5. Strictness Levels
- Smoke: key module load + primary end-to-end path.
- Standard: full test suite + targeted seam tests + invariant checks.
- Strict: standard plus stress/edge inputs and repeated runs for nondeterminism checks.

6. Report and Gate
- Fill [INTEGRATION_VALIDATION_TEMPLATE.md](./assets/INTEGRATION_VALIDATION_TEMPLATE.md).
- List blockers, risks, and required follow-up actions.
- Mark release readiness as pass/conditional/fail.

## Quality Gates
- All critical module seams validated.
- All required tests for chosen confidence level executed.
- Numerical robustness documented with explicit tolerances.
- Integration report completed with pass/fail rationale.

## Example Prompts
- Run strict global correctness and integration validation for the current branch and produce a release readiness report.
- Validate cross-module contracts after data-structure refactor in IF and MTk modules.
- Generate a standard integration report focused on module interface changes since last tag.
