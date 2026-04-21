---
name: julia-correctness-testing
description: 'Write and design Julia tests that assess mathematical correctness and coding correctness, including invariants, edge cases, regressions, and integration behavior across modules.'
argument-hint: 'Target module/function, expected math properties, and test depth (smoke, standard, strict)'
user-invocable: true
disable-model-invocation: false
---

# Julia Correctness Testing

## What This Skill Produces
- New and updated tests under `test/` for changed code.
- Mathematical correctness tests (invariants, tolerance-based equivalence, edge cases).
- Coding correctness tests (types, errors, API contracts, regressions, and integration paths).
- A completed test planning artifact in [TEST_STRATEGY_TEMPLATE.md](./assets/TEST_STRATEGY_TEMPLATE.md).

## When To Use
- After refactors, optimizations, or interface changes.
- When adding features or fixing defects.
- Before release validation.

## Inputs
- Target scope: file/module/function(s).
- Mathematical properties that must hold.
- Confidence level: smoke, standard, strict.
- Numerical tolerance policy where applicable.

## Workflow
1. Understand Correctness Requirements
- Extract explicit contracts from code, docs, and existing tests.
- Separate math properties from software-behavior contracts.

2. Build Test Matrix
- Use [TEST_STRATEGY_TEMPLATE.md](./assets/TEST_STRATEGY_TEMPLATE.md).
- Cover both categories:
  - Mathematical correctness: invariants, equivalence, monotonicity, conservation, bounds.
  - Coding correctness: return shapes/types, error conditions, argument validation, stable public behavior.

3. Write Tests
- Place tests in `test/` with module-aligned organization.
- Add deterministic fixtures and representative datasets.
- Include edge and corner cases, not only happy paths.

4. Add New Test Ideas Beyond Existing Coverage
- Add at least one new adversarial/edge test per changed hot path.
- Add at least one integration test for changed module seams.
- Add regression tests for each fixed bug.

5. Execute and Harden
- Run `test/runtests.jl`.
- Remove flaky behavior and stabilize randomness with seeded inputs.
- Ensure tests fail meaningfully when correctness breaks.

6. Finalize
- Summarize new tests and gaps.
- Record residual risks and deferred tests.

## Strictness Levels
- Smoke: core path tests + critical invariant checks.
- Standard: smoke plus edge/error cases and key integration seams.
- Strict: standard plus stress/property-style coverage and broader integration matrix.

## Quality Gates
- Mathematical properties tested with explicit tolerances where needed.
- Coding correctness tested for API contracts and error behavior.
- At least one new edge/adversarial test added for each changed area.
- Integration seams touched by changes have explicit tests.

## Example Prompts
- Write strict tests for src/IF/build.jl covering mathematical invariants and input validation errors.
- Add new integration tests for IF and MTk interfaces and include regression tests for recent fixes.
- Expand interpolation tests with numerical tolerance checks and adversarial edge inputs.
