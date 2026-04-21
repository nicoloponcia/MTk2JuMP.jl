---
name: julia-release-migration
description: 'Generate release changelog and downstream upgrade notes from Julia API diffs, including renamed symbols, signature changes, and data structure migrations.'
argument-hint: 'Base commit/tag, head commit/tag, and release version'
user-invocable: true
disable-model-invocation: false
---

# Julia Release Migration

## What This Skill Produces
- API diff summary across a release range.
- Release changelog using [CHANGELOG_TEMPLATE.md](./assets/CHANGELOG_TEMPLATE.md).
- Downstream consumer migration notes using [DOWNSTREAM_UPGRADE_NOTES_TEMPLATE.md](./assets/DOWNSTREAM_UPGRADE_NOTES_TEMPLATE.md).
- Updated interface migration record in [INTERFACE_MIGRATION_SUMMARY.txt](../julia-performance-refactor/INTERFACE_MIGRATION_SUMMARY.txt).

## When To Use
- Before tagging a release.
- After refactors that changed exports, names, signatures, or types.
- When dependent projects need concrete upgrade instructions.

## Inputs
- Release diff range (`base..head` or previous tag to current tag).
- Semver target and release date.
- Risk tolerance for breaking changes.

## Workflow
1. Collect API Diff Signals
- Compare exported names across the release range.
- Identify added/removed/renamed functions, changed signatures, and type layout changes.
- Classify each change as breaking, behavior-changing, or additive.

2. Validate Mathematical Robustness Notes
- Link each behavior-sensitive change to invariant checks or test evidence.
- Record tolerance changes explicitly for numerical paths.

3. Generate Changelog
- Fill [CHANGELOG_TEMPLATE.md](./assets/CHANGELOG_TEMPLATE.md).
- Include sections: Added, Changed, Deprecated, Removed, Fixed, Performance.
- Highlight breaking changes with migration pointers.

4. Generate Downstream Upgrade Notes
- Fill [DOWNSTREAM_UPGRADE_NOTES_TEMPLATE.md](./assets/DOWNSTREAM_UPGRADE_NOTES_TEMPLATE.md).
- Provide old-to-new call examples and type migration recipes.
- Include temporary compatibility shims and removal timeline when present.

5. Sync Interface Migration Record
- Update [INTERFACE_MIGRATION_SUMMARY.txt](../julia-performance-refactor/INTERFACE_MIGRATION_SUMMARY.txt).
- Ensure all public API and data layout changes are reflected.

6. Final Gate
- Ensure every breaking change has a concrete migration instruction.
- Ensure all examples are executable and aligned with current exports.

## Quality Gates
- Changelog and upgrade notes both generated for same commit range.
- All breaking API changes listed with explicit migration steps.
- Numerical behavior changes documented with robustness notes.
- Interface migration record synchronized.

## Example Prompts
- Build release migration docs from v0.4.0..HEAD and generate changelog plus downstream notes.
- Compare previous tag to current commit and produce upgrade notes for renamed solver APIs.
- Create release notes with API diff mapping and compatibility shim deprecation timeline.
