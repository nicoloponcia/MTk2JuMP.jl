---
name: julia-template-interface-sync
description: 'Update templates/Config.jl and templates/ScalesBounds.jl to track interface and functionality changes in src, including renamed symbols, new options, and removed behaviors.'
argument-hint: 'Change range or target module, plus whether to apply strict backward-compat notes'
user-invocable: true
disable-model-invocation: false
---

# Julia Template Interface Sync

## What This Skill Produces
- Updated template interface files:
  - [templates/Config.jl](../../../templates/Config.jl)
  - [templates/ScalesBounds.jl](../../../templates/ScalesBounds.jl)
- A synchronized change record in [TEMPLATE_SYNC_REPORT_TEMPLATE.md](./assets/TEMPLATE_SYNC_REPORT_TEMPLATE.md).
- Explicit mapping of old-to-new template fields, names, and callable options.

## When To Use
- After changing public interface behavior in `src/`.
- After adding/removing configuration options or scale/bounds behavior.
- After function renames or signature changes that affect template usage.

## Inputs
- Target diff range (`base..head`) or changed files list.
- Affected modules/functions.
- Migration strictness: soft (compat shims) or hard (breaking update only).

## Scope
- Always evaluate and update both template files together:
  - [templates/Config.jl](../../../templates/Config.jl)
  - [templates/ScalesBounds.jl](../../../templates/ScalesBounds.jl)
- Keep template examples aligned with current exported APIs and data expectations.

## Workflow
1. Detect Interface/Functionality Deltas
- Inspect changes in `src/` affecting configuration schema, bounds logic, scales logic, and construction workflow.
- Extract added/removed/renamed symbols and changed defaults.

2. Build Mapping Table
- Create old-to-new mapping for:
  - keys/fields
  - function calls
  - default values
  - expected data shape/type
- Identify behavior changes that require comments for downstream users.

3. Update Template Files
- Update [templates/Config.jl](../../../templates/Config.jl) for config schema and interface calls.
- Update [templates/ScalesBounds.jl](../../../templates/ScalesBounds.jl) for scale/bounds construction and semantics.
- Remove stale code paths unless marked as future extension points.

4. Validate Template Coherence
- Ensure template code references existing APIs only.
- Ensure options/fields in templates correspond to current runtime behavior.
- Keep deterministic defaults and clear inline guidance.

5. Record Migration Notes
- Fill [TEMPLATE_SYNC_REPORT_TEMPLATE.md](./assets/TEMPLATE_SYNC_REPORT_TEMPLATE.md).
- Document breaking vs non-breaking template changes.
- Provide downstream update snippets for projects that copy/use these templates.

6. Final Gate
- Both template files updated in same pass.
- No stale symbol references remain.
- Report includes migration mapping and verification checklist.

## Quality Gates
- `templates/Config.jl` reflects current interface/config behavior.
- `templates/ScalesBounds.jl` reflects current scales/bounds behavior.
- Old-to-new mapping documented for all renamed/removed items.
- Downstream usage notes included when template behavior changed.

## Example Prompts
- Sync templates with changes in src/IF/build.jl and src/IF/LGPoly.jl, then generate template migration notes.
- Update templates/Config.jl and templates/ScalesBounds.jl after interface renames and changed defaults.
- Track recent functionality changes and refresh both templates with a breaking-change summary.
