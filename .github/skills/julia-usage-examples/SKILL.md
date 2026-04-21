---
name: julia-usage-examples
description: 'Write practical Julia usage examples for package functionalities and place them in an examples folder separate from src, with runnable and documented workflows.'
argument-hint: 'Target functionality, intended audience, and example depth (quickstart, standard, advanced)'
user-invocable: true
disable-model-invocation: false
---

# Julia Usage Examples

## What This Skill Produces
- New runnable examples in `examples/` (never in `src/`).
- Task-oriented example scripts and/or markdown walkthroughs.
- API usage coverage for primary workflows and common edge usage patterns.
- A completed authoring checklist in [EXAMPLE_AUTHORING_TEMPLATE.md](./assets/EXAMPLE_AUTHORING_TEMPLATE.md).

## When To Use
- After adding or changing public functionalities.
- Before release to improve downstream adoption.
- When onboarding users to core package workflows.

## Folder Policy
- Put all usage examples in `examples/`.
- Keep implementation code in `src/` and tests in `test/`.
- If examples do not exist yet, create `examples/`.

## Inputs
- Target functionality or module.
- Audience level: quickstart, standard, advanced.
- Required assumptions: dependencies, data prep, and expected output behavior.

## Workflow
1. Select Representative Workflows
- Choose examples that map to the most common real usage.
- Include one minimal quickstart and one richer practical case.

2. Author Example Files in `examples/`
- Use clear filenames (`examples/quickstart_<feature>.jl`, `examples/advanced_<feature>.jl`).
- Keep examples executable and deterministic.
- Prefer concise scripts with explicit imports and setup.

3. Explain How To Run
- Include run commands and expected output characteristics.
- Note input assumptions and numerical tolerance notes when relevant.

4. Validate Against Current API
- Ensure examples match current exported names and signatures.
- Update examples whenever interface changes occur.

5. Coverage Check
- Confirm major public functionalities have at least one example path.
- Add edge usage examples for common user pitfalls.

6. Finalize
- Complete [EXAMPLE_AUTHORING_TEMPLATE.md](./assets/EXAMPLE_AUTHORING_TEMPLATE.md).
- Summarize what examples were added and what remains uncovered.

## Quality Gates
- All examples are in `examples/` and not in `src/`.
- Examples are runnable and aligned with current package API.
- Quickstart and practical workflows are both covered.
- Known assumptions and expected behavior are documented.

## Example Prompts
- Create quickstart and advanced examples for IF build pipeline and place them in examples/.
- Add runnable examples for MTk interpolation APIs, including one edge-case workflow.
- Refresh all existing examples after API renaming and verify they still run.
