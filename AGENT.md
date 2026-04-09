# AGENT Guide for This MFC Repository

This file is a working guide for coding agents operating in this repository root.
It is tailored to this MFC checkout and summarizes the expectations from the MFC contributing guide, with extra emphasis on the layout and risks that matter in this folder.

## Scope

This guide applies to `/home/user/Documents/GitHub/MFC-JRChreim` and everything under it.

## Repo Overview

MFC is organized as a three-phase pipeline:

- `pre_process`: reads case parameters, generates initial conditions, and writes the initial state.
- `simulation`: advances the solution in time. This is the GPU-accelerated phase.
- `post_process`: computes derived quantities and visualization output.

Shared code lives in `src/common/`. Changes there affect all three executables and should be treated as high-blast-radius edits.

Important top-level directories in this folder:

- `src/common/`: shared Fortran modules, derived types, global parameters, MPI, helpers, and shared physics logic.
- `src/pre_process/`: initial-condition generation and pre-processing logic.
- `src/simulation/`: main solver, reconstruction, Riemann solves, source terms, GPU-capable code paths.
- `src/post_process/`: derived output and visualization preparation.
- `toolchain/`: Python CLI, parameter definitions, validation, and build orchestration.
- `tests/`: regression coverage.
- `examples/`: runnable sample cases.

## Default Workflow

When making changes in this repository:

1. Read the relevant local modules first and identify whether the edit touches `src/common/`, a single phase, or both Fortran and Python toolchain code.
2. Keep edits narrowly scoped and consistent with surrounding style.
3. Prefer the project entrypoints over ad hoc commands:
   - `./mfc.sh build -j $(nproc)`
   - `./mfc.sh test -j $(nproc)`
   - `./mfc.sh format`
   - `./mfc.sh lint`
4. If a change affects `src/common/`, assume all three executables may be impacted.
5. If a change adds or modifies a simulation parameter, update both the Python toolchain and the relevant Fortran global-parameter modules.

## Fortran and Fypp Rules

Follow these repository conventions:

- Use modern Fortran style with `implicit none`.
- Give every dummy argument an explicit `intent`.
- Use lowercase Fortran keywords and intrinsics.
- Keep indentation at 2 spaces and align continuation lines under `&`.
- Module names should follow `m_<feature>`.
- Public subroutines should follow `s_<verb>_<noun>`.
- Public functions should follow `f_<verb>_<noun>`.
- Do not use `goto`, `COMMON`, or `save`-style global state unless already required by existing design.
- Use `s_mpi_abort(...)` for fatal errors. Do not use `stop` or `error stop`.
- Do not introduce raw OpenACC or OpenMP pragmas. Use the repository's Fypp GPU macros instead.

## GPU and Performance Notes

Only `simulation` and compatible shared dependencies are GPU-accelerated.

- Keep Fypp GPU macros simple and readable.
- Ensure macros expand correctly for both CPU and GPU builds.
- Avoid unnecessary host/device transfers in hot paths.
- When changing loop structure in GPU-capable files, think about `collapse(...)`, `private(...)`, and device coherence, but use the project's macro system rather than direct pragmas.

## Precision and Safety

MFC uses both storage and working precision concepts:

- Be deliberate about `stp` vs `wp`.
- Use generic intrinsics like `log`, `exp`, `sqrt`, `abs` with the appropriate kind.
- Do not introduce `real(8)`, `real(4)`, `dlog`, `dexp`, `dsqrt`, `dabs`, or `dble`.
- Check MPI type consistency when editing communication or I/O paths.

## Common Pitfalls in This Codebase

- Array bounds are often non-unity and may include ghost-cell ranges. Verify loop bounds against declarations.
- `src/common/` edits have wide blast radius.
- Public signature changes must be propagated carefully and called out clearly.
- Volume fractions, densities, and model-specific thermodynamic paths are sensitive to branch errors.
- If you allocate through the project's allocation macros or patterns, ensure matching deallocation under the same conditions.
- Compiler portability matters. Code should remain compatible with GNU, NVIDIA, Cray, and Intel toolchains.

## When Editing `src/common/`

Use extra care in shared modules such as `src/common/m_phase_change.fpp`.

- Avoid unnecessary interface churn.
- Preserve existing naming and error-handling patterns unless the surrounding code is already being refactored.
- Keep helper state local to the module instead of introducing new global state.
- Recheck all call sites in `pre_process`, `simulation`, and `post_process` if a public routine changes.
- Prefer small, reviewable refactors over sweeping rewrites.

## When Editing Parameters or Validation

If a new parameter or feature is added:

- Register it in `toolchain/mfc/params/definitions.py`.
- Add constraints and dependencies there if needed.
- Add cross-parameter or physics validation in `toolchain/mfc/case_validator.py` when appropriate.
- Declare the parameter in the relevant `m_global_parameters.fpp` file(s).
- Add GPU declaration support if the parameter is used in device code.

## Testing Expectations

Choose verification proportional to the change:

- Formatting-only or comment-only changes: run the minimal relevant check if possible.
- Localized Fortran changes: prefer at least formatting/linting and a relevant build or preprocess path.
- `src/common/` or parameter changes: prefer broader testing because they affect multiple executables.
- If the environment blocks a full build or test run, state what was attempted and where it failed.

## Commit and Review Hygiene

- Prefer small, atomic commits.
- Use concise imperative commit summaries.
- Bundle related fixes together to avoid unnecessary CI noise.
- Document behavioral changes, signature changes, and new parameters clearly.

## Agent-Specific Guidance

- Start by understanding the local module and nearby call sites before editing.
- Match the existing style of the file you are changing.
- Do not replace repository workflows with custom one-off build logic unless necessary for debugging.
- Favor minimal diffs that are easy for MFC maintainers to review.
- If instructions conflict, prefer local repository conventions and the MFC contributing guidance over generic habits.
