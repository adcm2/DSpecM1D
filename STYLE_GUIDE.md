# DSpecM1D C++ Style Guide (Eigen-Aligned)

This guide standardizes naming and structure to follow Eigen-style conventions
as closely as practical for this codebase.

## 1) Naming Rules

### Types and templates
- `Class/struct/enum` names: `PascalCase`
  - Example: `SEM`, `MeshModel`, `EarthquakeCMT`
- Type aliases (`using` / `typedef`): `PascalCase`
  - Example: `Scalar`, `RealScalar`, `MatrixType`, `Complex`
- Template type parameters: `PascalCase` with optional trailing `_` in class
  templates when needed for disambiguation.
  - Example: `template <typename MatrixType_>`

### Functions and methods
- Public and private functions: `lowerCamelCase`
  - Example: `sourceElement`, `receiverElements`, `calculateForce`
- Acronyms inside names use normal camel casing.
  - Prefer `ltgS` over `LtG_S` for new code.

### Variables
- Member variables: `m_` prefix + `lowerCamelCase`
  - Example: `m_mesh`, `m_lengthNorm`, `m_hasFluid`
- Local variables and parameters: `lowerCamelCase`
  - Example: `recElems`, `idxSource`, `maxIters`
- Loop indices may stay short (`i`, `j`, `k`, `idx`) where clear.

### Constants and macros
- `constexpr` and constant variables: `kPascalCase`
  - Example: `kPi`, `kTwoPi`
- Preprocessor macros and header guards: `UPPER_SNAKE_CASE` with project prefix.
  - Example: `DSPECM1D_FULL_SPEC_SINGLE_SEM_H`

### Namespaces
- External/public namespaces may keep current names for API stability.
- Internal implementation details should live under `internal`.
  - Example: `namespace SPARSESPEC { namespace internal { ... } }`

## 2) File and Header Conventions

- Class-centric headers: `PascalCase.h`
  - Example: `SEM.h`, `MeshModel.h`
- Module headers: choose one style and keep it consistent.
  - Recommended for this repo: keep existing `PascalCase` trend for core headers.
- Include paths are case-sensitive and must match filesystem names exactly.
- Keep one primary class per header when practical.

## 3) Code Organization Conventions

- Prefer `using` aliases near top of scope, with readable names.
- Keep public API names stable; use wrappers during transitions.
- Prefer `constexpr` over macros for constants.
- Keep comments focused on intent, not mechanics.

## 4) Migration Mapping for Current Codebase

Use these renames as the canonical pattern:

- `Receiver_Elements` -> `receiverElements`
- `Source_Element` -> `sourceElement`
- `CalculateForce_RED_Coefficients` -> `calculateForceRedCoefficients`
- `LtG_S` -> `ltgS`
- `_mesh` -> `m_mesh`
- `_length_norm` -> `m_lengthNorm`
- `_NQ` -> `m_nq`
- `MATRIX` (local alias) -> `MatrixC`
- `SMATRIX` (local alias) -> `SparseMatrixC`

## 5) Rollout Plan (Low-Risk)

1. Apply this style to all new code immediately.
2. Refactor one subsystem at a time (`SEM/*` first, then `FullSpec*`).
3. For public API renames, add temporary forwarding wrappers.
4. Remove wrappers only after benchmarks/tests and downstream code are updated.

## 6) Formatting

- `.clang-format` remains the source of formatting truth.
- Run formatter before commit; do not mix formatting-only edits with behavior
  changes unless requested.
