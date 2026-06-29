# ClimaAnalysis.jl — Typo & Bug Audit

## Context

The user asked for a **comprehensive deep dive** to find *all* typos and bugs across
`src`, `docs`, and `tests` (plus the `ext/` extension sources, which are also source code).
Two seed examples set the bar: a copy-paste bug (`push_preview = push_preview = all(...)`
in `docs/make.jl`) and a broken/no-op test (`@test_logs (:warn, "...")` with no expression
after it in `test/test_Var.jl:658`). So "typo" here spans spelling/grammar/markdown **and**
"bug" spans logic errors, copy-paste mistakes, dead code, and no-op tests.

This file is the deliverable: a categorized list of every finding, by **location**
(src / ext / docs / test) and **label** (bug / doc-error / typo / grammar / style),
followed by a fix plan and verification steps.

## Methodology (9 kinds of checks, cross-verified)

1. LLM per-file scan over **56 file-chunks** (whole repo) + an adversarial verification pass
   that re-read each candidate against the code (8 candidates rejected as false positives).
2. Direct code-reading of **all 43 bug findings** (every `src` bug confirmed by me; tests/docs spot-checked).
3. Duplicate-assignment regex (PCRE backref) → found `docs/make.jl:48`.
4. Duplicate-word regex → found `src/Var.jl:2163` ("the the").
5. `@test_logs` no-op detection → confirmed `test_Var.jl:658` is broken, `:994` is correct.
6. `aspell` dictionary spell-check over docs prose + source comments/docstrings (confirmed the
   single-word misspellings the LLM scan also caught — high overlap = thorough).
7. Bug-pattern greps: `@test x = y`, `== nothing`, repeated kwargs, anon-fn arrow typo `(x) - 1000x`.
8. Markdown link/anchor/image validation → confirmed `NEWS.md:158` bad anchor; image at
   `visualize_rmse_var.md:106` is generated at build time (false positive).
9. TODO/FIXME scan (informational only — all intentional).

## Summary

**Pass 1** findings by location/label:

| Location | bug | doc-error | typo | grammar | style | total |
|---|---|---|---|---|---|---|
| SRC | 12 | 27 | 18 | 53 | 1 | **111** |
| EXT | 0 | 7 | 4 | 9 | 0 | **20** |
| DOCS | 17 | 17 | 15 | 51 | 1 | **101** |
| TEST | 14 | 7 | 10 | 13 | 4 | **48** |
| **Total** | **43** | **58** | **47** | **126** | **6** | **280** |

**Cumulative across all three passes** (Pass 2 detail in the "Pass 2 — additional findings" section;
Pass 3 detail + per-location table in "Pass 3 — novel-lens hunt"):

| Pass | bug | doc-error | typo | grammar | style | total |
|---|---|---|---|---|---|---|
| Pass 1 | 43 | 58 | 47 | 126 | 6 | **280** |
| Pass 2 (≈) | 22 | 34 | 3 | 13 | 8 | **≈80** |
| Pass 3 | 9 | 16 | 0 | 17 | 5 | **47** |
| Pass 4 | 4 | 2 | 4 | 4 | 0 | **14** |
| Pass 5 | 3 | 3 | 0 | 0 | 1 | **7** |
| **Grand total (≈)** | **81** | **113** | **54** | **160** | **20** | **≈428** |

> Pass 2's split is approximate: the workflow contributed 75 (21 bug / 32 doc-error / 1 typo / 13 grammar
> / 8 style) plus ~5 manual sibling/deep-dive finds whose exact per-label count is fuzzy in the source notes.
> Pass 1, Pass 3, and Pass 4 counts are exact.

> Labels: **bug** = logic/correctness; **doc-error** = wrong docstring/comment/example/API name;
> **typo** = misspelling/markdown; **grammar** = subject-verb/article/tense; **style** = dead/unused code.

---

## ⭐ Highest-impact bugs (fix first — real runtime/logic defects)

- **`src/Var.jl:1472`** — `arecompatible` builds `y_dims_conventional_names` by iterating
  `x_dims` instead of `y_dims`, so the dimension-name mismatch check on line 1484 can *never*
  fire. **Fix:** `for name in y_dims`.
- **`src/Leaderboard.jl:351`** — the function arg `keys::AbstractVector` shadows `Base.keys`;
  the error message calls `keys(key2index)` → `MethodError` (vector called as a function) on the
  missing-key path. **Fix:** `Base.keys(...)` or rename the arg.
- **`src/Leaderboard.jl:96`** — `"...($length(units))..."` interpolates the `length` *function*,
  not the count (prints `length(units)`). **Fix:** `$(length(units))` / `$(length(model_names))`.
- **`src/Template.jl:226` & `:233`** — latitude is given `"degrees_east"` and longitude
  `"degrees_north"` — **swapped** (CF convention is lat=north, lon=east; cf. `test_Atmos.jl:224-225`).
- **`src/Utils.jl:375`** — empty-input return `Vector{Vector{eltype(dates)}}[]` is one nesting
  level too deep vs the non-empty path. **Fix:** `Vector{eltype(dates)}[]`.
- **`test/test_Var.jl:658`**, **`:2253-2255`**, **`:2141`**, **`:3581`** — broken/no-op tests
  (missing expression, missing `@test`, wrong variable passed). Details below.
- **`docs/make.jl:48`** — `push_preview = push_preview = all(...)` duplicated assignment.
- **`test/aqua.jl:10`** *(Pass 3)* — `Aqua.detect_ambiguities(ClimaAnalysis; recursive = true)` registers
  **no** `@test` (resolves to the non-asserting `Test.detect_ambiguities`, which only *returns* the list);
  every sibling line uses `Aqua.test_*`. **Method-ambiguity checking is silently disabled in CI.**
  **Fix:** `Aqua.test_ambiguities(ClimaAnalysis; recursive = true)`.
- **Aliasing quartet (missing `deepcopy`)** *(Pass 3+4)* — same class as the known `outvar_operators.jl:49`:
  **`src/outvar_selectors.jl:236/241`** (`window` shallow-`copy`s `dims`/`dim_attributes` → returned var
  shares mutable dim arrays + attr dicts with its input; verified `===`), **`src/Atmos.jl:90-93`/`:99-102`**
  (`to_pressure_coordinates` reuses original `v` for non-altitude dims/dim_attributes), and *(Pass 4)*
  **`ext/ClimaAnalysisUnitfulExt.jl:79`** (`convert_units` passes `var.dims`/`var.dim_attributes` straight into
  the constructor, bypassing `remake`, while its sibling `convert_dim_units`+`remake` both `deepcopy`).
  **Fix:** `deepcopy` the reused values in all four.
- **`src/Leaderboard.jl:84-98` `RMSEVariable` 5-arg constructor** *(Pass 4)* — two bugs: (1) **L90-91**
  `for key in keys(units) … delete!(units, key)` mutates the `Dict` while iterating its live `KeySet`
  (undefined behavior in Julia — can error "dict was modified during iteration" or skip entries); (2) **L84-98**
  the constructor mutates the **caller's** `units` dict in place (adds at L85, deletes at L91) and then stores
  that same reference (L98), violating the `deepcopy` convention (`add_model`/`add_category` use `units |> deepcopy`).
  **Fix:** `units = copy(units)` at the top — resolves both.
- **`test/test_Var.jl:2276-2277`** *(Pass 3)* — second occurrence of the bare-comparison no-op:
  `JJA/DJF.attributes["start_date"] == "2024-1-1"` with no `@test` → wrap in `@test`.
- **`docs/src/var.md:357`** & **`NEWS.md:331`** *(Pass 3)* — broken doc examples: `global_mse(sim, obs)`
  uses undefined `sim`/`obs` (block defines `sim_var`/`obs_var`); `cat(...; dim = "time")` uses singular
  `dim` but the signature is `Base.cat(vars::OutputVar...; dims::String)` → would error.
- **`src/Atmos.jl:118`** *(Pass 5)* — `to_pressure_coordinates` indexes the reduced `CartesianIndex` by the
  **original** axis number, so it `BoundsError`s / reads wrong columns whenever **altitude is not the last
  dimension** (e.g. `("z","lon","lat")` or 4-D `(lon,lat,z,time)`). Masked because all tests put `z` last; refutes
  Pass 4's "no dim-order assumption" note. **Fix:** consume the index by its own running position.
- **`src/Var.jl:839-840`** *(Pass 5)* — `_reduce_over` shallow-`copy`s `var.dims`/`var.dim_attributes`; the
  constructor's `OrderedDict(...)` keeps inner refs, so **every** reduction (`average_*`, `variance_*`,
  `integrate`, `slice`) returns a var aliasing its input's un-reduced dim arrays + attr dicts. Broadest member of
  the aliasing family; refutes Pass 3's "operator families aliasing-safe" note. **Fix:** `deepcopy` both.

---

## SRC

### `src/Atmos.jl`
- **L182** `[typo]` — `# Check if presure is a dimension` → `pressure`
- **L59** `[doc-error]` — `> :important:` pseudo-admonition → use a real Documenter block `!!! warning "Extrapolation"`
- **L144** `[grammar]` — "can chosen by supplying" → "can be chosen by supplying"

### `src/Catalog.jl`
- **L49** `[grammar]` — "pass in a short name or a pair of short name to alias" → reword to "a short-name-to-alias pair"
- **L166** `[grammar]` — `error("Cannot find variable with $short_name in catalog")` → "...with short name $short_name..."

### `src/Numerics.jl`
- **L60** `[bug]` — generic `_integrate_dim` errors `"Cannot integrate when z is a single point"` (hardcoded "z") → generic phrasing
- **L51** `[doc-error]` — "Integrate out a dimension `data`." (missing word) → "...a dimension from `data`."
- **L84** `[doc-error]` — docstring signature lists `int_weights,` with type that doesn't match the code
- **L157** `[doc-error]` — comment "we use the assumption that units are degrees" but generic path uses no `deg2rad` (copy-paste)
- **L13, L112** `[grammar]` — "each point correspond" → "corresponds"; "it make no contribution" → "makes"

### `src/masks.jl`
- **L93** `[bug]` — `error("Threshold ($threshold)")` echoes value only → "...must be between 0 and 1"
- **L273** `[bug]` — `Base.depwarn(..., :make_lonlat_mask!)` but the function is `make_lonlat_mask` (no `!`)
- **L116** `[typo]` — "multipled" → "multiplied"
- **L46, L48, L89, L252** `[grammar]` — "a `OutputVar`" → "an `OutputVar`" (several)
- **L60, L63** `[grammar]` — "Replace zeros/ones **to** this value" → "**with** this value"
- **L165** `[grammar]` — "if you need do any broadcasting" → "need to do"
- **L254** `[grammar]` — "takes in an element...and return" → "and returns"

### `src/Template.jl`
- **L226** `[bug]` — latitude units `"degrees_east"` → `"degrees_north"` (swapped; updates `test_Template.jl`)
- **L233** `[bug]` — longitude units `"degrees_north"` → `"degrees_east"` (swapped; updates `test_Template.jl`)
- **L61** `[typo]` — "Intialize a `TemplateVar`." → "Initialize"
- **L350, L351** `[typo]` — double space "then  `collect`"; missing space "then`collect`"
- **L257, L270** `[grammar]` — "which uses default values" → "which use default values"
- **L182** `[style]` — generated `add_*_dim` and `add_*_dim!` use different default `name`

### `src/Utils.jl`
- **L288** `[bug]` — docstring signature `seasons = (...)))` has an extra `)` (unbalanced)
- **L375** `[bug]` — empty case `Vector{Vector{eltype(dates)}}[]` over-nested → `Vector{eltype(dates)}[]`
- **L237** `[typo]` — `warp_string` → `wrap_string` (**breaking-ish rename**, see Decisions)
- **L411, L414** `[doc-error]` — docstring says "season" where it should say "month"
- **L159, L239, L464** `[grammar]` — wording fixes

### `src/Leaderboard.jl`
- **L96** `[bug]` — `$length(units)` interpolation bug → `$(length(units))` (also "number of unit" → "units")
- **L351** `[bug]` — `keys` arg shadows `Base.keys` → `MethodError` on error path
- **L550** `[bug]` — `rmses = rmse_var.RMSEs |> copy` is dead (overwritten on L552)
- **L653** `[bug]` — `error("Annual does not exist...")` hardcodes "Annual" but category is `category_name` (cf. L632)
- **L484** `[doc-error]` — docstring signature args `rmse_var1/rmse_var2` don't match real `rmse_var_src/_dest`
- **L40** `[typo]` — "Arrray" → "Array"
- **L266** `[typo]` — "Intialize" → "Initialize"
- **L83, L126, L146, L486** `[grammar]` — "has an unit"→"a unit"; "will has"→"will have"; "matches"→"match"

### `src/Var.jl`
- **L1472** `[bug]` — ⭐ `y_dims_conventional_names` iterates `x_dims` → `arecompatible` never detects dim-name mismatch. Use `y_dims`.
- **L1098/1099** `[doc-error]` — variance docstring describes "dividing the sample mean by n" (wrong formula)
- **L1300** `[doc-error]` — `!!! warn "Deprecated"` → `!!! warning` (invalid admonition keyword)
- **L1307** `[doc-error]` — deprecation points to `center_longitude` → should reference `shift_longitude`
- **L1424** `[doc-error]` — "Extrapolation is **now** allowed" → "is **not** allowed" (opposite of behavior)
- **L1991** `[doc-error]` — "(JDF)" → "(DJF)"
- **L2110** `[doc-error]` — docstring header `check_time_dim` → `_check_time_dim`
- **L2734/2736** `[doc-error]` — header `replace!(...)` but function is non-mutating `replace(...)`; also unbalanced backtick
- **L131, L297, L446, L2084, L2184** `[doc-error]` — name/signature/wording mismatches
- **L1361** `[typo]` — comment "represent the same latitude" → "same longitude"
- **L2041, L2138** `[typo]` — "Janauary" → "January"
- **L2163** `[typo]` — "promote the the type" → "the type"
- **L2190** `[typo]` — "Additonally, there is no checks...in the values" → "Additionally, there are no checks...on the values"
- **L2492** `[typo]` — "DateTime.Dates" → "Dates.DateTime"
- **L2698, L2717** `[typo]` — "occurences" → "occurrences"
- **L1572** `[grammar→couples to tests]` — `"...is not consistent..."` → `"are not consistent"` (also update `test_Var.jl:1278/1388/1467`)
- **L1918, L1937** `[grammar→couples to tests]` — `"var does not has longitude/latitude"` → `"does not have"` (update `test_Var.jl:2024/2043`)
- **L124, L135, L172, L223, L296, L391, L1463, L1949, L2068, L2536, L2751, L2768, L2835, L2837** `[grammar]` — agreement/article fixes

### `src/flat.jl`
- **L4, L150, L185, L189, L232** `[grammar]` — "is not the same"→"are not the same" (error msgs); "Representing...contain"→"Represents...contains"

### `src/outvar_dimensions.jl`
- **L83, L159, L274** `[grammar]` — "a `altitude`"→"an"; "dimension"→"dimensions"

### `src/outvar_operators.jl`
- **L135** `[doc-error]` — docstring names args `x`/`y` but signature uses `x_var`/`y_var`
- **L184** `[doc-error]` — comment "Start by copying all attributes from x" contradicts the code (empty set)
- **L4** `[typo]` — "`OutputVars`" → "`OutputVar`s" (move plural outside backticks)

### `src/outvar_selectors.jl`
- **L327** `[doc-error]` — error text "selectdim_as_reduction only supports..." → "slicing only supports..."
- **L410** `[doc-error]` — truncated/garbled comment
- **L49** `[typo]` — "dimesional" → "dimensional"
- **L28, L151, L264** `[grammar]` — "always take"→"takes"; "approximately match"→"matches"

### `src/split_apply_combine.jl`
- **L39** `[doc-error]` — docstring signature has stray trailing ` end`
- **L66, L132, L133, L148** `[grammar]` — "an [`GroupAll`]"→"a"; "cannot applies"→"apply"; "applied on"→"applied to"; "combine"→"combines"

### `src/Sim.jl`
- **L92** `[doc-error]` — example comment mixes `:`/`=>` dict syntax inconsistently

---

## EXT

### `ext/ClimaAnalysisGeoMakieExt.jl`
- **L163** `[doc-error]` — `GeoMakie.coastline` → `GeoMakie.coastlines`
- **L204** `[doc-error]` — documented signature lists `plot_contours = true,` (not a real kwarg)
- **L305** `[typo]` — "gloal" → "global"; and docstring `sim.data - var.data` likely → `sim.data - obs.data`
- **L148** `[grammar]` — "ClimaAnalysis do not support" → "does not support"

### `ext/ClimaAnalysisMakieExt.jl`
- **L219, L411, L528, L533** `[doc-error]` — "Syntactic sugar for `sliced_*`" should reference the `!` (mutating) variants; `line_plot!` → `plot!`
- **L361** `[doc-error]` — example claims a `lat-long` heatmap for a 1D-slice helper
- **L77** `[typo]` — "It we left it here" → "If we left it here"
- **L694** `[typo]` — `not_nan_idices` → `not_nan_indices` (identifier, used L694-696)
- **L76, L78, L270, L803, L805, L840, L942** `[grammar]` — agreement/article fixes

### `ext/ClimaAnalysisUnitfulExt.jl`
- **L14** `[typo]` — "Uniftul" → "Unitful"
- **L22** `[grammar]` — "This function in inherently" → "is inherently"

---

## DOCS

### `docs/make.jl`
- **L48** `[bug]` — `push_preview = push_preview = all(` → single assignment

### `docs/src/var.md`
- **L316, L324, L327** `[bug]` — `global_bias/mse/rmse(sim, obs)` use undefined `sim`/`obs` → `sim_var`/`obs_var`
- **L78** `[doc-error]` — `dim_set_units!` → real function `set_dim_units!`
- **L275** `[doc-error]` — season months "January, March, June, **August**, December" → **September**
- **L383** `[doc-error]` — "`apply_oceanmask` or `apply_oceanmask`" (dup) → "`apply_landmask` or `apply_oceanmask`"
- **L43, L109, L117, L118** `[typo]` — "set up` var`"; "which mean"; "Date.DateTime"→"Dates.DateTime"; "each dates...reference with"
- **L160, L375, L401** `[grammar]` — agreement fixes

### `docs/src/visualize.md`
- **L87, L121** `[bug]` — `plot_bias_on_globe!(fig, var, ...)` uses undefined `var` / one positional arg → `(fig, sim_var, obs_var, ...)`
- **L88** `[bug]` — uses `ca_kwargs(...)` but block never imports `kwargs as ca_kwargs`
- **L97** `[doc-error]` — references deprecated `make_lonlat_mask` → `generate_lonlat_mask`
- **L102** `[doc-error]` — wrong guidance on passing `nan_color`/`true_val`
- **L20** `[typo]` — "suppose `var` it the variable...with a ocean mask" → "is the variable...with an ocean mask"
- **L3, L100** `[grammar]` — punctuation; "do not support"→"does not support"

### `docs/src/flat.md`
- **L86** `[bug]` — unbalanced backtick "on the `OutputVar." → "`OutputVar`."
- **L96** `[bug]` — masking example measures `flat_nan_var.data` → should be `flat_masked_var.data`
- **L66** `[doc-error]` — `flatvar.metadata` → `flat_var.metadata`
- **L85** `[doc-error]` — keyword `mask_var` → `mask`
- **L132** `[typo]` — `Outputvar` → `OutputVar`
- **L50, L55, L89, L101** `[grammar]` — "a `OutputVar`"→"an"; "drop"→"drops"; "This flatten"→"flattens"

### `docs/src/howdoi.md`
- **L74, L75** `[bug]` — `add_dim("lon"/"lat", time, ...)` uses the `time` array → should use `lon`/`lat` (prose says len 15/20)
- **L80** `[doc-error]` — prose dimension lengths don't match the collected ranges
- **L329** `[typo]` — "Feburary" → "February"
- **L113, L158, L169** `[grammar]` — wording; "If this not the case"→"is not"; "others dimensions"→"other"

### `docs/src/index.md`
- **L85** `[doc-error]` — `ts_average = get(...reduction = "max"...)` misnamed → `ts_max`
- **L111** `[doc-error]` — prose says `.attributes` but code/field is `.dim_attributes`
- **L141** `[doc-error]` — `slice(...)` example uses already-sliced var / wrong dim keyword
- **L68, L93, L97** `[typo]` — "` OutputVar`" leading space inside backticks; "a"→"an"
- **L78, L122, L163** `[grammar]` — "stitch"→"stitches"; "These function use"→"functions"; "a `OutputVar`"→"an"

### `docs/src/read_files.md`
- **L54** `[doc-error]` — `rlut_var = get(catalog, "rsut", ...)` short-name doesn't match var name (`"rlut"`)

### `docs/src/rmse_var.md`
- **L103** `[doc-error]` — duplicated sentence (L103/L104); delete L103

### `docs/src/split_apply_combine.md`
- **L6, L18, L25, L73, L75, L78, L83, L84, L85, L94** `[grammar]` — "subjected to changes"→"subject to change"; "a `OutputVar`"→"an"; "does modifies"→"modify"; "have"→"has", etc.

### `docs/src/developer.md`
- **L37** `[doc-error]` — describes wrong keyword arguments for the dim-adding functions
- **L56, L57** `[grammar]` — "determine"→"determines"; "is not used"→"are not used"
- **L39** `[style]` — redundant "there are also functions"

### `docs/src/api.md`
- **L200** `[grammar]` — "For development and not" → "...internal utilities for development and not part of the public API."

### `docs/src/visualize_rmse_var.md`
- **L18** `[grammar]` — "dividing over the median model's RMSEs" → "dividing by..."

### `README.md`
- **L147** `[typo]` — "capital latter" → "capital letter"
- **L51, L78, L82, L95, L154** `[grammar]` — "lot of"→"a lot of"; "consists in"→"consists of"; "In either cases"→"case", etc.

### `NEWS.md`
- **L158** `[bug]` — TOC anchor `#...-flatVar-...` (uppercase) won't resolve (GitHub lowercases) → `flatvar`
- **L539** `[bug]` — `dim_names = ["longitude, "latitude"]` missing closing quote → `["longitude", "latitude"]`
- **L965** `[bug]` — `conversion_function = (x) - 1000x` (subtraction) → `(x) -> 1000x` (arrow)
- **L1109, L1117, L1120** `[bug]` — `global_bias/mse/rmse(sim, obs)` → `(sim_var, obs_var)`
- **L6, L803, L1231** `[doc-error]` — heading vs body term inconsistency; "end"→"start"; `add_units!`→`add_unit!`
- **L8, L280, L1236, L1243** `[typo]` — "implemention"→"implementation"; "In this example."→","; "mode than one"→"more than"; "Comparsion"→"Comparison"
- **L35, L614, L662, L677, L680, L785, L791, L823, L1001, L1081, L1084, L1107, L1246, L1247, L1248, L1268, L1302** `[grammar]` — agreement/article fixes (changelog prose)

---

## TEST

### `test/test_Var.jl`
- **L658** `[bug]` — `@test_logs (:warn, "...")` with **no expression** → no-op; `wrong_var` unused. Append `match_mode = :any ClimaAnalysis.average_lat(wrong_var; weighted = true)`
- **L2253, L2254, L2255** `[bug]` — bare `MAM/JJA/DJF.attributes["start_date"] == "2024-1-1"` with no `@test` → wrap in `@test`
- **L2141** `[bug]` — builds `minute_var` but asserts on `no_time_var` (dup of L2136); use `minute_var`
- **L3581** `[bug]` — builds `lat_var` (unused), passes `mask_var` instead → rename to `mask_var = make_template_var("lat")`
- **L217, L1873, L1999, L2324, L2441, L3806** `[typo]` — "condtions"; "with out"→"without"; "parsaeble"→"parseable"; "Cheak"→"Check" (×2); "dimenisions"
- **L236, L3557** `[grammar]` — "for and lat"→"for lon and lat"; "Permute the OutputVar"→"Permuting"
- **L1278, L1388, L1467, L2024, L2043** `[grammar]` — assert the buggy src messages; update together with `src/Var.jl:1572/1918/1937`

### `test/test_flat.jl`
- **L157** `[bug]` — compares `flat_var.metadata.*` (leftover from L134) instead of `flat_var_mask.*` (L157/158/160)
- **L190-210** `[bug]` — exact duplicate of L126-146 (no added coverage) → remove or repurpose
- **L212-230** `[bug]` — exact duplicate of L170-188 → remove or repurpose
- **L362** `[bug]` — `("lat","lon","time")` duplicated (L361/L362); missing perm `("lon","time","lat")`

### `test/test_outvar_selectors.jl`
- **L139** `[bug]` — `@test_throws ErrorException get_index(...) == 1` — the `== 1` is dead (throws first) → remove
- **L630-634** `[bug]` — identical to L613-617 (both test the missing-dimension path) → make L613 a real invalid-index test
- **L612** `[doc-error]` — "# Invalid indices" block doesn't test invalid indices (uses a missing dim)
- **L252** `[style]` — `fill(3)` (Int) vs Float64 data → `fill(3.0)`

### `test/test_Atmos.jl`
- **L233** `[bug]` — `pfull_attribs = Dict("short_name" => "pfull")` dead (overwritten on L239)
- **L13** `[grammar]` — "Let start" → "Let's start"
- **L17** `[style]` — unused `zvar`/`data`

### `test/test_Catalog.jl`
- **L65** `[bug]` — `@test precip.dims == precip.dims` (self-comparison, always true) → `precip.dims == pr.dims`
- **L96** `[grammar]` — "NCCatalog use the alias" → "uses"

### `test/test_split_apply_combine.jl`
- **L112, L113** `[doc-error]` — `# JJA`/`# MAM` comments are **swapped** vs the indices

### `test/test_Numerics.jl`
- **L35, L45, L55** `[doc-error]` — comments say "(not equispaced)" but the data is equispaced

### `test/test_Utils.jl`
- **L163** `[doc-error]` — comment `# 12/1` should list both dates (`# 12/1, 12/31`)
- **L264-266** `[style]` — duplicated `dict1`/`dict2` block

### `test/test_Sim.jl`
- **L156** `[typo]` — "empty tring" → "string"
- **L245, L273** `[grammar]` — "these variable mean"→"variables"; "does not exists"→"exist"

### `test/test_GeoMakieExt.jl` / `test/test_MakieExt.jl`
- **L39 / L41** `[typo]` — "Intialize" → "Initialize"
- **`test_MakieExt.jl` L33** `[style]` — unused `path = "a/b/c"` (also L54, L102)

### `test/test_Template.jl`
- **L109** `[typo]` — "TemplateVar with ones to n data" → "one to n"

### `test/test_Leaderboard.jl`
- **L371, L379** `[grammar]` — "Order are the same/different" → "Order is..."

---

## Cross-cutting notes (important for fixing)

- **src↔test message coupling:** `src/Var.jl:1572` ("is not consistent"→"are"), `:1918`/`:1937`
  ("does not has"→"have") are asserted verbatim by `test_Var.jl:1278/1388/1467/2024/2043` —
  must change source message **and** test assertion together or tests break.
- **`Template.jl` lat/lon units swap (bugs):** the helper currently produces lat=`degrees_east`,
  lon=`degrees_north`. `test_Template.jl` (~L145) asserts the swapped values, so it masks the bug —
  fixing the source requires updating that test. `test_Atmos.jl:224-225` already uses the correct convention.
- **`warp_string` → `wrap_string` rename:** touches the export (`Utils.jl:4`), definition (`:265`),
  docstring + 6 doctests (`:237-261`), `api.md:208`, both extensions (`MakieExt:70`, `GeoMakieExt:90`,
  `MakieExt:312`), and `test_Utils.jl:77-90`. It renames a (semi-internal) `Utils` symbol → mildly breaking.

## Fix plan (phased) — agreed scope: **bugs + doc-errors + typos**

1. **Bugs** — src logic first (`Var.jl:1472`, `Leaderboard.jl:96/351/550/653`, `Utils.jl:288/375`,
   `masks.jl:93/273`, `Numerics.jl:60`), then test no-ops, then docs/doctests/`make.jl`.
   (Hold `Template.jl:226/233` until the deferred test-coupled renames are approved.)
2. **Doc-errors** — wrong API names, signatures, opposite-meaning docstrings, invalid admonitions, broken examples.
3. **Typos** — spelling/markdown (one-liners). `warp_string` rename held (deferred set).
4. **Deferred (not now):** grammar + style (132); src↔test message edits; `warp_string` rename;
   `Template` units swap. Cataloged above so they're not lost.

Each finding's exact line/quote/fix is in `TYPO_BUG_AUDIT.md` (repo) and the session
`scratchpad/findings_byfile.json`. Pass 2's new findings will be merged before fixing begins.

## Verification

- `julia --project -e 'using Pkg; Pkg.test()'` — runs the suite; confirms the test-no-op fixes now
  assert (and that the `src↔test` message changes stay in sync).
- `julia --project=docs docs/make.jl` (or the doctest runner `test/doctest.jl`) — confirms the doctest
  fixes (`var.md`, `howdoi.md`, `NEWS.md` examples) build and that doctest outputs match.
- `test/format.jl` (JuliaFormatter) and `test/aqua.jl` — confirm formatting/quality after edits.
- Spot-check the highest-impact behavior change: add/verify a test that `arecompatible` returns
  `false` for two vars whose dimension *names* differ (this path was dead before the `Var.jl:1472` fix).

---

## Pass 2 — additional findings (in progress)

A deeper second pass (smaller chunks + logic/signature/doctest/message lenses, deduped vs pass 1)
plus a repeated-token sweep. New issues pass 1 missed:

### Repeated instances pass 1 caught only once
- **`src/Var.jl:1118-1119` & `:1144-1145`** `[doc-error]` — the variance docstring error
  ("dividing the **sample mean** by n-1" → should be "sum of squared deviations from the mean")
  also appears in `variance_lon` and `variance_lat`, not just `variance_time` (`:1098-1099`).
  All three should be reworded identically. *(Computation itself is correct — pure docstring fix.)*
- **`test/test_MakieExt.jl:106`** `[typo]` — third "Intialize" → "Initialize" in this file (pass 1 caught `:41`).
- **`test/test_Var.jl:2387`** `[typo]` — third "Cheak season" → "Check season" (pass 1 caught `:2324`, `:2441`).

### New bugs found in manual deep-dive
- **`src/Var.jl:864` & `:1002`** `[bug]` — the small-latitude/radians guard uses
  `abs(maximum(latitudes(var))) >= 0.5π`. That is `abs` of the maximum, not the maximum
  *magnitude*. For a southern-hemisphere or equator-ending domain in degrees (e.g. lats `-90:0`,
  `maximum = 0`), the guard is `0 >= 1.57` → **false** → it emits a spurious
  *"Detected latitudes are small. If units are radians, results will be wrong"* warning even
  though the units are degrees. **Fix (both lines):** `maximum(abs.(latitudes(var))) >= 0.5π`.
  Appears in `average_lat` (weighted branch) and `average_lonlat` (weighted branch).
- **`src/Leaderboard.jl:655`** `[doc-error]` — in `find_worst_single_model` the comment says
  "Replace all NaN with **Inf**" but line 658 correctly uses `NaN => -Inf` (mirrors `find_best`'s
  `findmax`). Comment copy-pasted from `find_best`; should say **-Inf**.

### Pass-2 workflow results — 75 new verified findings (21 bug, 32 doc-error, 1 typo, 13 grammar, 8 style)

*(After dropping 189 pass-1 duplicates + 20 false positives; then a quote-match re-dedup removed
20 cross-file restatements of pass-1 test findings.)* Highlighted new **bugs**:
`outvar_operators.jl:49` (no `deepcopy` → result aliases input dims), `Var.jl:1326`
(`circshift(data, shift)` shifts only dim 1; the unused `shift_tuple` shows lon dim was intended),
`Var.jl:2680` (`set_reference_date!` stores raw `reference_date` not parsed `ref_date`),
`GeoMakieExt.jl:85` (white mask default never applies — entry points pass `:mask => Dict()`),
`MakieExt.jl:864` (leaderboard `xticks`: 2 positions vs `num_models+1` labels → breaks for >1 model),
`Utils.jl:445` (`split_by_month` empty case over-nested type), `test_Var.jl:98` (aliased/tautological
`center_longitude!` assertion), `index.md:139` (`slice(var, 8_000)` — no positional value arg),
`index.md:245` (`import ... : x` with a space before `:` — verified to fail in Julia 1.12).

### SRC (pass 2)

**`src/Catalog.jl`** — L149 `[doc-error]` docstring sig `var_kwargs` → `var_kwargs = ()`; L151 `[grammar]` "a `OutputVar`"→"an".
**`src/Leaderboard.jl`** — L552 `[bug]` dead `rmses = RMSEs |> copy` (L550); L667 `[doc-error]` median docstring wrong; L402/469/502 `[grammar]` "Support"→"Supports"/"is not"→"are not"; L630 `[style]` `ann_idx`→`category_idx`.
**`src/Numerics.jl`** — L110 `[bug]` second `@test_throws` should test `_integrate_lon`; L72 `[doc-error]` docstring `int_weights::AbstractVector` type mismatch.
**`src/Sim.jl`** — L48 `[doc-error]` `variable_paths["short_name"]` should drop quotes; L192 `[doc-error]` `summary` docstring inaccurate.
**`src/Template.jl`** — L202 `[bug]` `add_dim!` default `name = $dim_name` vs curried `$default_dim_name` (inconsistent dim name); L428 `[grammar]` "is not"→"are not".
**`src/Var.jl`** — L1326 `[bug]` unused `shift_tuple`/`circshift` shifts wrong dim; L2680 `[bug]` `string(reference_date)`→`string(ref_date)`; L985 `[doc-error]` `average_lonlat` sig missing `weighted`; L1118 `[doc-error]` variance "sample mean" (variance_lon); L1738 `[doc-error]` `_resampled_as_partial` sig; L2373 `[doc-error]` "Add units back for bias"→"squared error"; L2466 `[doc-error]` malformed `shift_by` sentence; L107 `[typo]` "maps name array index"; L2152 `[grammar]` "squeeze"→"squeezes".
**`src/flat.jl`** — L168 `[style]` redundant `all(isapprox(...))` wrapper.
**`src/masks.jl`** — L3582 `[bug]` (in test_Var.jl) pass `lat_var` not `mask_var`; L71 `[doc-error]` docstring `var`→`mask_var`.
**`src/outvar_dimensions.jl`** — L142 `[doc-error]` "`dates`"→"`date`"; L145 `[grammar]` reword.
**`src/outvar_operators.jl`** — L49 `[bug]` OutputVar–OutputVar branch needs `deepcopy(x.dims)`/`deepcopy(x.dim_attributes)`.
**`src/outvar_selectors.jl`** — L455 `[bug]` `dim_attribs` positional `zip` of dims vs dim_attributes (should key by name); L124 `[doc-error]` sig needs `;` before `by`; L188/189 `[doc-error]` single→double quotes; L227 `[grammar]` message vs `>=`.

### EXT (pass 2)
**`ext/ClimaAnalysisGeoMakieExt.jl`** — L85 `[bug]` white mask default never applies; L143 `[doc-error]` "when computing the bias".
**`ext/ClimaAnalysisMakieExt.jl`** — L864 `[bug]` `xticks` positions vs labels length mismatch; L817 `[doc-error]` comment `rmse_model_vars`→`rmse_vars`.

### DOCS (pass 2)
**`README.md`** — L62 `[doc-error]` `> :note:` → `> [!NOTE]` GitHub alert.
**`docs/src/api.md`** — L125 `[doc-error]` `Base.cat(...; dim::String)` → `dims::String`.
**`docs/src/index.md`** — L139 `[bug]` `slice(var, 8_000)` → `slice(var, z = 8_000)`; L245 `[bug]` `import ... : x` space before colon; L54 `[doc-error]` orog tree/summary mismatch; L65 `[doc-error]` `period = "3.0h"`→`"1.0h"`; L131 `[doc-error]` bare `.dims =` blocks.
**`docs/src/rmse_var.md`** — L74 `[grammar]` "and rest"→"and the rest"; L135 `[grammar]` "for find"→"for finding".
**`docs/src/visualize.md`** — L85 `[bug]` `plot_bias_on_globe!(fig, var, ...)` → `(fig, sim_var, obs_var, ...)`; L156 `[doc-error]` `var`→`var2d`.

### TEST (pass 2)
**`test/test_GeoMakieExt.jl`** — L18 `[doc-error]` `@testset "MakieExt"`→"GeoMakieExt"; L113 `[doc-error]` filename says contours but is heatmap; L232 `[style]` double `#`.
**`test/test_Leaderboard.jl`** — L15 `[style]` extra space after `#`.
**`test/test_MakieExt.jl`** — L148 `[bug]` `:vertical => :false` (Symbol) → `false` (Bool); L160 `[bug]` duplicate output filename; L460 `[style]` stale "doesn't work?" comment.
**`test/test_Numerics.jl`** — L110 `[bug]` duplicate `@test_throws` should cover the longitude path.
**`test/test_Sim.jl`** — L126 `[doc-error]` comment `is_z_1d`→`is_z_1D`.
**`test/test_Template.jl`** — L144 `[bug]` encodes the swapped lat/lon units (masks `Template.jl` bug).
**`test/test_Utils.jl`** — L76 `[style]` `@testset "format_title"` should be "warp_string".
**`test/test_Var.jl`** — L98 `[bug]` aliased/tautological assertion (use `copy(data)`); L1623 `[bug]` `drop_lon_var`→`drop_both_var`; L96 `[doc-error]` "shifting by 91"→"-91"; L161 `[doc-error]` "Center"→"Shift"; L1221/1597/2526 `[doc-error]` wrong comments; L3666/3674 `[grammar]` "exists"→"exist"; L1715 `[style]` ongrid/oncenter vs oncell.
**`test/test_flat.jl`** — L160 `[bug]` `flat_var`→`flat_var_mask`; L348 `[style]` shadowed loop var `perm`.
**`test/test_outvar_selectors.jl`** — L174 `[bug]` duplicated slice block; L172 `[doc-error]` misleading comment; L552 `[grammar]` "a OutputVar"→"an".

---

### Pass 3 — novel-lens hunt (complete)

Seven file-group auditors fanned out across the whole repo (src split 4 ways, ext, all tests, all docs),
each applying the **sibling-asymmetry / aliasing-mutation / numerical-edge / kwarg-threading /
doctest-recompute / test-correctness** lenses and deduping against all 360+ prior-pass lines. Each agent
verified every candidate against the actual code (two agents ran the package to confirm aliasing via `===`).
**47 new verified findings** survived dedup: **9 bug, 16 doc-error, 0 typo, 17 grammar, 5 style.**

| Location | bug | doc-error | typo | grammar | style | total |
|---|---|---|---|---|---|---|
| SRC | 3 | 6 | 0 | 11 | 2 | **22** |
| EXT | 1 | 7 | 0 | 6 | 0 | **14** |
| DOCS | 2 | 3 | 0 | 0 | 1 | **6** |
| TEST | 3 | 0 | 0 | 0 | 2 | **5** |
| **Total** | **9** | **16** | **0** | **17** | **5** | **47** |

#### 🔑 Cross-cutting theme — aliasing (non-mutating funcs returning shared mutable state)
The codebase convention is `deepcopy(var.dims)` when building a derived `OutputVar` (used in ~8 places in
`Var.jl` and in `outvar_operators.jl:81/205`). Pass 1/2 found one violation (`outvar_operators.jl:49`).
Pass 3 found **three more** sibling functions that break it — they should all be fixed together:
- **`src/outvar_selectors.jl:236`** `[bug]` — `window`: `dims = copy(var.dims)` (shallow) → non-windowed
  dim arrays are shared with the input (verified `w.dims["lon"] === var.dims["lon"]`). Use per-array copy / `deepcopy`.
- **`src/outvar_selectors.jl:241`** `[bug]` — `window`: `dim_attributes = copy(var.dim_attributes)` (shallow) →
  inner per-dim attr `Dict`s shared (verified mutation leaks back to `var`). `select`/`view_select` (L456) copy per-inner-dict; `window` is the outlier.
- **`src/Atmos.jl:90-93` (+ `:99-102`)** `[bug]` — `to_pressure_coordinates`: the dims/dim_attributes
  comprehensions reuse the original `v` for non-altitude dims (lon/lat) → returned var aliases the input's arrays. `deepcopy` the reused `v`.

#### SRC (pass 3)

**`src/Var.jl`** —
- L88 `[grammar]` — struct docstring "Representing an output variable" → "Represents…" (distinct from `flat.jl:4`).
- L904 `[style]` — stray blank line between `"""` and the signature in the `average_lon` docstring; breaks the uniform signature-block rendering.
- L1489 `[doc-error]` — `_check_dims_consistent(x, y, dim_names = nothing)` header uses a comma, but the real arg is a **keyword** (`; dim_names`).
- L1648 `[doc-error]` — `resampled_as(src_var, dest_var, dim_names = nothing)` header comma → `; dim_names` (keyword in the real signature).
- L1778 `[doc-error]` — `resampled_as(src_var::OutputVar, kwargs...)` header comma → `; kwargs...`.
- L2983 `[grammar]` — `Base.show` zero-length branch prints "element" (singular) → "elements" (cf. the multi-element branch).

**`src/Leaderboard.jl`** — L390/412/424 `[grammar]` — three more "`Support` indexing"/"Do not support" → "`Supports`"/"Does not support" (new instances of the L402-class error, in `getindex`/`setindex!`/`keys` docstrings).

**`src/masks.jl`** —
- L158 `[doc-error]` — `_generate_binary_mask` docstring header (and body) names the first arg `mask`, but the real signature (L167) is `mask_var`.
- L252 `[grammar]` — "takes in a `OutputVar` and `mask` the data" → "masks" (subject-verb; distinct from the L252 article fix already listed).

**`src/Numerics.jl`** —
- L34/55 `[grammar]` — two more "each point `correspond`" → "corresponds" (`_integrate_lat`, `_integrate_dim`; pass 1 had only L13).
- L126/160 `[grammar]` — two more "it `make` no contribution" → "makes" (pass 1 had only L112).
- L172 `[doc-error]` — `_integration_weights_generic_equispaced` carries the comment "where we use the assumption that units are degrees", but the code (L174) does a bare `arr[begin+1]-arr[begin]` with **no** `deg2rad` → comment is false here (distinct from the L157 instance).

**`src/Atmos.jl`** — L192 `[doc-error]` — comment "Compute global mse and global rmse and store it as an attribute" is copy-pasted from `Var.jl:2387`, but `global_rmse_pfull` stores nothing — it just returns the number. Drop "and store it as an attribute".

**`src/Visualize.jl`** — L83 `[style]` — `print(io, "..."; )` has a stray trailing `;` before `)` (no-op keyword separator); same on L80-81. Cosmetic.

**`src/flat.jl`** — L28 `[grammar]` — `FlatVar` docstring "Representing a flat output variable" → "Represents…" (distinct line from the L4 fix).

#### EXT (pass 3)

**`ext/ClimaAnalysisMakieExt.jl`** —
- L903 `[bug]` — the leaderboard `Colorbar` uses `limits = extrema(rmse_no_nan_vec)` while the heatmap uses
  `colorrange = (1e-10, maximum(rmse_no_nan_vec))` with `lowclip = :white` (L894). The colorbar's lower bound
  disagrees with the heatmap's clipped range → colorbar can mislabel displayed colors. Make `limits` match `colorrange`.
- L915/931/957/967 `[grammar]` — four "a `OutputVar`"/"a OutputVar" → "an" (two `convert_arguments` docstrings, the `plottype` docstring, and the user-facing "Plotting a OutputVar with $N dimensions" error string).

**`ext/ClimaAnalysisGeoMakieExt.jl`** —
- L212 `[doc-error]` — the `contour2D_on_globe!` **GridLayout** signature block lists the phantom `plot_contours = true,` kwarg (no such keyword exists; same class as the L204 finding on the other signature).
- L232 `[doc-error]` — `contour2D_on_globe!` docstring says "applied to the `OutputVar`s when computing the bias", but it plots a single var and computes no bias (copy-paste from `plot_bias_on_globe!`; pass 1 caught the identical L143 in `heatmap2D_on_globe!`).
- L252/340 `[doc-error]` — two more "plotted from `GeoMakie.coastline`" → `coastlines` (plural; actual call at L98 is `coastlines()`; pass 1 caught only L163).
- L145/234/321 `[doc-error]` — three docstrings steer users to the **deprecated** `make_lonlat_mask` via `[…](@ref)` (and note text) → should reference `generate_lonlat_mask` (same nature as the `visualize.md:97` finding, in ext docstrings not previously listed).
- L237/325 `[grammar]` — two more "ClimaAnalysis do not support mask keyword arguments" → "does not support" (pass 1 caught only L148).

#### DOCS (pass 3)

**`docs/src/var.md`** —
- L357 `[bug]` — `julia> ClimaAnalysis.global_mse(sim, obs)` uses undefined `sim`/`obs`; the block (L345-353) defines `sim_var`/`obs_var` → `global_mse(sim_var, obs_var)` (fourth occurrence; pass 1 caught L316/324/327).
- L347 `[style]` — comment "# Load in 3D temperature variable…" is an exact duplicate of L344's comment but is attached to `sim_var`; should describe simulation (not obs) data.

**`docs/src/api.md`** — L120 `[doc-error]` *(borderline — verify vs doc policy)* — the `@docs` block lists the deprecated `Var.make_lonlat_mask` alongside `generate_lonlat_mask` (L117); documenting the deprecated alias may mislead users.

**`README.md`** — L90/123 `[doc-error]` — two more `> :note:` blocks → `> [!NOTE]` (GitHub renders `:note:` as literal text; pass 2 caught only L62).

**`NEWS.md`** — L331 `[bug]` — changelog example `DJF = cat(seasons[begin:4:end]..., dim = "time")` uses singular `dim`, but `Base.cat(vars::OutputVar...; dims::String)` (src/Var.jl:2821) requires `dims` → would error with an unsupported keyword (correct usage is at `howdoi.md:361`).

#### TEST (pass 3)

**`test/aqua.jl`** — L10 `[bug]` — ⭐ `Aqua.detect_ambiguities(ClimaAnalysis; recursive = true)` registers no `@test` (resolves to non-asserting `Test.detect_ambiguities`) → ambiguity checking is a silent no-op. Use `Aqua.test_ambiguities(...)`.

**`test/test_Var.jl`** — L2276/2277 `[bug]` — bare `JJA/DJF.attributes["start_date"] == "2024-1-1"` with no `@test` (no-op); second occurrence of the 2253-2255 defect, in the "Split by season with seasons kwarg" block.

**`test/test_GeoMakieExt.jl`** — L32 `[style]` — `path = "a/b/c"` assigned but never used (same dead-assignment class as `test_MakieExt.jl:33/54/102`, not previously listed here).

**`test/test_Leaderboard.jl`** — L29 `[style]` — double space after `#` in the constructor-test comment (second instance; pass 2 caught L15).

#### Verified non-issues worth recording (checked, deliberately NOT reported)
- `masks.jl:151` — `_integration_weights_lat_equispaced` uses `deg2rad.(scalar)` (0-d array) vs the lon version's `fill(deg2rad(...))`; numerically equivalent under broadcasting, pure style.
- `GeoMakieExt.jl:349` — `cmap_extrema` default ignores the `mask` while the plotted `bias_var` applies it; plausibly intentional design, left as borderline.
- `Var.jl:623` `convert_units` doctest uses `var` not `var_bob`, but the printed `(0.0, 500.0)` is still correct (non-parseable units take the `conversion_function` path) — not a defect.
- `test_Numerics.jl:5/15/25` "(not equispaced)" comments are **correct** (data really is non-equispaced); only L35/45/55 are mislabeled (already in plan).
- All Utils.jl / outvar_dimensions.jl jldoctests recompute correctly; `select`/`view_select` and the operator families are aliasing-safe — `window` and `to_pressure_coordinates` are the lone outliers. **[SUPERSEDED by Pass 5: `_reduce_over` (`Var.jl:839-840`) also aliases — see Pass 5 bugs.]**

> **Pass 3 fix-scope note:** the 9 bugs + 16 doc-errors fold into the existing phased plan (bugs first,
> then doc-errors, then typos). The aliasing trio should be fixed in one commit alongside `outvar_operators.jl:49`
> and covered by a regression test asserting `window`/`to_pressure_coordinates` outputs do **not** `===` their inputs'
> dim arrays. `aqua.jl:10` is the highest-leverage single fix (re-enables a whole CI check class). The 17 grammar +
> 5 style items join the deferred set.

---

### Pass 4 — static-analysis deep dive (complete, no code executed)

At the user's instruction this pass used **pure static analysis** (file reading + grep only — no julia, tests,
formatter, or doctests run). Six auditors applied lenses not used in passes 1-3: **math-formula verification /
control-flow & type / `@ref`+`@docs` integrity / test-semantics / IO+edge-cases / docs manual-examples**, each
deduping against all 400+ prior-pass lines. **14 new verified findings: 4 bug, 2 doc-error, 4 typo, 4 grammar, 0 style.**

| Location | bug | doc-error | typo | grammar | style | total |
|---|---|---|---|---|---|---|
| SRC | 2 | 0 | 1 | 4 | 0 | **7** |
| EXT | 1 | 0 | 0 | 0 | 0 | **1** |
| DOCS | 0 | 2 | 3 | 0 | 0 | **5** |
| TEST | 1 | 0 | 0 | 0 | 0 | **1** |
| **Total** | **4** | **2** | **4** | **4** | **0** | **14** |

#### Bugs (pass 4)
- **`ext/ClimaAnalysisUnitfulExt.jl:79`** `[bug]` — `convert_units` returns `OutputVar(new_attribs, var.dims,
  var.dim_attributes, new_data)` with **no `deepcopy`** → result aliases the input's dim arrays + dim-attribute dicts.
  4th member of the aliasing quartet (see Highest-impact section); fix alongside the others.
- **`src/Leaderboard.jl:90-91`** `[bug]` — `for key in keys(units) … delete!(units, key)` deletes from the `Dict`
  while iterating its live `KeySet` (undefined behavior in Julia). **Fix:** `for key in collect(keys(units))`.
- **`src/Leaderboard.jl:84-98`** `[bug]` — the `RMSEVariable` 5-arg constructor mutates the **caller's** `units`
  dict (L85 add, L91 delete) and stores that same reference (L98) — caller side-effect + aliasing, against the
  package's `deepcopy` convention. **Fix:** `units = copy(units)` at the top (also resolves the L90-91 bug).
- **`test/test_Template.jl:167-168`** `[bug]` — `"lon" => Dict("units" => "degrees_north")`, `"lat" => …"degrees_east"`
  encodes the **swapped** lat/lon units a *second* time (pass 1/2 caught L144-145), masking the `Template.jl:226/233`
  source swap. Fix together with the source swap + the L144-145 assertion.

#### Doc-errors (pass 4)
- **`docs/src/visualize.md:3`** `[doc-error]` — "consult `[`Visualize`](@ref)`" is an unresolvable cross-reference:
  the `Visualize` module has no docstring and appears in no `@docs` block, so Documenter can't resolve it. **Fix:**
  add `Visualize` to a `@docs` block (with a module docstring) or drop the `(@ref)`. (Prior passes noted only a
  grammar nit on this line.) *Verified: every other `@ref`/`@docs` target in the repo resolves.*
- **`NEWS.md:8`** `[doc-error]` — "an initial implemention of the **split-combine-apply** pattern" → "split-**apply**-combine"
  (correct order everywhere else: same release's L23, the page title, `split_apply_combine.md`). Distinct from the
  already-listed "implemention"→"implementation" typo on this line.

#### Typos (pass 4)
- **`src/masks.jl:140`** `[typo]` — "ocean are **multipled** by `ocean`" → "multiplied" (2nd occurrence; pass 1 caught only L116 in `generate_land_mask`).
- **`NEWS.md:785`** `[typo]` — "rely on the **interpolats**" → "interpolants" (codebase uses "interpolant" 17× in src).
- **`NEWS.md:816`** `[typo]` — "an **interpolat** will not be made" → "interpolant".
- **`docs/src/var.md:109`** `[typo]` — "the **interpolats**, such as `resampled_as`" → "interpolants" (distinct from the adjacent "which mean" grammar nit already listed).

#### Grammar (pass 4)
- **`src/masks.jl:46`** `[grammar]` — `LonLatMask` docstring "**Representing** a mask…" → "Represents…" (same class as `Var.jl:88`/`flat.jl:4/28`; this line was only listed for its "a `OutputVar`"→"an" article fix).
- **`src/Var.jl:1887/1910/1929`** `[grammar]` — three integration docstrings (`integrate_lonlat`/`integrate_lon`/`integrate_lat`): "each point **correspond** to the midpoint" → "corresponds" (pass 1/3 tracked this verb only in `Numerics.jl`).

#### Strong negative results (verified correct — recorded so they aren't re-investigated)
- **Core math is correct.** Hand-verified: trapezoidal/midpoint integration weights, the `deg2rad`/`cos(lat)`
  area-weighting (right axis, normalized `Σ data·cos / Σ cos`), variance n−1 (`corrected`), `bias = sim−obs`,
  `mse = ∫(sim−obs)²/∫1`, `rmse = √mse` with consistent masked normalization, the 0/1 mask multiply, and the
  land/ocean inversion. No new computational defect survived verification.
- **No resource leaks.** Every `NCDataset(...)` open in `Catalog.jl`/`Var.jl`/`masks.jl` uses the `do…end` form;
  path handling uses `joinpath` and correct `.nc` extension checks.
- **No fixed dim-order assumption** in `to_pressure_coordinates` (uses `dim2index[z_name]` + per-index `ntuple`). **[WRONG — REFUTED by Pass 5: `Atmos.jl:118` assumes altitude is the last dim and breaks otherwise.]**

#### Near-misses worth recording (verified, deliberately NOT reported)
- **`src/Utils.jl:504-505`** `_isequispaced` indexes `arr[begin+1]` with no length check → latent `BoundsError`
  on a single-element/empty array; same unguarded pattern at `src/Var.jl:187`/`:221`. Currently unreachable (all
  callers guard single-point dims earlier), so a latent edge, not a live bug — flagged here in case a future caller removes the guard.
- **`ext/ClimaAnalysisMakieExt.jl:855`** trailing adjoint `'` on `vcat(...)'` is superfluous (lengths match) — cosmetic.
- **`test/test_Var.jl:3344`** a var named `MAM` actually concatenates `seasons[2]` (MAM) + `seasons[5]` (Dec/DJF) —
  misleading name only; the expected data/dims/attributes are all self-consistent, so no assertion is wrong.
- **`src/Leaderboard.jl:616`** `first(units)` on an empty dict — unreachable (an `RMSEVariable` always has ≥1 model).

> **Pass 4 fix-scope note:** the 4 bugs slot into the bugs-first phase — the `convert_units` aliasing joins the
> quartet commit; the two `Leaderboard` constructor bugs share one `copy(units)` fix; `test_Template.jl:167-168`
> rides along with the already-deferred `Template.jl` units-swap decision. The 2 doc-errors + 4 typos join their
> phases; the 4 grammar items join the deferred set. Static-only: none of these were confirmed by execution —
> recommend the suite/doctests be run once before committing (the user asked not to run code this pass).

---

### Pass 5 — new-dimension static-analysis hunt (complete, no code executed)

Per the user's instruction this pass stayed **pure static analysis** (read + grep only — no julia, tests,
formatter, or doctests) and deliberately opened **dimensions/lenses no prior pass had used**. Five auditors
fanned out, one per new lens, each deduping against all ~421 prior-pass lines:

1. **Configuration / build / CI files** — an entirely unexamined *location*: `Project.toml`, `Manifest.toml`,
   `Artifacts.toml`, `docs/Project.toml`, the five `.github/workflows/*.yml`, `dependabot.yml`,
   `.JuliaFormatter.toml`, `.gitignore`, `NOTICE`, `LICENSE`.
2. **Comparison & boolean operators + off-by-one / bounds / slicing** (`>=`/`>`, `&&`/`||`, negation,
   `begin`/`end`, range arithmetic, float `==` vs `≈`).
3. **Mutation `!`-convention & in-place semantics** (name-vs-behavior; non-`!` siblings that secretly mutate).
4. **Export ↔ definition ↔ `@docs` ↔ NEWS/deprecation integrity** (export-list/definition/`@docs` surface,
   not the `@ref`-resolution angle Pass 4 used).
5. **Collection / iteration & `nothing`-handling** (`findfirst`/`get`/`zip`/`enumerate`/`only`,
   dict-mutate-during-iteration, `CartesianIndex` arithmetic).

Plus a deterministic, repo-wide **string-interpolation sweep** (`$ident(` = interpolating a function instead
of its result): the only hit besides metaprogramming (`outvar_operators.jl`'s `$op(` inside `quote`, benign)
was the already-known `Leaderboard.jl:96`. Clean negative result.

**7 new verified findings: 3 bug, 3 doc-error, 0 typo, 0 grammar, 1 style.** Two of the three bugs **refute
explicit prior-pass "verified-correct" notes** (see "Corrections" below) — the highest-value outcome of this pass.

| Location | bug | doc-error | typo | grammar | style | total |
|---|---|---|---|---|---|---|
| SRC | 2 | 1 | 0 | 0 | 0 | **3** |
| DOCS (NEWS) | 0 | 2 | 0 | 0 | 0 | **2** |
| CONFIG/CI | 1 | 0 | 0 | 0 | 1 | **2** |
| **Total** | **3** | **3** | **0** | **0** | **1** | **7** |

#### ⚠️ Corrections to prior-pass "verified non-issue" notes (both wrong)
- Pass 3's note (this file, "the operator families are aliasing-safe — `window` and `to_pressure_coordinates`
  are the lone outliers") is **incomplete**: `Var.jl:839-840` (`_reduce_over`) aliases too, and it underlies
  *every* reduction (`average_*`, `variance_*`, `integrate`, `slice`). See bug below.
- Pass 4's note ("No fixed dim-order assumption in `to_pressure_coordinates`") is **wrong**: `Atmos.jl:118`
  bakes in an implicit "altitude is the last dimension" assumption and breaks otherwise. See bug below.

#### Bugs (pass 5)
- **`src/Atmos.jl:118`** `[bug]` — ⭐ *(found independently by the comparison **and** the collection auditors)* —
  `to_pressure_coordinates`'s column loop mis-indexes the reduced `CartesianIndex` when **altitude is not the
  last dimension.** `ranges = [1:size(var.data)[i] for i in 1:num_dims if i != z_index]` (L114) has
  `num_dims-1` entries, so `idx` is a `CartesianIndex{num_dims-1}` whose *k*-th component maps to the *k*-th
  **non-z** axis. But L118 reads `idx[i]` using the **original** full-array axis number `i`:
  `indices = ntuple(i -> (i == z_index ? Colon() : idx[i]), num_dims)`. Only correct when `z_index == num_dims`.
  When `z_index < num_dims`: every axis after z reads the wrong component, and the top axis `i = num_dims`
  evaluates `idx[num_dims]` on a length-`(num_dims-1)` index → **`BoundsError`** (e.g. `dims = ("z","lon","lat")`,
  or any 4-D `(lon, lat, z, time)` var). **Reachable:** `to_pressure_coordinates` is public, performs **no**
  dim reorder before the loop, and `global_rmse_pfull` forwards vars unchanged; the bug is masked only because
  every multi-dim test puts `z` last. **Fix:** consume `idx` by its own running position, e.g.
  `let j = 0; ntuple(i -> i == z_index ? Colon() : (j += 1; idx[j]), num_dims) end` (or precompute
  `non_z_axes = filter(!=(z_index), 1:num_dims)` and map `idx[k] → non_z_axes[k]`). Add a z-not-last
  regression test (e.g. dims `("z","lon","lat")`).
- **`src/Var.jl:839-840`** `[bug]` — ⭐ **broadest member of the aliasing family** (same class as the
  `outvar_operators.jl:49` / `window` / `to_pressure_coordinates` / `convert_units` quartet). `_reduce_over` does
  `dims_dict = copy(var.dims)` / `dim_attributes = copy(var.dim_attributes)` — *shallow* container copies. The
  4-arg `OutputVar` constructor then rebuilds the containers with `OrderedDict(dims)` (L267-268), which **keeps
  the inner array/dict references**, so the returned var's un-reduced dimensions satisfy
  `result.dims["lon"] === var.dims["lon"]` (and likewise per-dim attribute `Dict`s). Because `_reduce_over`
  backs **all** reductions — `average_lon`/`average_x`/`average_y`/`average_time`, `variance_*`,
  generic `_average_dims`, `integrate`, and `slice` (via `_slice_general`) — every non-mutating reduction returns
  a var that shares mutable dim state with its input. (The `average_lat`/`average_lonlat` *weighted* branches are
  incidentally safe only because they `copy(var)` first; the non-weighted paths are not.) **Fix:**
  `dims_dict = deepcopy(var.dims)` and `dim_attributes = deepcopy(var.dim_attributes)` at L839-840 — fixes all
  paths uniformly. Fold into the aliasing commit + the "result must not `===` input dim arrays" regression test.
- **`.github/workflows/TagBot.yml:8`** `[bug]` — `if: github.event_name == 'workflow_dispatch' || github.actor
  == 'JuliaTagBot'` references a `workflow_dispatch` event, but `on:` (L2-5) declares only `issue_comment` →
  the `workflow_dispatch` half is **dead**, so a maintainer can never manually re-run TagBot (e.g. to retry a
  missed tag). The canonical JuliaRegistries template includes the trigger. **Fix:** add `workflow_dispatch:`
  under `on:`. *(Config files were untouched by passes 1-4 — entirely new location.)*

#### Doc-errors (pass 5)
- **`src/Var.jl:2748-2752`** `[doc-error]` — the docstring heading the **in-place** `Base.replace!(new, var; count)`
  (L2754, mutates `var.data`, returns `nothing`) reads *"Return each value of `var.data` by `new(x)`."* — wrong
  mutation semantics (an in-place mutator described as "Return … value") and grammatically broken (missing verb).
  Distinct from the already-reported L2734/2736, which heads the **non-mutating** `replace`. **Fix:** e.g.
  "Replace each value of `var.data` with `new(x)` in place. … See [`replace`](@ref) for the non-mutating version."
- **`NEWS.md:98` & `NEWS.md:191`** `[doc-error]` — both reference a nonexistent **`select_view`** ("you can now
  use `select` or `select_view` …"); the real exported function is **`view_select`** (`src/outvar_selectors.jl:6`
  export, def at `:408`). `select_view` has 0 occurrences anywhere else in the repo. **Fix:** `select_view` →
  `view_select` on both lines.

#### Style (pass 5)
- **`Project.toml:37`** `[style]` — `NaNStatistics = "=0.6.8, 0.6.8 - 0.6.50, 0.6.53 - 0.6.56, 0.6.58"` — the
  leading exact-pin `=0.6.8` is redundant (the next clause `0.6.8 - 0.6.50` already includes 0.6.8); resolution is
  unaffected. The gap structure correctly excludes the bad 0.6.51/0.6.52 (NEWS:402) and 0.6.57 (NEWS:51).
  **Fix:** drop the first clause → `"0.6.8 - 0.6.50, 0.6.53 - 0.6.56, 0.6.58"`.

#### Refinement of an existing entry (not counted as new)
- **`src/Var.jl:1306-1307`** — Pass 1 logged "deprecation points to `center_longitude` → should reference
  `shift_longitude`" as a *docstring* issue. The docstring (`!!! warn` block, L1300-1302) is actually **fine**;
  the defect is the **runtime `Base.depwarn` message** (L1306-1307), and it has a *second* error beyond the name:
  it invents a `; shift_by = 0.0` kwarg that `shift_longitude(var, lower_lon, upper_lon)` (L1353) does not have.
  **Fix:** message → `` `shift_longitude(var, lower_lon, upper_lon)` `` (no `shift_by`).

#### Strong negative results (verified, recorded so they aren't re-investigated)
- **CONFIG lens is otherwise clean:** every `[deps]`/`[weakdeps]`/`[extras]` package + `julia` has a `[compat]`
  entry; `[extensions]`↔`[weakdeps]`↔`ext/` filenames all consistent; the `test` target covers every test-only
  `using`; no malformed/duplicate UUIDs; CI matrix (`1.9/1.10/1.11`) matches the `julia = "1.9"` lower bound;
  action pins (`checkout@v7`, `codecov@v7`, `cache@v3`) are intentional & consistent; `Downgrade.yml` skip-list is
  correct; `dependabot.yml`/`.JuliaFormatter.toml`/`.gitignore` keys all valid. `ClimaAnalysisUnitfulExt.jl` having
  no `[extensions]` entry is **intentional** (its header says Unitful is loaded directly for now).
- **API surface fully self-consistent:** all 137 export entries resolve to a definition; all ~120 `@docs` entries
  resolve; exported set == documented set; every `src/*.jl` is `include`d; no undeclared `import`s.
- **Comparison/boolean lens otherwise clean:** `wrap_longitude` half-open interval, `arecompatible` strict `==`
  on grids, season-month partitions, `find_best`/`find_worst` Inf/-Inf pairing, `_slice_over` reshape arithmetic,
  and `_constrained_cmap`'s correct `maximum(abs, …)` form all verified correct.
- **`nothing`-handling clean:** every `findfirst`/`indexin`/`match`/`only`/`first`/`pop!` site is either
  `isnothing`-guarded or reached only after a `length(...) == 1` / validation check; the `zip`/`enumerate`
  pairings (`Var.jl:2020/1731/2852`) are length-aligned by dispatch/construction.
- **Latent (not live) edges** worth a glance only if a future caller changes: `_transform_dates`'s `@. date_arr =
  …` could in principle mutate `var.dims["date"]` from the non-mutating `transform_dates`, but is unreachable (the
  function requires seconds-`time` and errors on `DateTime` time arrays); `MakieExt:894/903` `maximum`/`extrema`
  throw on an all-NaN leaderboard column (degenerate input).

> **Pass 5 fix-scope note:** the 3 bugs are bugs-first — `Atmos.jl:118` is a standalone correctness fix (+ a
> z-not-last regression test); `Var.jl:839-840` folds into the aliasing-family commit (it is the most impactful
> member); `TagBot.yml:8` is a one-line CI fix. The 2 doc-errors + 1 style + the `Var.jl:1306-1307` message
> refinement join their phases. Static-only, as instructed: none confirmed by execution — run the suite/doctests
> once before committing, and in particular add the `to_pressure_coordinates` z-not-last test, which would have
> caught `Atmos.jl:118` immediately.
