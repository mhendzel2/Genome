R script placeholder convention

This project uses a simple placeholder convention to pass input file paths to R scripts that are executed via `utils.r_integration.RIntegration._execute_r_script`.

Placeholder format

- Use the exact pattern `{file:<name>}` inside R script templates.
- The RIntegration layer will write provided file contents to a temporary directory and then replace each `{file:<name>}` with the corresponding absolute path (using forward slashes) before running R.

Why this exists

- Avoids embedding large file contents or quoting issues inside generated R script strings.
- Works for binary or large text files (e.g., gzipped BED files) because data is written to disk first.
- Keeps R script templates simple and readable: the template treats `{file:counts}` like a normal filename argument to `read.table()`.

Recommended placeholder names and examples

- `{file:input}` — generic single input file
  - Example R usage: `data <- read.table("{file:input}", sep="\t", header=FALSE)`

- `{file:counts}` and `{file:samples}` — differential expression pipelines
  - Example R usage: `count_matrix <- read.table("{file:counts}", sep='\t', header=TRUE, row.names=1)`
  - Example R usage: `sample_info <- read.table("{file:samples}", sep='\t', header=TRUE)`

- `{file:tissue1}` and `{file:tissue2}` — tissue/two-sample comparisons
  - Example R usage: `t1 <- read.table("{file:tissue1}", sep='\t', header=FALSE)`

- `{file:hic}` — Hi-C BEDPE-like input for Hi-C analysis
  - Example R usage: `hic <- read.table("{file:hic}", sep='\t', header=FALSE)`

Notes for R script authors

- Use `jsonlite::toJSON(...)` in R to return structured results to Python via stdout.
- Avoid calling system-level functions in templates unless you trust the environment; RIntegration filters certain dangerous patterns but runtime privilege separation is recommended.
- Treat placeholders as plain file paths; do not try to embed raw file content inside the R script.

How Python caller provides files

- The Python caller passes a mapping of names to file contents to `RIntegration._execute_r_script(script, data_files=payload)` where `payload` is like `{'counts': <string content>, 'samples': <string content>}`.
- RIntegration will create temp files with sanitized filenames and substitute `{file:counts}` -> `/tmp/.../counts` before invoking R.

Example end-to-end (Python)

1. Prepare payload:

```python
payload = {
    'counts': counts_file_content,  # str or bytes
    'samples': samples_file_content
}

r_script = get_r_script('differential_expression')
ri._execute_r_script(r_script, data_files=payload)
```

2. In the R template, use the matching placeholders (as shown above).

Security notes

- The placeholder mechanism sanitizes filenames to avoid path traversal.
- R is still executed with the system R binary; ensure running environment is trusted and consider chrooting or running in an isolated container for untrusted uploads.
