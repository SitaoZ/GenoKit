#!/usr/bin/env bash
# Run from any directory. Optional argument: an alternate HTML output directory.
set -euo pipefail
repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
python_cmd="${PYTHON:-python}"
output_dir="${1:-$repo_dir/docs/build/html}"
"$python_cmd" -m sphinx -b html -W --keep-going "$repo_dir/docs/source" "$output_dir"
