#!/usr/bin/env bash

set -euo pipefail

branch=${1:-origin/codex/raw-lc-daily-rms}

echo "Fetching latest GitHub refs..."
git fetch origin

echo "Restoring code files from ${branch}"
git restore --source="$branch" -- \
  '*.py' \
  '*.sh' \
  '*.sbatch' \
  '.gitignore' \
  '.gitattributes' \
  'flow_readme.txt'

echo "Done. HPCC run products were left in place."
echo "Review code changes with: git status --short"
