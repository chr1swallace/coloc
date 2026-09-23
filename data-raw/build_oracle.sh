#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Build a 6.0.3 oracle install, so the goldens stay regenerable once R/ has
# moved on. The goldens themselves are the truth; this is only the machine
# that produces them.
#
#   ./data-raw/build_oracle.sh [commit]
#   COLOC_ORACLE_LIB=$(./data-raw/build_oracle.sh) \
#       Rscript data-raw/make_golden.R
#
# PIN BY COMMIT, NEVER BY VERSION: three commits report Version: 6.0.3 and two
# of them differ in behaviour.
# ---------------------------------------------------------------------------
set -euo pipefail

COMMIT="${1:-8f20f0bc5e60ffc99e4c2f787bd55fd30cfe7c45}"
REPO="$(git rev-parse --show-toplevel)"
LIB="${COLOC_ORACLE_LIB:-${REPO}/.oracle-lib}"
SCRATCH="$(mktemp -d "${TMPDIR:-/tmp}/coloc-oracle-XXXXXX")"

cleanup() {
    git -C "$REPO" worktree remove --force "$SCRATCH" 2>/dev/null || true
    rm -rf "$SCRATCH"
}
trap cleanup EXIT

mkdir -p "$LIB"
git -C "$REPO" worktree add --detach "$SCRATCH" "$COMMIT" >&2
R CMD INSTALL --no-multiarch --with-keep.source -l "$LIB" "$SCRATCH" >&2

echo "oracle installed at $LIB (commit $COMMIT)" >&2
echo "now run:" >&2
echo "  COLOC_ORACLE_LIB=$LIB COLOC_ORACLE_SHA=$COMMIT Rscript data-raw/make_golden.R" >&2
echo "$LIB"
