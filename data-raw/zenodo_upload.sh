#!/usr/bin/env bash
## Upload one release of split eQTL Catalogue files (from data-raw/eqtl_cat_split.R) to a new Zenodo draft record.
##
## Usage:  ZENODO_TOKEN=xxxx bash zenodo_upload.sh <release_dir> <release> [sandbox]
##   eg    ZENODO_TOKEN=xxxx bash zenodo_upload.sh ~/scratch/eqtl_cat_split/r7 r7
##         ZENODO_TOKEN=xxxx bash zenodo_upload.sh ~/scratch/eqtl_cat_split/r8-beta r8-beta
##   add "sandbox" as a third argument to try everything on https://sandbox.zenodo.org first (needs a sandbox token)
##
## The token needs scopes deposit:write and deposit:actions (Zenodo > Applications > Personal access tokens).
## The script creates a draft, uploads every file in <release_dir>, and checks each file's md5 against Zenodo.
## It does NOT publish: review the draft in your browser, then press Publish there (publishing cannot be undone).
## Needs curl and python3. Uploads are large, so run it inside screen/tmux on a node with internet access.

set -euo pipefail
dir=${1:?release directory}; release=${2:?release name, eg r7 or r8-beta}
host=https://zenodo.org; [[ ${3:-} == sandbox ]] && host=https://sandbox.zenodo.org
: "${ZENODO_TOKEN:?set ZENODO_TOKEN}"
auth=(-H "Authorization: Bearer $ZENODO_TOKEN")
json() { python3 -c "import json,sys; d=json.load(sys.stdin); print($1)"; } # extract a field from JSON on stdin
## curl wrapper that works with old curl versions: prints the response body, or the HTTP status and Zenodo's
## error message on failure (returning non-zero)
api() {
  local out code
  out=$(curl -sS -w '\n%{http_code}' "$@") || return 1
  code=${out##*$'\n'}; out=${out%$'\n'*}
  if (( code >= 400 )); then echo "HTTP $code: $out" >&2; return 1; fi
  printf '%s' "$out"
}

case $release in
  r7) source_doi=10.5281/zenodo.23350544 ;;
  r8-beta) source_doi=10.5281/zenodo.22846842 ;;
  *) echo "unknown release $release" >&2; exit 1 ;;
esac

nfiles=$(find "$dir" -maxdepth 1 -type f ! -name '*.part*' | wc -l)
(( nfiles <= 100 )) || { echo "$nfiles files: over Zenodo's 100 files per record" >&2; exit 1; }

## 1. create a draft with its metadata. Edit creators/description as you wish; you can also edit them in the browser later.
metadata=$(python3 - "$release" "$source_doi" <<'EOF'
import json, sys
release, doi = sys.argv[1], sys.argv[2]
print(json.dumps({"metadata": {
  "upload_type": "dataset",
  "title": f"eQTL Catalogue {release} fine-mapped credible sets, split by genomic position for coloc",
  "creators": [{"name": "Wallace, Chris", "affiliation": "University of Cambridge"}],
  "description": (f"<p>The eQTL Catalogue {release} fine-mapped credible sets (<a href='https://doi.org/{doi}'>doi:{doi}</a>), "
                  "reorganised for fast regional queries by the R package coloc (query_eqtl_cat()). Each quantification method's "
                  "file is split into chunks sorted by chromosome and credible set start position; no credible set is split "
                  "across chunks. The index file lists, for each chunk and chromosome, the range of variant positions, "
                  "with file sizes and md5 checksums. The data are unchanged from the original release.</p>"
                  "<p>Please cite the eQTL Catalogue: Kerimov N et al. (2021) A compendium of uniformly processed human gene "
                  "expression and splicing quantitative trait loci. Nature Genetics 53:1290-1299.</p>"),
  "access_right": "open",
  "license": "cc-by-4.0",
  "keywords": ["eQTL Catalogue", "fine mapping", "credible sets", "colocalisation", "coloc"],
  "related_identifiers": [{"identifier": doi, "relation": "isDerivedFrom", "resource_type": "dataset"}],
  "version": release}}))
EOF
)
resp=$(api "${auth[@]}" -H "Content-Type: application/json" -X POST "$host/api/deposit/depositions" -d "$metadata")
id=$(json "d['id']" <<<"$resp"); bucket=$(json "d['links']['bucket']" <<<"$resp")
echo "Created draft $id"

## 2. upload every file (retrying network errors). Rerunning the script creates a new draft, so if it fails part way,
##    discard that draft in the browser, or resume by setting id/bucket by hand and rerunning this loop.
i=0
for f in "$dir"/*; do
  [[ -f $f && $f != *.part* ]] || continue
  i=$((i + 1)); echo "[$i/$nfiles] $(basename "$f")"
  for try in 1 2 3 4 5; do # retry network errors and server errors
    api "${auth[@]}" --upload-file "$f" "$bucket/$(basename "$f")" > /dev/null && break
    (( try == 5 )) && { echo "Giving up on $(basename "$f")" >&2; exit 1; }
    echo "  retrying in $((try * 10))s" >&2; sleep $((try * 10))
  done
done

## 3. check every file arrived intact: compare md5 on Zenodo with the local file
remote=$(api "${auth[@]}" "$host/api/deposit/depositions/$id/files")
bad=0
while read -r name md5; do
  [[ $(json "{f['filename']: f['checksum'] for f in d}.get('$name', 'missing').replace('md5:', '')" <<<"$remote") == "$md5" ]] || { echo "MISMATCH or missing: $name" >&2; bad=1; }
done < <(cd "$dir" && for f in *; do [[ -f $f && $f != *.part* ]] && echo "$f $(md5sum "$f" | cut -d' ' -f1)"; done)
(( bad == 0 )) && echo "All $nfiles files uploaded and verified."

echo "Review the draft, then press Publish: $host/uploads/$id"
echo "After publishing, the record id ($id) goes in the coloc package."
