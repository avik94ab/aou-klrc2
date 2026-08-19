#!/bin/bash
# Re-emit the paste-able transfer for nkc_extract.py.
#
# The Workbench is a controlled environment with no route back to Wynton, and a
# single file does not justify a bucket round trip: 20 kB of base64 pastes into
# a terminal in one go and carries its own checksum.
#
#   bash make_paste.sh > nkc_extract.paste.sh
set -euo pipefail
cd "$(dirname "$0")"
src="${1:-nkc_extract.py}"
md5=$(md5sum "$src" | cut -d' ' -f1)
printf '# Generated from %s (md5 %s).\n' "$src" "$md5"
printf '# Paste this whole file into a Workbench terminal to recreate the script.\n'
printf '# Regenerate after editing %s:\n' "$src"
printf '#   bash make_paste.sh > %s.paste.sh\n' "${src%.py}"
printf 'base64 -d <<%s | gunzip > %s\n' "'PAYLOAD'" "$src"
gzip -c9 "$src" | base64 -w 100
printf 'PAYLOAD\n'
printf 'md5sum %s   # expect %s\n' "$src" "$md5"
