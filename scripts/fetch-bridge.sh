#!/bin/bash
# Download one BridgeDb .bridge file from its (Figshare) URL.
#
#   scripts/fetch-bridge.sh <downloadURL> <target path>
#
# The URL comes from the bridgedb/data gene.json `downloadURL` field, which is
# an opaque Figshare file id -- it cannot be constructed from the species name.
#
# Two Figshare behaviours to know about:
#
#   - It redirects to a presigned S3 URL that expires in ten seconds, so the
#     redirect has to be followed immediately (-L) rather than re-requested.
#
#   - It answers "HTTP 202 Accepted" with an empty body when it decides the
#     caller is an interactive browser that should poll while the file is
#     prepared. Sending a full Chrome user-agent string reliably triggers this;
#     curl's own user-agent gets a plain 302 -> 200. So this script deliberately
#     does NOT spoof a browser. 202 is still treated as a retryable failure
#     below, because it is a success code: a client that only checks for errors
#     accepts it and silently writes a zero-byte file. (That is precisely how
#     the previous wget-based download failed -- wget exits 0 on 202, and
#     --retry-on-http-error cannot match a 2xx code.)
#
# A download is accepted only on an explicit 200 that also passes a zip test.
# Exits non-zero if no valid file could be fetched.
set -euo pipefail

url="${1:?usage: fetch-bridge.sh <downloadURL> <target>}"
target="${2:?usage: fetch-bridge.sh <downloadURL> <target>}"

attempts=10

# figshare.com/ndownloader redirects; the ndownloader host serves directly.
norm="${url/https:\/\/figshare.com\/ndownloader/https:\/\/ndownloader.figshare.com}"

if [ -s "$target" ] && unzip -t "$target" >/dev/null 2>&1; then
  echo "OK: $target already present and valid"
  exit 0
fi

mkdir -p "$(dirname "$target")"
tmp="$(mktemp)"
trap 'rm -f "$tmp"' EXIT

delay=2
i=1
while :; do
  echo "fetching $norm -> $target (attempt $i/$attempts)"
  # --speed-limit/--speed-time abandon a transfer that drops below 10 kB/s for
  # a minute, so a stalled connection is retried in ~1 minute rather than
  # sitting until --max-time. The cap is generous because these files reach
  # ~840 MB and Figshare is sometimes genuinely slow rather than stuck.
  http=$(curl -sS -L --max-time 1800 \
              --speed-limit 10240 --speed-time 60 \
              --retry 3 --retry-delay 2 --retry-connrefused \
              -o "$tmp" -w '%{http_code}' "$norm") || http="000"

  if [ "$http" = "200" ] && unzip -t "$tmp" >/dev/null 2>&1; then
    mv -f "$tmp" "$target"
    echo "saved $target ($(du -h "$target" | cut -f1))"
    exit 0
  fi

  if [ "$http" = "202" ]; then
    echo "  HTTP 202 with no body: Figshare wants the caller to poll"
  else
    echo "  HTTP $http, or the download failed its zip test"
  fi

  if [ "$i" -ge "$attempts" ]; then
    echo "ERROR: no valid $target after $attempts attempts" >&2
    exit 1
  fi
  sleep "$delay"
  delay=$((delay * 2))
  i=$((i + 1))
done
