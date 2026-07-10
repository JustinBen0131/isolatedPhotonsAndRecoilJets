#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
if [[ -f "${script_dir}/../../runtime/io/recoiljets_io_paths.sh" ]]; then
  source "${script_dir}/../../runtime/io/recoiljets_io_paths.sh"
elif [[ -f scripts/sdcc/runtime/io/recoiljets_io_paths.sh ]]; then
  source scripts/sdcc/runtime/io/recoiljets_io_paths.sh
else
  echo "[ERROR] Cannot locate scripts/sdcc/runtime/io/recoiljets_io_paths.sh" >&2
  exit 2
fi

repo_root="$(rj_find_repo_root "$script_dir")"
cd "$repo_root"

say() { printf '[THE96-HB] %s\n' "$*"; }
die() { printf '[THE96-HB][ERROR] %s\n' "$*" >&2; exit "${2:-2}"; }

mode="${1:-status}"
campaign_tag="${2:-${RJ_THE96_CAMPAIGN_TAG:-}}"
[[ -n "$campaign_tag" ]] || die "campaign tag required: $0 <status|watch> <campaign_tag>"
rj_validate_campaign_tag "$campaign_tag" || exit 2

config_dir="${RJ_THE96_CONFIG_DIR:-${repo_root}/condor_generated_configs/${campaign_tag}}"
manifest="${RJ_THE96_MANIFEST:-${config_dir}/the96_campaign_manifest.tsv}"
interval="${RJ_THE96_HEARTBEAT_SECONDS:-900}"
ready_marker="${config_dir}/the96_ready_to_plot.txt"

[[ -f "$manifest" ]] || die "manifest not found: ${manifest}"

queue_snapshot() {
  condor_q "${USER:-patsfan753}" -af ClusterId JobStatus Args 2>/dev/null | grep -F "$campaign_tag" || true
}

queue_counts() {
  queue_snapshot | awk '
    BEGIN { total=0; idle=0; running=0; held=0; other=0; }
    NF >= 2 {
      total++;
      if ($2 == 1) idle++;
      else if ($2 == 2) running++;
      else if ($2 == 5) held++;
      else other++;
    }
    END { printf "total=%d idle=%d running=%d held=%d other=%d\n", total, idle, running, held, other; }
  '
}

held_count() {
  queue_snapshot | awk 'NF >= 2 && $2 == 5 { n++ } END { print n+0 }'
}

total_count() {
  queue_snapshot | awk 'NF >= 2 { n++ } END { print n+0 }'
}

expected_final_roots() {
  tail -n +2 "$manifest" | awk -F '\t' 'NF >= 8 { n++ } END { print n+0 }'
}

final_roots() {
  awk -F '\t' 'NR > 1 && NF >= 8 { print $6 }' "$manifest" | sort -u | while IFS= read -r merge_base; do
    [[ -n "$merge_base" && -d "$merge_base" ]] || continue
    find "$merge_base" -type f -name 'RecoilJets_*_MERGED.root' -size +50000c -print 2>/dev/null
  done | sort -u
}

final_root_count() {
  final_roots | awk 'END { print NR+0 }'
}

print_status() {
  local now counts held total expected found
  now="$(date '+%Y-%m-%dT%H:%M:%S%z')"
  counts="$(queue_counts)"
  held="$(held_count)"
  total="$(total_count)"
  expected="$(expected_final_roots)"
  found="$(final_root_count)"
  say "timestamp=${now}"
  say "campaign_tag=${campaign_tag}"
  say "manifest=${manifest}"
  say "queue ${counts}"
  say "final_roots=${found}/${expected}"
  if (( found > 0 )); then
    final_roots | sed 's/^/[THE96-HB] final_root=/'
  fi
  if (( held > 0 )); then
    say "held job rows detected; heartbeat will not run job-control automatically"
    queue_snapshot | awk 'NF >= 2 && $2 == 5 { print "[THE96-HB] held_row="$0 }' | head -20
  fi
  if (( total == 0 && found >= expected && expected > 0 )); then
    {
      printf 'THE96_READY_TO_PLOT_V1\n'
      printf 'timestamp=%s\n' "$now"
      printf 'campaign_tag=%s\n' "$campaign_tag"
      printf 'manifest=%s\n' "$manifest"
      printf 'final_roots=%s/%s\n' "$found" "$expected"
    } > "$ready_marker"
    say "ready_marker=${ready_marker}"
    return 0
  fi
  return 1
}

case "$mode" in
  status)
    print_status || true
    ;;
  watch)
    say "watching campaign ${campaign_tag} every ${interval}s"
    while true; do
      if print_status; then
        say "READY_TO_PLOT"
        exit 0
      fi
      if (( "$(held_count)" > 0 )); then
        die "held jobs present; stop for read-only diagnosis and explicit recovery approval" 20
      fi
      sleep "$interval"
    done
    ;;
  *)
    die "unknown mode: ${mode}; expected status or watch"
    ;;
esac
