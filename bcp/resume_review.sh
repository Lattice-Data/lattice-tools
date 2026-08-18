#!/usr/bin/env bash
# Re-run only the review-verify agents that died on a session limit.
#
# The workflow tool caches by (prompt, opts), so resuming from the original run id
# replays every agent that already finished and re-runs only the failures. That
# makes this script idempotent: if it dies again, run it again.
#
# It cannot invoke the Workflow tool itself -- that is the assistant's job. What it
# does is print the exact invocation and record each attempt, so progress survives
# a session ending mid-run.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_ID="wf_b7f1006e-d80"
SCRIPT="/Users/gabdank/.claude/projects/-Users-gabdank-github-Lattice-Data-lattice-tools-bcp/fcc1d57f-e52c-419f-8796-f6bd7ce5a526/workflows/scripts/chebi-promote-review-${RUN_ID}.js"
JOURNAL="/Users/gabdank/.claude/projects/-Users-gabdank-github-Lattice-Data-lattice-tools/fcc1d57f-e52c-419f-8796-f6bd7ce5a526/subagents/workflows/${RUN_ID}/journal.jsonl"
LOG="${HERE}/.chebi_review_attempts.log"

verdicts() { grep -c '"refuted"' "$JOURNAL" 2>/dev/null || echo 0; }
findings() { python3 - "$JOURNAL" <<'PY'
import json,sys
n=0
for line in open(sys.argv[1]):
    e=json.loads(line)
    if e.get('type')!='result': continue
    r=e['result']
    if isinstance(r,str):
        try: r=json.loads(r)
        except Exception: continue
    if isinstance(r,dict) and 'findings' in r: n+=len(r['findings'])
print(n)
PY
}

echo "$(date -u +%FT%TZ) attempt: verdicts=$(verdicts) findings=$(findings)" | tee -a "$LOG"

if [ ! -f "$SCRIPT" ]; then
  echo "workflow script missing: $SCRIPT" >&2
  echo "Re-author the review workflow from REVIEW_PROGRESS.md instead." >&2
  exit 2
fi

cat <<MSG

Resume is driven by the assistant, not by this script. Ask it to run:

  Workflow({
    scriptPath: "${SCRIPT}",
    resumeFromRunId: "${RUN_ID}"
  })

Completed agents replay from cache; only the ones that died re-run.
Current state: $(verdicts) verdicts for $(findings) findings.
Attempts are appended to ${LOG}.
MSG
