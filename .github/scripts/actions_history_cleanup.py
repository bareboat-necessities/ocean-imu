#!/usr/bin/env python3
"""Repository-wide Actions retention without relying on capped status search."""
import datetime as dt
import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request

TOKEN = os.environ["GH_TOKEN"]
REPO = os.environ["REPO"]
DRY = os.getenv("DRY_RUN", "false").lower() == "true"
MAX = int(os.getenv("MAX_DELETIONS", "20000"))
BATCH = max(1, int(os.getenv("DELETE_BATCH_SIZE", "250")))
KEEP_FAILED = 15
MIN_AGE = int(os.getenv("MIN_AGE_MINUTES", "60"))
BASE = f"https://api.github.com/repos/{REPO}"
HEAD = {
    "Accept": "application/vnd.github+json",
    "Authorization": f"Bearer {TOKEN}",
    "X-GitHub-Api-Version": "2022-11-28",
}


def api(path, params=None, method="GET"):
    url = BASE + path
    if params:
        url += "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers=HEAD, method=method)
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req, timeout=60) as response:
                return json.load(response) if method == "GET" else None
        except urllib.error.HTTPError as exc:
            if method == "DELETE" and exc.code == 404:
                return None
            if exc.code not in (403, 429, 500, 502, 503, 504) or attempt == 4:
                raise
            time.sleep(3 * (attempt + 1))
    return None


def total(status, created=None):
    params = {"per_page": 1, "status": status}
    if created:
        params["created"] = created
    return int(api("/actions/runs", params)["total_count"])


def page_all(status, created=None):
    out = []
    page = 1
    while True:
        params = {"per_page": 100, "page": page, "status": status}
        if created:
            params["created"] = created
        rows = api("/actions/runs", params)["workflow_runs"]
        out.extend(rows)
        if len(rows) < 100:
            return out
        page += 1
        if page > 10:
            raise RuntimeError(f"capped shard for {status} {created}")


def iso(value):
    return value.strftime("%Y-%m-%dT%H:%M:%SZ")


def parse(value):
    return dt.datetime.fromisoformat(value.replace("Z", "+00:00"))


def shard(status, start, end):
    query = f"{iso(start)}..{iso(end)}"
    count = total(status, query)
    if count <= 1000:
        return page_all(status, query)
    span = end - start
    unit = dt.timedelta(days=1) if span > dt.timedelta(days=1) else dt.timedelta(hours=1)
    if span <= dt.timedelta(hours=1):
        raise RuntimeError(f"more than 1000 {status} runs in one hour: {query}")
    out = []
    cur = start
    while cur < end:
        nxt = min(cur + unit, end)
        out += shard(status, cur, nxt - dt.timedelta(seconds=1))
        cur = nxt
    return out


def enumerate_status(status):
    created = parse(api("")["created_at"])
    now = dt.datetime.now(dt.timezone.utc)
    cur = created.replace(day=1, hour=0, minute=0, second=0, microsecond=0)
    out = {}
    while cur <= now:
        year, month = cur.year, cur.month
        next_year, next_month = (year + 1, 1) if month == 12 else (year, month + 1)
        nxt = cur.replace(year=next_year, month=next_month, day=1)
        end = min(nxt - dt.timedelta(seconds=1), now)
        for run in shard(status, cur, end):
            out[int(run["id"])] = run
        cur = nxt
    return list(out.values())


def delete(ids, label):
    if DRY:
        return 0
    done = 0
    ids = list(ids)
    for offset in range(0, len(ids), BATCH):
        batch = ids[offset : offset + BATCH]
        print(
            f"{label}: deleting batch {offset // BATCH + 1} "
            f"size={len(batch)} progress={offset}/{len(ids)}"
        )
        for run_id in batch:
            api(f"/actions/runs/{run_id}", method="DELETE")
            done += 1
        # Give the Actions index and abuse/rate limiter a short boundary between
        # batches. Individual DELETE retries are handled by api().
        if offset + BATCH < len(ids):
            time.sleep(2)
    return done


def summary(lines):
    path = os.getenv("GITHUB_STEP_SUMMARY")
    if path:
        with open(path, "a") as handle:
            handle.write("\n".join(lines) + "\n")


def failed_snapshot():
    rows = [
        run
        for run in enumerate_status("failure")
        if run.get("status") == "completed" and run.get("conclusion") == "failure"
    ]
    return rows, total("failure")


# Failed-run retention is deliberately separate because we keep the newest 15
# repository-wide. Fresh sharded snapshots handle the eventually-consistent
# Actions index after DELETE.
initial = total("failure")
budget = MAX
deleted = 0
passes = 0
last_count = None
stable = 0
while True:
    passes += 1
    all_failed, reported = failed_snapshot()
    actual = len(all_failed)
    if actual != reported:
        if not DRY and deleted and passes <= 12:
            print(f"Failed pass {passes}: index settling API={reported} enumerated={actual}")
            time.sleep(min(2 * passes, 15))
            continue
        raise RuntimeError(
            f"failed shard enumeration mismatch: API={reported} enumerated={actual}"
        )
    all_failed.sort(
        key=lambda run: (run["updated_at"], run["created_at"], int(run["id"])),
        reverse=True,
    )
    excess = max(0, actual - KEEP_FAILED)
    print(
        f"Failed pass {passes}: API={reported} enumerated={actual} "
        f"excess={excess} budget={budget}"
    )
    if excess == 0 or DRY or budget == 0:
        break
    keep = {int(run["id"]) for run in all_failed[:KEEP_FAILED]}
    candidates = [run for run in reversed(all_failed) if int(run["id"]) not in keep]
    selected = [int(run["id"]) for run in candidates[: min(excess, budget)]]
    if not selected:
        raise RuntimeError(f"failed cleanup has {actual} runs but no selectable candidate")
    deleted += delete(selected, "failed")
    budget -= len(selected)
    time.sleep(3)
    if last_count == actual:
        stable += 1
    else:
        stable = 0
    last_count = actual
    if stable >= 8:
        raise RuntimeError(
            f"failed cleanup index made no progress after repeated DELETEs: {actual} remain"
        )

final_rows, final_reported = failed_snapshot()
final = len(final_rows)
if final != final_reported and not DRY:
    for wait in (2, 4, 8, 12, 15):
        time.sleep(wait)
        final_rows, final_reported = failed_snapshot()
        final = len(final_rows)
        if final == final_reported:
            break
summary(
    [
        "### Failed-run retention",
        "",
        f"- Failed runs before cleanup: `{initial}`",
        f"- Whole failed runs deleted: `{deleted}`",
        f"- Cleanup passes: `{passes}`",
        f"- Failed runs after cleanup: `{final}`",
        f"- Required maximum: `{KEEP_FAILED}`",
    ]
)
if not DRY and budget > 0 and final > KEEP_FAILED:
    raise RuntimeError(
        f"retention invariant not met after settle/requery: {final} failures remain"
    )

# Cancelled/skipped retention. The previous implementation passed
# exclude_pull_requests=true, which excluded exactly the large population this
# repository accumulates when superseded PR checks are cancelled. Enumerate ALL
# cancelled/skipped runs, preserve the newest run for each workflow+branch pair,
# and drain the old remainder in explicit DELETE batches.
cutoff = dt.datetime.now(dt.timezone.utc) - dt.timedelta(minutes=MIN_AGE)
pool = []
for conclusion in ("cancelled", "skipped"):
    pool += [
        run
        for run in enumerate_status(conclusion)
        if run.get("status") == "completed"
        and run.get("conclusion") in ("cancelled", "skipped")
    ]
pool = {int(run["id"]): run for run in pool}.values()
ordered = sorted(
    pool,
    key=lambda run: (
        int(run["workflow_id"]),
        run.get("head_branch") or "-",
        run["updated_at"],
        int(run["id"]),
    ),
    reverse=True,
)
protected = set()
seen = set()
for run in ordered:
    key = (int(run["workflow_id"]), run.get("head_branch") or "-")
    if key not in seen:
        seen.add(key)
        protected.add(int(run["id"]))
self_id = int(os.getenv("SELF_RUN_ID", "0") or 0)
eligible = [
    run
    for run in ordered
    if int(run["id"]) not in protected
    and int(run["id"]) != self_id
    and parse(run["updated_at"]) < cutoff
]
eligible.sort(key=lambda run: (run["updated_at"], int(run["id"])))
selected_cancelled = [int(run["id"]) for run in eligible[:MAX]]
deleted_cancelled = delete(selected_cancelled, "cancelled/skipped")
print(
    f"Cancelled/skipped eligible={len(eligible)} selected={len(selected_cancelled)} "
    f"deleted={deleted_cancelled} protected={len(protected)}"
)
summary(
    [
        "### Cancelled/skipped cleanup",
        "",
        f"- Eligible old non-retained runs: `{len(eligible)}`",
        f"- Protected newest workflow/branch runs: `{len(protected)}`",
        f"- Delete batch size: `{BATCH}`",
        f"- Whole runs deleted: `{deleted_cancelled}`",
    ]
)
# No immediate count assertion: DELETE is eventually consistent. A subsequent
# scheduled/manual cleanup is idempotent and will drain anything still visible.
