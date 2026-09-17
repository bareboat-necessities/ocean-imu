#!/usr/bin/env python3
"""Repository-wide Actions retention without relying on GitHub's capped status search."""
import datetime as dt, json, os, time, urllib.error, urllib.parse, urllib.request

TOKEN=os.environ['GH_TOKEN']; REPO=os.environ['REPO']; DRY=os.getenv('DRY_RUN','false').lower()=='true'
MAX=int(os.getenv('MAX_DELETIONS','5000')); KEEP_FAILED=15; MIN_AGE=int(os.getenv('MIN_AGE_MINUTES','60'))
BASE=f'https://api.github.com/repos/{REPO}'
HEAD={'Accept':'application/vnd.github+json','Authorization':f'Bearer {TOKEN}','X-GitHub-Api-Version':'2022-11-28'}

def api(path, params=None, method='GET'):
    url=BASE+path
    if params: url+='?'+urllib.parse.urlencode(params)
    req=urllib.request.Request(url,headers=HEAD,method=method)
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req,timeout=60) as r:
                return json.load(r) if method=='GET' else None
        except urllib.error.HTTPError as e:
            if method=='DELETE' and e.code==404:return None
            if e.code not in (403,429,500,502,503,504) or attempt==4:raise
            time.sleep(3*(attempt+1))

def total(status,created=None,exclude_pr=False):
    p={'per_page':1,'status':status}
    if created:p['created']=created
    if exclude_pr:p['exclude_pull_requests']='true'
    return int(api('/actions/runs',p)['total_count'])

def page_all(status,created=None,exclude_pr=False):
    out=[];page=1
    while True:
        p={'per_page':100,'page':page,'status':status}
        if created:p['created']=created
        if exclude_pr:p['exclude_pull_requests']='true'
        rows=api('/actions/runs',p)['workflow_runs'];out.extend(rows)
        if len(rows)<100:return out
        page+=1
        if page>10:raise RuntimeError(f'capped shard for {status} {created}')

def iso(x):return x.strftime('%Y-%m-%dT%H:%M:%SZ')
def parse(x):return dt.datetime.fromisoformat(x.replace('Z','+00:00'))

def shard(status,start,end,exclude_pr=False):
    query=f'{iso(start)}..{iso(end)}';n=total(status,query,exclude_pr)
    if n<=1000:return page_all(status,query,exclude_pr)
    span=end-start
    unit=dt.timedelta(days=1) if span>dt.timedelta(days=1) else dt.timedelta(hours=1)
    if span<=dt.timedelta(hours=1):raise RuntimeError(f'more than 1000 {status} runs in one hour: {query}')
    out=[];cur=start
    while cur<end:
        nxt=min(cur+unit,end);out+=shard(status,cur,nxt-dt.timedelta(seconds=1),exclude_pr);cur=nxt
    return out

def enumerate_status(status,exclude_pr=False):
    created=parse(api('')['created_at']);now=dt.datetime.now(dt.timezone.utc);cur=created.replace(day=1,hour=0,minute=0,second=0,microsecond=0);out={}
    while cur<=now:
        y,m=cur.year,cur.month;ny,nm=(y+1,1) if m==12 else (y,m+1);nxt=cur.replace(year=ny,month=nm,day=1);end=min(nxt-dt.timedelta(seconds=1),now)
        for r in shard(status,cur,end,exclude_pr):out[int(r['id'])]=r
        cur=nxt
    return list(out.values())

def delete(ids):
    if DRY:return 0
    done=0
    for rid in ids:api(f'/actions/runs/{rid}',method='DELETE');done+=1
    return done

def summary(lines):
    p=os.getenv('GITHUB_STEP_SUMMARY')
    if p:
        with open(p,'a') as f:f.write('\n'.join(lines)+'\n')

def failed_snapshot():
    rows=[r for r in enumerate_status('failure') if r.get('status')=='completed' and r.get('conclusion')=='failure']
    return rows,total('failure')

# GitHub's DELETE endpoint is asynchronous with respect to the Actions run index:
# a successful 204 can remain in status=failure total_count for several seconds.
# Therefore drive retention from fresh sharded snapshots, never from arithmetic
# initial-delete_count, and tolerate stale post-delete reads while making progress.
initial=total('failure');budget=MAX;deleted=0;passes=0;last_count=None;stable=0
while True:
    passes+=1
    all_failed,reported=failed_snapshot()
    actual=len(all_failed)
    if actual!=reported:
        # The search index can be briefly inconsistent immediately after DELETE.
        # Retry rather than failing a valid cleanup.
        if not DRY and deleted and passes<=12:
            print(f'Failed pass {passes}: index settling API={reported} enumerated={actual}')
            time.sleep(min(2*passes,15));continue
        raise RuntimeError(f'failed shard enumeration mismatch: API={reported} enumerated={actual}')
    all_failed.sort(key=lambda r:(r['updated_at'],r['created_at'],int(r['id'])),reverse=True)
    excess=max(0,actual-KEEP_FAILED)
    print(f'Failed pass {passes}: API={reported} enumerated={actual} excess={excess} budget={budget}')
    if excess==0 or DRY:break
    if budget==0:break
    keep={int(r['id']) for r in all_failed[:KEEP_FAILED]}
    candidates=[r for r in reversed(all_failed) if int(r['id']) not in keep]
    selected=[int(r['id']) for r in candidates[:min(excess,budget)]]
    if not selected:raise RuntimeError(f'failed cleanup has {actual} runs but no selectable candidate')
    deleted+=delete(selected);budget-=len(selected)
    # Wait for the Actions search index to reflect successful DELETEs. A stale
    # count is not itself a failure; repeated fresh snapshots decide progress.
    time.sleep(3)
    if last_count==actual:stable+=1
    else:stable=0
    last_count=actual
    if stable>=8:raise RuntimeError(f'failed cleanup index made no progress after repeated DELETEs: {actual} remain')

final_rows,final_reported=failed_snapshot()
final=len(final_rows)
if final!=final_reported and not DRY:
    # Final bounded settle loop avoids the exact 27 -> delete 12 -> stale 19 race.
    for wait in (2,4,8,12,15):
        time.sleep(wait);final_rows,final_reported=failed_snapshot();final=len(final_rows)
        if final==final_reported:break
summary(['### Failed-run retention','',f'- Failed runs before cleanup: `{initial}`',f'- Whole failed runs deleted: `{deleted}`',f'- Cleanup passes: `{passes}`',f'- Failed runs after cleanup: `{final}`',f'- Required maximum: `{KEEP_FAILED}`'])
if not DRY and budget>0 and final>KEEP_FAILED:raise RuntimeError(f'retention invariant not met after settle/requery: {final} failures remain')

# Cancelled/skipped retention. Preserve newest run for every workflow+branch pair.
cutoff=dt.datetime.now(dt.timezone.utc)-dt.timedelta(minutes=MIN_AGE);pool=[]
for status in ('cancelled','skipped'):
    pool += [r for r in enumerate_status(status,True) if r.get('status')=='completed' and r.get('conclusion') in ('cancelled','skipped')]
pool={int(r['id']):r for r in pool}.values();ordered=sorted(pool,key=lambda r:(int(r['workflow_id']),r.get('head_branch') or '-',r['updated_at'],int(r['id'])),reverse=True)
protected=set();seen=set()
for r in ordered:
    key=(int(r['workflow_id']),r.get('head_branch') or '-')
    if key not in seen:seen.add(key);protected.add(int(r['id']))
self_id=int(os.getenv('SELF_RUN_ID','0') or 0)
eligible=[r for r in ordered if int(r['id']) not in protected and int(r['id'])!=self_id and parse(r['updated_at'])<cutoff]
eligible.sort(key=lambda r:(r['updated_at'],int(r['id'])))
selected2=[int(r['id']) for r in eligible[:MAX]];deleted2=delete(selected2)
print(f'Cancelled/skipped eligible={len(eligible)} selected={len(selected2)} deleted={deleted2}')
summary(['### Cancelled/skipped cleanup','',f'- Eligible old non-retained runs: `{len(eligible)}`',f'- Whole runs deleted: `{deleted2}`'])
# Do not immediately assert against the eventually-consistent search index here;
# the next scheduled run is idempotent and will drain any entries still visible.
