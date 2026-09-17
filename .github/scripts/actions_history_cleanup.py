#!/usr/bin/env python3
"""Repository-wide Actions retention without relying on GitHub's capped status search."""
import calendar, datetime as dt, json, os, sys, time, urllib.error, urllib.parse, urllib.request

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
            if method=='DELETE' and e.code==404: return None
            if e.code not in (403,429,500,502,503,504) or attempt==4: raise
            time.sleep(3*(attempt+1))

def total(status, created=None, exclude_pr=False):
    p={'per_page':1,'status':status}
    if created:p['created']=created
    if exclude_pr:p['exclude_pull_requests']='true'
    return int(api('/actions/runs',p)['total_count'])

def page_all(status, created=None, exclude_pr=False):
    out=[]; page=1
    while True:
        p={'per_page':100,'page':page,'status':status}
        if created:p['created']=created
        if exclude_pr:p['exclude_pull_requests']='true'
        rows=api('/actions/runs',p)['workflow_runs']; out.extend(rows)
        if len(rows)<100:return out
        page+=1
        if page>10: raise RuntimeError(f'capped shard for {status} {created}')

def iso(x): return x.strftime('%Y-%m-%dT%H:%M:%SZ')
def parse(x): return dt.datetime.fromisoformat(x.replace('Z','+00:00'))

def shard(status,start,end,exclude_pr=False):
    query=f'{iso(start)}..{iso(end)}'; n=total(status,query,exclude_pr)
    if n<=1000:return page_all(status,query,exclude_pr)
    span=end-start
    if span>dt.timedelta(days=1):
        out=[]; cur=start
        while cur<end:
            nxt=min(cur+dt.timedelta(days=1),end); out+=shard(status,cur,nxt-dt.timedelta(seconds=1),exclude_pr); cur=nxt
        return out
    if span>dt.timedelta(hours=1):
        out=[]; cur=start
        while cur<end:
            nxt=min(cur+dt.timedelta(hours=1),end); out+=shard(status,cur,nxt-dt.timedelta(seconds=1),exclude_pr); cur=nxt
        return out
    raise RuntimeError(f'more than 1000 {status} runs in one hour: {query}')

def enumerate_status(status,exclude_pr=False):
    created=parse(api('')['created_at']); now=dt.datetime.now(dt.timezone.utc); cur=created.replace(day=1,hour=0,minute=0,second=0,microsecond=0); out={}
    while cur<=now:
        y,m=cur.year,cur.month; ny,nm=(y+1,1) if m==12 else (y,m+1); nxt=cur.replace(year=ny,month=nm,day=1)
        end=min(nxt-dt.timedelta(seconds=1),now)
        for r in shard(status,cur,end,exclude_pr): out[int(r['id'])]=r
        cur=nxt
    return list(out.values())

def delete(ids):
    if DRY:return 0
    done=0
    for rid in ids:
        api(f'/actions/runs/{rid}',method='DELETE'); done+=1
    return done

def summary(lines):
    p=os.getenv('GITHUB_STEP_SUMMARY')
    if p:
        with open(p,'a') as f:f.write('\n'.join(lines)+'\n')

# Failed retention: page 1 is sufficient to identify the newest 15, but not to
# enumerate old candidates. Enumerate by creation-time shards and verify against
# the API total before deleting anything.
initial=total('failure'); newest=page_all('failure')[:100]
newest=[r for r in newest if r.get('status')=='completed' and r.get('conclusion')=='failure']
newest.sort(key=lambda r:(r['updated_at'],r['created_at'],int(r['id'])),reverse=True)
keep={int(r['id']) for r in newest[:KEEP_FAILED]}
all_failed=[r for r in enumerate_status('failure') if r.get('status')=='completed' and r.get('conclusion')=='failure']
if len(all_failed)!=initial: raise RuntimeError(f'failed shard enumeration mismatch: API={initial} enumerated={len(all_failed)}')
candidates=[r for r in all_failed if int(r['id']) not in keep]
candidates.sort(key=lambda r:(r['updated_at'],r['created_at'],int(r['id'])))
expected=min(max(0,initial-KEEP_FAILED),MAX)
if len(candidates)<expected: raise RuntimeError(f'failed candidate mismatch: expected={expected} available={len(candidates)}')
selected=[int(r['id']) for r in candidates[:expected]]
print(f'Failed workflow runs before cleanup: {initial}; enumerated={len(all_failed)}; selected={len(selected)}')
deleted=delete(selected); final=initial if DRY else total('failure')
summary(['### Failed-run retention','',f'- Failed runs before cleanup: `{initial}`',f'- Shard-enumerated failed runs: `{len(all_failed)}`',f'- Whole failed runs deleted: `{deleted}`',f'- Failed runs after cleanup: `{final}`',f'- Required maximum: `{KEEP_FAILED}`'])
if not DRY and initial<=KEEP_FAILED+MAX and final>KEEP_FAILED: raise RuntimeError(f'retention invariant not met: {final} failures remain')

# Cancelled/skipped retention. Enumerating by shards also removes the old 1,000
# result blind spot here. Preserve newest run for every workflow+branch pair.
cutoff=dt.datetime.now(dt.timezone.utc)-dt.timedelta(minutes=MIN_AGE); pool=[]
for status in ('cancelled','skipped'):
    pool += [r for r in enumerate_status(status,True) if r.get('status')=='completed' and r.get('conclusion') in ('cancelled','skipped')]
pool={int(r['id']):r for r in pool}.values(); ordered=sorted(pool,key=lambda r:(int(r['workflow_id']),r.get('head_branch') or '-',r['updated_at'],int(r['id'])),reverse=True)
protected=set(); seen=set()
for r in ordered:
    key=(int(r['workflow_id']),r.get('head_branch') or '-')
    if key not in seen: seen.add(key); protected.add(int(r['id']))
self_id=int(os.getenv('SELF_RUN_ID','0') or 0)
eligible=[r for r in ordered if int(r['id']) not in protected and int(r['id'])!=self_id and parse(r['updated_at'])<cutoff]
eligible.sort(key=lambda r:(r['updated_at'],int(r['id'])))
selected2=[int(r['id']) for r in eligible[:MAX]]; deleted2=delete(selected2)
print(f'Cancelled/skipped eligible={len(eligible)} selected={len(selected2)} deleted={deleted2}')
summary(['### Cancelled/skipped cleanup','',f'- Eligible old non-retained runs: `{len(eligible)}`',f'- Whole runs deleted: `{deleted2}`'])
if not DRY and len(eligible)<=MAX:
    # Re-enumerate to prove no old non-retained run remains when the budget was sufficient.
    pool2=[]
    for status in ('cancelled','skipped'): pool2 += enumerate_status(status,True)
    pool2={int(r['id']):r for r in pool2 if r.get('status')=='completed' and r.get('conclusion') in ('cancelled','skipped')}.values()
    ordered2=sorted(pool2,key=lambda r:(int(r['workflow_id']),r.get('head_branch') or '-',r['updated_at'],int(r['id'])),reverse=True)
    seen=set(); protected2=set()
    for r in ordered2:
        key=(int(r['workflow_id']),r.get('head_branch') or '-')
        if key not in seen:seen.add(key);protected2.add(int(r['id']))
    remain=[r for r in ordered2 if int(r['id']) not in protected2 and int(r['id'])!=self_id and parse(r['updated_at'])<cutoff]
    if remain: raise RuntimeError(f'cancelled/skipped invariant not met: {len(remain)} eligible runs remain')
