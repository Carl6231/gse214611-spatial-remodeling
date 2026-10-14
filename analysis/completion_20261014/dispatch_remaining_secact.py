"""Resume only missing sample computations after the previous process has exited."""
from pathlib import Path
import csv, os, shutil, subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
root=Path('/home/zhenhanli/myocardial/scop/GSE214611_spatial')
base=root/'22_JTM_completion_20261014'; prior=base/'P06_SECACT_compat2'; out=base/'P06_SECACT_parallel'
def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
status=read(prior/'functions_final_status.tsv')
pilot=read(prior/'tables/pilot_output_audit.tsv')
assert len(pilot)==17 and all(x['pass']=='TRUE' for x in pilot)
# Do not replace or terminate a live analysis.
for proc in Path('/proc').iterdir():
 if not proc.name.isdigit():continue
 try:cmd=(proc/'cmdline').read_bytes().split(b'\0')
 except (OSError,PermissionError):continue
 if any(b'/usr/lib/R/bin/exec/R'==a for a in cmd) and any(b'P06_run_resumable.R' in a for a in cmd):
  raise RuntimeError('Original R process is still running; no dispatch allowed')
mem={x.split(':')[0]:int(x.split()[1]) for x in Path('/proc/meminfo').read_text().splitlines()}
assert mem['MemAvailable']>80*1024*1024, 'Need at least 80 GiB available before a four-process batch'
ids=['SRR21777596','SRR21777597','SRR21777598','SRR21777581','SRR21777582','SRR21777599']+[f'MI_D3_rep{i}' for i in range(1,4)]+[f'MI_D7_rep{i}' for i in range(1,4)]
complete=[];missing=[]
for sid in ids:
 rows=[r for r in status if r['sample_id']==sid]
 expected=2 if sid.startswith('SRR') else 15
 if len(rows)==expected and all(r['status'] in ['PASS','NEGATIVE_COMPLETE','NOT_APPLICABLE'] for r in rows):complete.append(sid)
 else:missing.append(sid)
print('Reusing completed samples:',','.join(complete),flush=True)
print('Missing samples:',','.join(missing),flush=True)
assert set(missing)=={'SRR21777599','MI_D3_rep2','MI_D3_rep3','MI_D7_rep1','MI_D7_rep2','MI_D7_rep3'}
# The sole prior incomplete CCC sample is handled first to validate the empty-table fix.
def run(sid):
 job=out/'jobs'/sid
 for name in ['script','logs','tables','objects','figures']:(job/name).mkdir(parents=True,exist_ok=True)
 assert not (job/'STATUS.tsv').exists(),'Refuse to silently overwrite an attempted job'
 for name in ['P06_run_sample.R','secact_compat.R']:shutil.copy2(out/'script'/name,job/'script'/name)
 # Reuse only producer outputs explicitly marked PASS in the prior saved status.
 for r in status:
  if r['sample_id']==sid and r['status']=='PASS' and r['evidence_path'].endswith('.rds'):
   p=Path(r['evidence_path'])
   if p.is_file():shutil.copy2(p,job/'objects'/p.name)
 env=dict(os.environ,SECACT_JOB_ID=sid,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',NUMEXPR_NUM_THREADS='1')
 with (job/'logs/run.log').open('wb') as log:
  rc=subprocess.call(['Rscript',str(job/'script/P06_run_sample.R'),str(root),str(job),str(root/'21_scop_expansion/_runtime/spacet/R-4.5/site-library')],env=env,stdout=log,stderr=subprocess.STDOUT)
 assert rc==0 and (job/'_SUCCESS').exists(),f'{sid} failed; see its retained log'
 print('Completed:',sid,flush=True)
run('SRR21777599')
with ThreadPoolExecutor(max_workers=4) as pool:
 jobs={pool.submit(run,s):s for s in missing if s!='SRR21777599'}
 for future in as_completed(jobs):future.result()
print('ALL_REMAINING_JOBS_COMPLETE',flush=True)
