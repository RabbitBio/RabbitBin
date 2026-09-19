"""Compare full-file and byte-sharded BAM depth, graph and bin outputs.

Usage: test_bam_depth_equivalence.py CURRENT WORK_DIR SAMTOOLS [BASELINE]
The optional separate baseline can check changes across RabbitBin versions.
"""

import csv, hashlib, json, os, random, subprocess, sys
from pathlib import Path
CURRENT = Path(sys.argv[1]).resolve()
WORK = Path(sys.argv[2]).resolve()
SAMTOOLS = sys.argv[3]
BASELINE = Path(sys.argv[4]).resolve() if len(sys.argv) > 4 else CURRENT
WORK.mkdir(parents=True, exist_ok=True)
rng = random.Random(20260919)
refs = [(f'ctg_{i}', (700, 1400, 3200, 4600)[i % 4]) for i in range(96)]
seqs = {name: ''.join(rng.choice('ACGT') for _ in range(length)) for name,length in refs}
def fasta(path, rows):
    path.write_text(''.join(f'>{name}\n{sequence}\n' for name,sequence in rows))
ordered = WORK/'ordered.fa'
fasta(ordered, [(n,seqs[n]) for n,_ in refs])
reverse = WORK/'reverse.fa'
fasta(reverse, [(n,seqs[n]) for n,_ in refs[::-1]])
subset = WORK/'subset.fa'
fasta(subset, [(n,seqs[n]) for i,(n,_) in enumerate(refs) if i%11])
duplicate = WORK/'duplicate.fa'
fasta(duplicate, [(n,seqs[n]) for n,_ in refs] + [(refs[-1][0], seqs[refs[-1][0]])])
cross = WORK/'cross-class.fa'
fasta(cross, [(n,seqs[n]) for n,_ in refs] + [(refs[-1][0], seqs[refs[-1][0]][:1400])])
bams=[]
patterns=[('100M',0),('10S90M',1),('40M2I58M',3),('50M2D50M',3),('40=2X58=',2),('50M10N50M',0),('100M',8)]
for sample in range(3):
    sam=WORK/f's{sample}.sam'; bam=WORK/f's{sample}.bam'
    with sam.open('w') as f:
        f.write('@HD\tVN:1.6\tSO:coordinate\n')
        for name,length in refs: f.write(f'@SQ\tSN:{name}\tLN:{length}\n')
        for ci,(name,length) in enumerate(refs):
            if ci%13==0: continue
            reads=120+((ci*(sample+3))%13)*30
            positions=sorted(rng.randrange(1,length-115) for _ in range(reads))
            for ri,pos in enumerate(positions):
                cigar,nm=patterns[(ri+sample)%len(patterns)]
                mq=(0,3,10,30,60)[(ri+ci+sample)%5]
                seq=''.join(rng.choice('ACGT') for _ in range(100))
                f.write(f's{sample}_{ci}_{ri}\t0\t{name}\t{pos}\t{mq}\t{cigar}\t*\t0\t0\t{seq}\t'+('I'*100)+f'\tNM:i:{nm}\n')
            # A repaired right-edge overrun exercises invalidation of cached ends.
            f.write(f'edge_{sample}_{ci}\t0\t{name}\t{length-40}\t60\t100M\t*\t0\t0\t'+('A'*100)+'\t'+('I'*100)+'\tNM:i:0\n')
        for ri in range(6):
            f.write(f'unmapped_{sample}_{ri}\t4\t*\t0\t0\t*\t*\t0\t0\t'+('A'*100)+'\t'+('I'*100)+'\n')
    subprocess.run([SAMTOOLS,'view','-b','-o',str(bam),str(sam)],check=True)
    bams.append(bam)
    sam.unlink()
listing=WORK/'bams.list'; listing.write_text(''.join(str(p)+'\n' for p in bams))
env={k:v for k,v in os.environ.items() if not k.startswith(('RB_','RABBIT_','OMP_'))}
cases=[('one',ordered,1,5,1),('boundaries',ordered,31,5,1),('compact',ordered,31,5,1000),('single-channel',ordered,31,0,1000),('reordered',reverse,17,5,1000),('extra-bam-names',subset,17,5,1000),('duplicate-fasta',duplicate,17,5,1000),('cross-class-name',cross,17,5,1000)]
results=[]
for label,assembly,shards,dual,minlen in cases:
    outputs=[]
    for version,binary,route in [('baseline',BASELINE,'0'),('sharded',CURRENT,'1')]:
        prefix=WORK/f'{label}-{version}'; dump=Path(str(prefix)+'.depth.tsv')
        runenv=dict(env,RABBIT_DEPTH_B2=route,RABBIT_DEPTH_B2_SHARDS=str(shards),RABBIT_FUSE_DUMP_DEPTH=str(dump))
        cmd=[str(binary),'bin','--fasta',str(assembly),'--bam-list',str(listing),'--output',str(prefix),'--threads','4','--seed','42','--min-bin-size','0','--dual-depth',str(dual),'--min-contig-length',str(minlen)]
        with Path(str(prefix)+'.log').open('w') as log:
            subprocess.run(cmd,env=runenv,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=90)
        with Path(str(prefix)+'.bins.tsv').open() as f:
            bins=[{k:v for k,v in row.items() if k!='FileName'} for row in csv.DictReader(f,delimiter='\t')]
        outputs.append((dump.read_bytes(),Path(str(prefix)+'.members.tsv').read_bytes(),bins))
    if outputs[0]!=outputs[1]:
        differences=[kind for kind,a,b in zip(('depth','members','bins'),outputs[0],outputs[1]) if a!=b]
        raise AssertionError((label,differences))
    record=dict(case=label,depth_identical=True,members_identical=True,bin_stats_identical=True,shards_requested=shards,dual=dual,min_length=minlen,depth_sha256=hashlib.sha256(outputs[1][0]).hexdigest())
    results.append(record)
    print('PASS',label,flush=True)
(WORK/'results.json').write_text(json.dumps(results,indent=2)+'\n')
print(f'{len(results)} BAM matrix/member/bin comparisons passed',flush=True)

# A resynchronization failure must discard partial writes before rescanning.
prefix=WORK/'forced-tail-fallback'
env={k:v for k,v in os.environ.items() if not k.startswith(('RB_','RABBIT_','OMP_'))}
env.update(RABBIT_DEPTH_B2='1',RABBIT_DEPTH_B2_SHARDS='4096',RABBIT_FUSE_DUMP_DEPTH=str(prefix)+'.depth.tsv')
cmd=[str(CURRENT),'bin','--fasta',str(WORK/'ordered.fa'),'--bam-list',str(WORK/'bams.list'),'--output',str(prefix),'--threads','4','--seed','42','--min-bin-size','0','--dual-depth','5','--min-contig-length','1000']
with Path(str(prefix)+'.log').open('w') as f:
 subprocess.run(cmd,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=90,check=True)
log=Path(str(prefix)+'.log').read_text()
if 'B2 re-sync failed on' not in log: raise AssertionError('Fallback was not exercised')
base=WORK/'compact-baseline'
for ext in ('.depth.tsv','.members.tsv'):
 if Path(str(base)+ext).read_bytes()!=Path(str(prefix)+ext).read_bytes(): raise AssertionError(ext)
def bins(p):
 with Path(str(p)+'.bins.tsv').open() as f:
  return [{k:v for k,v in r.items() if k!='FileName'} for r in csv.DictReader(f,delimiter='\t')]
if bins(base)!=bins(prefix): raise AssertionError('bin statistics')
result=dict(fallback_exercised=True,depth_identical=True,members_identical=True,bin_stats_identical=True)
Path(str(prefix)+'.json').write_text(json.dumps(result,indent=2)+'\n')
print('PASS fallback after partial interior writes: depth/member/bin identical')
