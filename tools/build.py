"""Build the reviewed source list with Intel Fortran and oneMKL; no calculations start here."""
from pathlib import Path
import argparse,subprocess,os,shutil,json,re,sys
ROOT=Path(__file__).resolve().parents[1]
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mode',choices=['release','debug','smoke'],default='release')
    p.add_argument('--target',choices=['runtime','all','pardiso','stokes','preprocess','tests'],default='runtime')
    a=p.parse_args()
    fc=shutil.which(os.environ.get('FC','ifx')) or shutil.which('ifort')
    if not fc:raise SystemExit('Intel Fortran not found. Load the oneAPI compiler and MKL environment first.')
    windows=os.name=='nt';suffix='.exe' if windows else ''
    dest=ROOT/'build'/a.mode;bindir=dest/'bin';bindir.mkdir(parents=True,exist_ok=True)
    common=ROOT/'src/common/sparse_triplet_builder.f90'
    pardiso=ROOT/'src/common/marquardt_pardiso_solver.f90'
    targets=[]
    if a.target in ['runtime','all','pardiso']:
        targets.append(('unsteady_ns_marquardt_pardiso',[common,pardiso,*sorted((ROOT/'src/solvers').glob('*.f90'))]))
    if a.target in ['all','stokes']:
        targets.append(('stokes_initial',[common,*sorted((ROOT/'src/initial').glob('*.f90'))]))
    if a.target in ['runtime','all','preprocess']:
        targets.extend((s.stem,[s]) for s in sorted((ROOT/'src/preprocess').glob('*.f90')))
    if a.target in ['all','tests']:
        for s in sorted((ROOT/'tests').glob('test_*.f90')):
            targets.append((s.stem,[common,*([pardiso] if 'pardiso' in s.stem else []),s]))
    report=[]
    for name,sources in targets:
        obj=dest/'obj'/name;obj.mkdir(parents=True,exist_ok=True)
        local=[]
        for src in sources:
            if a.mode=='smoke' and src.parent.name in ['solvers','initial']:
                text=src.read_text(encoding='utf8')
                text,count=re.subn(r'\bKNAA\s*=\s*10000\b','KNAA=256',text,flags=re.I)
                if count==0:raise RuntimeError('Small-test allocation adaptation did not match '+str(src))
                out=obj/src.name;out.write_text(text,encoding='utf8');local.append(out)
            else:local.append(src)
        exe=bindir/(name+suffix)
        if windows:
            flags=['/O2'] if a.mode=='release' else ['/Od','/check:bounds']
            cmd=[fc,*flags,'/fp:strict','/Qprotect-parens','/Qfma-','/Qopenmp','/Qmkl','/heap-arrays','/F536870912','/traceback','/module:'+str(obj),'/exe:'+str(exe),*map(str,local)]
        else:
            flags=['-O2'] if a.mode=='release' else ['-O0','-check','bounds']
            cmd=[fc,*flags,'-fp-model','strict','-assume','protect_parens','-no-fma','-qopenmp','-qmkl','-heap-arrays','-traceback','-module',str(obj),'-o',str(exe),*map(str,local)]
        log=dest/(name+'_build.log')
        with log.open('w',encoding='utf8') as f:
            f.write('Compiler: '+fc+'\nMode: '+a.mode+'\nSources:\n'+'\n'.join(map(str,sources))+'\n');f.flush()
            result=subprocess.run(cmd,cwd=obj,stdout=f,stderr=subprocess.STDOUT)
        report.append(dict(target=name,exit_code=result.returncode,executable=exe.relative_to(ROOT).as_posix(),log=log.relative_to(ROOT).as_posix()))
        print(('BUILT ' if result.returncode==0 else 'FAILED ')+name,flush=True)
        if result.returncode:
            print(log.read_text(encoding='utf8',errors='replace')[-8000:]);break
    (dest/'build_report.json').write_text(json.dumps(dict(mode=a.mode,compiler=fc,results=report),indent=2),encoding='utf8')
    if not report or any(r['exit_code'] for r in report):raise SystemExit(1)
if __name__=='__main__':main()
