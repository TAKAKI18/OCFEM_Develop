"""Small linear-module and preprocessor checks; no production NS calculation."""
from pathlib import Path
import subprocess,json,shutil,os,time,math,uuid
ROOT=Path(__file__).resolve().parents[1]
BIN=ROOT/'build/smoke/bin';RUN=ROOT/'build/smoke'/('checks_'+uuid.uuid4().hex[:8])
SUFFIX='.exe' if os.name=='nt' else ''
def main():
    RUN.mkdir(parents=True)
    env=os.environ.copy();env.setdefault('OMP_NUM_THREADS','2');env.setdefault('MKL_NUM_THREADS','2')
    results=[]
    def run(name,stdin=''):
        exe=BIN/(name+SUFFIX);start=time.perf_counter()
        try:
            proc=subprocess.run([str(exe)],input=stdin,cwd=RUN,capture_output=True,text=True,encoding='utf8',errors='replace',env=env,timeout=60)
            (RUN/(name+'.log')).write_text(proc.stdout+proc.stderr,encoding='utf8')
            results.append(dict(name=name,exit_code=proc.returncode,elapsed_s=time.perf_counter()-start))
            print(name,proc.returncode,flush=True)
            if proc.returncode:raise RuntimeError(proc.stdout[-1000:]+proc.stderr[-1000:])
        except subprocess.TimeoutExpired:
            results.append(dict(name=name,timeout_s=60));raise
    try:
        subprocess.run([os.sys.executable,str(ROOT/'tools/verify_sources.py')],check=True)
        # Collect independent failures without relabeling a failed test as a pass.
        failures=[]
        for name in ['test_sparse_triplet_builder','test_marquardt_pardiso']:
            try:run(name)
            except Exception as e:failures.append(dict(name=name,error=str(e)))
        run('calcollQ_TAKAKI_newton_q_xy','2 2\n0 0\ntime_points.dat\n')
        run('calcollQ_TAKAKI_newton_q_xyz','2 2 2\n0 0\nspace_points.dat\n')
        run('COMV1T_TAKAKI_q','time_points.dat\ntime_operator.dat\n')
        run('COMV3d_TAKAKI_q_Mxyz','space_points.dat\nspace_operator.dat\n')
        run('calweight_TAKAKI_2D_q','2 2\n0 0\nweights.dat\n')
        def numbers(s):return [float(v.replace('D','E').replace('d','e')) for v in s.split()]
        lines=(RUN/'time_operator.dat').read_text().splitlines();n=int(lines[0])+2
        nodes=[numbers(v)[0] for v in lines[1:n+1]]
        entries=[numbers(v) for v in lines[n+1:]]
        assert len(entries)==n*n
        A=[[entries[i*n+j][0] for j in range(n)] for i in range(n)]
        C=[[entries[i*n+j][1] for j in range(n)] for i in range(n)]
        e0=max(abs(sum(row)) for row in A)
        e1=max(abs(sum(A[i][j]*nodes[j] for j in range(n))-1) for i in range(n))
        e2=max(abs(sum(C[i][j]*nodes[j]**2 for j in range(n))-2) for i in range(n))
        assert max(e0,e1,e2)<1e-10,(e0,e1,e2)
        results.append(dict(name='time_operator_polynomial_identities',constant=e0,linear=e1,quadratic=e2,passed=True))
        lines=(RUN/'space_operator.dat').read_text().splitlines();counts=list(map(int,lines[0].split()));offset=1+sum(v+2 for v in counts)
        size=math.prod(v+2 for v in counts);rows=[numbers(v) for v in lines[offset:]]
        assert len(rows)==size*size and all(len(row)==9 for row in rows)
        constant=max(abs(sum(rows[i*size+j][k] for j in range(size))) for i in range(size) for k in range(9))
        assert constant<1e-9,constant
        results.append(dict(name='space_operators_annihilate_constant',nodes=size,max_abs=constant,passed=True))
        lines=(RUN/'weights.dat').read_text().splitlines();n1,n2=map(int,lines[0].split());v=[numbers(row) for row in lines[1:]]
        sums=[sum(row[1] for row in v[:n1]),sum(row[1] for row in v[n1:n1+n2])]
        assert max(abs(value-1) for value in sums)<1e-12,sums
        results.append(dict(name='section_weights_integrate_constant',sums=sums,passed=True))
        report=dict(passed=not failures,scope='linear modules and preprocessing only; NS production not run',run_directory=RUN.relative_to(ROOT).as_posix(),results=results,failures=failures)
    except Exception as e:
        report=dict(passed=False,error=str(e),run_directory=RUN.relative_to(ROOT).as_posix(),results=results)
    (ROOT/'build/smoke/runtime_report.json').write_text(json.dumps(report,indent=2),encoding='utf8')
    print(json.dumps(report),flush=True)
    if not report['passed']:raise SystemExit(1)
if __name__=='__main__':main()
