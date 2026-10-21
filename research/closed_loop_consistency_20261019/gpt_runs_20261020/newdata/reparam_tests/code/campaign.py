"""Resume-safe sequential campaign stages; workers write each endpoint separately."""
from profile_core import *
import subprocess

def run(name,args):
    with (ROOT/'logs'/f'{name}.log').open('a',encoding='utf8') as out, (ROOT/'logs'/f'{name}.err').open('a',encoding='utf8') as err:
        p=subprocess.run([sys.executable,str(ROOT/'code'/args[0]),*args[1:]],stdout=out,stderr=err,cwd=ROOT)
    js(ROOT/'analysis'/f'stage_{name}.json',dict(returncode=p.returncode,command=args,finished=time.strftime('%Y-%m-%d %H:%M:%S')))
    return p.returncode

if __name__=='__main__':
    # The initial full P6 pool may already be active. Do not overlap two fitting pools.
    import ctypes
    pid=int(sys.argv[1]) if len(sys.argv)>1 else None
    if pid:
        k32=ctypes.windll.kernel32
        k32.OpenProcess.restype=ctypes.c_void_p
        handle=k32.OpenProcess(0x100000,False,pid)
        if handle:
            while k32.WaitForSingleObject(ctypes.c_void_p(handle),1000)==258:pass
            k32.CloseHandle(ctypes.c_void_p(handle))
    run('p6_full_resume',['run_p6.py','--workers','4'])
    run('p6_metrics',['metrics.py'])
    run('p7_real',['run_p7_real.py'])
    run('p7_sim',['run_p7_sim.py','--workers','4'])
    run('p6_jackknife',['run_p6.py','--jackknife','--methods','M3','--workers','4'])
    run('p6_fixed_low',['run_fixed_low.py'])
    # Additional analyses/report are separate commands, so partial work is inspectable.
