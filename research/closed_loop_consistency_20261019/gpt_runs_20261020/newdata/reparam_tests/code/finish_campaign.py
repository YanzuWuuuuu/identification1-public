from campaign import *
if __name__=='__main__':
    import ctypes
    pid=int(sys.argv[1])
    k32=ctypes.windll.kernel32;k32.OpenProcess.restype=ctypes.c_void_p
    h=k32.OpenProcess(0x100000,False,pid)
    if h:
        while k32.WaitForSingleObject(ctypes.c_void_p(h),1000)==258:pass
        k32.CloseHandle(ctypes.c_void_p(h))
    run('p6_M1_jackknife',['run_p6.py','--jackknife','--methods','M1','--warm-hybrid','--workers','4'])
    run('final_analysis',['analyze.py'])
    run('final_report',['report.py'])
    run('final_audit',['audit.py'])
