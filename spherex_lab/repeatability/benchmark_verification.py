"""Scientific equivalence and process-memory helpers for the full workload benchmark."""
import math
import os
from pathlib import Path
import json


def peak_memory_bytes():
    """Peak resident process memory, including native arrays (not just Python heap)."""
    if os.name == 'nt':
        import ctypes
        from ctypes import wintypes
        class Counters(ctypes.Structure):
            _fields_ = [('cb', wintypes.DWORD), ('PageFaultCount', wintypes.DWORD)] + [(name, ctypes.c_size_t) for name in ('PeakWorkingSetSize','WorkingSetSize','QuotaPeakPagedPoolUsage','QuotaPagedPoolUsage','QuotaPeakNonPagedPoolUsage','QuotaNonPagedPoolUsage','PagefileUsage','PeakPagefileUsage')]
        counters = Counters(); counters.cb = ctypes.sizeof(counters)
        kernel = ctypes.WinDLL('kernel32', use_last_error=True)
        kernel.GetCurrentProcess.restype = wintypes.HANDLE
        psapi = ctypes.WinDLL('psapi', use_last_error=True)
        psapi.GetProcessMemoryInfo.argtypes = [wintypes.HANDLE, ctypes.POINTER(Counters), wintypes.DWORD]
        if not psapi.GetProcessMemoryInfo(kernel.GetCurrentProcess(), ctypes.byref(counters), counters.cb):
            raise ctypes.WinError(ctypes.get_last_error())
        return int(counters.PeakWorkingSetSize)
    import resource, sys
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (1 if sys.platform == 'darwin' else 1024)


def compare_reports(baseline, actual):
    """Compare every accepted source, uncertainty, flag and model/confirmation field."""
    differences = []; maxima = {}; checked = 0
    def walk(a, b, name):
        if isinstance(a, bool) or a is None or isinstance(a, str):
            if a != b: differences.append(name)
        elif isinstance(a, (int,float)) and isinstance(b,(int,float)):
            if not (math.isfinite(a) and math.isfinite(b)):
                if a != b: differences.append(name)
                return
            absolute = abs(a-b); relative = absolute/max(abs(a),1e-300)
            parts=name.split('/')
            group = parts[-2] if parts[-1].isdigit() else parts[-1]
            current = maxima.setdefault(group, {'max_absolute':0.,'max_relative':0.})
            current['max_absolute'] = max(current['max_absolute'],absolute)
            current['max_relative'] = max(current['max_relative'],relative)
            if not math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-14): differences.append(name)
        elif isinstance(a, dict) and isinstance(b,dict):
            if set(a)!=set(b):differences.append(name+'/keys')
            for k in a.keys() & b.keys():walk(a[k],b[k],name+'/'+k)
        elif isinstance(a,list) and isinstance(b,list):
            if len(a)!=len(b):differences.append(name+'/length')
            for i,(x,y) in enumerate(zip(a,b)):walk(x,y,name+'/'+str(i))
        elif a != b: differences.append(name)
    old={o['publisher_id']:o for o in baseline['observations']}
    new={o['publisher_id']:o for o in actual['observations']}
    if set(old)!=set(new):differences.append('observation_ids')
    for ident in old.keys() & new.keys():
        a=json.loads(Path(old[ident]['measurement_file']).read_text())
        b=json.loads(Path(new[ident]['measurement_file']).read_text())
        walk(a['measurements'],b['measurements'],'sources/'+ident)
        walk(a['rejected'],b['rejected'],'rejected/'+ident)
        checked+=len(a['measurements'])
    walk(baseline.get('confirmation'),actual.get('confirmation'),'confirmation')
    walk(baseline['pair_evaluations'],actual['pair_evaluations'],'pair_evaluations')
    walk(baseline['status'],actual['status'],'status')
    return {'passed':not differences,'accepted_sources_checked':checked,
            'selected_pair_identical':baseline.get('confirmation',{}).get('selected_pair')==actual.get('confirmation',{}).get('selected_pair'),
            'confirmation_exactly_identical':baseline.get('confirmation')==actual.get('confirmation'),
            'tolerance':{'relative':1e-10,'absolute':1e-14},'differences':differences,'numeric_maxima_by_field':maxima}
