"""Windows worker-tree working-set sampling and scoped task-tree termination.

Standard library only. This module does not start scientific computation.
"""
import ctypes
from ctypes import wintypes
import subprocess


class PROCESSENTRY32W(ctypes.Structure):
    _fields_ = [
        ('dwSize', wintypes.DWORD), ('cntUsage', wintypes.DWORD),
        ('th32ProcessID', wintypes.DWORD), ('th32DefaultHeapID', ctypes.c_size_t),
        ('th32ModuleID', wintypes.DWORD), ('cntThreads', wintypes.DWORD),
        ('th32ParentProcessID', wintypes.DWORD), ('pcPriClassBase', wintypes.LONG),
        ('dwFlags', wintypes.DWORD), ('szExeFile', wintypes.WCHAR * 260),
    ]


class PROCESS_MEMORY_COUNTERS(ctypes.Structure):
    _fields_ = [('cb', wintypes.DWORD), ('PageFaultCount', wintypes.DWORD)] + [
        (name, ctypes.c_size_t) for name in (
            'PeakWorkingSetSize', 'WorkingSetSize', 'QuotaPeakPagedPoolUsage',
            'QuotaPagedPoolUsage', 'QuotaPeakNonPagedPoolUsage',
            'QuotaNonPagedPoolUsage', 'PagefileUsage', 'PeakPagefileUsage')
    ]


_kernel = ctypes.WinDLL('kernel32', use_last_error=True)
_psapi = ctypes.WinDLL('psapi', use_last_error=True)
_kernel.CreateToolhelp32Snapshot.argtypes = [wintypes.DWORD, wintypes.DWORD]
_kernel.CreateToolhelp32Snapshot.restype = wintypes.HANDLE
_kernel.Process32FirstW.argtypes = [wintypes.HANDLE, ctypes.POINTER(PROCESSENTRY32W)]
_kernel.Process32FirstW.restype = wintypes.BOOL
_kernel.Process32NextW.argtypes = [wintypes.HANDLE, ctypes.POINTER(PROCESSENTRY32W)]
_kernel.Process32NextW.restype = wintypes.BOOL
_kernel.OpenProcess.argtypes = [wintypes.DWORD, wintypes.BOOL, wintypes.DWORD]
_kernel.OpenProcess.restype = wintypes.HANDLE
_kernel.CloseHandle.argtypes = [wintypes.HANDLE]
_kernel.CloseHandle.restype = wintypes.BOOL
_kernel.GetExitCodeProcess.argtypes = [wintypes.HANDLE, ctypes.POINTER(wintypes.DWORD)]
_kernel.GetExitCodeProcess.restype = wintypes.BOOL
_psapi.GetProcessMemoryInfo.argtypes = [
    wintypes.HANDLE, ctypes.POINTER(PROCESS_MEMORY_COUNTERS), wintypes.DWORD]
_psapi.GetProcessMemoryInfo.restype = wintypes.BOOL


def _processes():
    handle = _kernel.CreateToolhelp32Snapshot(0x00000002, 0)  # TH32CS_SNAPPROCESS
    if handle == ctypes.c_void_p(-1).value:
        raise ctypes.WinError(ctypes.get_last_error())
    result = {}
    entry = PROCESSENTRY32W()
    entry.dwSize = ctypes.sizeof(entry)
    try:
        more = _kernel.Process32FirstW(handle, ctypes.byref(entry))
        while more:
            result[entry.th32ProcessID] = {
                'pid': int(entry.th32ProcessID),
                'parent_pid': int(entry.th32ParentProcessID),
                'executable': entry.szExeFile,
            }
            more = _kernel.Process32NextW(handle, ctypes.byref(entry))
        error = ctypes.get_last_error()
        if error != 18:  # ERROR_NO_MORE_FILES is normal snapshot exhaustion.
            raise ctypes.WinError(error)
        return result
    finally:
        _kernel.CloseHandle(handle)


def _working_set(pid):
    handle = _kernel.OpenProcess(0x1000, False, pid)  # PROCESS_QUERY_LIMITED_INFORMATION
    if not handle:
        error = ctypes.get_last_error()
        if error == 87:  # Process exited between snapshot and opening its handle.
            return 0
        raise ctypes.WinError(error)
    counters = PROCESS_MEMORY_COUNTERS()
    counters.cb = ctypes.sizeof(counters)
    try:
        if not _psapi.GetProcessMemoryInfo(handle, ctypes.byref(counters), counters.cb):
            error = ctypes.get_last_error()
            exit_code = wintypes.DWORD()
            if _kernel.GetExitCodeProcess(handle, ctypes.byref(exit_code)) and exit_code.value != 259:
                return 0
            raise ctypes.WinError(error)
        return int(counters.WorkingSetSize)
    finally:
        _kernel.CloseHandle(handle)


def process_tree_snapshot(pid):
    """Current sum of WorkingSetSize over pid and all snapshot descendants."""
    processes = _processes()
    children = {}
    for row in processes.values():
        children.setdefault(row['parent_pid'], []).append(row['pid'])
    pending = [int(pid)]
    visited = set()
    rows = []
    while pending:
        current = pending.pop()
        if current in visited:
            continue
        visited.add(current)
        pending.extend(children.get(current, ()))
        if current in processes:
            row = dict(processes[current])
            row['working_set_bytes'] = _working_set(current)
            rows.append(row)
    return {'root_pid': int(pid), 'rss_bytes': sum(row['working_set_bytes'] for row in rows),
            'processes': sorted(rows, key=lambda row: row['pid'])}


def process_tree_rss(pid):
    """Return process-tree current RSS/working-set sum in bytes."""
    return process_tree_snapshot(pid)['rss_bytes']


def terminate_tree(pid):
    """Force-stop this owned worker PID and descendants via taskkill /PID /T /F.

    Return the actual command result for the caller's budget/stop log. Call while
    the owned worker is still alive; no global process-name kill is used.
    """
    command = ['taskkill', '/PID', str(int(pid)), '/T', '/F']
    result = subprocess.run(command, capture_output=True, text=True,
                            creationflags=subprocess.CREATE_NO_WINDOW)
    return {'command': command, 'returncode': result.returncode,
            'stdout': result.stdout, 'stderr': result.stderr}
