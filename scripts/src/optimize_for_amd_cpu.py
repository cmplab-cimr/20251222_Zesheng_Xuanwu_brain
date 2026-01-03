import os
import platform


def is_amd_cpu() -> bool:
    """Detect if the current CPU is AMD."""
    try:
        if platform.system() == "Linux":
            with open("/proc/cpuinfo", "r") as f:
                cpuinfo = f.read().lower()
                return 'amd' in cpuinfo or 'authenticamd' in cpuinfo
        elif platform.system() == "Windows":
            import winreg
            key = winreg.OpenKey(
                winreg.HKEY_LOCAL_MACHINE,
                r"HARDWARE\DESCRIPTION\System\CentralProcessor\0"
            )
            vendor, _ = winreg.QueryValueEx(key, "VendorIdentifier")
            winreg.CloseKey(key)
            return 'amd' in vendor.lower()
        elif platform.system() == "Darwin":  # macOS
            import subprocess
            result = subprocess.run(
                ["sysctl", "-n", "machdep.cpu.vendor"],
                capture_output=True, text=True
            )
            return 'amd' in result.stdout.lower()
    except Exception:
        pass
    return False


def set_thread_limits(num_threads: int = 4) -> bool:
    """
    Set BLAS/OpenMP thread limits for AMD CPUs.

    Must be called BEFORE importing numpy, scipy, or other numerical libraries.

    Args:
        num_threads: Maximum number of threads per process.

    Returns:
        True if limits were set (AMD CPU detected), False otherwise.
    """
    if not is_amd_cpu():
        return False

    thread_str = str(num_threads)

    # OpenMP settings (used by many BLAS implementations)
    os.environ.setdefault("OMP_NUM_THREADS", thread_str)
    os.environ.setdefault("OMP_DYNAMIC", "FALSE")  # Fixed typo: was "." instead of ","

    # OpenBLAS (common on AMD systems)
    os.environ.setdefault("OPENBLAS_NUM_THREADS", thread_str)

    # BLIS (AMD's optimized BLAS library)
    os.environ.setdefault("BLIS_NUM_THREADS", thread_str)

    # AMD AOCL (AMD Optimizing CPU Libraries)
    os.environ.setdefault("AOCL_NUM_THREADS", thread_str)

    # NumExpr (used by pandas and others)
    os.environ.setdefault("NUMEXPR_NUM_THREADS", thread_str)
    os.environ.setdefault("NUMEXPR_MAX_THREADS", thread_str)

    # MKL settings (in case MKL is still used on AMD)
    os.environ.setdefault("MKL_NUM_THREADS", thread_str)
    os.environ.setdefault("MKL_DYNAMIC", "FALSE")

    # macOS Accelerate framework (unlikely on AMD, but included for completeness)
    os.environ.setdefault("VECLIB_MAXIMUM_THREADS", thread_str)

    return True


def optimize_for_amd_cpu(thread_num_constraint=None):
    if thread_num_constraint is not None:
        if set_thread_limits(num_threads=thread_num_constraint):
            print("Thread limits set for AMD CPU")
        else:
            print("Non-AMD CPU detected, no limits set")

