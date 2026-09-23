"""Isolated stdlib-only launcher; TensorFlow is imported after kernel confinement."""

from __future__ import annotations

import ctypes
import errno
import hashlib
import importlib.metadata
import json
import os
import platform
import sys
from pathlib import Path

MODEL_SHA256 = "6210a674a0e7690b7b6c03184f5732a65037cf00ff4383e1892f26e90fdcb217"


def confine_network() -> None:
    """Deny network syscalls for this process and every descendant (fail closed)."""
    if sys.platform != "linux":
        raise RuntimeError("native MARLIN requires Linux seccomp confinement")
    lib = ctypes.CDLL("libseccomp.so.2", use_errno=True)
    lib.seccomp_init.argtypes = [ctypes.c_uint32]
    lib.seccomp_init.restype = ctypes.c_void_p
    lib.seccomp_syscall_resolve_name.argtypes = [ctypes.c_char_p]
    lib.seccomp_syscall_resolve_name.restype = ctypes.c_int
    lib.seccomp_rule_add.argtypes = [ctypes.c_void_p, ctypes.c_uint32, ctypes.c_int, ctypes.c_uint]
    lib.seccomp_load.argtypes = [ctypes.c_void_p]
    lib.seccomp_release.argtypes = [ctypes.c_void_p]
    ctx = lib.seccomp_init(0x7FFF0000)  # SCMP_ACT_ALLOW
    if not ctx:
        raise RuntimeError("seccomp initialization failed")
    try:
        for name in (
            "socket",
            "socketpair",
            "connect",
            "bind",
            "listen",
            "accept",
            "accept4",
            "sendto",
            "sendmsg",
            "sendmmsg",
            "recvfrom",
            "recvmsg",
            "recvmmsg",
            "socketcall",
            "io_uring_setup",
        ):
            number = lib.seccomp_syscall_resolve_name(name.encode("ascii"))
            if number >= 0 and lib.seccomp_rule_add(ctx, 0x00050000 | errno.EPERM, number, 0):
                raise RuntimeError("seccomp network rule failed")
        if lib.seccomp_load(ctx):
            raise RuntimeError("seccomp network confinement failed")
    finally:
        lib.seccomp_release(ctx)
    # Exercise the installed kernel rule rather than trusting a library return alone.
    import socket

    try:
        sock = socket.socket()
    except PermissionError:
        pass
    else:
        sock.close()
        raise RuntimeError("network confinement self-test failed")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    config = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
    threads = config["threads"]
    if not isinstance(threads, int) or not 1 <= threads <= 8:
        raise ValueError("invalid bounded thread plan")
    for key in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "TF_NUM_INTRAOP_THREADS",
    ):
        os.environ[key] = str(threads)
    os.environ["TF_NUM_INTEROP_THREADS"] = "1"
    os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
    os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"
    confine_network()
    if platform.python_version() != config["python_version"]:
        raise ValueError("Python runtime version mismatch")
    if importlib.metadata.version("tensorflow-cpu") != "2.13.1":
        raise ValueError("TensorFlow runtime version mismatch")
    if importlib.metadata.version("keras") != "2.13.1":
        raise ValueError("Keras runtime version mismatch")
    if config.get("preflight"):
        print(json.dumps({"ready": True, "confinement": "linux-seccomp-network-deny-v1"}))
        return
    model_path = Path(config["model_path"])
    if config["model_sha256"] != MODEL_SHA256 or sha256(model_path) != MODEL_SHA256:
        raise ValueError("unapproved or changed MARLIN model")
    tensor_path = Path(config["tensor_path"])
    if sha256(tensor_path) != config["tensor_sha256"]:
        raise ValueError("MARLIN tensor checksum mismatch")
    import numpy as np

    tf = importlib.import_module("tensorflow")

    if tf.__version__ != "2.13.1":
        raise ValueError("TensorFlow import identity mismatch")
    tf.config.threading.set_intra_op_parallelism_threads(threads)
    tf.config.threading.set_inter_op_parallelism_threads(1)
    vector = np.frombuffer(tensor_path.read_bytes(), dtype="<f4").reshape(1, -1)
    if vector.shape != (1, 357340) or not np.isin(vector, (-1, 0, 1)).all():
        raise ValueError("invalid MARLIN feature tensor")
    model = tf.keras.models.load_model(str(model_path), compile=False)
    if model.input_shape != (None, 357340) or model.output_shape != (None, 42):
        raise ValueError("MARLIN model tensor shapes differ")
    scores = model(vector, training=False).numpy()
    if (
        scores.shape != (1, 42)
        or not np.isfinite(scores).all()
        or (scores < 0).any()
        or (scores > 1).any()
        or not np.isclose(scores.sum(), 1, atol=1e-5)
    ):
        raise ValueError("invalid MARLIN model output scores")
    payload = {
        "scores": scores[0].tolist(),
        "tensorflow": tf.__version__,
        "keras": importlib.metadata.version("keras"),
        "python": platform.python_version(),
        "tensor_sha256": config["tensor_sha256"],
        "model_sha256": MODEL_SHA256,
        "confinement": "linux-seccomp-network-deny-v1",
        "threads": threads,
        "score_vector_sha256": hashlib.sha256(scores.astype("<f4").tobytes()).hexdigest(),
    }
    with Path(config["output_path"]).open("x", encoding="utf-8") as handle:
        json.dump(payload, handle, allow_nan=False)
        handle.flush()
        os.fsync(handle.fileno())


if __name__ == "__main__":
    main()
