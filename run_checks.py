from pathlib import Path
import hashlib
import json
import os
import platform
import re
import subprocess
import sys
import time
import tomllib

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "verification"
OUT.mkdir(exist_ok=True)
commands = {
    "compgate_tests": [sys.executable, str(ROOT / "tests" / "test_recap.py")],
    "joint_model_tests": [sys.executable, str(ROOT / "tests" / "test_joint_model.py")],
    "quickstart": [sys.executable, str(ROOT / "scripts" / "quickstart.py"), "--output", str(OUT / "quickstart.json")],
    "figure4": [sys.executable, str(ROOT / "scripts" / "rebuild_figure4.py"), "--output", str(OUT / "figure4")],
    "all_figures": [sys.executable, str(ROOT / "scripts" / "rebuild_all_figures.py"), "--output", str(OUT / "all_figures")],
}
env = os.environ.copy()
env["PYTHONPATH"] = str(ROOT) + os.pathsep + env.get("PYTHONPATH", "")
results = {}
for name, command in commands.items():
    started = time.perf_counter()
    proc = subprocess.run(command, cwd=ROOT, env=env, capture_output=True, text=True, timeout=300)
    seconds = time.perf_counter() - started
    log = proc.stdout + ("\n" + proc.stderr if proc.stderr else "")
    (OUT / f"{name}.log").write_text(log, encoding="utf-8")
    results[name] = {"exit_code": proc.returncode, "seconds": seconds, "log": f"{name}.log"}

pyproject = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
init = (ROOT / "compgate" / "__init__.py").read_text(encoding="utf-8")
init_version = re.search(r'__version__\s*=\s*["\']([^"\']+)', init).group(1)
required = ["LICENSE", "README.md", "CHANGELOG.md", "CITATION.cff", "requirements.lock", "docs/INPUT_OUTPUT_SCHEMA.md", "docs/INSTALL.md", "docs/TROUBLESHOOTING.md", "docs/REPRODUCTION_RECORD.md"]
missing = [name for name in required if not (ROOT / name).is_file()]
absolute_path_hits = []
drive_marker = "E:" + "\\桌面"
workspace_marker = "林煌伟" + "要求要做的"
for pattern in ("*.py", "*.md", "*.toml"):
    for path in ROOT.rglob(pattern):
        if "verification" in path.parts or "manuscript_verification" in path.parts:
            continue
        text = path.read_text(encoding="utf-8", errors="ignore")
        if drive_marker in text or workspace_marker in text:
            absolute_path_hits.append(path.relative_to(ROOT).as_posix())

quick = json.loads((OUT / "quickstart.json").read_text(encoding="utf-8")) if (OUT / "quickstart.json").is_file() else {}
anchors = json.loads((OUT / "figure4" / "figure4_anchors.json").read_text(encoding="utf-8")) if (OUT / "figure4" / "figure4_anchors.json").is_file() else {}
compgate_passed = (OUT / "compgate_tests.log").read_text(encoding="utf-8").count("PASS test_") if (OUT / "compgate_tests.log").is_file() else 0
joint_passed = (OUT / "joint_model_tests.log").read_text(encoding="utf-8").count("PASS test_") if (OUT / "joint_model_tests.log").is_file() else 0
expected_figure_hashes = {
    "Figure_1.png": "68e36fc9ddde901389d70e344636776f559ff28854cea64008cf407ef03427c8",
    "Figure_2.png": "b547eb25e1d21970a87990c554d76a3aa760c4e9bf0e0311a4b67f6b81c41b39",
    "Figure_3.png": "d368fdbd2987bf6d7af8973c06215e5c81c052c970fc5de7e28b6c8cef51e2c5",
    "Figure_4.png": "b5088b840e4a96325eba5799d4138a278860fe5451cbcf3fb3e7b35657bfc2a6",
    "Figure_5.png": "ba06dced1358409a74dee697600735a7b596465e3504f176c228a96e259ec2e2",
}
actual_figure_hashes = {}
for name in expected_figure_hashes:
    path = OUT / "all_figures" / "figures" / name
    actual_figure_hashes[name] = hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None
figure_hash_match = actual_figure_hashes == expected_figure_hashes
source_data_csv_count = len(list((ROOT / "source_data").glob("*.csv")))
report = {
    "status": "PASS",
    "scope": "Local candidate verification; not Linux or independent-team validation",
    "platform": platform.platform(),
    "python": sys.version,
    "pyproject_version": pyproject["project"]["version"],
    "package_version": init_version,
    "compgate_tests_passed": compgate_passed,
    "joint_model_tests_passed": joint_passed,
    "quickstart_failures": quick.get("failures"),
    "quickstart_repetitions": quick.get("repetitions"),
    "quickstart_runtime_ms_median": quick.get("runtime_ms_median"),
    "quickstart_peak_memory_kib": quick.get("peak_memory_kib_single_run"),
    "figure4_anchors": anchors,
    "all_figure_png_sha256_match": figure_hash_match,
    "all_figure_png_sha256": actual_figure_hashes,
    "source_data_csv_count": source_data_csv_count,
    "required_files_missing": missing,
    "local_absolute_path_hits": absolute_path_hits,
    "commands": results,
}
failures = []
if pyproject["project"]["version"] != "0.4.0" or init_version != "0.4.0":
    failures.append("version mismatch")
if compgate_passed != 11 or joint_passed != 7:
    failures.append("test count mismatch")
if quick.get("failures") != 0 or quick.get("repetitions") != 200:
    failures.append("quickstart mismatch")
if anchors.get("paired_rows") != 79 or anchors.get("common_six_arm_rows") != 50 or anchors.get("bootstrap_replicates") != 200:
    failures.append("Figure 4 anchor mismatch")
if not figure_hash_match:
    failures.append("one or more rebuilt figure PNG hashes differ from the frozen delivery")
if source_data_csv_count != 50:
    failures.append("source data CSV count mismatch")
if missing or absolute_path_hits:
    failures.append("package completeness or portability")
if any(item["exit_code"] for item in results.values()):
    failures.append("subprocess failure")
if failures:
    report["status"] = "FAIL"
    report["failures"] = failures
(OUT / "local_verification.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps(report, ensure_ascii=False, indent=2))
if failures:
    raise SystemExit(1)
