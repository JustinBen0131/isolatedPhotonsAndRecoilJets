#!/usr/bin/env python3
"""Build and invoke sealed context packs for the Claude slide-policy worker."""

from __future__ import annotations

import argparse
import json
import os
import re
import shlex
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
CLAUDE_BIN = Path.home() / ".local" / "bin" / "claude"
WORKER_ID = "claude_slide_worker"
AGENT_NAME = "thesis-slide-policy-worker"
AGENT_FILE = REPO_ROOT / ".claude" / "agents" / f"{AGENT_NAME}.md"
POLICY_FILE = REPO_ROOT / "agent_context" / "policies" / "DELEGATION_KERNEL.md"
TEMPLATE_FILE = (
    REPO_ROOT
    / "agent_context"
    / "templates"
    / "delegation"
    / "CLAUDE_SLIDE_WORKER_CONTEXT_PACK.md"
)
DEFAULT_OUTPUT_ROOT = (
    REPO_ROOT / "agent_context" / "local" / "delegations" / WORKER_ID
)
DEFAULT_POLICY_FILES = [
    REPO_ROOT / "agent_context" / "policies" / "DELEGATION_KERNEL.md",
    REPO_ROOT / "agent_context" / "policies" / "SLIDES_WORKFLOW.md",
    REPO_ROOT / "agent_context" / "policies" / "PLOTTING.md",
    REPO_ROOT / "agent_context" / "SLIDE_STYLE_MAP.md",
]
FISH_CANDIDATES = [
    Path("/opt/homebrew/bin/fish"),
    Path("/usr/local/bin/fish"),
    Path("/bin/fish"),
]


def local_now() -> datetime:
    return datetime.now().astimezone()


def iso_now() -> str:
    return local_now().isoformat(timespec="seconds")


def slugify(value: str, limit: int = 56) -> str:
    slug = re.sub(r"[^A-Za-z0-9]+", "-", value.strip().lower()).strip("-")
    return (slug or "slide-worker-pack")[:limit].strip("-")


def repo_path(value: str | Path) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = REPO_ROOT / path
    return path.resolve()


def display_path(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def run_cmd(argv: list[str], timeout: int = 20) -> dict[str, Any]:
    try:
        completed = subprocess.run(
            argv,
            cwd=REPO_ROOT,
            text=True,
            capture_output=True,
            timeout=timeout,
            check=False,
        )
        return {
            "argv": argv,
            "returncode": completed.returncode,
            "stdout": completed.stdout.strip(),
            "stderr": completed.stderr.strip(),
            "timeout": False,
        }
    except FileNotFoundError as exc:
        return {
            "argv": argv,
            "returncode": 127,
            "stdout": "",
            "stderr": str(exc),
            "timeout": False,
        }
    except subprocess.TimeoutExpired as exc:
        return {
            "argv": argv,
            "returncode": 124,
            "stdout": (exc.stdout or "").strip()
            if isinstance(exc.stdout, str)
            else "",
            "stderr": (exc.stderr or "").strip()
            if isinstance(exc.stderr, str)
            else "",
            "timeout": True,
        }


def auth_status(argv: list[str]) -> tuple[bool, str]:
    result = run_cmd(argv, timeout=20)
    detail = result["stdout"] or result["stderr"]
    ok = False
    if result["stdout"]:
        try:
            payload = json.loads(result["stdout"])
            ok = bool(payload.get("loggedIn"))
        except json.JSONDecodeError:
            ok = False
    return ok, detail


def fish_binary() -> Path | None:
    return next((path for path in FISH_CANDIDATES if path.exists()), None)


def shell_quote_command(argv: list[str]) -> str:
    return " ".join(shlex.quote(arg) for arg in argv)


def fish_wrapped_command(argv: list[str], fish_bin: Path) -> list[str]:
    return [str(fish_bin), "-lc", "exec " + shell_quote_command(argv)]


def select_claude_runner(prompt: str, max_budget_usd: str) -> tuple[list[str], str]:
    direct_argv = claude_command(prompt, max_budget_usd)
    direct_ok, _ = auth_status([str(CLAUDE_BIN), "auth", "status", "--json"])
    if direct_ok:
        return direct_argv, "direct"
    fish_bin = fish_binary()
    if fish_bin:
        fish_ok, _ = auth_status(
            [str(fish_bin), "-lc", "claude auth status --json"]
        )
        if fish_ok:
            return fish_wrapped_command(direct_argv, fish_bin), "fish-login-shell"
    return direct_argv, "direct-unauthenticated"


def read_excerpt(path: Path, max_chars: int = 9000) -> str:
    if not path.exists():
        return f"[missing: {display_path(path)}]"
    if path.is_dir():
        return f"[directory: {display_path(path)}]"
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        return f"[unreadable: {display_path(path)}: {exc}]"
    if len(text) <= max_chars:
        return text
    return text[:max_chars] + "\n\n[excerpt truncated by claude_slide_worker.py]\n"


def format_artifact_list(paths: list[Path]) -> str:
    if not paths:
        return "- None."
    lines = []
    for path in paths:
        status = "exists" if path.exists() else "missing"
        kind = "dir" if path.is_dir() else "file"
        lines.append(f"- `{display_path(path)}` ({kind}, {status})")
    return "\n".join(lines)


def build_policy_excerpts(paths: list[Path]) -> str:
    blocks = []
    seen: set[Path] = set()
    for path in paths:
        if path in seen:
            continue
        seen.add(path)
        blocks.append(
            "\n".join(
                [
                    f"### {display_path(path)}",
                    "",
                    "```text",
                    read_excerpt(path),
                    "```",
                ]
            )
        )
    return "\n\n".join(blocks) if blocks else "- None."


def unique_run_dir(output_root: Path, objective: str) -> Path:
    stamp = local_now().strftime("%Y%m%dT%H%M%S%z")
    base = output_root / f"{stamp}-{slugify(objective)}"
    candidate = base
    suffix = 2
    while candidate.exists():
        candidate = output_root / f"{base.name}-{suffix}"
        suffix += 1
    return candidate


def create_context_pack(args: argparse.Namespace) -> Path:
    output_root = repo_path(args.out_root or DEFAULT_OUTPUT_ROOT)
    output_root.mkdir(parents=True, exist_ok=True)
    run_dir = unique_run_dir(output_root, args.objective)
    run_dir.mkdir(parents=True, exist_ok=False)

    artifact_paths: list[Path] = []
    for value in args.png or []:
        artifact_paths.append(repo_path(value))
    for value in args.script or []:
        artifact_paths.append(repo_path(value))
    for value in args.include_file or []:
        artifact_paths.append(repo_path(value))

    policy_paths = [repo_path(value) for value in args.policy or []]
    if not policy_paths:
        policy_paths = DEFAULT_POLICY_FILES.copy()

    notes_chunks: list[str] = []
    if args.notes:
        notes_chunks.append(args.notes.strip())
    if args.notes_file:
        notes_chunks.append(read_excerpt(repo_path(args.notes_file), max_chars=12000))
    if not notes_chunks:
        notes_chunks.append("No additional notes.")

    artifact_detail = [format_artifact_list(artifact_paths)]
    for path in artifact_paths:
        if path.suffix.lower() in {
            ".md",
            ".txt",
            ".py",
            ".yaml",
            ".yml",
            ".json",
            ".tex",
            ".js",
            ".ts",
            ".css",
            ".html",
        }:
            artifact_detail.append(
                "\n".join(
                    [
                        "",
                        f"### Excerpt: {display_path(path)}",
                        "",
                        "```text",
                        read_excerpt(path, max_chars=12000),
                        "```",
                    ]
                )
            )

    template = TEMPLATE_FILE.read_text(encoding="utf-8")
    pack_id = run_dir.name
    context_pack = (
        template.replace("{{generated_at}}", iso_now())
        .replace("{{pack_id}}", pack_id)
        .replace("{{objective}}", args.objective.strip())
        .replace("{{referenced_artifacts}}", "\n".join(artifact_detail))
        .replace("{{notes}}", "\n\n".join(notes_chunks))
        .replace("{{policy_excerpts}}", build_policy_excerpts(policy_paths))
    )

    prompt = "\n".join(
        [
            "You are the ThesisAnalysis Claude slide-policy worker under Codex control.",
            "Use only the sealed context pack below and the named local files.",
            "Return exactly one Markdown Slide Worker Report matching the output contract.",
            "Do not edit files, mutate external state, browse the web, install software,",
            "or ask for broader authority. If a task exceeds scope, refuse that portion.",
            "",
            "<sealed_context_pack>",
            context_pack,
            "</sealed_context_pack>",
            "",
        ]
    )

    manifest = {
        "schema_version": 1,
        "worker_id": WORKER_ID,
        "agent_name": AGENT_NAME,
        "created_at": iso_now(),
        "repo_root": str(REPO_ROOT),
        "run_dir": str(run_dir),
        "objective": args.objective.strip(),
        "artifacts": [str(path) for path in artifact_paths],
        "policy_files": [str(path) for path in policy_paths],
        "notes_file": str(repo_path(args.notes_file)) if args.notes_file else None,
        "context_pack": str(run_dir / "context_pack.md"),
        "prompt": str(run_dir / "prompt.txt"),
    }

    (run_dir / "context_pack.md").write_text(context_pack, encoding="utf-8")
    (run_dir / "prompt.txt").write_text(prompt, encoding="utf-8")
    (run_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return run_dir


def cmd_doctor(args: argparse.Namespace) -> int:
    output_root = DEFAULT_OUTPUT_ROOT
    output_root.mkdir(parents=True, exist_ok=True)
    fish_bin = fish_binary()
    checks: list[dict[str, Any]] = []

    def add_check(name: str, ok: bool, detail: str, critical: bool = True) -> None:
        checks.append(
            {
                "name": name,
                "ok": ok,
                "critical": critical,
                "detail": detail,
            }
        )

    add_check(
        "claude_binary",
        CLAUDE_BIN.exists() and os.access(CLAUDE_BIN, os.X_OK),
        str(CLAUDE_BIN),
    )
    version = run_cmd([str(CLAUDE_BIN), "--version"], timeout=20)
    add_check(
        "claude_version",
        version["returncode"] == 0 and "Claude Code" in version["stdout"],
        version["stdout"] or version["stderr"],
    )
    auth_ok, auth_detail = auth_status(
        [str(CLAUDE_BIN), "auth", "status", "--json"]
    )
    fish_auth_ok = False
    fish_auth_detail = "fish not found"
    if fish_bin:
        fish_auth_ok, fish_auth_detail = auth_status(
            [str(fish_bin), "-lc", "claude auth status --json"]
        )
    add_check(
        "claude_auth_direct",
        auth_ok,
        auth_detail
        or "run /Users/patsfan753/.local/bin/claude auth login --claudeai",
        critical=False,
    )
    add_check(
        "claude_auth_fish",
        fish_auth_ok,
        fish_auth_detail,
        critical=False,
    )
    add_check(
        "claude_invocation_auth",
        auth_ok or fish_auth_ok,
        "direct" if auth_ok else "fish-login-shell" if fish_auth_ok else "none",
    )
    add_check("agent_file", AGENT_FILE.exists(), display_path(AGENT_FILE))
    add_check("delegation_policy", POLICY_FILE.exists(), display_path(POLICY_FILE))
    add_check("context_template", TEMPLATE_FILE.exists(), display_path(TEMPLATE_FILE))
    add_check(
        "output_root",
        output_root.exists() and output_root.is_dir(),
        display_path(output_root),
    )

    zsh = run_cmd(["zsh", "-lc", "command -v claude && claude --version"], timeout=20)
    add_check(
        "zsh_login_path",
        zsh["returncode"] == 0,
        zsh["stdout"] or zsh["stderr"],
        critical=False,
    )
    if fish_bin:
        fish = run_cmd(
            [str(fish_bin), "-lc", "command -v claude; claude --version"],
            timeout=20,
        )
        add_check(
            "fish_login_path",
            fish["returncode"] == 0,
            fish["stdout"] or fish["stderr"],
            critical=False,
        )
    else:
        add_check("fish_login_path", False, "fish not found", critical=False)

    critical_ok = all(check["ok"] for check in checks if check["critical"])
    status = "ok" if critical_ok else "fail"
    payload = {
        "status": status,
        "checked_at": iso_now(),
        "repo_root": str(REPO_ROOT),
        "checks": checks,
    }
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print(f"claude_slide_worker doctor: {status}")
        for check in checks:
            marker = "ok" if check["ok"] else "FAIL" if check["critical"] else "warn"
            print(f"- {marker}: {check['name']}: {check['detail']}")
    return 0 if critical_ok else 1


def claude_command(prompt: str, max_budget_usd: str) -> list[str]:
    return [
        str(CLAUDE_BIN),
        "--agent",
        AGENT_NAME,
        "--permission-mode",
        "dontAsk",
        "--output-format",
        "json",
        "--no-session-persistence",
        "--setting-sources",
        "project,local",
        "--max-budget-usd",
        max_budget_usd,
        "--tools",
        "Read,Grep,Glob",
        "-p",
        prompt,
    ]


def extract_report(stdout: str) -> str | None:
    stripped = stdout.strip()
    if not stripped:
        return None
    try:
        payload = json.loads(stripped)
    except json.JSONDecodeError:
        return stripped if "# Slide Worker Report" in stripped else None
    for key in ("result", "content", "text", "message"):
        value = payload.get(key)
        if isinstance(value, str) and value.strip():
            return value.strip()
    return None


def invoke_pack(
    pack_dir: Path,
    max_budget_usd: str,
    timeout: int,
    print_summary: bool = True,
) -> int:
    pack_dir = pack_dir.resolve()
    prompt_file = pack_dir / "prompt.txt"
    if not prompt_file.exists():
        print(f"missing prompt file: {prompt_file}", file=sys.stderr)
        return 2
    if not CLAUDE_BIN.exists():
        print(f"missing Claude CLI: {CLAUDE_BIN}", file=sys.stderr)
        return 2

    prompt = prompt_file.read_text(encoding="utf-8")
    report_file = pack_dir / "worker_report.md"
    error_file = pack_dir / "worker_error.md"
    for stale_file in (report_file, error_file):
        if stale_file.exists():
            stale_file.unlink()
    argv, runner = select_claude_runner(prompt, max_budget_usd)
    started_at = iso_now()
    try:
        completed = subprocess.run(
            argv,
            cwd=REPO_ROOT,
            text=True,
            capture_output=True,
            timeout=timeout,
            check=False,
        )
        returncode = completed.returncode
        raw_returncode = completed.returncode
        stdout = completed.stdout
        stderr = completed.stderr
        timed_out = False
    except subprocess.TimeoutExpired as exc:
        returncode = 124
        raw_returncode = 124
        stdout = exc.stdout if isinstance(exc.stdout, str) else ""
        stderr = exc.stderr if isinstance(exc.stderr, str) else ""
        timed_out = True

    finished_at = iso_now()
    (pack_dir / "claude_stdout.txt").write_text(stdout or "", encoding="utf-8")
    (pack_dir / "claude_stderr.txt").write_text(stderr or "", encoding="utf-8")
    report = extract_report(stdout or "")
    contract_ok = bool(report and "# Slide Worker Report" in report)
    if returncode == 0 and not contract_ok:
        returncode = 3
        report = None
    elif returncode != 0 and report and "# Slide Worker Report" not in report:
        report = None
    invocation = {
        "schema_version": 1,
        "worker_id": WORKER_ID,
        "agent_name": AGENT_NAME,
        "started_at": started_at,
        "finished_at": finished_at,
        "returncode": returncode,
        "raw_returncode": raw_returncode,
        "report_contract_ok": contract_ok,
        "timed_out": timed_out,
        "timeout_seconds": timeout,
        "max_budget_usd": max_budget_usd,
        "runner": runner,
        "cwd": str(REPO_ROOT),
        "claude_bin": str(CLAUDE_BIN),
        "argv_without_prompt": claude_command("@prompt.txt", max_budget_usd),
    }
    (pack_dir / "claude_invocation.json").write_text(
        json.dumps(invocation, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    if report:
        report_file.write_text(report + "\n", encoding="utf-8")
    elif returncode != 0:
        error_file.write_text(
            "\n".join(
                [
                    "# Claude Slide Worker Error",
                    "",
                    f"- returncode: {returncode}",
                    f"- raw_returncode: {raw_returncode}",
                    f"- timed_out: {timed_out}",
                    f"- report_contract_ok: {contract_ok}",
                    f"- stdout: {(stdout or '').strip()[:1000]}",
                    f"- stderr: {(stderr or '').strip()[:1000]}",
                    "",
                ]
            ),
            encoding="utf-8",
        )

    if print_summary:
        print(f"pack_dir={display_path(pack_dir)}")
        print(f"returncode={returncode}")
        print(f"timed_out={timed_out}")
        print(f"runner={runner}")
        if report:
            print(f"worker_report={display_path(report_file)}")
        elif returncode != 0:
            print(f"worker_error={display_path(error_file)}")
        if stderr:
            print(f"stderr={stderr.strip()[:500]}")
    return returncode


def cmd_pack(args: argparse.Namespace) -> int:
    run_dir = create_context_pack(args)
    print(display_path(run_dir))
    print(
        "next: python3 scripts/os/delegation/claude_slide_worker.py invoke "
        f"--pack-dir {display_path(run_dir)}"
    )
    return 0


def cmd_invoke(args: argparse.Namespace) -> int:
    return invoke_pack(
        repo_path(args.pack_dir),
        max_budget_usd=args.max_budget_usd,
        timeout=args.timeout,
    )


def cmd_smoke(args: argparse.Namespace) -> int:
    pack_args = argparse.Namespace(
        objective=(
            "Smoke test the Claude slide worker contract using only the sealed "
            "context pack and local delegation policy."
        ),
        png=[],
        script=[],
        include_file=[],
        policy=[str(POLICY_FILE), str(AGENT_FILE)],
        notes=(
            "This is a smoke test. Confirm the worker boundaries, identify any "
            "missing setup evidence, and return the required report shape. Do "
            "not inspect unrelated project files."
        ),
        notes_file=None,
        out_root=args.out_root,
    )
    run_dir = create_context_pack(pack_args)
    print(display_path(run_dir))
    if not args.run:
        print(
            "next: python3 scripts/os/delegation/claude_slide_worker.py invoke "
            f"--pack-dir {display_path(run_dir)} --max-budget-usd {args.max_budget_usd}"
        )
        return 0
    return invoke_pack(
        run_dir,
        max_budget_usd=args.max_budget_usd,
        timeout=args.timeout,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Create and run sealed context packs for the ThesisAnalysis Claude "
            "slide-policy worker."
        )
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    doctor = subparsers.add_parser("doctor", help="Check local Claude worker setup.")
    doctor.add_argument("--json", action="store_true", help="Emit JSON.")
    doctor.set_defaults(func=cmd_doctor)

    pack = subparsers.add_parser("pack", help="Create a sealed context pack.")
    pack.add_argument("--objective", required=True)
    pack.add_argument("--png", action="append", help="Slide PNG path.")
    pack.add_argument("--script", action="append", help="Slide generator/script path.")
    pack.add_argument("--include-file", action="append", help="Additional local file.")
    pack.add_argument("--policy", action="append", help="Policy/style file to excerpt.")
    pack.add_argument("--notes", help="Inline notes from Codex.")
    pack.add_argument("--notes-file", help="File containing Codex notes.")
    pack.add_argument("--out-root", help="Override output root.")
    pack.set_defaults(func=cmd_pack)

    invoke = subparsers.add_parser("invoke", help="Invoke Claude on a pack.")
    invoke.add_argument("--pack-dir", required=True)
    invoke.add_argument("--max-budget-usd", default="0.50")
    invoke.add_argument("--timeout", type=int, default=240)
    invoke.set_defaults(func=cmd_invoke)

    smoke = subparsers.add_parser("smoke", help="Create a minimal smoke pack.")
    smoke.add_argument("--run", action="store_true", help="Also invoke Claude.")
    smoke.add_argument("--max-budget-usd", default="0.25")
    smoke.add_argument("--timeout", type=int, default=180)
    smoke.add_argument("--out-root", help="Override output root.")
    smoke.set_defaults(func=cmd_smoke)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
