import subprocess
import re

out = subprocess.check_output(["cargo", "run", "-q", "--bin", "zani", "--", "--help"]).decode("utf-8")
subcommands = []
in_cmds = False
for line in out.splitlines():
    if line.startswith("Commands:"):
        in_cmds = True
        continue
    if in_cmds and line.strip() == "":
        break
    if in_cmds and line.startswith("  "):
        m = re.match(r"^  (\w+)\s+(.+)", line)
        if m and m.group(1) != "help":
            subcommands.append((m.group(1), m.group(2).strip()))

with open("docs/cli.md", "w") as f:
    f.write("# CLI Reference\n\n")
    f.write("```text\n")
    f.write(out)
    f.write("```\n\n")
    for cmd, desc in subcommands:
        f.write(f"## `zani {cmd}`\n\n")
        f.write(f"{desc}\n\n")
        f.write("```text\n")
        cmd_out = subprocess.check_output(["cargo", "run", "-q", "--bin", "zani", "--", cmd, "--help"]).decode("utf-8")
        f.write(cmd_out)
        f.write("```\n\n")
