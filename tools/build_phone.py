#!/usr/bin/env python3
"""Build a portable .pyz with the standard library; no pip or compiler needed."""
from pathlib import Path
import argparse
import shutil
import tempfile
import zipapp


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("dist/explorative-math.pyz"))
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="explorative-math-build-") as temporary:
        source = Path(temporary)
        shutil.copytree(root / "explorative_math", source / "explorative_math",
                        ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
        shutil.copy2(root / "LICENSE", source / "LICENSE")
        (source / "__main__.py").write_text(
            "import runpy\nrunpy.run_module('explorative_math', run_name='__main__')\n",
            encoding="utf-8")
        zipapp.create_archive(source, target=output, compressed=True)
    print(output)


if __name__ == "__main__":
    main()
