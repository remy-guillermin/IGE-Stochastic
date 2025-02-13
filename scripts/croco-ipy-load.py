#!/usr/bin/env python3
import argparse
from textwrap import dedent

def start_ipython_hello():
    """Load a CROCO simulation directory"""
    
    parser = argparse.ArgumentParser(
        prog="croco-ipy-load",
        description="Start IPython and load croco_plot in a folder."
    )
    parser.add_argument("path_dir", nargs="?", default=None)
    parser.add_argument("--clear", action="store_true", help="Clear the terminal before starting IPython")
    args = parser.parse_args()
    
    from IPython import start_ipython
    
    argv = ["--matplotlib", "-i", "-c"]
    
    code = dedent(
        """
        import os
        import croco_plot as cplot
        import glob 
        print("Listing files in", os.getcwd())
        print(os.listdir())
    """
    )
    
    lines = code.strip().split("\n")
    
    if args.path_dir is not None:
        lines.insert(1, f"os.chdir('{args.path_dir}')")
    
    if args.clear:
        import os
        os.system('clear')
    
    argv.append("; ".join(lines))
    start_ipython(argv=argv)

if __name__ == "__main__":
    start_ipython_hello()
