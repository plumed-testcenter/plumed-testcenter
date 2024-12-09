import os
import getopt
import subprocess
from zipfile import ZipFile


def zip_and_remove(myzipfile, paths):
    """
    Compresses the specified files into a zip archive and removes the original files.

    Args:
        myzipfile (str): The path to the zip file to be created.
        paths (list): A list of file paths to be compressed and removed.

    The function compresses the files at the specified paths into a zip archive
    located at 'myzipfile'. Each file is added to the archive without its directory path,
    resulting in all files being located at the root of the zip archive. After each file
    is added to the archive, it is removed from the file system.
    """
    with ZipFile(myzipfile, "w") as f_out:
        for path in paths:
            # a text file can ve compressed with maximum compression effor without problems
            nopathfile = path.split("/")[-1]
            # with removing the path of the file the files will be compressed in the "home" of the zip (with less clicks to check)
            f_out.write(path, arcname=nopathfile, compresslevel=9)
            os.remove(path)


def buildInstallPage(code):
    stable_version = (
        subprocess.check_output("plumed info --version", shell=True)
        .decode("utf-8")
        .strip()
    )
    # Save the stable version of plumed to a file to pass to the update.py script (not classy Gareth)
    with open("stable_version.md", "w+") as vf:
        vf.write(stable_version)
    # Zip all the logs
    zip_and_remove(
        f"tests/{code}/stable_output.zip",
        [f"tests/{code}/stdout.txt", f"tests/{code}/stderr.txt"],
    )
    zip_and_remove(
        f"tests/{code}/master_output.zip",
        [f"tests/{code}/stdout_master.txt", f"tests/{code}/stderr_master.txt"],
    )

    with open("templates/install.md", "r") as templatefile:
        template = templatefile.read()
    with open(f"tests/{code}/install.sh", "r") as sf:
        script = sf.read()

    with open(f"tests/{code}/install.md", "w+") as of:
        of.write(template.format(code=code, script=script))


if __name__ == "__main__":
    import sys

    code, argv = "", sys.argv[1:]
    try:
        opts, args = getopt.getopt(argv, "hc:", ["code="])
    except Exception as _:
        print("setupInstall.py -c <code>")

    for opt, arg in opts:
        if opt in ["-h"]:
            print("setupInstall.py -c <code>")
            sys.exit()
        elif opt in ["-c", "--code"]:
            code = arg
    # Setup compile page
    buildInstallPage(code)
