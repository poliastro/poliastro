import os
import subprocess


"""
This python script runs all the GMAT scripts provided in the
```GMAT_Tests/Scripts``` folder. Each script saves a report wih the same name
in the output folder in the GMAT folder. The new scripts then replace the
old reports in ```poliastro/src/poliastro/GMAT_Tests/Reports/``` folder.

The script assumes that both the poliastro and GMAT folder are present
in the root directory.

The script has to be run in the folder that contains the GMAT console.
ex. In the folder GMAT/R2018a/bin in Linux
"""


def generate_scripts():
    # Make the file paths
    base_dir = os.path.expanduser("~/poliastro/src/poliastro/GMAT_Tests/")
    script_dir = os.path.join(base_dir, "Scripts/")
    report_dir = os.path.join(base_dir, "Reports/")

    try:
        # Create target Directory to save the generated reports
        os.mkdir(report_dir)
    except FileExistsError:
        pass

    try:
        # Create target Directory in the GMAT folder to temporarily store the new reports
        os.mkdir("../output/Reports")
        print("Done")
    except FileExistsError:
        pass

    dirs = os.listdir(script_dir)
    out = []
    for file in dirs:
        b = subprocess.check_output(["./GmatConsole", script_dir + file]).decode(
            "utf-8"
        )
        b = [i for i in b.split("\n") if i != ""]
        out.append([file, " : ", b[-2]])

    # Prints result of each script run
    for i in out:
        print(" ".join(i))

    # Replaces the old reports with the new ones
    command_run = subprocess.call(["rsync", "-r", "../output/Reports/", report_dir])
    if command_run != 0:
        print(
            "The reports have been generated in the output/Reports in GMAT directory folder. But, the old reports were not replaced."
        )


if __name__ == "__main__":
    generate_scripts()
