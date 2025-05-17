import subprocess
import os
import sys

def check_testRun_execution():
    test_command = [
        "./bin/testRun",
        "--numticks", "10",  # small number for minimal run
        "--inputfile", "configFiles/config_Scaffold_GH2.txt",
        "--wxw", "0.1",
        "--wyw", "0.1",
        "--wzw", "0.1"
    ]

    try:
        print("Running testRun with minimal parameters...")
        result = subprocess.run(test_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        if result.returncode != 0:
            print("testRun failed with return code:", result.returncode)
            print("stderr:", result.stderr.decode())
            return False

        output_file = "output/Output_Biomarkers.csv"
        if not os.path.exists(output_file):
            print("Output file was not created:", output_file)
            return False

        print("testRun executed successfully and output file was created.")
        return True

    except Exception as e:
        print("Exception occurred during testRun execution:", str(e))
        return False


if __name__ == "__main__":
    success = check_testRun_execution()
    sys.exit(0 if success else 1)
