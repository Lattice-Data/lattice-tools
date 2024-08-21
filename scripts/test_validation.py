import os
import multiprocessing
import subprocess


files = [f for f in os.listdir() if ".h5ad" in f]
workers = len(files)


def validate(file_name):
    validate = subprocess.run(["cellxgene-schema", "validate", file_name], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print(f"Validation for {file_name}...")
    for line in validate.stdout.decode('utf-8').split('\n'):
        # results[file]["stdout"].append(line)
        print(line)
    for line in validate.stderr.decode('utf-8').split('\n'):
        # results[file]["stderr"].append(line)
        print(line)
    print("=" * 113 + "\n")


def validate_all_files(files):
    with multiprocessing.Pool(processes=workers) as pool:
        pool.map(validate, files)


if __name__ == "__main__":
    validate_all_files(files)
