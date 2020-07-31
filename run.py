import argparse
import os
import subprocess
import sys
import urllib.request
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor
from multiprocessing import RLock

lock = RLock()


# Download instance from university server
# (sequencing by hybridization with position information and positive only errors)
def get_instance_from_server(sequence_length, oligo_length, positive_errors_percent, confidence_interval_width):
    url = f'http://www.piotr.e.wawrzyniak.doctorate.put.poznan.pl/bio.php?' \
          f'n={sequence_length}&' \
          f'k={oligo_length}&' \
          f'mode=basic&intensity=0&position=1&' \
          f'sqpep={positive_errors_percent}&' \
          f'sqnep=0&' \
          f'pose={confidence_interval_width}'

    with urllib.request.urlopen(url) as r:
        xml = r.read().decode('utf-8')
        instance_id = ET.fromstring(xml).attrib['key']

    output_dir = './instances'
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    output_filename = f'{sequence_length}_{oligo_length}_{positive_errors_percent}_' \
                      f'{confidence_interval_width}_{instance_id}.xml'
    with open(f'{output_dir}/{output_filename}', "w") as f:
        f.write(xml)


def initialize_results_file(filename):
    with open(filename, "w") as f:
        f.write("INSTANCE,ALGORITHM,SEED,STATUS,OBJECTIVE,TIME\n")


def write_results(filename, results):
    content = ",".join([results["instance"], results["algorithm"], results["seed"],
                        results["status"], results["objective"], results["nodes"]])

    with open(filename, "a") as f:
        f.write(content + "\n")


def run_entry(instance, algorithm, data, seed):
    # Keep the result
    results = dict()
    results["instance"] = instance
    results["algorithm"] = algorithm
    results["seed"] = str(seed)
    results["status"] = "Error"
    results["objective"] = ""
    results["time"] = ""

    try:

        # Create the command to run
        command = ([data["command"], "--seed", str(seed)] + data["algorithms"][algorithm] +
                   ["--file", data["instances"][instance]])

        # Run the command
        output = subprocess.run(command, stdout=subprocess.PIPE, universal_newlines=True)

        # Check if the command run without errors
        if output.returncode == 0:

            # Get the output as a string
            str_output = str(output.stdout).strip('\n\t ').split(" ")

            # Process the output
            status = str_output[0]
            results["status"] = status

            if status != "Error":
                results["nodes"] = str_output[3]
                results["time"] = str_output[4]

            if status == "Feasible" or status == "Optimal":
                results["objective"] = str_output[1]

    except:
        results["status"] = "Error"
        results["objective"] = ""
        results["nodes"] = ""
        results["time"] = ""

    # Write results in the output file and log the progress
    with lock:

        # Write results
        write_results(data["output-file"], results)

        # Progress log
        data["progress"] += 1
        total_entries = len(data["instances"]) * len(data["algorithms"]) * len(data["seeds"])
        progress = (data["progress"] / total_entries) * 100

        print(f"[{data['progress']:3} of {total_entries:3} ({progress:6.2f}%) completed] "
              f"{algorithm:15} -> {instance:16} -> {seed:4} -> {results['status']:8}")

        sys.stdout.flush()


def main():
    # Configure argument parser
    parser = argparse.ArgumentParser(description="Perform the computational experiments.")
    parser.add_argument("-c", "--continue", dest="continue_previous", action="store_true",
                        help="Resume the experiment from a previous state.")

    # Parse input arguments
    args = parser.parse_args()

    # Map used to store data necessary to run the experiments
    data = dict()
    data["command"] = "./build/xseq"  # program
    data["threads"] = 2  # entries to solve simultaneously
    data["output-file"] = "results.csv"  # file to write the results
    data["instances-path"] = "./instances"  # base path to instances
    data["instances"] = dict()  # path to each instance
    data["algorithms"] = dict()  # settings of each algorithm
    data["seeds"] = None  # seeds
    data["progress"] = 0  # num. of entries finished

    # Set the seeds used at each repetition of (instance, algorithm)
    data["seeds"] = [13, 20, 38, 45, 52, 66, 71, 89, 94, 107]

    # Load the list of instances
    positive_error_percents = (5, 10, 15, 20, 25)
    for pep in positive_error_percents:
        get_instance_from_server(500, 10, pep, 20)

    for instance in os.listdir('./instances'):
        data["instances"][instance] = data["instances-path"] + "/" + instance + ".xml"

    # Algorithms' settings
    params_exact = ("--algorithm exact".split(" ") +
                    # "--verbose" +
                    "--use-vertex-lists true".split(" "))

    params_genetic = ("--algorithm genetic".split(" ") +
                      "--use-binary-encoding true".split(" ") +
                      "--population-size 50".split(" ") +
                      "--crossover-probability 0.9".split(" ") +
                      "--mutation-probability 0.001".split(" ") +
                      "--selection srm".split(" ") +
                      "--keep-best true".split(" "))

    data["algorithms"]["exact"] = params_exact
    data["algorithms"]["genetic"] = params_genetic
    # data["algorithms"]["hybrid"] = params_hybrid

    # Create and initialize the output file
    if not os.path.exists(data["output-file"]) or not args.continue_previous:
        initialize_results_file(data["output-file"])

    # Check entries already solved
    entries_completed = []
    if args.continue_previous:
        if os.path.exists(data["output-file"]):
            with open(data["output-file"], "r") as f:
                lines = f.readlines()[1:]
                for line in lines:
                    values = line.strip().split(",")
                    print(values)
                    if len(values) > 0:
                        entries_completed.append((values[0], values[1], int(values[2])))

    data["progress"] = len(entries_completed)

    # Run each entry (instance, algorithm, seed)
    # with ThreadPoolExecutor(max_workers=data["threads"]) as executor:
    #     for instance in list(data["instances"].keys()):
    #         for algorithm in list(data["algorithms"].keys()):
    #             for seed in data["seeds"]:
    #                 if entries_completed.count((instance, algorithm, seed)) == 0:
    #                     executor.submit(run_entry, instance, algorithm, data, seed)


###################################################################################################
# Main statements
#
if __name__ == "__main__":
    main()
