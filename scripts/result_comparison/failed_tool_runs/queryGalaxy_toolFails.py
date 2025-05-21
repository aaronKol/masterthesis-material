from bioblend.galaxy import GalaxyInstance
from bioblend.galaxy.histories import HistoryClient
from bioblend.galaxy.datasets import DatasetClient
import re

GALAXY_URL = "https://usegalaxy.eu/"
API_KEY = "-"

gi = GalaxyInstance(url=GALAXY_URL, key=API_KEY)

# script for studying fail reasons, that contains 4 functions:
# 1. query VAPORs stdOut&Err in a dict
# 2. query IRMAs stdOut&Err in a dict
# 3. write IRMA dict to tsv file
# 4. write VAPOR dict to tsv file
# querying and writing differs for the tools. Still with effort, parts could be generalized into one function for each task in the future.

def queryVaporJobStdOutStdErr(history_client, dataset_client, history_id, tool_id = None):
    datasets = gi.histories.show_history(history_id, contents=True)
    resultDict = {}
    for dataset in datasets:
        if 'state' in dataset and dataset['state'] == "error":
            #print(dataset)
            if 'id' in dataset:
                dataset_id = dataset['id']
                provenance_data = history_client.show_dataset_provenance(history_id, dataset_id)
                if tool_id:
                    if provenance_data["tool_id"] != tool_id:
                        print("Skipping")
                        continue
                std_err = provenance_data["stderr"]
                std_out = provenance_data["stdout"]
                
                for param_name, param_value in provenance_data['parameters'].items():
                    if isinstance(param_value, dict) and 'id' in param_value:
                        input_dataset_id = param_value['id']
                        try:
                            input_dataset = dataset_client.show_dataset(input_dataset_id)
                            print(input_dataset)
                            if '-' in input_dataset['name']:
                                continue
                            print("Got the right set")

                            
                            job = gi.jobs.show_job(input_dataset['creating_job']) #job that craeated the fastp file

                            print("JOB", job)

                            command_line = job['command_line'] # commando zeile, weil da die samplenamen zu finden sind, mit folgender regex:
                            pattern = r"(\S+)(?=\.fastq\.gz)"
                            match = re.search(pattern, command_line)


                            if match:
                                sample_name = match.group(1)
                            else:
                                sample_name = "ERROR"
                                
                            print("=============")
                            print(sample_name)
                            print("=============")
                            break
                        except Exception as e:
                            print(param_name, ": Error retrieving input dataset info: ", e)

                resultDict[sample_name] = {'std_out': std_out, 'std_err': std_err}
                print(sample_name, std_out, std_err)
                print("=================================")
    return resultDict

def queryIrmaJobStdOutStdErr(history_client, dataset_client, history_id, tool_id = None):
    datasets = gi.histories.show_history(history_id, contents=True)
    resultDict = {}
    for dataset in datasets:
        if 'state' in dataset and dataset['state'] == "error":
            if 'id' in dataset:
                dataset_id = dataset['id']
                provenance_data = history_client.show_dataset_provenance(history_id, dataset_id)
                if tool_id:
                    if provenance_data["tool_id"] != tool_id:
                        print("Skipping")
                        continue
                std_err = provenance_data["stderr"]
                std_out = provenance_data["stdout"]
                
                for param_name, param_value in provenance_data['parameters'].items():
                    if isinstance(param_value, dict) and 'id' in param_value:
                        input_dataset_id = param_value['id']
                        try:
                            input_dataset = dataset_client.show_dataset(input_dataset_id)
                            sample_name = input_dataset['name'].split(':')[0]
                            break
                        except Exception as e:
                            print(param_name, ": Error retrieving input dataset info: ", e)

                resultDict[sample_name] = {'std_out': std_out, 'std_err': std_err}
    return resultDict


def writeIRMA(history_client, dataset_client, history_id):
    tool_id = "toolshed.g2.bx.psu.edu/repos/iuc/irma/irma/1.2.0+galaxy3"
    result = queryIrmaJobStdOutStdErr(history_client, dataset_client, history_id, tool_id)
    print("DONE", len(result.keys()))
    with open("stdOutOfFailedIrmaOfHistory_" + history_id, 'w') as f:
        f.write("Sample\tErrorType\n")
        for sample, value in result.items():
            if "found no QC'd data" in value['std_out']:
                f.write(sample + "\tNo Read survived QC\n")
            elif "R1 aborted" in value['std_out']:
                f.write(sample + "\tToo few survived QC\n")
            else:
                f.write(sample + "UNKOWNn")
                print(sample, value)


def writeVapor(history_client, dataset_client, history_ids):
    tool_id = "toolshed.g2.bx.psu.edu/repos/iuc/vapor/vapor/1.0.2+galaxy3"
    resultDict = {}
    for history_id in history_ids:
        print(history_id)
        result = queryVaporJobStdOutStdErr(history_client, dataset_client, history_id, tool_id)
        resultDict = resultDict | result
    
    print("DONE", len(resultDict.keys()))
    with open("stdOutOfFailedVaporOfMultiHistories.tsv", 'w') as f:
        f.write("Sample\tErrorType\n")
        for sample, value in result.items():
            if "lower -m thresh" in value['std_err']:
                f.write(sample + "\tVapor_m_thresh\n")
            elif "No virus found in" in value['std_err']:
                f.write(sample + "\tNoVirusFound\n")
            elif "Try a lower coverage cutoff -c" in value['std_err']:
                f.write(sample + "\tVapor_c_thresh\n")
            else:
                f.write(sample + "UNKOWNn")
                print(sample, value)


if __name__ == '__main__':
    history_client = HistoryClient(gi)
    dataset_client = DatasetClient(gi)
    history_id_irma = "80a0ec8f7d91eca5"
    history_ids_vapor = ['3ed8ba57b430d1fd'] # liste von den neuen vapor history ids, wenn fertig

    writeIRMA(history_client, dataset_client, history_id_irma)
    writeVapor(history_client,dataset_client, history_ids_vapor)
