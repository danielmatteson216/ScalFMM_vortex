#!/usr/bin/python
import getopt
import sys
import math
import copy
import os
import socket
import subprocess
import re
import types
import glob

class ScalFMMConfig(object):
    num_threads    = 1
    num_nodes      = 1
    algorithm      = "implicit"
    model          = "cube"
    num_particules = 10000
    height         = 4
    bloc_size      = 100
    order          = 5

    def show(self):
        print ("=== Simulation parameters ===")
        print ("Number of nodes: " + str(self.num_nodes))
        print ("Number of threads: " + str(self.num_threads))
        print ("Model: " + str(self.model))
        print ("Number of particules: " + str(self.num_particules))
        print ("Height: " + str(self.height))
        print ("Bloc size: " + str(self.bloc_size))
        print ("Order: " + str(self.order))

    def gen_header(self):
        columns = [
            "model",
            "algo",
            "nnode",
            "nthreads",
            "npart",
            "height",
            "bsize",
            "global_time",
            "runtime_time",
            "task_time",
            "idle_time",
            "scheduling_time",
            "communication_time",
            "communication_vol",
            "rmem",
        ]
        header = ""
        for i in range(len(columns)):
            if not i == 0:
                header += ","
            header += "\"" + columns[i] + "\""
        header += "\n"
        return header

    def gen_header_gantt(self):
        columns = [
            "model",
            "algo",
            "nnode",
            "nthreads",
            "npart",
            "height",
            "bsize",
            "filename",
            "start_profiling",
            "stop_profiling",
        ]
        header = ""
        for i in range(len(columns)):
            if not i == 0:
                header += ","
            header += "\"" + columns[i] + "\""
        header += "\n"
        return header

    def gen_record(self, global_time, runtime_time, task_time, idle_time, scheduling_time, communication_time, communication_vol, rmem):
        columns = [
            self.model,
            self.algorithm,
            self.num_nodes,
            self.num_threads,
            self.num_particules,
            self.height,
            self.bloc_size,
            global_time,
            runtime_time,
            task_time,
            idle_time,
            scheduling_time,
            communication_time,
            communication_vol,
            rmem,
        ]
        record = ""
        for i in range(len(columns)):
            if not i == 0:
                record += ","
            if (type(columns[i]) is bool or
                type(columns[i]) == str):
                record += "\""
            record += str(columns[i])
            if (type(columns[i]) == bool or
                type(columns[i]) == str):
                record += "\""
        record += "\n"
        return record

    def gen_record_gantt(self, filename, start_profiling, stop_profiling):
        columns = [
            self.model,
            self.algorithm,
            self.num_nodes,
            self.num_threads,
            self.num_particules,
            self.height,
            self.bloc_size,
            filename,
            start_profiling,
            stop_profiling,
        ]
        record = ""
        for i in range(len(columns)):
            if not i == 0:
                record += ","
            if (type(columns[i]) is bool or
                type(columns[i]) == str):
                record += "\""
            record += str(columns[i])
            if (type(columns[i]) == bool or
                type(columns[i]) == str):
                record += "\""
        record += "\n"
        return record

def get_times_from_trace_file(filename):
    return 1.0, 1.0, 1.0, 1.0, 1.0
    cmd = "starpu_trace_state_stats.py " + filename
    proc = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    if not proc.returncode == 0:
        sys.exit("FATAL: Failed to parse trace.rec!")
        return proc.returncode
    task_time = 0.0
    idle_time = 0.0
    runtime_time = 0.0
    scheduling_time = 0.0
    communication_time = 0.0
    for line in stdout.decode().splitlines():
        arr = line.replace("\"", "").split(",")
        if arr[0] == "Name":
            continue
        if len(arr) >= 4:
            value = float(arr[3])
            if value < 0:
                print("\033[31mError negative " + arr[2] + ":" + arr[0] + " -> " + arr[3] + "\033[39m")
            value = math.fabs(value)
            if arr[2] == "Runtime":
                if arr[0] == "Scheduling":
                    scheduling_time += value
                else:
                    runtime_time += value
            elif arr[2] == "Task":
                task_time += value
            elif arr[2] == "Other":
                idle_time += value
            elif arr[2] == "MPI":
                communication_time += value
            elif arr[2] == "User":
                if arr[0] == "Decoding task for MPI":
                    communication_time += value
                elif arr[0] == "Preparing task for MPI":
                    communication_time += value
                elif arr[0] == "Post-processing task for MPI":
                    communication_time += value
                else:
                    runtime_time += value
            else:
                print("Error type " + arr[2])
            # sys.exit("Invalid time!")
    return runtime_time, task_time, idle_time, scheduling_time, communication_time

def generate_gantt(config, gantt_filename, initial_dir):
    if (os.path.isfile(gantt_filename)):
        gantt_file = open(gantt_filename, "a")
    else:
        gantt_file = open(gantt_filename, "w")
        gantt_file.write(config.gen_header_gantt())
    csv_paje_filename = initial_dir + "/paje.csv"
    simple_paje_filename = initial_dir + "/paje.trace"
    if not os.path.exists(csv_paje_filename):
        return

    proc = subprocess.Popen(['wc', '-l', csv_paje_filename], stdout=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    line = stdout.decode().splitlines()[0].split(" ")[0]
    if int(line) <= 2:
        print("There is too m few line in " + csv_paje_filename)
        return

    #Get start profiling time
    cmd = "grep start_profiling " + simple_paje_filename
    proc = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    if not proc.returncode == 0:
        sys.exit("FATAL: Failed to parse " + simple_paje_filename + "!")
        return proc.returncode
    start_profiling = float("inf")
    for line in stdout.decode().splitlines():
        arr = line.split()
        start_marker = float(arr[1])
        if(start_marker < start_profiling):
            start_profiling = start_marker
    
    #Get stop profiling time
    cmd = "grep stop_profiling " + simple_paje_filename
    proc = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    if not proc.returncode == 0:
        sys.exit("FATAL: Failed to parse " + simple_paje_filename + "!")
        return proc.returncode
    stop_profiling = 0.0
    for line in stdout.decode().splitlines():
        arr = line.split()
        stop_marker = float(arr[1])
        if(stop_marker > stop_profiling):
            stop_profiling = stop_marker

    gantt_file.write(config.gen_record_gantt(csv_paje_filename, start_profiling, stop_profiling))

def main():
    output_trace_file=""
    trace_filename="trace.rec"
    output_filename="loutre.db"
    gantt_database=""
    only_gantt=False
    long_opts = ["help",
                 "trace-file=",
                 "output-trace-file=",
                 "gantt-database=",
                 "output-file=",
                 "only-gantt"]

    opts, args = getopt.getopt(sys.argv[1:], "ht:i:o:g:", long_opts)
    for o, a in opts:
        if o in ("-h", "--help"):
            # usage()
            print("No help")
            sys.exit()
        elif o in ("-t", "--trace-file"):
            trace_filename = str(a)
        elif o in ("-i", "--output-trace-file"):
            output_trace_file = str(a)
        elif o in ("-o", "--output-file"):
            output_filename = str(a)
        elif o in ("-g", "--gantt-database"):
            gantt_database = str(a)
        elif o in ("--only-gantt"):
            only_gantt = True
        else:
            assert False, "unhandled option"

    config=ScalFMMConfig()
    rmem = 0
    global_time = -1.0
    runtime_time = 0.0
    task_time = 0.0
    idle_time = 0.0
    scheduling_time = 0.0
    communication_time = 0.0
    communication_vol = 0.0

    if (os.path.isfile(output_filename)): #Time in milli
        output_file = open(output_filename, "a")
    else:
        output_file = open(output_filename, "w")
        output_file.write(config.gen_header())

    with open(output_trace_file, "r") as ins:
        for line in ins:
            if re.search("Average", line):
                a = re.findall("[-+]?\d*\.\d+|\d+", line)
                if len(a) == 1:
                    global_time = float(a[0])*1000 # Else it is in sec
            elif re.search("Total particles", line):
                a = re.findall("[-+]?\d*\.\d+|\d+", line)
                if len(a) == 1:
                    config.num_particules = int(a[0])
            elif re.search("Group size", line):
                a = re.findall("[-+]?\d*\.\d+|\d+", line)
                if len(a) == 1:
                    config.bloc_size = int(a[0])
            elif re.search("Nb node", line):
                a = re.findall("[-+]?\d*\.\d+|\d+", line)
                if len(a) == 1:
                    config.num_nodes = int(a[0])
            elif re.search("Tree height", line):
                a = re.findall("[-+]?\d*\.\d+|\d+", line)
                if len(a) == 1:
                    config.height = int(a[0])
            elif re.search("Nb thread", line):
                a = re.findall("[-+]?\d*\.\d+|\d+", line)
                if len(a) == 1:
                    config.num_threads = int(a[0])
            elif re.search("Model", line):
                config.model = line[line.index(":")+1:].strip()
            elif re.search("Algorithm:", line):
                config.algorithm = line[line.index(":")+1:].strip()
            elif re.search("TOTAL", line) and re.search("starpu_comm_stats", line):
                a = re.findall("(\d*\.\d+|\d+).MB", line)
                if len(a) == 1:
                    communication_vol += float(a[0])
                
    if config.num_nodes > 1 and communication_vol == 0:
        communication_vol = -1

    if(gantt_database != ""):
        generate_gantt(config, gantt_database, os.path.dirname(trace_filename))

    if not only_gantt:
        print("Generating time ...")
        if (os.path.isfile(trace_filename)): #Time in milli
            runtime_time, task_time, idle_time, scheduling_time, communication_time = get_times_from_trace_file(trace_filename)
        else:
            print("File doesn't exist " + trace_filename)

        # sum_time = (runtime_time + task_time + scheduling_time + communication_time + idle_time)/(config.num_nodes*config.num_threads)
        # diff_time = float('%.2f'%(abs(global_time-sum_time)/global_time))

        # if diff_time > 0.01:   
            # print('\033[31m/!\\Timing Error of ' + str(diff_time) + '\033[39m')
            # print('\033[31m Global ' + str(global_time) + ' Sum ' + str(sum_time) + '\033[39m')
            # print('\033[31m Nodes number ' + str(config.num_nodes) + ' CPU ' + str(config.num_threads) + '\033[39m')

        # Write a record to the output file.
        output_file.write(config.gen_record(global_time,
                          float(runtime_time),
                          float(task_time),
                          float(idle_time),
                          float(scheduling_time),
                          float(communication_time),
                          float(communication_vol),
                          int(rmem)))

main()
