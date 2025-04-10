import os
import re
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

try:
    import pandas as pd
except Exception:
    pass


csv_filepath = "stats.csv"
log_dir = "./logs"

lfn_pattern = r"simplex_d(\d+)_n(\d+)_([a-zA-Z0-9]+)_l(\d+)_L(\d+)\.log"
line_pattern = r".*Tital time for\s*([a-zA-Z])\s*:\s(\d+)"

sections_local = [
        "Total Time",
        "Iteration Step",
        "owned_query_reply",
        "sort_peers",
        "assign_element_nodes",
        "collect_local_tets",
        "populate_sharers",

        # "Iter Volume",
        # "Iter Face",
        # "Iter Edge",
        # "Iter Corner",
        ]

sections_comm = [
        "post_query_reply",
        "wait_query_reply",
        "allgather_counts",

]

sections_peers = [
        "sort_peers",
]



def list_files(dir_path):
    for f in os.listdir(dir_path):
        f_path = os.path.join(dir_path, f)
        if os.path.isfile(f_path) and f.endswith(".log"):
            yield f

def parse_lfile(lf):
    values = {}
    with open(os.path.join(log_dir, lf)) as f:
        for l in f:
            if "global triangles" in l:
                m = re.match(r".*global triangles (\d+).*nodes (\d+)", l)
                vals = m.groups()
                values['global tets'] = vals[0]
                values['global nodes'] = vals[1]
                # print(l, m.groups())
            if "Total time for" in l:
                m = re.match(r".* for \[([\w\s]+)\]: ([\d\.]+)", l)
                if not m:
                    continue
                # print(l, m.groups())
                values[m.groups()[0]] = float(m.groups()[1])
            if "Done p8est_balance" in l:
                m = re.match(r".* with ([\d\.]+) total quadrants", l)
                if not m:
                    continue
                # print(l, m.groups())
                values['global quads'] = float(m.groups()[0])


    return values
            # match = re.match(line_pattern, l)
            # print(match)

            # if match:
            #     print(l)

def extract_from_log_files():
    keys = set()
    data = []
    for lf in list_files(log_dir):
        # simplex_d%d_n%d_%s_l%d_L%d.log
        match = re.match(lfn_pattern, lf)
        if match:
            d, n, conn, l, L = match.groups()
            values = {
                "logfilepath": lf,
                "dim": d, "procn": n, "conn": conn,
                "l": l,  "L": L }
            values |= parse_lfile(lf)
            keys |= values.keys()
            data.append(values)

            # print(lf, values)
    return keys, data

def write_to_csv(keys, data):
    # write_to_csv(keys, data);
    with open(csv_filepath, "w") as f:
        wrtier = csv.DictWriter(f, fieldnames=keys)
        wrtier.writeheader()
        wrtier.writerows(data)

def plot_against_nodes(df):
    df = df.sort_values(by="global nodes")

    nodes = df["global nodes"]
    tets = df["global tets"]
    ttime = df["Total Time"]
    istime = df["Iteration Step"]

    # Plotting Total Time against global tets
    plt.figure(figsize=(8, 5))  # Optional: Adjust the figure size

    plt.plot(nodes, ttime, marker="o", linestyle="-", color="b", label="Total Time")
    plt.plot(nodes, istime, marker="o", linestyle="-", color="b", label="Iteration Step")

    # Adding labels and title
    plt.xlabel("Global Node count", fontsize=12)
    # plt.ylabel("Total Time (seconds)", fontsize=12)
    # plt.title("Total Time vs. Global Tets", fontsize=14)
    plt.legend()  # Show legend
    # plt.grid(True)  # Optional: Add grid lines
    plt.grid(False)  # Optional: Add grid lines
    plt.tight_layout()  # Optional: Adjust layout to fit everything nicely

    plt.show()


def plot_loglog(df, xs, sections):
    uniform_cond = df["l"] == df["L"] # Select uniform

    # Create the figure
    fig, ax = plt.subplots()

    vmin=10**5
    vmax=10**-5

    # Plot all points
    for section in sections:
        ax.loglog(xs, df[section], marker="o", label=section)
        vmin = min(vmin, df[section].min())
        vmax = max(vmax, df[section].max())

    # Mark uniform points
    # ax.loglog(n[uniform_cond], ttime[uniform_cond], marker="o", linestyle="None", color="r", label="Uniform")

    # Plot log lines
    ax.set_clip_on(True)
    n_lines = 5;
    for k in range(-15, 10):
        y_background = 10**k * xs # * np.log(xs)
        ax.plot(
            xs, y_background, '--', color="grey", alpha=0.5,
            zorder=-1, clip_on=True  # Ensure lines are in the background
        )


    ax.set_xlim(xs.min(), xs.max())
    ax.set_ylim(vmin / 2, vmax * 2)

    return plt



# def plot_against_procn():


def plot_time_cpx(df):
    quads = df["global quads"]

    plot_loglog(df, quads, sections_comm)
    # plot_loglog(df, quads, sections_local)

    plt.xlabel("Global Quads (log_10 count)", fontsize=12)
    plt.ylabel("Total Time (log_10 seconds)", fontsize=12)
    plt.title("Log-Log Plot: Total Time vs. Global Quads", fontsize=14)

    # plt.ylabel("Global Node count", fontsize=12)
    # plt.title("Total Time vs. Global Tets", fontsize=14)
    plt.legend()  # Show legend
    plt.grid(False)  # Optional: Add grid lines
    plt.tight_layout()  # Optional: Adjust layout to fit everything nicely

    plt.show()

def plot_ncproc(df, sections=None):

    df = df.sort_values(by="procn")

    sections = sections_local + sections_comm


    xs = df["procn"]

    print(xs)

    plt = plot_loglog(df, xs, sections)

    # for section in sections:
    #     ax.loglog(xs, df[section], marker="o", label=section)

    plt.legend()
    plt.xlabel("Process count (log_10 #)", fontsize=12)
    plt.ylabel("Total Time (log_10 seconds)", fontsize=12)
    plt.title("Plot: Total Time vs. Process count", fontsize=14)

    plt.show()

def main(options, args):

    if options.parse_logs:
        keys, data = extract_from_log_files()
        print(keys)
        df = pd.DataFrame(data)
        df.to_csv(csv_filepath, index=False)
    df = pd.read_csv(csv_filepath)

    df = df[df["conn"] == "sphere"]

    selection = (df["l"] == df["L"])
    df = df[ selection ]
    plot_time_cpx(df)

    # selection = (df["l"] == 5) * (df["L"] == 8)
    # df = df[ selection ]
    # plot_ncproc(df)


import sys, os, optparse

__description = '''
    A simple script for parsing the statistics in the log files
'''

if __name__ == "__main__":
    usage = "%prog [options] <arg1>"
    parser = optparse.OptionParser(usage = usage, description=__description)
    parser.add_option('-l', dest='parse_logs', default=False, action='store_true', help='enable debug output')

    (options, args) = parser.parse_args()

    main(options, args)
