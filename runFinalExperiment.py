import subprocess
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams.update({'font.size': 14})
import os
from tqdm import tqdm
from glob import glob
import numpy as np

SAVE_FIGURES = False

def get_diversity(output):
    o = str(output).split("\\n")[:-3]
    return o

def get_score(output):
    o = str(output).split("\\n")[-3] # get the score field
    o = float(o.split(" ")[1])
    return o


def save_plot(lines, directory, run):
    lines = [l for l in lines if len(l) > 0 and l[0].isdigit()]
    lines = [l.split(",") for l in lines]

    fig, ax1 = plt.subplots()

    t = [int(l[0]) for l in lines]
    d = [float(l[1]) for l in lines]
    s = [float(l[2]) for l in lines]

    s = [max(s[:i]) for i in range(1,len(s)+1)]

    ax1.semilogy(t, d, 'b-')
    ax1.set_xlabel('iterations')
    ax1.set_ylabel('diversity', color='b')
    ax1.tick_params('y', colors='b')

    ax2 = ax1.twinx()
    ax2.plot(t, s, 'r')
    ax2.set_ylabel('score', color='r')
    ax2.tick_params('y', colors='r')

    fig.tight_layout()
    fig.savefig(directory+str(run)+".png",bbox_inches='tight')
    plt.clf()


def read_txt(file):
    with open(file, "r") as f:
        doc = f.read()

    scores = [x for x in doc.split("\n")]
    scores = [float(x) for x in scores if len(x) > 0]
    return scores

def show_all_stats():
    directory = "results/"

    files = [y for x in os.walk(directory) for y in glob(os.path.join(x[0], '*.txt'))]

    for file in files:
        params = file.split("/")[-2]
        problem = file.split("/")[1]
        scores = read_txt(file)

        print(params)
        print(problem)
        print("MEAN: %.5f" % np.mean(scores))
        print("STD: %.5f" % np.std(scores))


def runExperiment(mut_args, args_dict, f):
    args = []
    for key in args_dict.keys():
        args.append(key + "=" + str(args_dict[key]))


    num_iterations = 100
    plot_limit = 20

    directory = "results/%s/%s/" % (f, "-".join(mut_args+args))

    if not os.path.exists(directory):
        os.makedirs(directory)

    scores = []

    for i in tqdm(range(num_iterations-1)):
        test = subprocess.Popen(["java"] + mut_args + args + ["-jar", "testrun.jar", "-submission=player23", "-evaluation=%s" % f, "-seed=%d" % i], stdout=subprocess.PIPE)
        output = test.communicate()[0]

        if i < plot_limit and SAVE_FIGURES:
            save_plot(get_diversity(output), directory, i)

        scores += [get_score(output)]

    with open(directory+"results.txt", "w") as f:
        for score in scores:
            f.write("%s\n" % score)


def run_experiment(function, alg):
    if alg == "LDM":
        if function == "BentCigarFunction":
            params = {'-Dmu': 80, '-Dlambda': 169, '-Dvariance': 0.003817733496358963, '-Dk': 6, '-Dparsize': 64}
            mut_args = ["-DmutationType=mutateStandard", "-DmutationFunction=linearDecreasing", "-DmutationDistribution=gaussian"]
        elif function == "SchaffersEvaluation":
            params = {'-Dmu': 67, '-Dlambda': 186, '-Dvariance': 0.02931366086141189, '-Dk': 8, '-Dparsize': 67}
            mut_args = ["-DmutationType=mutateStandard", "-DmutationFunction=linearDecreasing", "-DmutationDistribution=gaussian"]
        else: #KatsuuraEvaluation
            params = {'-Dmu': 17, '-Dlambda': 152, '-Dvariance': 0.0010615728075645916, '-Dk': 4, '-Dparsize': 4}
            mut_args = ["-DmutationType=mutateStandard", "-DmutationFunction=linearDecreasing", "-DmutationDistribution=gaussian"]
    elif alg == "DGEA":
        if function == "BentCigarFunction":
            params = {'-Dmu': 80, '-Dlambda': 114, '-DdLow': 0.0020283390188316447, '-DdHigh': 0.3126058325379971, '-Dk': 9, '-Dparsize': 23}
            mut_args = ['-DmutationType=mutateDG', '-DmutationFunction=static', '-DmutationDistribution=gaussian']
        elif function == "SchaffersEvaluation":
            params = {'-Dmu': 128, '-Dlambda': 128, '-DdLow': 0.0022876039949272374, '-DdHigh': 0.23176941855081462, '-Dk': 8, '-Dparsize': 35}
            mut_args = ['-DmutationType=mutateDG', '-DmutationFunction=static', '-DmutationDistribution=gaussian']
        else: #KatsuuraEvaluation
            params = {'-Dmu': 93, '-Dlambda': 217, '-DdLow': 0.016926841434103237, '-DdHigh': 0.23935788204849787, '-Dk': 9, '-Dparsize': 46}
            mut_args = ['-DmutationType=mutateDG', '-DmutationFunction=static', '-DmutationDistribution=gaussian']
    else: #CMEAS
        params = {'-Dmu': 40}
        mut_args = ["-DmutationType=cmaes"]

    runExperiment(mut_args, params, function)


if __name__ == "__main__":
    f = "KatsuuraEvaluation" #KatsuuraEvaluation SchaffersEvaluation BentCigarFunction
    a = "CMEAS" #LDM DGEA CMEAS

    #run_experiment(f,a)

    show_all_stats()
