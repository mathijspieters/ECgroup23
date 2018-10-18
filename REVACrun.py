from REVAC import run

if __name__ == "__main__":
    """ This bit of code allows for automated runs """


    #External REVAC parameters
    population_size = 100 #Number of parameter vectors.
    generations = 100 #Number of new vectors tested.

    params = {"-Dmu":[10, 100], "-Dlambda":[1, 200], "-Dvariance":[0.0001, 1], "-Dk":[2, 10], "-Dparsize":[2, 100]} #"-DdLow":[0.00005, 0.04], "-DdHigh":[0.07, 0.5], "-Dvariance":[0.0001, 1]
    mut_args = ["-DmutationType=mutateStandard", "-DmutationFunction=linearDecreasing", "-DmutationDistribution=gaussian"]
    f = "KatsuuraEvaluation" #KatsuuraEvaluation #SchaffersEvaluation #BentCigarFunction
    player = "player23Tweak"



    for i in range(10):
        run(population_size, generations, params, mut_args, f, player)
