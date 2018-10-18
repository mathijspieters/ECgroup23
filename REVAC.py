import numpy as np
import subprocess
from bisect import bisect
from tqdm import tqdm
import time
import os


class x():
    """ This class stores candidate parameters """
    def __init__(self, params):
        self.params = params
        self.dic = {}
        self.fitness = 0
        self.scores = []

        #Used to determine necessary evaluations.
        self.max_evals = 1
        self.eval_count = 0

    def initialise(self):
        """ Initialises a parameter vector between specified limits """

        for key in self.params.keys():
            if type(self.params[key][0]) == type(0):
                self.dic[key] = np.random.randint(self.params[key][0], self.params[key][1])
            elif type(self.params[key][0]) == type(1.0):
                self.dic[key] = (np.random.random()*(self.params[key][1] - self.params[key][0]) + self.params[key][0])

        #Check conditions among variables.
        if self.dic["-Dmu"] < self.dic["-Dparsize"]:
            self.dic["-Dparsize"] = np.random.randint(self.params["-Dparsize"][0] + 1, self.dic["-Dmu"])
        if self.dic["-Dmu"] > self.dic["-Dlambda"]:
            self.dic["-Dlambda"] = np.random.randint(self.dic["-Dmu"] + 1, self.params["-Dlambda"][1])


    def __str__(self):
        """ String representation """
        return str(self.dic)

    def get(self):
        """ Returns personal parameter dictionary """
        return self.dic

    def get_std(self):
        return np.std(self.scores)

    def get_xi(self):
        """ Xi value of hoeffding inequality """
        return np.sqrt(1*np.log(2./0.2)/(2.*self.eval_count))

    def set_fitness(self):
        """ Set new fitness value """
        self.fitness = np.mean(self.scores)
        self.eval_count += 1

    def more_to_go(self):
        """ Determine whether to continue """
        return self.eval_count < self.max_evals

class Population():
    """ This class stores a population of parameter vectors. It also creates
        new candidates and evaluates them. It might also calculate D and Shannon
        entropy, I'm not sure yet. """

    def __init__(self, m, n, params, mut_args, f, player):
        self.m = m #Population size.
        self.n = n #Parent size.
        self.mut_args = mut_args
        self.params = params
        self.f = f
        self.player = player
        self.pop = []
        self.parents = []
        self.children = []
        self.max_evals = 3
        self.best_fitness = 0
        self.best = -100000


        #Create write file.
        year, month, day, hour, minute = \
            time.strftime("%Y,%m,%d,%H,%M").split(',')
        self.file = open("REVAC/" + self.f + " || " + str(self.mut_args) + year + "-" + month +\
            "-" + day + " " + hour + "-" + hour + "-" + minute + ".txt", "a+")

        #Create population.
        for i in range(self.m):
            eve = x(self.params)
            eve.initialise()
            self.pop.append(eve)

        self.evaluate(self.pop)
        self.rank(self.pop)

        print("POPULATION INITIALISED [highscore:%0.6f]" % self.pop[0].fitness)

    def __str__(self):
        """ String representation """
        output = []
        for i in range(self.m):
            output.append(self.pop[i].get())
        return str(output)

    def rank(self, pop):
        """ Sort by fitness decreasingly """

        pop.sort(key=lambda x: x.fitness, reverse=True)

    def sort(self, pop, key):
        """ Sort by parameter increasingly """

        pop.sort(key=lambda x: x.dic[key])

    def next(self, offspring, h):
        self.breed(offspring)
        self.mutate(h)


    def breed(self, offspring):
        """ Create offspring by uniform scanning """

        #Parents are the best n individuals.
        for i in range(offspring):
            child = x(self.params)
            for key in self.params.keys():
                parent = self.pop[np.random.randint(0, self.n)]
                val = parent.dic[key]
                child.dic[key] = val
            self.children.append(child)

    def mutate(self, h):
        """ Mutate according to h-neighbour uniform """
        self.parents = self.pop[:self.n]

        if 2*h > self.n:
            print("You have a problem with your h. Better fix it quick.")

        #Lots of loops to determine mutation domain.
        for child in self.children:
            for key in self.params.keys():
                #List to hold xi's.
                vals = []
                for parent in self.parents:
                    vals.append(parent.dic[key])
                vals.sort()
                index = vals.index(child.dic[key])

                #Determine domain.
                if index < h:
                    #Low end mirroring.
                    high = vals[index + h]
                    low = 2*child.dic[key] - high
                elif index + h > (self.n - 1):
                    #High end mirroring.
                    low = vals[index - h]
                    high = 2*child.dic[key] - low
                else:
                    low = vals[index - h]
                    high = vals[index + h]

                if low < 0:
                    if type(low) == type(0):
                        low = 1
                    elif type(low) == type(1.0):
                        low = 0

                #Mutate child uniformly.
                if type(child.dic[key]) == type(0):
                    child.dic[key] = np.random.randint(low, high+1)
                elif type(child.dic[key]) == type(1.):
                    child.dic[key] = np.random.random()*(high - low) + low

            #Sanity checks.
            if child.dic["-Dmu"] < child.dic["-Dparsize"]:
                child.dic["-Dparsize"] = max([child.dic["-Dmu"] - 1, 2])
            if child.dic["-Dmu"] > child.dic["-Dlambda"]:
                child.dic["-Dlambda"] = child.dic["-Dmu"]




        #Test children.
        self.evaluate(self.children)
        self.rank(self.children)

        #Deterministically select child.
        self.pop[-1] = self.children[0]
        self.children = []

        self.rank(self.pop)
        self.write()


    def write(self):
        self.file.write(str(self.pop[0]) + ", " + str(self.pop[0].fitness) + ", " +str(self.pop[0].get_std()) + "\n")

    def more_to_go(self, group):
        """ Determine whether anyone has the right to go on """
        for individual in group:
            if individual.more_to_go():
                return True
        return False

    def evaluate(self, group):
        while self.more_to_go(group):
            for individual in group:
                if not individual.more_to_go():
                    continue
                args = []
                for key in individual.dic.keys():
                    args.append(key + "=" + str(individual.dic[key]))
                    #Evaluate once more.

                #Evaluate once more.
                test = subprocess.Popen(["java"] + self.mut_args + args + ["-jar", "testrun.jar", "-submission=%s" %self.player, "-evaluation=%s" %self.f, "-seed=%d" % individual.eval_count], stdout=subprocess.PIPE)
                output = test.communicate()[0]
                try:
                    individual.scores.append(get_score(output))
                except:
                    print("TROUBLEMAKER: ", str(individual))
                individual.set_fitness()

                #Update best fitness for Sharpening.
                if individual.fitness > self.best_fitness:
                    self.best_fitness = individual.fitness
                    #SHARPENING:
                    self.max_evals = 2 + int(0.1*self.best_fitness**2)

                #Update best minimal score
                if individual.fitness - individual.get_xi() > self.best:
                    self.best = individual.fitness - individual.get_xi()

                #RACING: Increase evaluations if there's a good chance this one's a winner!
                if individual.fitness + individual.get_xi() > self.best:
                    if individual.max_evals < self.max_evals:
                        individual.max_evals += 1
                else:
                    individual.eval_count = individual.max_evals
                    #print("eliminated", individual.fitness)


        # print("%s --- %s" % (f,str(mut_args + args)))
        # print("Min:  %5f" % np.min(scores))
        # print("Max:  %5f" % np.max(scores))
        # print("Mean: %5f" % np.mean(scores))
        # print("Std:  %5f" % np.std(scores))

def get_score(output):
	o = str(output).split("\\n")[-3] # get the score field
	o = float(o.split(" ")[1])
	return o

def run(population_size, generations, params, mut_args, f, player):
    """ Automated runner function """
    #Create data directory.
    cwd = os.getcwd()
    if not os.path.exists(cwd+"/REVAC"):
        os.makedirs(cwd+"/REVAC")

    print("Thank you for choosing REVAC. We will be testing " + str(mut_args) + " on the " + f +". Please stand by.")
    print()

    Henk = Population(population_size, int(population_size/2.), params, mut_args, f, player)
    for i in tqdm(range(generations)):
        Henk.next(1, int(population_size/10.))

    Henk.file.close()
    print("Best fitness: ", Henk.pop[0].fitness, Henk.pop[0].get_std())
    print(Henk.pop[0].dic)


if __name__ == "__main__":
    #External REVAC parameters
    population_size = 100 #Number of parameter vectors.
    generations = 100 #Number of new vectors tested.

    params = {"-Dmu":[10, 100], "-Dlambda":[1, 200], "-Dvariance":[0.0001, 1.], "-Dk":[2, 10], "-Dparsize":[2, 100]}
    mut_args = ["-DmutationType=mutateUncorrelated", "-DmutationFunction=Ten", "-DmutationDistribution=uniform"]
    f = "BentCigarFunction" #KatsuuraEvaluation #SchaffersEvaluation #BentCigarFunction
    player = "player23Tweak"

    run(population_size, generations, params, mut_args, f, player)
