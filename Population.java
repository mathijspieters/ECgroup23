import org.vu.contest.ContestEvaluation;
import java.util.Random;
import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Map;
import java.lang.Math;

public class Population{
    // Population Class with all evolutionary functions, population size, total fitness, diversity and
    // dictionary that describes the mutation function including necessary parameters we want to use

    int popSize_; // Population size
    double totFitness = 0; // Total Fitness of Population
    List<Individual> population_ = new ArrayList<Individual>();
    Random rnd_;
    GenUtils util_;
    String mutationType_;
    Map<String, Object> dict_;
    double diversity_;
    double best_fitness_;
    boolean DGEA = false;

    // Variables for correlated and uncorrelated mutation
    double tauPrime = 1/Math.sqrt(20);
    double tau = 1/Math.sqrt(2*Math.sqrt(10));
    double tauOne = 1/Math.sqrt(10);


    // Initialize Population
    public Population(int populationSize, Random rnd, GenUtils util, Map<String, Object> dict){
        popSize_ = populationSize;
    	rnd_ = rnd;
        util_ = util;
        mutationType_ = (String) dict.get("mutationType"); // Mutationfunction
        dict_ = dict; // Input for mutation function
        best_fitness_ = 0;
    }

    // Fill population with individuals
    public void init(){
        for (int i = 0; i < popSize_; i++){
            Individual ind = new Individual(util_);
            ind.setFitness();
            population_.add(ind);

            // Add fitness of individual to total fitness
            totFitness += ind.getFitness();
    	}

        // Set diversity of total population
        setDiversity();
    }

    public void setDiversity(){
        // Calculates diversity of population according to paper Diversity Guided EA's (Ursem)

        double diagonal = Math.sqrt(1000);
        double[] average = new double[10];
        double diversity = 0;

        for(int i = 0; i < 10; i++){

            average[i] = 0;

            for(int j = 0; j < popSize_; j++){
                average[i] += population_.get(j).getDna()[i];
            }

            average[i] = average[i]/popSize_;
        }

        for(int j = 0; j < popSize_; j++){
            double temp_div = 0;

            for(int i = 0; i < 10; i++){
                temp_div += Math.pow(population_.get(j).getDna()[i] - average[i],2);
            }

            diversity += Math.sqrt(temp_div);
        }

        diversity_ = diversity/(popSize_*diagonal);
    }

    public double getDiversity(){
        return diversity_;
    }

    public double getBestFitness(){
        return best_fitness_;
    }

    public int rouletteWheel(double tempFitness, List<Integer> indices){
        // Returns index from population list using roulette wheel selection
        // Remark: for the first generation this function doesn't work correctly due to rounding
        double partialSum = 0;

        // Randomly select value between 0 and total fitness
        double rand = util_.randDouble(0, tempFitness);

        // Loop through indivuduals, stop when partialSum exceeds rand
        for (int i=0; i < population_.size(); i++){

            // Ensure that index isn't already in index array
            if (indices.contains(i) == false){

                partialSum += population_.get(i).getFitness();

                // Index found
                if (partialSum >= rand){
                    return i;
                }
            }
        }

        // Default return
        return population_.size()-1;
    }

    public void selectParents(int numParents){
        // Calls roulettewheel selection numParents times, changes population into selected parents
        List<Integer> parentIndices = new ArrayList<Integer>();
        List<Individual> newParents = new ArrayList<Individual>();

        double tempFitness = totFitness;

        for (int i = 0; i < numParents; i++){

            // Select parent and add to indices
            parentIndices.add(rouletteWheel(tempFitness, parentIndices));

            newParents.add(population_.get(parentIndices.get(i)));

            // Update tempFitness
            tempFitness -= population_.get(parentIndices.get(i)).getFitness();
        }


        // Only selected parents survive
        population_ = newParents;
    }

    public void selectTournament(int numParents){
        //Use tournament selection to select parents.
        // Change to: with replacement.
        List<Individual> newParents = new ArrayList<Individual>();
        List<Integer> players = new ArrayList<Integer>();
        List<Integer> indices = new ArrayList<Integer>();
        //4-player rounds.
        int k = (int) (Integer) dict_.get("k");
        int index;
        double bestFitness;
        int bestIndex = 0;

        while (newParents.size()<numParents){
            bestFitness = 0;
            players.clear();
            for (int i = 0; i < popSize_; i++){
                indices.add(i);
            }
            while(players.size() < k){
                int indexRandom = rnd_.nextInt(indices.size());
                index = indices.get(indexRandom);
                indices.remove(indexRandom);
                players.add(index);
                if (population_.get(index).fitness_ > bestFitness){
                    bestIndex = index;
                    bestFitness = population_.get(index).fitness_;

                  }
            }
            newParents.add(population_.get(bestIndex));

        }
        population_ = newParents;
    }


    public void crossover(int kidSize, int parSize){
        // Per child, per parameter, randomly select dna of a parent
        for(int i = 0; i < kidSize; i++){
            Individual child = new Individual(util_);
            double[] dna = new double[10];
            double[] sigma = new double[10];
            double[] alpha = new double[45];

            for(int j = 0; j < 10; j++){
                // Randomly select dna of a single parent
                int randIndex = util_.randInt(0,parSize-1);
                dna[j] = population_.get(randIndex).getDna()[j];
                sigma[j] = population_.get(randIndex).getSigma()[j];
            }

            for(int j = 0; j < 45; j++){
                int randIndex = util_.randInt(0,parSize-1);
                alpha[j] = population_.get(randIndex).getAlpha()[j];
            }
            // Add child to population
            child.setDna(dna);
            child.setSigma(sigma);
            child.setAlpha(alpha);
            population_.add(child);

        }
        population_ = population_.subList(parSize, parSize +kidSize);


    }

    // sets all fitness values for population
    public void setFitnessPopulation(){
        for(int i = 0; i < popSize_; i++){
          population_.get(i).setFitness();

          if(population_.get(i).getFitness() > best_fitness_){
              best_fitness_ = population_.get(i).getFitness();
          }
        }
    }

    public void mutate(){
        // Calls mutation that is given as input, using switch

        switch( mutationType_) {
            case "mutateStandard":
                mutateStandard();
                break;
            case "mutateUncorrelated":
                mutateUncorrelated();
                break;
            case "mutateCorrelated":
                 mutateCorrelated();
                 break;
            default:
                // System.out.println("ERROR, undefined mutation type");
                break;
        }
    }

    public List<Individual> rank(List<Individual> pop){
        // Ranks population according to fitness

        pop.sort((o1, o2) -> Double.compare(o1.getFitness(),o2.getFitness()));
        Collections.reverse(pop);
        return pop;
    }

    public void selectSurvivors(){
        // Selects best individuals of generation as survivors using rank

        population_ = rank(population_);
        population_ = population_.subList(0, popSize_);

        totFitness = 0;
        for(int i = 0; i < popSize_; i++){
            // Update total fitness
            totFitness += population_.get(i).getFitness();
        }
    }

    public void mutateStandard(){
        // Mutates every parameter of every individual with a uniform sampled chance

        double step = (double) util_.getEvals() ;

        double variance = (double) (Double) dict_.get("variance");

        // Get function according to which mutation step size will decrease
        String mutationFunction = (String) dict_.get("mutationFunction");
        // Get distribution according to which mutation step size will be made
        String mutationDistribution = (String) dict_.get("mutationDistribution");
        // Set mutation according to function
        double mutation = 0;
        switch( mutationFunction){
            case "static":
                mutation = variance;
                break;

            case "linearDecreasing":
                mutation = variance - variance * step / (double) util_.getEvaluationLimit();
                break;

            case "sigmoidDecreasing":
                // Do something
                mutation = variance * (1 / (1 + Math.exp(step * (20 / util_.getEvaluationLimit()) - 10)));
                break;

            case "exponentialDecreasing":
                // Do something
                mutation = variance * Math.exp(- step / (double) util_.getEvaluationLimit());
                break;

            default:
                System.out.println("ERROR, undefined mutation function");
                break;
        }

        for (int j = 0; j < population_.size(); j++){

            // Mutate every parameter per individual
            for(int i = 0; i < 10; i++){

                // Sample delta gaussian, with mean zero and corresponding variance
                double delta =  0;
                switch(mutationDistribution){
                    case "uniform":
                        delta = util_.randDouble(-mutation, mutation);
                        break;

                    case "gaussian":
                        delta = util_.sampleGaussian(0, mutation);
                        break;

                    default:
                        System.out.println("ERROR, undefined mutation distribution");
                        break;
                }
                // Update dna
                population_.get(j).getDna()[i] = population_.get(j).getDna()[i] + delta;
            }
          if(!DGEA){
          // Update fitness
          population_.get(j).setFitness();

          if(population_.get(j).getFitness() > best_fitness_){
              best_fitness_ = population_.get(j).getFitness();
          }
        }
        }
    }
    public void mutateDG(double dLow, double dHigh, int par_size, int kid_size){
      String mode = "Exploit";
      DGEA = true;
      setDiversity();
      double div = getDiversity();
      if(div < dLow){
        mode = "Explore";

      }
      else if (div > dHigh) {
        mode = "Exploit";
      }

      if (mode == "Exploit"){

        // Select Parents
        selectTournament(par_size);

        // Do crossover
        crossover(kid_size, par_size);
        setFitnessPopulation();
        // Select survivors with best fitness
        selectSurvivors();

      }
      else{

      // Do mutation for all parents and generated children
      mutateStandard();
      }

        // Mutates every parameter of every individual with a uniform sampled chance

    }

    public void mutateUncorrelated(){
        //Apply uncorrelated 10-or-1-dimensional mutation.
        double[] sigma;
        double[] dna;
        double N;
        double Ni;
        String mutationFunction = (String) dict_.get("mutationFunction");

        for(int i=0;i<population_.size();i++){
            dna = population_.get(i).getDna();
            sigma = population_.get(i).getSigma();

            N = util_.sampleGaussian(0, 1);
            int dim = 0;

            switch( mutationFunction){
                case "One":
                    dim = 1;
                    break;
                case "Ten":
                    dim = 10;
                    break;

                default:
                    System.out.println("ERROR, undefined mutation function");
                    break;
            }

            for(int j=0;j<dim;j++){
                //Update sigma.
                Ni = util_.sampleGaussian(0, 1);

                switch( mutationFunction){
                    case "One":
                        sigma[j] = sigma[j]*Math.exp(tauOne*N);
                        break;
                    case "Ten":
                        sigma[j] = sigma[j]*Math.exp(tauPrime*N + tau*Ni);
                        break;
                }

                //Update dna.
                Ni = util_.sampleGaussian(0, 1);
                dna[j] += sigma[j]*Ni;

                //Check boundaries.
                if (dna[j] < -5){
                    dna[j] = -5;
                }
                else if(dna[j] > 5){
                    dna[j] = 5;
                }
            }
            //Apply the changes.
            population_.get(i).setDna(dna);
            population_.get(i).setSigma(sigma);
            population_.get(i).setFitness();

            if(population_.get(i).getFitness() > best_fitness_){
                best_fitness_ = population_.get(i).getFitness();
            }
        }
    }

        public double[][] setCovariance(int k){
        // Initializes covariance matrix for an individual k as described on page 61 and 62 of the book
        double[] sigma = population_.get(k).getSigma();
        double[] dna = population_.get(k).getDna();
        double[] alpha = population_.get(k).getAlpha();

        double N; double Ni; double Nj;

        double beta = 5;

        double[][] covarianceMatrix = new double[10][10];


        N = util_.sampleGaussian(0,1);
        // Loop through covariancematrix rows

        int count = 0;

        for (int i = 0; i < 10; i++){
            Ni = util_.sampleGaussian(0,1);

            // Update sigma, according to page 62 book
            sigma[i] = sigma[i] * Math.exp(tauPrime *N + tau * Ni);

            // Loop through covariancematrix columns
            for (int j = i; j < 10; j++){

                Nj = util_.sampleGaussian(0,1);

                // Fill covariancematrix
                if (i == j){
                    // If diagonal, only variance (sigma_j^2)
                    covarianceMatrix[i][j] = Math.pow(sigma[j], 2);
                }
                else {
                    // Else, covariance between i and j

                    alpha[count] = alpha[count] + beta * Nj;

                    if(Math.abs(alpha[count]) >  Math.PI){
                        if(alpha[count] > 0){
                            alpha[count] = alpha[count] - 2 * Math.PI;
                        }else{
                            alpha[count] = alpha[count] + 2 * Math.PI;
                        }
                    }

                    covarianceMatrix[i][j] = (1 / 2) * (Math.pow(sigma[i],2) - Math.pow(sigma[j],2)) * Math.tan(Math.toRadians(2*alpha[count]));
                    covarianceMatrix[j][i] = covarianceMatrix[i][j];
                    count = count + 1;
                }


            }
        }

        // Update sigma's of individual
        population_.get(k).setSigma(sigma);

        population_.get(k).setAlpha(alpha);

        // Return
        return covarianceMatrix;
    }

    public void mutateCorrelated(){
        // Update every individual with a step drawn from a multivariate distribution
        double[][] covarianceMatrix;
        double[] mean = new double[10];
        double[] dna;


        // Loop through population
        for(int k = 0; k < popSize_; k++){
            // Get Covariance Matrix
            covarianceMatrix = setCovariance(k);

            // Create multivariate normal distribution object
            double[] mutationstep = MatrixFunctions.multivariateGaussian(mean, covarianceMatrix, util_);

            dna = population_.get(k).getDna();

            // Mutate individual parameters bij adding sampled vector from multivariate distribution
            for(int i=0; i<10; i++){
                dna[i] = dna[i] + mutationstep[i];
                if (dna[i] < -5){
                    dna[i] = -5;
                }
                else if(dna[i] > 5){
                    dna[i] = 5;
                }
            }

            // Update indidivuals
            population_.get(k).setDna(dna);
            population_.get(k).setFitness();

            if(population_.get(k).getFitness() > best_fitness_){
                best_fitness_ = population_.get(k).getFitness();
            }


        }
    }



}
