import org.vu.contest.ContestEvaluation;
import java.util.Random;
import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Map;
import java.lang.Math;

public class PopulationCmaes{
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

    // Variables for CMA-ES, mean update
    int lambda_; int mu_;
    double[] mean_ = new double[10];
    double[][] covarianceMatrix_ = new double[10][10];
    double sigma_;
    double[] weights_;
    double muEff_; double muEffMinus_; double alphaMu_;
    double alphaMuEff_; double alphaPosDef_;
    double cM_;

    double mindna_; double maxdna_; // Specify search interval

    // Variables for CMA-ES, step-size control
    double cSigma_; double dSigma_;

    // Variables for CMA-ES, covariance matrix adaption
    double cCovariance_; double cOne_; double cMu_;

    // Psigma and pc
    double[] pSigma_ = new double[10];
    double[] pC_ = new double[10];

    double generation_;

    double sigma_bound_;

    double bound_used_ = 0;

    // Initialize Population
    public PopulationCmaes(int populationSize, Random rnd, GenUtils util, Map<String, Object> dict, double mindna, double maxdna){
        popSize_ = populationSize;
        lambda_ = populationSize;
        mu_ = populationSize/2;
    	rnd_ = rnd;
        util_ = util;
        mutationType_ = (String) dict.get("mutationType"); // Mutationfunction
        dict_ = dict; // Input for mutation function
        best_fitness_ = 0;
        generation_ = 0;
        mindna_ = (double) -mindna;
        maxdna_ = (double) maxdna;
        sigma_bound_ = (double) dict.get("sigma_bound"); // Mutationfunction
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
        setParamsCMAES();
        setMean();
        setSigma();
        setCovariance();
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


    public List<Individual> rank(List<Individual> pop){
        // Ranks population according to fitness

        pop.sort((o1, o2) -> Double.compare(o1.getFitness(),o2.getFitness()));
        Collections.reverse(pop);
        return pop;
    }


    public double sumArray(double[] array, int start, int end){
        double sum = 0;
        for (int i = start; i < end; i++){
            sum += array[i];
        }
        return sum;
    }

    public double sumArrayPos(double[] array, int start, int end){
        double sum = 0;
        for (int i = start; i < end; i++){
            if(array[i] > 0){
                sum += array[i];
            }
        }
        return sum;
    }

    public double sumArrayMin(double[] array, int start, int end){
        double sum = 0;
        for (int i = start; i < end; i++){
            if(array[i] < 0){
                sum += array[i];
            }
        }
        return sum;
    }

    public double sumArraySquared(double[] array, int start, int end){
        double sumSquared = 0;
        for (int i = start; i < end; i++){
            sumSquared += array[i]*array[i];
        }
        return sumSquared;
    }

    public void setParamsCMAES(){
        // Initializes default strategy params for CMA-ES according to table 1

        weights_ = new double[lambda_];

        // Set Weights Prime (eq 49, table 1)
        double[] weightsPrime = new double[lambda_];
        for (int i = 0; i < lambda_; i++){
            weightsPrime[i] = Math.log((lambda_+1)/2) - Math.log(i+1);
        }

        // Initialize muEff and muEffMinus according to table 1
        muEff_ = Math.pow(sumArray(weightsPrime, 0, mu_),2) / sumArraySquared(weightsPrime, 0, mu_);
        muEffMinus_ = Math.pow(sumArray(weightsPrime, mu_, lambda_),2) / sumArraySquared(weightsPrime, mu_, lambda_);

        // Initialize Covariance matrix adaptation parameters (eq 56-58, table 1)
        cCovariance_ = (4 + muEff_ / 10) / (10 + 4 + 2 * muEff_ / 10);
        cOne_ = 2 / (Math.pow((10 + 1.3) , 2) + muEff_);
        cMu_ = Math.min(1 - cOne_, 2 * (muEff_ - 2 + 1 / muEff_) / (Math.pow((10 + 2),2) + 2 * muEff_/ 2));

        // Initialize step-size control parameters (eq. 55, table 1)
        cSigma_ = (muEff_ + 2) / (10 + muEff_ + 5);
        dSigma_ = 1 + 2 * Math.max(0, Math.sqrt((muEff_ - 1) / (10 + 1)) - 1) + cSigma_;

        // Initialize selection and recombination params alphaMu, alphaMuEff (eq 50 - 52 ,table 1)
        alphaMu_ = 1 + cOne_ / cMu_;
        alphaMuEff_ = 1 + (2 * muEff_) / (muEff_  + 2);
        alphaPosDef_ = (1 - cOne_ - cMu_) / (10 * cMu_);

        // Initialize weights (eq 53, table 1)
        for (int i = 0; i < lambda_; i++){
            if (weightsPrime[i] >= 0){
                weights_[i] = (1 / (sumArrayPos(weightsPrime, 0, lambda_))) * weightsPrime[i];
            }
            else{
                weights_[i] = (Math.min(Math.min(alphaMu_, alphaMuEff_), alphaPosDef_) / Math.abs(sumArrayMin(weightsPrime, 0, lambda_))) * weightsPrime[i];
            }
        }
        cM_ = 1;

    }

    public void setMean(){
        // Initialize mean randomly, uniform (-5, 5)
        for(int i = 0; i < mean_.length; i++){
            mean_[i] = util_.randDouble(-4, 4);
        }
    }

    public void setSigma(){
        // Set sigma according to subscript page 29
        sigma_ = 0.5;
    }

    public void setCovariance(){
        // initial matrix is identity matrix
        for (int i = 0; i < 10; i++){
            for (int j = 0; j < 10; j++){
                if (i==j){
                    covarianceMatrix_[i][j] = 1;
                }
                else{
                    covarianceMatrix_[i][j] = 0;
                }
            }
        }
    }


    public boolean runCMAES(){
        // Set initial mean, sigma and covariance

        double yk[] = new double[10]; double xk[] = new double[10];

        double m[] = {0,0,0,0,0,0,0,0,0,0};

        // Loop through individuals and update dna according to eq 40 in figure 6
        for (int i = 0; i < lambda_; i++){
            yk = MatrixFunctions.multivariateGaussian(m, covarianceMatrix_, util_);
            xk = MatrixFunctions.add(mean_, MatrixFunctions.multiply(yk, sigma_));

            for(int j=0; j<10; j++){
                if(xk[j] > 5.0){
                    xk[j] = 5.0;
                }else if(xk[j] < -5.0){
                    xk[j] = -5.0;
                }
            }
            // Set dna, new fitness and yk of individual
            population_.get(i).setDna(xk);
            population_.get(i).setY(yk);
            population_.get(i).setFitness();

            if(population_.get(i).getFitness() > best_fitness_){
                best_fitness_ = population_.get(i).getFitness();
            }

        }

        population_ = rank(population_);

        double yw[] = new double[10];

        for(int i=0; i<mu_; i++){
            yw = MatrixFunctions.add(yw, MatrixFunctions.multiply(population_.get(i).getY(),weights_[i]));
        }

        mean_ = MatrixFunctions.add(mean_, MatrixFunctions.multiply(yw,sigma_));

        //MatrixFunctions.printMatrix(covarianceMatrix_);

        double[][] cRoot = Matrix.matrix_inversion(covarianceMatrix_);

        //MatrixFunctions.printMatrix(cRoot);

        double[] pSigma_temp1 = MatrixFunctions.multiply(pSigma_, (1-cSigma_));
        double[] pSigma_temp2 = MatrixFunctions.multiply(MatrixFunctions.multiply(cRoot, yw), (Math.sqrt(cSigma_*(2-cSigma_)*muEff_)));


        pSigma_ = MatrixFunctions.add(pSigma_temp1,pSigma_temp2);

        double pSigma_norm = MatrixFunctions.normVec(pSigma_);

        double E = 3.0847265652;

        sigma_ = sigma_ * Math.exp( (cSigma_/dSigma_) * ((pSigma_norm/E) -1));

        sigma_ = Math.max(0.00001, Math.min(sigma_, sigma_bound_));

        bound_used_ = sigma_;

        double hSigma = 0;

        if(pSigma_norm / Math.sqrt(1-Math.pow(1-cSigma_, 2*(generation_+1))) < (1.4+2/(10+1)) * E){
            hSigma = 1;
        }

        double[] pC_temp1 = MatrixFunctions.multiply(pC_, (1-cCovariance_));
        double[] pC_temp2 = MatrixFunctions.multiply(yw, hSigma * Math.sqrt(cCovariance_*(2-cCovariance_) * muEff_) );


        pC_ = MatrixFunctions.add(pC_temp1, pC_temp2);

        double[] weightsDot = new double[lambda_];

        for(int i=0; i<lambda_; i++){
            if(weights_[i] >= 0){
                weightsDot[i] = weights_[i];
            }else{
                double norm = MatrixFunctions.normVec(MatrixFunctions.multiply(cRoot,population_.get(i).getY()));
                weightsDot[i] = weights_[i] * (10 / Math.pow(norm,2));
            }
        }

        double deltaHSigma = (1-hSigma)*cCovariance_*(2-cCovariance_);

        double[][] cTemp1 = MatrixFunctions.multiply(covarianceMatrix_, 1+cOne_*deltaHSigma - cOne_ - cMu_);

        double[][] cTemp2 = MatrixFunctions.multiply(MatrixFunctions.multiply(pC_,pC_) ,cOne_);

        double[][] cTemp3 = new double[10][10];

        for(int i=0; i<lambda_; i++){
            double[] yI = population_.get(i).getY();
            cTemp3 = MatrixFunctions.add(cTemp3, MatrixFunctions.multiply(MatrixFunctions.multiply(yI,yI), weightsDot[i]));
        }

        covarianceMatrix_ = MatrixFunctions.add(cTemp1, MatrixFunctions.add(cTemp2, MatrixFunctions.multiply(cTemp3, cMu_)));

        if(!MatrixFunctions.isSymmetric(covarianceMatrix_) || !MatrixFunctions.isPositiveDefinite(covarianceMatrix_)){
            return false;
        }

        generation_ = generation_ + 1;

        return true;

    }

}
