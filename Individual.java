import org.vu.contest.ContestEvaluation;
import java.util.Random;
import java.util.Arrays;

public class Individual{
    double fitness_ = -1;
    double[] dna_ = new double[10];
    double[] sigma_ = new double[10];
    double[] y_ = new double[10]; //for CMA-ES
    double[] alpha_ = new double[45];
    GenUtils util_;

    // Constructor class to initialise an individual.
    public Individual(GenUtils util){
        util_ = util;
        init();
    }

    // Initialises an individual.
    public void init(){
    	for(int i = 0; i < 10; i++){
    		dna_[i] = util_.randDouble(-5.0,5.0);
            sigma_[i] = util_.randDouble(0.01, 0.5);
    	}
        for(int i = 0; i < 45; i++){
            alpha_[i] = 0;
        }
    }

    // set the dna of the child
    public void setDna(double[] dnaKid){
      dna_ = dnaKid;
    }

    public void setSigma(double[] sigmaKid){
        sigma_ = sigmaKid;
    }

    public void setAlpha(double[] alpha){
        alpha_ = alpha;
    }

    public void setY(double[] y){
        y_ = y;
    }

    public double[] getDna(){
        return dna_;
    }

    public double[] getSigma(){
        return sigma_;
    }

    public double[] getAlpha(){
        return alpha_;
    }

    public double[] getY(){
        return y_;
    }

    // Set fitness value.
    public void setFitness(){
        fitness_ = util_.determineFitness(dna_);
    }

    // Get fitness value.
    public double getFitness(){
        return fitness_;
    }

}
