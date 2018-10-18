import java.util.Random;
import java.util.List;
import java.util.Collections;
import java.util.Comparator;
import org.vu.contest.ContestEvaluation;
import java.lang.Math;
import java.util.ArrayList;
import java.util.Arrays;

public class GenUtils{
	Random rnd_;
	ContestEvaluation eval_;
	int evals_;
	int evaluations_limit_;

	public GenUtils(Random rnd, ContestEvaluation eval, int evaluations_limit){
		rnd_ = rnd;
		eval_ = eval;
		evaluations_limit_ = evaluations_limit;
		evals_ = 0;
	}

	public double randDouble(double min, double max) {
		// Returns randomly selected double within (min,max)

	    double result = rnd_.nextDouble() * (max - min) + min;

	    return result;
	}

	public double sampleGaussian(double mean, double variance) {
		// Samples item from a normal distribution with input mean and variance
		double result = rnd_.nextGaussian() * Math.sqrt(variance) + mean;

		return result;
	}

	public int randInt(int min, int max){
		int result = rnd_.nextInt(max - min + 1) + min;

        return result;
	}

	public double determineFitness(double[] dna){
		evals_ += 1;
		return (double) eval_.evaluate(dna);
	}

	public int getEvals(){
		return evals_;
	}

	public int getEvaluationLimit(){
		return evaluations_limit_;
	}

}
