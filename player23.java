import org.vu.contest.ContestSubmission;
import org.vu.contest.ContestEvaluation;

import java.util.Random;
import java.util.Properties;
import java.util.Comparator;
import java.lang.Math;
import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.stream.IntStream;
import java.util.Map;
import java.util.HashMap;

public class player23 implements ContestSubmission
{
	Random rnd_;
	ContestEvaluation evaluation_;
  private int evaluations_limit_;
  double variance;
  String mutationType_;
  String mutationFunction_;
	String mutationDistribution_;

	int parSize;
	int kidSize;
	int popSize;
	int k;
	double dLow;
	double dHigh;
	boolean DGEA = false;
	boolean CMAES = false;
	boolean continueCmeas;
	double sigma_bound;

	public player23()
	{
		rnd_ = new Random();
	}

	public void setSeed(long seed)
	{
		// Set seed of algortihms random process
		rnd_.setSeed(seed);
	}

	public void setEvaluation(ContestEvaluation evaluation)
	{
		// Set evaluation problem used in the run
		evaluation_ = evaluation;

		// Get evaluation properties
		Properties props = evaluation.getProperties();
    // Get evaluation limit
    evaluations_limit_ = Integer.parseInt(props.getProperty("Evaluations"));
		// Property keys depend on specific evaluation
		// E.g. double param = Double.parseDouble(props.getProperty("property_name"));
    boolean isMultimodal = Boolean.parseBoolean(props.getProperty("Multimodal"));
    boolean hasStructure = Boolean.parseBoolean(props.getProperty("Regular"));
    boolean isSeparable = Boolean.parseBoolean(props.getProperty("Separable"));

		if (System.getProperty("variance") == null){
			variance = 0.1;
		}
		else{
			variance = Double.parseDouble(System.getProperty("variance"));
		}
		if (System.getProperty("dLow") == null){
			dLow = 0.1;
		}
		else{
			dLow = Double.parseDouble(System.getProperty("dLow"));
		}
		if (System.getProperty("dHigh") == null){
			dHigh = 0.1;
		}
		else{
			dHigh = Double.parseDouble(System.getProperty("dHigh"));
		}


		mutationType_ = System.getProperty("mutationType", "mutateStandard");
    mutationFunction_ = System.getProperty("mutationFunction", "linearDecreasing");
		mutationDistribution_ = System.getProperty("mutationDistribution", "gaussian");
		parSize = Integer.parseInt(System.getProperty("parSize", "10"));
		kidSize = Integer.parseInt(System.getProperty("lambda", "70"));
		popSize = Integer.parseInt(System.getProperty("mu", "30"));
		k = Integer.parseInt(System.getProperty("k", "4"));

		if(mutationType_.equals("cmaes")){
			CMAES = true;
		}
		if(mutationType_.equals("mutateDG")){
			DGEA = true;
		}
		// Do sth with property values, e.g. specify relevant settings of your algorithm
        if(isMultimodal){
            if(hasStructure){
                // Schaffers
									sigma_bound = 1.5;
            }else{
							// Katsuura
							sigma_bound = 0.5;
            }
        }else{
            // Bentcigar
						sigma_bound = 2;
        }
    }

    // hier begint het
	public void run(){
			boolean print_diversity = false;

   		GenUtils util = new GenUtils(rnd_, evaluation_, evaluations_limit_);

   		// Create dictionary with all items necessary for mutationfunction
   		Map<String, Object> dict = new HashMap<String, Object>();
   		dict.put("mutationType", mutationType_);
   		dict.put("mutationFunction", mutationFunction_);
		  dict.put("mutationDistribution", mutationDistribution_);
   		dict.put("variance", variance);
		  dict.put("k", k);
			dict.put("sigma_bound", sigma_bound);

			// Initialize new population randomly
   		Population pop = new Population(popSize, rnd_, util, dict);
   		pop.init();

   		PopulationCmaes popCmaes = new PopulationCmaes(popSize, rnd_, util, dict, 5, 5);
   		popCmaes.init();

   		int evalsNext = parSize + kidSize;


   		if(CMAES){
   			evalsNext = popSize;
   		}

   		// Evolve while evaluations is smaller than limit
			while (util.getEvals() + evalsNext < evaluations_limit_){

				if(DGEA){
					pop.mutateDG(dLow, dHigh, parSize, kidSize);
				}
				else if (CMAES){
					continueCmeas = popCmaes.runCMAES();
					if(!continueCmeas){
						popCmaes = new PopulationCmaes(popSize, rnd_, util, dict, 5, 5);
	   					popCmaes.init();
					}
				}
				else{


					// Select Parents
			 	pop.selectTournament(parSize);
				//pop.selectTournament(parSize);

			 	// Do crossover
			 	pop.crossover(kidSize, parSize);

			 	// Do mutation for all parents and generated children
			 	pop.mutate();

			 	// Select survivors with best fitness
			 	pop.selectSurvivors();

			}

			 	if(print_diversity){
			 		if(CMAES){
			 			// Update diversity
			 			popCmaes.setDiversity();
			 			System.out.print(util.getEvals());
			 			System.out.print(",");
			 			System.out.print(popCmaes.getDiversity());
			 			System.out.print(",");
			 			System.out.println(popCmaes.getBestFitness());
			 		}else{
			 			// Update diversity
			 			pop.setDiversity();
			 			System.out.print(util.getEvals());
			 			System.out.print(",");
			 			System.out.print(pop.getDiversity());
			 			System.out.print(",");
			 			System.out.println(pop.getBestFitness());
			 		}

			 	}

		}
	}
}
