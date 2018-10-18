# Evolutionary-Computing

'''
javac -cp contest.jar player23.java Individual.java Population.java GenUtils.java Matrix.java PopulationCmaes.java MatrixFunctions.java EigenValueDecomposition.java
jar cmf MainClass.txt submission.jar player23.class Individual.class Population.class GenUtils.class Matrix.class MatrixFunctions.class EigenValueDecomposition.class PopulationCmaes.class
'''

To run an experiment, execute:
'''
java -jar testrun.jar -submission=player23 -evaluation=KatsuuraEvaluation -seed=1
'''
