# Evolutionary-Computing

```
javac -cp contest.jar player23.java Individual.java Population.java GenUtils.java Matrix.java PopulationCmaes.java MatrixFunctions.java EigenValueDecomposition.java
jar cmf MainClass.txt submission.jar player23.class Individual.class Population.class GenUtils.class Matrix.class MatrixFunctions.class EigenValueDecomposition.class PopulationCmaes.class
```

To run an experiment, execute:
```
java -jar testrun.jar -submission=player23 -evaluation=KatsuuraEvaluation -seed=1
```


## LDM
```
java -DmutationType=mutateStandard -DmutationFunction=linearDecreasing -DmutationDistribution=gaussian -Dmu=80 -Dlambda=169 -Dvariance=0.003817733496358963 -Dk=6 -Dparsize=64 -jar testrun.jar -submission=player23 -evaluation=BentCigarFunction -seed=1
```
```
java -DmutationType=mutateStandard -DmutationFunction=linearDecreasing -DmutationDistribution=gaussian -Dmu=67 -Dlambda=186 -Dvariance=0.02931366086141189 -Dk=8 -Dparsize67 -jar testrun.jar -submission=player23 -evaluation=SchaffersEvaluation -seed=1
```
```
java -DmutationType=mutateStandard -DmutationFunction=linearDecreasing -DmutationDistribution=gaussian -Dmu=17 -Dlambda=152 -Dvariance=0.0010615728075645916 -Dk=4 -Dparsize=4 -jar testrun.jar -submission=player23 -evaluation=KatsuuraEvaluation -seed=1
```

## CMA-ES

```
java -DmutationType=cmaes -Dmu=40  -jar testrun.jar -submission=player23 -evaluation=BentCigarFunction -seed=1
```

```
java -DmutationType=cmaes -Dmu=40  -jar testrun.jar -submission=player23 -evaluation=SchaffersEvaluation -seed=1
```

```
java -DmutationType=cmaes -Dmu=40  -jar testrun.jar -submission=player23 -evaluation=KatsuuraEvaluation -seed=1
```
