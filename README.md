# Evolutionary-Computing

In order to compile the code, run:
```
javac -cp contest.jar player23.java Individual.java Population.java GenUtils.java Matrix.java PopulationCmaes.java MatrixFunctions.java EigenValueDecomposition.java
jar cmf MainClass.txt submission.jar player23.class Individual.class Population.class GenUtils.class Matrix.class MatrixFunctions.class EigenValueDecomposition.class PopulationCmaes.class
```


## LDM
### BentCigarFunction
```
java -DmutationType=mutateStandard -DmutationFunction=linearDecreasing -DmutationDistribution=gaussian -Dmu=80 -Dlambda=169 -Dvariance=0.003817733496358963 -Dk=6 -Dparsize=64 -jar testrun.jar -submission=player23 -evaluation=BentCigarFunction -seed=1
```
### SchaffersEvaluation
```
java -DmutationType=mutateStandard -DmutationFunction=linearDecreasing -DmutationDistribution=gaussian -Dmu=67 -Dlambda=186 -Dvariance=0.02931366086141189 -Dk=8 -Dparsize67 -jar testrun.jar -submission=player23 -evaluation=SchaffersEvaluation -seed=1
```
### KatsuuraEvaluation
```
java -DmutationType=mutateStandard -DmutationFunction=linearDecreasing -DmutationDistribution=gaussian -Dmu=17 -Dlambda=152 -Dvariance=0.0010615728075645916 -Dk=4 -Dparsize=4 -jar testrun.jar -submission=player23 -evaluation=KatsuuraEvaluation -seed=1
```

## DGEA
### BentCigarFunction
```
java -DmutationType=mutateDG -DmutationFunction=static -DmutationDistribution=gaussian  -Dmu=80 -Dlambda=114 -DdLow=0.0020283390188316447 -DdHigh=0.3126058325379971 -Dk=9  -Dparsize=23 -jar testrun.jar -submission=player23 -evaluation=BentCigarFunction -seed=1
```
### SchaffersEvaluation
```
java -DmutationType=mutateDG -DmutationFunction=static -DmutationDistribution=gaussian -Dmu=128 -Dlambda=128 -DdLow=0.0022876039949272374 -DdHigh=0.23176941855081462 -Dk=8 -Dparsize=35 -jar testrun.jar -submission=player23 -evaluation=SchaffersEvaluation -seed=1
```
### KatsuuraEvaluation
```
java -DmutationType=mutateDG -DmutationFunction=static -DmutationDistribution=gaussian -Dmu=93 -Dlambda=217 -DdLow=0.016926841434103237 -DdHigh=0.23935788204849787 -Dk=9 -Dparsize=46 -jar testrun.jar -submission=player23 -evaluation=KatsuuraEvaluation -seed=1
```

## CMA-ES
### BentCigarFunction
```
java -DmutationType=cmaes -Dmu=40  -jar testrun.jar -submission=player23 -evaluation=BentCigarFunction -seed=1
```
### SchaffersEvaluation
```
java -DmutationType=cmaes -Dmu=40  -jar testrun.jar -submission=player23 -evaluation=SchaffersEvaluation -seed=1
```
### KatsuuraEvaluation
```
java -DmutationType=cmaes -Dmu=40  -jar testrun.jar -submission=player23 -evaluation=KatsuuraEvaluation -seed=1
```
