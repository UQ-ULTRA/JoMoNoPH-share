1.read relevant simulation scenario from scenario.R
2.generate data by calling datagen.R
3.fit model by calling stan models
4.save results in /hpc/output_yyyymmddhhmm folder in a specific format (e.g. .rds) for later aggregation 
5.aggregate results by calling aggregate_hpc.R
6.produce tables
7.produce figures by calling plot_hpc.R