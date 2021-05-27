library(dplyr)
x = read.table("stats", header=TRUE)

means = x %>% group_by(N, p, recrate) %>% summarise(mean_pi = mean(div), var_pi = var(div), cv_pi = sd(div)/mean(div), .groups = "drop")

print(means)
