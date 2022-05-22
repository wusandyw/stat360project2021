source("mars.R")
set.seed(123)
n = 10
x = data.frame(x1=rnorm(n), x2=rnorm(n))
y = rnorm(n)
rp = fwd_stepwise(y,x,Mmax=2)


#########################################################

library(ISLR)
data(Wage)
# Inputs
y = Wage$wage
x = matrix(Wage$age, ncol=1)
Mmax=4

# Test fwd_stepwise function
fwd = fwd_stepwise(y,x,Mmax)
# fwd should be list of three objects

# Sample output
length(fwd$y)
# [1] 3000
head(fwd$B)
# B0 B1 B2 B3 B4
# [1,]  1  0 20  0 45
# [2,]  1  0 14  0 39
# [3,]  1  7  0  0 18
# [4,]  1  5  0  0 20
# [5,]  1 12  0  0 13
# [6,]  1 16  0  0  9
length(fwd$splits)
# [1] 5  



# Test bwd_stepwise function

bwd = bwd_stepwise(fwd)

# bwd should be list of three objects

# Sample output
length(bwd$y)
# [1] 3000
head(bwd$B)
# B0 B2 B3
# [1,]  1 20  0
# [2,]  1 14  0
# [3,]  1  0  0
# [4,]  1  0  0
# [5,]  1  0  0
# [6,]  1  0  0
length(bwd$splits)
# [1] 3