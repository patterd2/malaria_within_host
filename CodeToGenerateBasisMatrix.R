# to generate a basis matrix for use in Matlab simulations

library(splines)
h = 0.125;
# spline will be defined over 1000 days at 7.5 minute increments
timeVals = seq(0,1000, by = h/24)

# basis matrix (cubic spline), no knots, nonzero intercept allowed
basisMatrixNoKnots = bs(timeVals, intercept = T)

# to visualize:
matplot(x = timeVals, y = basisMatrixNoKnots, type = 'l')

# basis matrix (cubic spline), one knot at day 75, nonzero intercept allowed
basisMatrixKnot = bs(timeVals, knots = c(75), intercept = T)

# to visualize:
matplot(x = timeVals, y = basisMatrixKnot, type = 'l')
abline(v=75, col = 'purple', lwd = 2)

write.table(basisMatrixNoKnots, "basisMatrixNoKnots_1000_0.125.txt")
write.table(basisMatrixKnot, "basisMatrixKnot_1000_0.125.txt")

# not to generate output, but if you wanted a knot halfway through the time series:
basisMatrixKnotHalfway = bs(timeVals, df = 5, intercept = T)

# plot for comparison with previous two:
matplot(x = timeVals, y = basisMatrixKnotHalfway, type = 'l')

# code to generate conversion rate (transmission investment) strategy:
cparms = runif(length(basisMatrixNoKnots[1,]), -10, 10)
cVal = exp(-exp(basisMatrixNoKnots %*% cparms))

plot(timeVals, cVal, type = 'l')
