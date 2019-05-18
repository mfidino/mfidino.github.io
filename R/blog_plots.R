

x <- seq(-10,10, 0.01)
y <- dnorm(x, 0, 2)

jpeg("./blog/images/normal_plot1.jpg", quality = 100)
par(mar = c(5,5,0.5,0.5))
plot(y~x, type ='l', bty = 'l', las = 1, lwd = 4, xlab = 'Experimental outcome',
     ylab = 'Probability density of outcome', col = 'blue', cex.lab = 1.5)
polygon(c(x, rev(x)), c(y,rep(0, length(y))), col = scales::alpha('blue', 0.5),
        border = NA)
dev.off()


jpeg("./blog/images/normal_plot2.jpg", quality = 100)
par(mar = c(5,5,0.5,0.5))
x <- seq(-10,10, 0.01)
y1 <- dnorm(x, 0, 2)
y2 <- dnorm(x, -2, 1)
plot(y1~x, type ='l', bty = 'l', las = 1, lwd = 4, xlab = 'Experimental outcome',
     ylab = 'Probability density of outcome', col = 'blue', cex.lab = 1.5, ylim = c(0, 0.5))
polygon(c(x, rev(x)), c(y1,rep(0, length(y))), col = scales::alpha('blue', 0.5),
        border = NA)
polygon(c(x, rev(x)), c(y2,rep(0, length(y))), col = scales::alpha('gray30', 0.5),
        border = NA)
lines(y2 ~ x, col ='gray30', lwd = 4)
text(x = 2.75, y = 0.1, pos = 4, labels = 'Normal(0, 2)', col = 'blue', cex = 1.5)
text(x = -0.75, y = 0.3, pos = 4, labels = 'Normal(-2, 1)', col = 'gray30', cex=1.5)
dev.off()

jpeg("./blog/images/normal_plot3.jpg", quality = 100)
par(mar = c(5,5,0.5,0.5))
# generate 2000 random normal data points with mean = 500 and sd = 100.
our_normal_data <- rnorm(n = 2000, mean = 500, sd = 100)

# plot out the density of the data as a histogram
hist(our_normal_data,
     freq = FALSE,
     breaks = 50,
     col = 'lightblue',
     main = "",
     xlab = "Our normally distributed data",
     cex.lab = 1.5)

# Calculate the probability density of a range of points from
#  0 to 1000 from Normal(mean = 500, sd = 100)
range_of_values <- seq(0,1000,1)
normal_density <- dnorm(range_of_values, mean = 500, sd = 100)

# add this to the histogram
lines(normal_density ~ range_of_values, lwd = 2)
dev.off()
