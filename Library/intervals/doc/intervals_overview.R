### R code from vignette source 'intervals_overview.Rnw'

###################################################
### code chunk number 1: width
###################################################
options( width = 80 )


###################################################
### code chunk number 2: Intervals
###################################################
library( intervals )
x <- Intervals( matrix( 1:6, ncol = 2 ) )
x
x[2,2] <- NA
x[3,1] <- 6
x


###################################################
### code chunk number 3: Intervals_full
###################################################
y <- as( x, "Intervals_full" )
y <- c( y, Intervals_full( c(2,3,5,7) ) )
closed(y)[2:3,1] <- FALSE
closed(y)[4,2] <- FALSE
rownames(y) <- letters[1:5]
y


###################################################
### code chunk number 4: size
###################################################
size(x)
empty(x)
empty(y)


###################################################
### code chunk number 5: plotting
###################################################
plot( y )


###################################################
### code chunk number 6: set_operations
###################################################
reduce( y )
interval_intersection( x, x + 2 )
interval_complement( x )


###################################################
### code chunk number 7: set_operations
###################################################
interval_union( x, interval_complement( x ) )


###################################################
### code chunk number 8: distance
###################################################
B <- 100000
left <- runif( B, 0, 1e8 )
right <- left + rexp( B, rate = 1/10 )
v <- Intervals( cbind( left, right ) )
head( v )
mean( size( v ) )
dim( reduce( v ) )
system.time( d <- distance_to_nearest( sample( 1e8, B ), v ) )


###################################################
### code chunk number 9: distanceplot
###################################################
hist( d, main = "Distance to nearest interval" )


###################################################
### code chunk number 10: overlap
###################################################
rownames(v) <- sprintf( "%06i", 1:nrow(v) )
io <- interval_overlap( v, v )
head( io, n = 3 )
n <- sapply( io, length )
sum( n > 1 )
k <- which.max( n )
io[ k ]
v[ k, ]
v[ io[[ k ]], ]


###################################################
### code chunk number 11: clusters
###################################################
B <- 100
left <- runif( B, 0, 1e4 )
right <- left + rexp( B, rate = 1/10 )
y <- reduce( Intervals( cbind( left, right ) ) )
w <- 100
c2 <- clusters( y, w )
c2[1:3]


###################################################
### code chunk number 12: expand
###################################################
delta <- .Machine[[ "double.eps" ]]^0.5
y1 <- Intervals( c( .5, 1 - delta / 2 ) )
y2 <- Intervals( c( .25, 1, .75, 2 ) )
y1
y2
all.equal( y1[1,2], y2[2,1] )
interval_intersection( y1, y2 )


###################################################
### code chunk number 13: expand
###################################################
contract( y1, delta, "relative" )


###################################################
### code chunk number 14: expand
###################################################
inner <- interval_intersection(
                               contract( y1, delta, "relative" ),
                               contract( y2, delta, "relative" )
                               )
inner
outer <- interval_intersection(
                               expand( y1, delta, "relative" ),
                               expand( y2, delta, "relative" )
                               )
outer


###################################################
### code chunk number 15: expand
###################################################
interval_difference( outer, inner )


###################################################
### code chunk number 16: gaps
###################################################
x <- Intervals( c(1,10,100,8,50,200), type = "Z" )
x
w <- 2
close_intervals( contract( reduce( expand(x, w/2) ), w/2 ) )


###################################################
### code chunk number 17: integer_range
###################################################
.Machine$integer.max
numeric_max <- with( .Machine, double.base^double.digits )
options( digits = ceiling( log10( numeric_max ) ) )
numeric_max


###################################################
### code chunk number 18: sessionInfo
###################################################
si <- as.character( toLatex( sessionInfo() ) )
cat( si[ -grep( "Locale", si ) ], sep = "\n" )


