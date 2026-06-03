## Forecasting objects of class varest
args(vars:::predict.varest)
predictions <- predict(varsimest, n.ahead = 25,
                       ci = 0.95)
class(predictions)
args(vars:::plot.varprd)
## Plot of predictions for y1
plot(predictions, names = "y1")
## Fanchart for y2
args(fanchart)
fanchart(predictions, names = "y2")
