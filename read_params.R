readInputParams <- function(){
args <- commandArgs(trailingOnly = TRUE)
argsIter = 1 
# metrics
metrics = args[argsIter]
populationNumber = strtoi(args[argsIter<-argsIter+1])
#print(metrics)
#print(populationNumber)

if (is.na(populationNumber) == TRUE){
  print("Wrong pop number provided")
}

agnesType = strtoi(args[argsIter<-argsIter+1])


# metody grupowania dla GA
groupingNumberGA = strtoi(args[argsIter<-argsIter+1])
gaGroup = c()
if (groupingNumberGA == 3){
  # all 
}else {
  # metody muszą być wyspecyfikowane
  for (i in 1:groupingNumberGA)
    gaGroup <- union(gaGroup, c(args[argsIter<-argsIter+1]))
}

# parametry dla oceny GA: groupingNumberGA + gaGroup (wector z nazwami grup)
#print("GA")
#print(groupingNumberGA)
#print(gaGroup)


# metody grupowania dla RBGA
# groupingNumberRBGA = strtoi(args[argsIter<-argsIter+1])
# rbgaGroup = c()
# if (groupingNumberRBGA == 3){
#   # all 
# }else {
#   for (i in 1:groupingNumberRBGA)
#     rbgaGroup <- union(rbgaGroup, c(args[argsIter<-argsIter+1]))
# }

# parametry dla oceny RBGA: groupingNumberRBGA + rbgaGroup (wector z nazwami grup)
#print("RBGA")
#print(groupingNumberRBGA)
#print(rbgaGroup)
# metody grupowania dla DE
# groupingNumberDE= strtoi(args[argsIter<-argsIter+1])
# deGroup = c()
# if (groupingNumberDE == 3){
#   # all 
# }else {
#   for (i in 1:groupingNumberDE)
#     deGroup <- union(deGroup, c(args[argsIter<-argsIter+1]))
# }

# parametry dla oceny DE: groupingNumberDE + deGroup (wector z nazwami grup)
#print("DE")
#print(groupingNumberDE)
#print(deGroup)
# opcjonalne parametry GA na końcu
gaArgs = c()
gaProvided = strtoi(args[argsIter<-argsIter+1])
for (i in 1:gaProvided){
  gaArgs <- union(gaArgs, c(args[argsIter<-argsIter+1]))
  
}

#print("OPTIONAL GA ARGS: ")
#print(gaProvided)
#print(gaArgs)

result <- list("metrics"= metrics, "popNr"=populationNumber, 
               "GaGroupsNr" =  groupingNumberGA, "gaGroups" = gaGroup, "GAparamsNr" = gaProvided, "GAparams" = gaArgs)
# print(result)
return(result)
}

