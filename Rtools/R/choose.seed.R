#Randy Johnson
#Randomly choose a random seed that can be saved and reused...
#Created April 23, 2004

choose.seed = function(len=9)
{
    x = sample(c(0:9), len, replace=TRUE)
    return(sum(x*(rep(10,len)^c(0:(len-1)))))
}
