#Vector

data <- c(1,2,3,4,5)
data1 <- c("a","b","c")

#List

my_list <- list(22,"ab",TRUE) #Create a list
my_list <- list(data,"ab",TRUE) #Create a list
data
my_list

#Matrix
my_matrix <- matrix(1:6,nrow=3,ncol=2)
my_matrix
my_matrix <- matrix(1:6,nrow=3,ncol=2,byrow=T)
my_matrix

a <- 1:5
b <- 6:10
c <- 11:15

cbind(a,b,c) #column wise bind
rbind(a,b,c) #row wise bind

#Data-frame

df <- data.frame(name =c("pravin","ajay","mamta","pallavi"), score=c(87,56,91,67))

df[3:4,]
df[c(1,3),]
#Subsetting the R-objects
data[3] #Print 3rd position from vector data
my_list[[2]] #Print 2nd position from list my_list
my_list[[1]][3]
my_list
my_matrix[2,1] #Print row 2 col 1 value
my_matrix[3,] #Print all the values of row 3
my_matrix[,2] #Print all the values of column 2
df[3,1] #Print info stored in 3rd row first column
df[2,] #Print second row
df[2:3,] #Print all info  of 2nd and 3rd row 

df
#Subset the dataframe
df1 <- subset(df, score > 60)
df1
data[3] > 2
df$score > 60

#Exercice

a <- c(10,20,30,40)
b <- c('book','pen','textbook','pencil_case')
c <- c(TRUE,FALSE,TRUE,FALSE)
d <- c(2.5,8,10,7)
df <- data.frame(a,b,c,d)
colnames(df) <- c("ID","ITEMS","STORE","PRICE")
class(df)
str(df) #Print structure of dataframe

#1. Select row 1 in column 2
df[1,2]
#2. Select rows 1 to 3
df[1:3,]
#3. Select rows 1 to 3 and columns 3 to 4
df[1:3,3:4]
#4. Select items if price is more than 8
df1 <- subset(df,PRICE > 8)
df$PRICE > 8
#5. Select row info where item is book
subset(df,ITEMS == "book")
df

install.packages("dplyr",dependencies = T)

#How to
