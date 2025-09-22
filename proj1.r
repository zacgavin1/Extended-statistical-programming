#setwd("Extended-statistical-programming") ## comment out of submitted
a <- scan("shakespeare.txt",what="character",skip=83,nlines=196043-83,
          fileEncoding="UTF-8")


####### START OF PREPROCESSING ###############

## (a) identify and remove stage directions 

# First we manually fix all bracketing errors in the text
a[463343] <- "So"    # The original text was erroneously "[So"
a[322751] <- "Rest.]" # missing end of direction


direction_starts = grep("[", a, fixed = TRUE) # find where stage directions start

direction_ends <- c()
for (direction in direction_starts) {
   close_brackets <- grep("]", a[direction:(direction+100)]); # this gets ALL close brackets
   direction_ends <- append(direction_ends, close_brackets[1]) # and then we choose the first one after each [
}



direction_indexes <- c() # to hold positions of stage directions
# the following loop may be able to be done more nicely by some vector operation
for (i in 1:length(direction_starts)){ 
  direction_indexes <- append(direction_indexes, direction_starts[i]:(direction_starts[i]+direction_ends[i]-1)) ;
}

a <- a[-direction_indexes] # removes all stage directions 


## (d) isolating punctuation
## this section probably wants to go before the removing capital words, numerals etc, to avoid deleting "I, "

# this task is similar to the one from notes
split_punc <- function(v, marks){
  regex = paste(marks, collapse="|") # This contains . which will give problems 
  i_punc <- which(grepl(regex, v)) # using regex to see where elements of marks are
  vs <- rep("", length(v)+length(i_punc)) # initialise new vector of req length
  ips <- i_punc+ 1:length(i_punc) # i assume here all punctuation is at end of words
  vs[ips]<- substr(v[i_punc], nchar(v[i_punc]), nchar(v[i_punc])) 
  vs[ips-1] <- substr(v[i_punc], 1, nchar(v[i_punc])-1) # words before puncs
  vs[-(c(ips, ips-1))] <- v[-i_punc] # all other words
  vs
}


# (e) actually split up the puncs and the words
a <- split_punc(a, c(",","\\.",";","!",":","\\?"))    #need to use \\ before symbols that have regex meaning like .



## (b) Removing character names and arabic numerals #####

# note: some roman numeral I's are likely left in
a <- a[-which((toupper(a)==a) & (a!= "I") & (a != "A") &(a != "O")) ] 


## (c) removing _ and - from words
a <- gsub("-", "", a ) # it may be better to split the words rather than combining them
a <- gsub("_", "", a )




# (f) making a lower case 

a <- tolower(a)


########### END OF THE PREPROCESSING ############

######## Question 5 ############

words <- unique(a) # b is the set of all words used
e <- match(a,words) # where each element of a appears in b

occurences <- tabulate(e)
names(occurences) <- words  # list of occurences, with entry names the words

occurences <- sort(occurences, decreasing=TRUE)
b <- names(occurences[1:1000])  # now have the top 1000 words


######## Question 6 ###########
n    <- length(a) # note this is the same length as tknised
mlag <- 4
# (a)
tknised <- match(a,b) # now have a tokenised version of a. In instructions he calls this M1


# (b)
M <- tknised[1:(n-mlag)] 
for (col in 2:(mlag+1)){ # loop appends shifted copies of tknised to the right of M
  M <- cbind(M, tknised[col:(n-mlag-1 + col)])
}


####### Question 7 ###########

next.word <- function(key, M, M1, w=rep(1,ncol(M)-1)){
  # need to find where the row matches are for each length string 
  # need to make probabilites vector corresponding to els of b
  # need to add probs to it based on how many row matches there are (weighted by w_i)
  # repeat this for each length of word
  # then sample from this distbn
  
  nw_measure <- rep(0, 1000) # i would do length(b) but he doesn't seem to want us to use b in this function
   
  for (mc in 1:mlag){
    
    # the double loop is ugly but this works to identify where the matches are
    # its not based off his suggestion from the notes but i couldn't work that
    for (row in 1:nrow(M)){
      matches[row] <- isTRUE(all.equal(M[row,mc:mlag], c(1,2,28,30), check.attributes=FALSE))
    }
  
    #now we add to the distribution we've created 
    
  
  }
}









