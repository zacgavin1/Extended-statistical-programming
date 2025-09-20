#setwd("Extended-statistical-programming") ## comment out of submitted
a <- scan("Shakespeare_complete_works.txt",what="character",skip=83,nlines=196043-83,
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

## (b) Removing character names and arabic numerals #####

# note: some roman numeral I's are likely left in
a <- a[-which((toupper(a)==a) & (a!= "I") & (a != "A") &(a != "O")) ] 


## (c) removing _ and - from words
a <- gsub("-", "", a ) # it may be better to split the words rather than combining them
a <- gsub("_", "", a )

## (d) isolating punctuation

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
a <- split_punc(a, c(",",":"))   


# (f) making a lower case 

a <- tolower(a)


########### END OF THE PREPROCESSING ############

######## Question 5 ############

words <- unique(a) # b is the set of all words used
e <- match(a,words) # where each element of a appears in b

occurences <- tabulate(e)
names(occurences) <- words  # list of occurences, with entry names the words

occurences <- sort(occurences, decreasing=TRUE)
b <- occurences[1:1000]  # now have the top 1000 words


######## Question 6 ###########

mlag = 4


  



