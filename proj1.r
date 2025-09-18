#setwd("Extended-statistical-programming") ## comment out of submitted
a <- scan("Shakespeare_complete_works.txt",what="character",skip=83,nlines=196043-83,
          fileEncoding="UTF-8")

######## identify and remove stage directions #########

## First we manually fix all bracketing errors in the text
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

###### Removing character names and arabic numerals #####



