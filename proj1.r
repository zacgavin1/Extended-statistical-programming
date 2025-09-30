#setwd("Extended-statistical-programming") ## comment out of submitted
#setwd("C:/Users/shaeh/Desktop/edinburgh-notes/Extended-Statistical-Programming/project-ESP/Extended-statistical-programming") ## comment out of submitted
#setwd("C:/Users/Brandon Causing/Downloads/Extended Statistical Programming/Extended-statistical-programming")

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
   close_brackets <- grep("]", a[direction:(direction+100)]); # this gets ALL close brackets following an open bracket
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
a <- split_punc(a, c(",","\\.",";","!",":","\\?"))    #need to use \\ before symbols that have regex meaning like .

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
n    <- length(a) # note this is the same length as M1
mlag <- 4
# (a)
M1 <- match(a,b) # now have a tokenised version of a. In instructions he calls this M1


# (b)
M <- M1[1:(n-mlag)] 
for (col in 2:(mlag+1)){ # loop appends shifted copies of M1 to the right of M
  M <- cbind(M, M1[col:(n-mlag-1 + col)])
}


####### Question 7 ###########

# pick a next-word based on key (sequence of tokens), M, M1, w (weights)
next.word <- function(key, M, M1, w=rep(1,ncol(M)-1)){
  
  ###### 1 - Set-Up ######
  mlag <- ncol(M) - 1 # use first mlag cols to come up with next word
  
  if (length(key) > mlag) { # only use the last mlag tokens
    key <- key[(length(key) - mlag + 1):length(key)]   
  }
  
  all_next_words <- c() # initialize list of next-words to sample from
  length_i <- c() # no. of next-words in iteration i
  
  ##### 2 - Search for Next Words ######
  # find key match of length mlag, length mlag - 1, ..., 1
  for (i in 1:length(key)) {
    
    ### 2a - Pick columns to match key ###
    current_key <- key[i:length(key)] # use last mlag - i + 1 tokens of key
    context_len <- length(current_key)  
    cols_to_match <- (mlag - context_len + 1):mlag # which cols in M to match
    
    ### 2b - Find matches ###
    # return F if key[j] matches M[k, j] for some row k of M, T otherwise
    ii <- colSums(!(t(M[, cols_to_match, drop=FALSE])==current_key))
    # if sum of components in a row = 0, entire key matches
    matching_rows <- which(ii == 0)
    
    ### 2c - Store next-words and no. of next-words ###
    if (length(matching_rows) > 0) { # if there is a matching row
      
      # get last column of M (possible next-words) for matching rows
      next_words_found <- M[matching_rows, mlag + 1]
      # collect found next-words, omitting rare words (NA's)
      all_next_words <- c(all_next_words, na.omit(next_words_found))
      # store no. of (non-rare) next-words found in iteration i
      length_i <- c(length_i, length(na.omit(next_words_found)))
      
    } else { # if there is not a matching row
      length_i <- c(length_i, 0) # length_i = 0 if no matches
    }
  } # end for-loop
  
  #### 3 - Assign Weights ####
  # assign weights to next-words corresponding to length of matched string
  weights <- rep(w[1:length(key)]/length_i, length_i) #normalize within lengths
  next_words_table <- cbind(all_next_words, weights = weights)
  
  #### 4 - Pick a Word ####
  # pick a next-word
  if (nrow(next_words_table) > 0) { # if there are next-words found
    # sample from all possible next-words w/ assigned prob. weights
    next_token <- sample(next_words_table[ , "all_next_words"], 
                         1, 
                         prob = next_words_table[ , "weights"])
  } else { # if no next-words are found
    next_token <- sample(na.omit(M1), 1) # sample from all words in text
  }
  
  #### 5 - Output ####
  return(next_token)
}


######## Question 8 ###########

# remove punctuation to make list exclusively of words to sample from
M2 <- a[grepl("[A-Za-z]", a)] 

# generate token of random common word from text
generate_word <- function(word_list) {
  repeat {
    start_word <- sample(word_list, 1)
    start_token <- match(start_word, b)
    if (!is.na(start_token)) {
      break
    }
  }
  return(start_token)
}

start_token <- generate_word(M2)

######## Question 9 ###########

# sentence generator
# generate token based on previous token(s) until full stop is generated
sentence <- function(start_token, w = c(1, 1, 1, 1)) {
  
  # initialize loop with start token
  token.v <- start_token
  start_word <- b[start_token]
  output_list <- c(start_word)
  
  # generate first mlag words
  for(i in 1:mlag) {
    nw.token <- next.word(token.v, M, M1, w)
    token.v <- append(token.v, nw.token)
    nw <- b[nw.token]
    output_list <- append(output_list, nw)
    if (nw == ".") break
  }
  
  # generate word if start token > mlag words
  while(nw != ".") {
    token.v <- token.v[2:(1+mlag)] # use last mlag - 1 tokens
    nw.token <- next.word(token.v, M, M1, w)
    nw <- b[nw.token]
    token.v <- append(token.v, nw.token)
    output_list <- append(output_list, nw)
  }
  
  # generate sentence with any punctuation collapsed onto the previous word
  output <- c()
  for (i in 1:length(output_list)) {
    if (grepl("[A-Za-z]", output_list[i]) == TRUE) {
      output <- append(output, output_list[i])
      output <- paste(output, collapse = " ")
    } else {
      output <- append(output, output_list[i])
      output <- paste(output, collapse = "")
    }
  }
  cat(output)
}

sentence(start_token)
sentence(start_token, w = c(1000, 100, 10, 1))
