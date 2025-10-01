# Zachary Gavin s2222962, Shaeroz Khalid ..., Brandon Causing ...
# Zachary Gavin: Preprocessing, setting up M and b, tokenized etc (Q4-6), debugging and commenting
# Shaeroz Khalid: next.word function and code to find starting word, debuging and commenting
# Brandon Causing: sentence function, altering next.word to include weights, debugging and commenting
# division of work was close to 1/3 each 

#setwd("Extended-statistical-programming") ## comment out of submitted
#setwd("C:/Users/shaeh/Desktop/edinburgh-notes/Extended-Statistical-Programming/project-ESP/Extended-statistical-programming") ## comment out of submitted
#setwd("C:/Users/Brandon Causing/Downloads/Extended Statistical Programming/Extended-statistical-programming")

a <- scan("shakespeare.txt",what="character",skip=83,nlines=196043-83,
          fileEncoding="UTF-8")


## GENERAL DESCRIPTION
# The project aim is to make a "Shakespeare text generator" using a Markov chain idea.
# The first section preprocesses the text, mainly dealing with stage directions and punctuation
# The next token prediction is based off the idea of feeding in a key of some length, and finding 
# all matches of that key in the text, where we can then sample from all the words that follow the key
# in the text. When we predict words, we mix together matches from the previous 1, 2,...,n 
# words, for some chosen n, to ensure we have some matches. 
# Finally we feed in a starting word, and then repeatedly sample to produce a full sentence. 


####### ------- START OF PREPROCESSING ------- ###############

# {include description of pre-processing here}

#### --- Identify and remove stage directions --- ####

# First we manually fix all bracketing errors in the text
a[463343] <- "So"    # The original text was erroneously "[So"
a[322751] <- "Rest.]" # missing end of direction

direction_starts = grep("[", a, fixed = TRUE) # find start of stage directions

direction_ends <- c()
for (direction in direction_starts) {
  
   # get all closed brackets in the 100 words following an open bracket
   close_brackets <- grep("]", a[direction:(direction+100)])
   # and then we choose the first appearance after each [
   direction_ends <- append(direction_ends, close_brackets[1])

}

# {brief sentence of what this is doing}
direction_indexes <- c() # to hold positions of stage directions
for (i in 1:length(direction_starts)){ 
  
  direction_indexes <- append(direction_indexes, 
                              direction_starts[i]:(direction_starts[i]+direction_ends[i]-1))

}

a <- a[-direction_indexes] # removes all stage directions 


#### --- Removing character names and Arabic numerals --- #####

# note: some roman numeral I's are likely left in
a <- a[-which((toupper(a)==a) & (a!= "I") & (a != "A") &(a != "O")) ] 


#### --- Removing _ and - from words --- ####

# we combine words at hyphens rather than split them, as it made  sense in more 
# of the cases in the text
a <- gsub("-", "", a ) 
a <- gsub("_", "", a )


#### --- Isolating punctuation --- ####

# v is a string vector, marks is a string vector of characters (where symbols that have 
# regex meaning are preceded by \\). split_punc picks out all strings in v containing any punctuation,
# splits the words and the punctuations, and inserts all created strings into new vector vs
split_punc <- function(v, marks){
  regex = paste(marks, collapse="|") 
  i_punc <- which(grepl(regex, v)) # returns the indexes of a where there is punctuation
  vs <- rep("", length(v)+length(i_punc)) # initialise new vector of req length
  ips <- i_punc+ 1:length(i_punc) # assume here all punctuation is at end of words
  vs[ips]<- substr(v[i_punc], nchar(v[i_punc]), nchar(v[i_punc])) 
  vs[ips-1] <- substr(v[i_punc], 1, nchar(v[i_punc])-1) # words before puncs
  vs[-(c(ips, ips-1))] <- v[-i_punc] # all other words
  vs
}


#### --- Call function to split up the puncs and the words --- ####
a <- split_punc(a, marks = c(",","\\.",";","!",":","\\?")) # need to use \\ before symbols that have regex meaning like .


#### --- Making whole text lower case --- ####
a <- tolower(a)


########### ------ END OF PREPROCESSING ------ ############

### {summary of what goes on in this section}

######## ------ Find Most Common Words ------ ############

words <- unique(a) # words is the set of all words used
e <- match(a,words) # where each element of a occurs first in a (length(e)=length(a))

occurences <- tabulate(e) # length is length(words), with the number of times each word shows up counted
names(occurences) <- words  # list of occurrences, with entry names the words

occurences <- sort(occurences, decreasing=TRUE)
b <- names(occurences[1:1000])  # now have the top 1000 words


######## ------ Set Up Match Sequences ------ ###########

n    <- length(a) 
mlag <- 4
# {summary of what goes on here}
M1 <- match(a,b) # now have a tokenised version of a. In instructions he calls this M1


# {summary of what goes on here}
M <- M1[1:(n-mlag)] 
for (col in 2:(mlag+1)){ # loop appends shifted copies of M1 to the right of M
  M <- cbind(M, M1[col:(n-mlag-1 + col)])
}


####### ------ Pick the Next Word ------ ###########

# pick a next-word based on key (sequence of tokens), M, M1, w (weights)
next.word <- function(key, M, M1, w=rep(1,ncol(M)-1)){
  
  ###### -- 1 - Set-Up -- ######
  mlag <- ncol(M) - 1 # define mlag in terms of function argument M
  
  if (length(key) > mlag) { # only use the last mlag tokens
    key <- key[(length(key) - mlag + 1):length(key)]   
  }
  
  all_next_words <- c() # initialize list of next-words to sample from
  length_i <- c() # no. of next-words in iteration i
  
  ##### -- 2 - Search for Next Words -- ######
  # find key match of length mlag, length mlag - 1, ..., 1
  for (i in 1:length(key)) {
    
    ### - 2a - Pick columns to match key - ###
    current_key <- key[i:length(key)] # use last mlag - i + 1 tokens of key
    context_len <- length(current_key)  
    cols_to_match <- (mlag - context_len + 1):mlag # which cols in M to match
    
    ### - 2b - Find matches - ###
    # return F if key[j] matches M[k, j] for some row k of M, T otherwise
    ii <- colSums(!(t(M[, cols_to_match, drop=FALSE])==current_key))
    # if sum of components in a row = 0, entire key matches
    matching_rows <- which(ii == 0)
    
    ### - 2c - Store next-words and no. of next-words - ###
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
  
  #### -- 3 - Assign Weights -- ####
  # assign weights to next-words corresponding to length of matched string. vectorised form of rep is used here
  weights <- rep(w[1:length(key)]/length_i, length_i) # weights are w_i/n_i, where n_i is number of matches found with context length (mlag-i+1)
  next_words_table <- cbind(all_next_words, weights = weights)
  
  #### -- 4 - Pick a Word -- ####
  # pick a next-word
  if (nrow(next_words_table) > 0) { # if there are next-words found
    # sample from all possible next-words w/ assigned prob. weights
    next_token <- sample(next_words_table[ , "all_next_words"], 
                         1, 
                         prob = next_words_table[ , "weights"])
  } else { # if no next-words are found
    next_token <- sample(na.omit(M1), 1) # sample from all words in text
  }
  
  #### -- 5 - Output -- ####
  return(next_token)
}


######## ------ Generate a Start Word ------ ###########

# remove punctuation to make list exclusively of words to sample from
M2 <- a[grepl("[A-Za-z]", a)] 

# sample token of random (common) word from text
generate_word <- function(word_list, b_match = b) {
  
  start_tokens <- match(word_list, b_match)
  start_token <- sample(start_tokens[!is.na(start_tokens)], 1)
  return(start_token)
  
}

# generate the start word
start_token <- generate_word(M2)


######## ------ Generate a Sentence ------ ###########

# sentence generator
# generate token based on previous token(s) until full stop is generated
sentence <- function(start_token, 
                     M.seq = M, M1.tokens = M1, b_match = b,
                     w = c(1, 1, 1, 1)) {
  
  mlag <- ncol(M) - 1
  
  # initialize loop with start token
  token.v <- start_token
  start_word <- b_match[start_token]
  output_list <- c(start_word)
  
  # generate first mlag words
  for(i in 1:mlag) {
    nw.token <- next.word(token.v, M.seq, M1.tokens, w)
    token.v <- append(token.v, nw.token)
    nw <- b_match[nw.token]
    output_list <- append(output_list, nw)
    if (nw == ".") break
  }
  
  # generate word if start token > mlag words
  while(nw != ".") {
    token.v <- token.v[2:(1+mlag)] # use last mlag - 1 tokens
    nw.token <- next.word(token.v, M.seq, M1.tokens, w)
    nw <- b_match[nw.token]
    token.v <- append(token.v, nw.token)
    output_list <- append(output_list, nw)
  }
  
  # generate sentence with any punctuation collapsed onto the previous word
  output <- c()
  for (i in 1:length(output_list)) {
    if (grepl("[A-Za-z]", output_list[i]) == TRUE) { #if string contains letter
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
