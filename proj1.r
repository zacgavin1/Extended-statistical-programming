# Zachary Gavin s2222962, Shaehroz Khalid s2869421, Brandon Causing s2901457
# Zachary Gavin: Preprocessing, setting up M and b, tokenized etc (Q4-6), debugging and commenting
# Shaehroz Khalid: next.word function and code to find starting word, debuging and commenting
# Brandon Causing: sentence function, altering next.word to include weights, debugging and commenting
# division of work was close to 1/3 each 

#setwd("Extended-statistical-programming") ## comment out of submitted
#setwd("C:/Users/shaeh/Desktop/edinburgh-notes/Extended-Statistical-Programming/project-ESP/Extended-statistical-programming")
#setwd("C:/Users/Brandon Causing/Downloads/Extended Statistical Programming/Extended-statistical-programming")


## GENERAL DESCRIPTION
# The project aim is to make a "Shakespeare text generator" using a Markov chain idea.
# The first section pre-processes the text, mainly dealing with stage directions and punctuation.
# The next token prediction is based off the idea of feeding in a key of some length, and finding 
# all matches of that key in the text, where we can then sample from all the words that follow the key
# in the text. When we predict next words, we mix together matches from the previous 1, 2,...,n 
# words, for some chosen n, to ensure we have some matches.
# Finally we feed in a starting word, and then repeatedly sample to produce a full sentence. 


a <- scan("shakespeare.txt",what="character",skip=83,nlines=196043-83,
          fileEncoding="UTF-8")


###########################################################
######### ------- START OF PREPROCESSING ------- ##########
###########################################################

# We begin by removing parts of the text that are not part of the literature itself (i.e 
# stage directions, character names), and putting the text into a tokenisable form. We do 
# this to generate proper sentences based only on the dialogue and with accurate sampling 
# probabilities.

#### --- Identify and remove stage directions --- ####

# First we manually fix all bracketing errors in the text
a[463343] <- "So"    # The original text was erroneously "[So"
a[322751] <- "Rest.]" # missing end of direction

# Find starts of stage directions
direction_starts = grep("[", a, fixed = TRUE) 

# Find ends of stage directions (relative to the start)
direction_ends <- c()
for (direction in direction_starts) {
   # get all closed brackets in the 100 words following an open bracket
   close_brackets <- grep("]", a[direction:(direction+100)])
   # and then we choose the first appearance after each "["
   direction_ends <- append(direction_ends, close_brackets[1])
}

# Collect indices of all words within stage directions
direction_indexes <- c() # to hold positions of stage directions
for (i in 1:length(direction_starts)) {
  dir_start_to_end <- direction_starts[i]:(direction_starts[i]+direction_ends[i]-1)
  direction_indexes <- append(direction_indexes, dir_start_to_end)
}

# Remove all stage directions
a <- a[-direction_indexes] 


#### --- Removing character names and Arabic numerals --- #####

# note: some roman numeral I's are likely left in
a <- a[-which((toupper(a)==a) & (a!= "I") & (a != "A") &(a != "O")) ] 


#### --- Removing _ and - from words --- ####

# we combine words at hyphens rather than split them, as it made sense in more cases in the text
a <- gsub("-", "", a ) 
a <- gsub("_", "", a )


#### --- Isolating punctuation --- ####

# split_punc() splits words and punctuation into individual strings
# * v = string vector
# * marks = string vector of characters (symbols w/ regex meaning preceded by \\) 
split_punc <- function(v, marks){
  regex = paste(marks, collapse="|") 
  i_punc <- which(grepl(regex, v)) # indexes of a where there is punctuation
  vs <- rep("", length(v)+length(i_punc)) # initialise new vector of req length
  ips <- i_punc+ 1:length(i_punc) # assume here all punctuation is at end of words
  vs[ips]<- substr(v[i_punc], nchar(v[i_punc]), nchar(v[i_punc])) # puncs 
  vs[ips-1] <- substr(v[i_punc], 1, nchar(v[i_punc])-1) # words before puncs
  vs[-(c(ips, ips-1))] <- v[-i_punc] # all other words
  vs
}


#### --- Call function to split up the puncs and the words --- ####

# need to use \\ before symbols that have regex meaning like "."
a <- split_punc(a, marks = c(",","\\.",";","!",":","\\?"))


#### --- Make whole text lower case --- ####

# for simplicity in use
a <- tolower(a)

###########################################################
########### ------ END OF PREPROCESSING ------ ############
###########################################################


###########################################################
######### ------ SETTING UP THE GENERATOR ------- #########
###########################################################

# We proceed by ranking words by usage and contextualising strings of words to logically
# pick each following word, given a starting point.

######## ----- Find Most Common Words ----- ############

words <- unique(a) # set of all words used
e <- match(a,words) # where each element of a occurs first in a

occurences <- tabulate(e) # count number of times each word shows up
names(occurences) <- words  # list of occurrences, with entry names the words

occurences <- sort(occurences, decreasing=TRUE)
b <- names(occurences[1:1000])  # now have the top 1000 words


###### ----- Set Up Word Contextualisation ----- #######

n    <- length(a) 
mlag <- 4

# tokenize a; words not in b (most common words) are tokenenised as NA
M1 <- match(a,b) 

# -- rows of M are all (mlag+1) length sequences of (tokenised) words from text
# -- match input string to rows of M to find set of possible next-words
M <- M1[1:(n-mlag)] 
for (col in 2:(mlag+1)){ # loop appends shifted copies of M1 to the right of M
  M <- cbind(M, M1[col:(n-mlag-1 + col)])
}


####### ----- Pick the Next Word ----- ###########

# pick a next-word based on key (sequence of tokens), M, M1, w (weights)
next.word <- function(key, M, M1, w=rep(1,ncol(M)-1)){
  
  ###### -- 1 - Set-Up -- ######
  mlag <- ncol(M) - 1 # define mlag in terms of function argument M
  
  if (length(key) > mlag) { # only use the last mlag tokens
    key <- key[(length(key) - mlag + 1):length(key)]   
  }
  
  all_next_words <- c() # list of next-words to sample from
  length_i <- c() # no. of possible next-words in iteration i
  
  ##### -- 2 - Search for Next Words -- ######
  # find key match of length mlag, mlag - 1, ... , 1
  for (i in 1:length(key)) {
    
    ### - 2a - Pick columns to match key to - ###
    current_key <- key[i:length(key)] # use last mlag - i + 1 tokens of key
    context_len <- length(current_key)  
    cols_to_match <- (mlag - context_len + 1):mlag # which cols in M to match
    
    ### - 2b - Find matches - ###
    # return F if key[j] matches M[k, j] for some row k of M; T otherwise
    ii <- colSums(!(t(M[, cols_to_match, drop=FALSE])==current_key))
    # if sum of components in a row = 0, entire key matches
    matching_rows <- which(ii == 0)
    
    ### - 2c - Store next-words and no. of next-words - ###
    if (length(matching_rows) > 0) { # if there is a matching row
      
      # get last column of M (possible next-words) for all matching rows
      next_words_found <- M[matching_rows, mlag + 1]
      # collect found next-words, omitting rare words (NA's)
      all_next_words <- c(all_next_words, na.omit(next_words_found))
      # store no. of (non-rare) next-words found in iteration i (length_i[i])
      length_i <- c(length_i, length(na.omit(next_words_found)))
      
    } else { # if there is not a matching row, length_i[i] = 0
      length_i <- c(length_i, 0)
    }
  } # end for-loop
  
  #### -- 3 - Assign Weights -- ####
  # assign weights to next-words corresponding to length of matched string
  # -- vectorised form of rep is used here
  # -- weights = w_i/n_i 
  # -- n_i = no. of matches found with context length mlag-i+1
  weights <- rep(w[1:length(key)]/length_i, length_i)
  next_words_table <- cbind(all_next_words, weights = weights)
  
  #### -- 4 - Pick a Next-Word -- ####
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

###########################################################
############ ------ SENTENCE GENERATOR ------ #############
###########################################################

# Finally, we are able to randomly generate Shakespearean dialogue.

####### ----- Generate a Starting Word ----- ##########

# remove punctuation to make list exclusively of words to sample from
M2 <- a[grepl("[A-Za-z]", a)] 

# sample token of any (common) word from text
generate_word <- function(word_list, b_match = b) {
  
  start_tokens <- match(word_list, b_match)
  start_token <- sample(start_tokens[!is.na(start_tokens)], 1)
  return(start_token)
  
}

# generate the start word
start_token <- generate_word(M2)


######## ------ Generate a Sentence ------ ###########

# sentence generator; input either a single token or character
# -- generate token based on previous token(s) until full stop is generated
sentence <- function(start_token, 
                     M.seq = M, M1.tokens = M1, b_match = b,
                     w = rep(1, ncol(M) - 1)) {
  
  mlag <- ncol(M) - 1
  
  # initialize loop with start token
  token.v <- start_token
  start_word <- b_match[start_token]
  
  # flexibility to input characters instead of tokens, whether common or not
  if (is.character(start_token) == TRUE) { # if character
    start_word <- start_token
    if (start_token %in% b) { # if common (character)
      token.v <- which(b == start_token) 
    } 
  } 
  
  # initialize sentence
  output_list <- start_word
  
  # generate first mlag words
  for(i in 1:mlag) {
    nw.token <- next.word(token.v, M.seq, M1.tokens, w)
    token.v <- append(token.v, nw.token)
    nw <- b_match[nw.token]
    output_list <- append(output_list, nw)
    if (nw == ".") break
  }
  
  # generate word if start token longer than mlag words
  while(nw != ".") {
    token.v <- token.v[2:(1+mlag)] # use last mlag - 1 tokens
    nw.token <- next.word(token.v, M.seq, M1.tokens, w)
    nw <- b_match[nw.token]
    token.v <- append(token.v, nw.token)
    output_list <- append(output_list, nw)
  }
  
  # generate sentence (with punctuation collapsed onto the previous word)
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

#### ---- Make a sentence ---- ####
# using default weights
sentence(start_token)

# heavier weights on longer matches 
# -- more likely to produce more coherent sentences
# -- more likely to be shorter
sentence(start_token, w = c(1000, 100, 10, 1)) 
