#setwd("Extended-statistical-programming") ## comment out of submitted
#setwd("C:/Users/shaeh/Desktop/edinburgh-notes/Extended-Statistical-Programming/project-ESP/Extended-statistical-programming") ## comment out of submitted


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
#we define the function next.word below. key represents the sequence of word tokens we are using
#as context to predict the next word. M is the matrix of all 5-word sequences from Shakespeare
#M1 is the tokenised version of the whole text
#w is the mixture weights
next.word <- function(key, M, M1, w=rep(1,ncol(M)-1)){

#mlag is set to 4 as M has 5 columns.
  mlag <- ncol(M) - 1

#Below, we make sure the key is not too long. The model is built for only 4 words being mlag here.
#Sps key was c(10, 20, 30, 40, 50, 60) this would become c(30, 40, 50, 60)

    if (length(key) > mlag) {   
    key <- key[(length(key) - mlag + 1):length(key)]   
  }
#Here we create an empty vector called "all_next_words". This is all the possible words that can
#come next.This is our "raffle drum" as we drop a "ticket" for the word that shows up into the vector
  all_next_words <- c()
  
#This for loop below allows us to try out all different lengths. If the input key has 4 words
#the loop will run 4 times: The first time using the full 4-word context. The second time it
#will search using the last 3 words of the context. The third time using the last 2 words.
#Finally using the last word.

  for (i in 1:length(key)) {
    
    current_key <- key[i:length(key)]
    context_len <- length(current_key)
    cols_to_match <- (mlag - context_len + 1):mlag

#If length(key) = 4. An example for the above code will be, sps i = 2. Current_key <- key[2,4]
#Context_len = 3 which is the number of words remaining in key. cols_to_match will be (4-3+1):4
# which is just 2:4 as required because we are going from 2:4.

#For this command from the PDF below, check WA gc
#An example would be sps we want c(1, 2, 3) and we find that pattern in rows 5, 200 and 512
#Then matching_rows will be c(5, 200, 512)
    
    ii <- colSums(!(t(M[, cols_to_match, drop=FALSE])==current_key), na.rm = TRUE)

    matching_rows <- which(ii == 0)

#Here we collect the "prizes" for our spsd "raffle" This code only runs if we found any matches
#The first line in the if statement looks at the matching rows and grabs the value from the 5th column
#The second line takes the tokens and adds it to our "raffle drum", na.omit cleans out any NA's (rare words)
#i.e. if the sentence "to be," was followed by "alas" (token NA) once and "my" (token 15) twice,
#This step would add 15 and 15 into all_next_words vector and exclude "alas"
    
    if (length(matching_rows) > 0) {
      next_words_found <- M[matching_rows, mlag + 1]
      all_next_words <- c(all_next_words, na.omit(next_words_found))
    }
  }

#The first part of the if the statement samples one "ticket" from our "raffle drum", if "all_next_words"
#has any tickets in it. If a token was found more often -> has more tickets -> higher prob of being selected
#the else part of the if statement just gives up and chooses a random word from the entire book
#The entire book here is 'M1' which is the tokenised vector
  
  if (length(all_next_words) > 0) {
    next_token <- sample(all_next_words, 1)
  } else {
    next_token <- sample(na.omit(M1), 1)
  }

#Below we return the single token that was chosen which is the model's final answer for the word that
#comes next.
  return(next_token)
}

######## Question 8 ###########

#M2 here is the vector of all the words from the vector a excluding punctuation now as well

M2 <- a[grepl("[A-Za-z]", a)] 

#start_word is a random word chosen from M2
#We repeat finding this start_word until we don't have an N/A (low prob anyway)
repeat {
  start_word <- sample(M2, 1)
  
  start_token <- match(start_word, b)
  
  if (!is.na(start_token)) {
    break
  }
}

#here we get our start_word and the token for it
#we get rid of N/A's as we don't have those rare words in tokenised form
print(start_word)
print(start_token)

######## Question 9 ###########
#Now we will simulate from the model until a full stop is reached.Then we can convert the generated tokens back to words and print them nicely

########################

# sentence generator
sentence <- function(start_token) {
  
  #initiate loop
  token.v <- start_token
  start_word <- b[start_token]
  output <- c(start_word)
  
  #generate first four words
  for(i in 1:4) {
    nw.token <- next.word(token.v, M, M1)
    token.v <- append(token.v, nw.token)
    nw <- b[nw.token]
    output <- append(output, nw)
    if (nw == ".") break
    #print(token.v) #uncomment to see tokens inputted + token output
    #print(nw)
  }
  
  #generate if longer than four words
  while(nw != ".") {
    token.v <- token.v[2:5] #use last four tokens to generate next token
    nw.token <- next.word(token.v, M, M1)
    nw <- b[nw.token]
    token.v <- append(token.v, nw.token)
    output <- append(output, nw)
    #print(token.v)
    #print(nw)
  }
  
  cat(output)
  
}

# generate random word from text
generate_word <- function(word_list) {
  repeat {
    start_word <- sample(word_list, 1)
    start_token <- match(start_word, b)
    if (!is.na(start_token)) {
      break
    }
  }
  return(c(start_word, start_token))
}

start_token <- as.numeric(generate_word(M2)[2])

sentence(start_token)


