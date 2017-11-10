# Markovian_SIR_Deaths_Model

I noticed that traditional methods to predict a disease outbreak was by performing sentiment analysis on Twitter posts and Google Search terms.  Unfortunately, these methods were inadequate, as Twitter and Google is not popular in all countries. So, I created a system to model and predict outbreaks without the need for social media.  

The system was able to update the probabilities of a virus from spreading from A to B in real time, and I plan to release it to the public next year.  I also used Machine Learning and Deep Learning to predict larger long-term virus trends with Google Trends, and this acted as a validator for the MSIRD model.

# Research Highlights

1. Mosquito Distribution Formulation
   a. Temperature, Environmental Factors
   b. Living Standards
   
2. Basic Markov Modelling
   a. Growth and Decay
   b. New Route Derivation - Encoding and Decoding graphs
   
3. MSIR Model
   a. Matrix block concentation
   b. States and average states
   c. [S->I], [S->S], [S->R], [I->R], [I->I], [R->R] blocks
   d. Markov Updating
   
4. Demonstration
   a. 9 squares
   b. 49 squares
   
5. MSIRD Model
   a. Revised Formulas
   b. Deaths incoporation
      i. [S->D], [I->D], [R->D]
