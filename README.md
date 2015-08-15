# Expansion Format of Storing Floating point data

This project aims to determine if performance gains can be made by storing arbitrary precision floating point data in 
the expansion format, i.e., as $x=x_1+x_2+\cdots+x_n$, where the $x_i$ are all of the same precision and non overlapping. 
 It is further assumed that they are arranged in decreasing order of magnitude allowing straightforward truncation to 
 the desired precision.