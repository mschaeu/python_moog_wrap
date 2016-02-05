In this repository, I have placed a version of MOOG that will result in a .so file that can be called by any other language. To compile this .so file, simply run the ‘make’ command.


I’ve also included a python file that briefly shows how to call this moog.so library from python. Currently, this python code is only set up to work with the synth and abfind drivers. In the future, I will expand it to work with all of the drivers. 

To make this python code work correctly, please update line 16 to reflect the location of the moog.so file on your system. The rest of the .py file should provide enough guidance to tell you how to use the .so library and the python output.

NOTE: If the code should not work correctly after you have adjusted line 16, there might be an additional thing you need to change: the name of the actual, callable library, as specified in line 17. To confirm this name, please type:

nm moog.so

in your terminal and carefully look at the output (which will be quite long). Make sure the name _moogsilent_ appears. If it doesn’t, replace line 17 with the corresponding name. Always keep in mind to omit the first underscore, but not the 2nd one. 


