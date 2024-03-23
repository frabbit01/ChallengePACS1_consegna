# ChallengePACS1_consegna

# ChallengePACS1
First challenge for polimi PACS course, a.a 2023-2024. (Visalli Francesca cp 10697667)

The main file is called challenge_jason.cpp

# Parameters choices and data.json
The parameters can be input from a json file as the one that has been uploaded in the repository with name "data.json".
It would be better to just modify the provided file as wished, otherwise it may be necessary to change the names of the fields as they are called in the main.
The chosen values are then grouped in a struct by the main program. Everything is initialized to the proposed parameters for testing. For convergence sake, I would advice to keep alpha at 0.1. 
The fields in the json file are:
- tolerance for the norm of the gradient
- tolerance to control the step length
- the sigma parameter (if you choose to apply the Armijo rule) or mu parameter (for inverse decay or exponential decay)
- function to be optimized as a string
- first component of the gradient as a string
- second component of the gradient as a string (I considered functions with two variables)
- initial step alpha
- maximum number of iterations
- first component of the initial guess
- second component of the initial guess
- learning rate (for the nesterov method)

# Additional user choices
At runtime the user will be asked to input their preferred choices for methods and strategies or eventually finite differences methods. Further details are present in the code comments.

# Makefile modifications
In order to compile the main file (which is challenge_jason.cpp), two libraries are necessary: json and muparser. I included the header files "json.hpp" and "muparser.h" that are not in this repository, so the paths in the makefile stored in the variables "CPPFLAGS" and "LDFLAGS" may need to be changed depending on where the libraries are stored. 
It may also be required to execute the following command if the library is still not found after compiling:
$export LD_LIBRARY_PATH=path/of/the/directory/of/the/library/you/use:$LD_LIBRARY_PATH

# What you will find

The main function is called minimize, which will call one of three auxiliary functions, depending on the choice the user has input when asked. They will also be asked to choose their preferred strategy for choosing the stepsize and, if so desired, which finite difference method to be used by the program. Other auxiliary functions have been coded with this purpose. 
