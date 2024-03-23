#include <json.hpp>
#include <muParser.h>
#include "muparser_fun.hpp"
#include <iostream>
#include<vector>
#include<functional>
#include <cmath>
#include<fstream>

using json = nlohmann::json;

//PARAMETERS
//I define a struct to hold all the parameters that are required by the main function
//I put in some default values, which were the ones suggested in the exercise, since some fields would have not been default initialized and to speed up testing
//note: I did not put two different fields (for the mu and sigma parameters) but just one for sigma, as they are never used together, so the sigma parameter corresponds also to mu if the selected strategy is not Armijo
//after having tested different parameters, one of the better values for alpha to ensure convergence for all methods is 0.1
struct Parameters{
    double er=1.0e-6;
    double es=1.0e-6;
    double sigma=0.3;
    std::function<double (double, double)> f=[] (double x,double y){return x*y+4*pow(x,4)+y*y+3*x;};
    std::vector<std::function<double (double,double)>> grad={[] (double x,double y){return y+16*pow(x,3)+3;},[] (double x,double y){return x+2*y;}};
    double alpha=0.1;
    int maxiter=1000;
    std::vector<double> x0={0.0,0.0};
    double learning_rate=0.1;
};

//FUNCTION DECLARATIONS

//auxiliary functions

//vector operations
double norm(std::vector<double> const & vec);
std::vector<double> sub(std::vector<double> a, std::vector<double> b);
std::vector<double> sum(std::vector<double> a, std::vector<double> b);
std::vector<double> mul(std::vector<double> const & a, double c);

//different ways of selecting the learning rate
double Armijo(Parameters parameters,int k,std::vector<double> x);
double exponential_decay(Parameters parameters,int k,std::vector<double> x);
double inverse_decay(Parameters parameters,int k,std::vector<double> x);

//Different strategies to minimize
std::vector<double> gradient_descent(Parameters const & parameters);
std::vector<double> Nesterov(Parameters const & parameters);
std::vector<double> momentum(Parameters const & parameters);

//FINITE DIFFERENCES
std::vector<double> centred_differences(std::function<double (double, double)> f,std::vector<double> x); //centred finite difference method
std::vector<double> forward_differences(std::function<double (double, double)> f, std::vector<double> x);//forward finite difference method
std::vector<double> bkd_differences(std::function<double (double, double)> f, std::vector<double> x);//bkd finite difference method

//user's choices
//Theses will be called in the minimization strategies method, to ask the user to choose a strategy or a finite difference method if so desired
bool choice(std::function<double (Parameters,int,std::vector<double>)> & method, Parameters parameters);
bool choice_fd(std::function<std::vector<double> (std::function<double (double, double)>,std::vector<double>)>  & method, Parameters parameters);

//here I define my main function
std::vector<double> minimize(Parameters const & parameters);

int main(){
    std::ifstream file("data.json");
    json data = json::parse(file);
    Parameters parameters;
    //I initialize the parameters struct from the json file
    std::string funString = data.value("f","");
    std::vector<std::string> gradString=data["grad"].get<std::vector<std::string>>();
    parameters.grad={};
    for(auto x:gradString)
        parameters.grad.push_back(MuparserFun(x));
    parameters.f=MuparserFun(funString);
    parameters.er=data.value("er",1e-6);
    parameters.es=data.value("es",1e-6);
    parameters.alpha=data.value("alpha",0.1);
    parameters.learning_rate=data.value("learning_rate",0.1);
    parameters.maxiter=data.value("maxiter",1000);
    parameters.sigma=data.value("sigma",0.3);
    parameters.x0=data["x0"].get<std::vector<double>>(); 
    

    //Here I call the main function requested by the exercise
    auto x_min= minimize(parameters);
    //I print the results 
    std::cout<<x_min[0]<<","<<x_min[1]<<std::endl;
    std::cout<<"The value of the function is: "<<parameters.f(x_min[0],x_min[1])<<std::endl;
}
//euclidian norm
double norm(std::vector<double> const & vec){
    double sum=0;
    for(auto x:vec){
        sum+=x*x;
    }
    return sqrt(sum);
}
//I define sum and subrtraction between vectors and scalar-vector multiplication
//In the sum and subtraction functions there is also a check for the vectors dimentions
std::vector<double> sub(std::vector<double> a, std::vector<double> b){
    if(a.size()!=b.size()) {
        std::cerr << "Dimensions are not compatible"<<std::endl;
        return std::vector<double>{};
    }
    else{
        std::vector<double> res;
        for(std::size_t i=0;i<a.size();++i) {
            res.push_back(a[i] - b[i]);
        }
        return res;
    }
}

std::vector<double> sum(std::vector<double> a, std::vector<double> b){
    if(a.size()!=b.size()) {
        std::cerr << "Dimensions are not compatible"<<std::endl;
        return std::vector<double>{};
    }
    else{
        std::vector<double> res;
        for(std::size_t i=0;i<a.size();++i) {
            res.push_back(a[i] + b[i]);
        }
        return res;
    }
}

std::vector<double> mul(std::vector<double> const & a,double c){
    std::vector<double> res;
    for(auto x: a){
        res.push_back(c*x);
    }
    return res;
}

//step size related functions
//Armijo method
double Armijo(Parameters parameters,int k,std::vector<double> x){  
    double alpha=parameters.alpha; 
    //check on the allowed range of the sigma parameter
    if(parameters.sigma>0.5 || parameters.sigma<0){
        std::cerr<<"The value for sigma cannot be used, please input another value in the interval (0,0.5)"<<std::endl;
        return -1.0;
    }
    std::vector<double> grad;
    for(auto f:parameters.grad)
        grad.push_back(f(x[0],x[1]));
    while(parameters.f(x[0],x[1])-parameters.f(sub(x,mul(grad,alpha))[0],sub(x,mul(grad,alpha))[1])<parameters.sigma*alpha*norm(grad)*norm(grad))
        alpha/=2;
    return alpha;
}

//exponential decay
double exponential_decay(Parameters parameters,int k,std::vector<double> x) {
        return parameters.alpha*exp(-1*parameters.sigma*k);
}

//inverse decay
double inverse_decay(Parameters parameters,int k,std::vector<double> x){  
    return parameters.alpha/(1+parameters.sigma*k);
}

//minimizing strategies
//At the beginning of each of the minimizing strategies functions I ask the user to choose a method to compute alpha and call back the same minimizing
//function if a wrong choice is input

//gradient descent
std::vector<double> gradient_descent(Parameters const & parameters){
    std::vector<double> x_old;
    std::vector<double> x=parameters.x0;
    double alpha=parameters.alpha;
    //strategy choice for alpha
    bool check;
    std::function<double (Parameters,int,std::vector<double>)> method;
    check=choice(method,parameters);
    //check on wether the user has typed an allowed value
    if(!check) {
        std::cerr << "You have not input one of the three allowed values, please compile again";
        return gradient_descent(parameters);
    }
    //finite differences choice
    std::function<std::vector<double> (std::function<double (double, double)>,std::vector<double>)> fd_method;
    bool fd=choice_fd(fd_method,parameters);
    for(int k=0;k<parameters.maxiter;++k){
        x_old=x;
        alpha= method(parameters,k,x);
        std::vector<double> grad;
        if(!fd){
        for(auto f:parameters.grad)
            grad.push_back(f(x[0],x[1]));
        }
        else{
            grad=fd_method(parameters.f,x);
        }
        x=sub(x,mul(grad,alpha));
        if(norm(sub(x,x_old))<parameters.es){
            break;
        }
        if(norm(grad)<parameters.er)
            break;
    }
    
    return x;
}

//Nesterov method
std::vector<double> Nesterov(Parameters const & parameters){
    std::vector<double> x_old=parameters.x0;
    std::vector<double> grad;
    std::vector<double> x=sub(x_old,mul(grad,parameters.alpha));
    double alpha=parameters.alpha;
    //strategy choice for alpha
    bool check;
    std::function<double (Parameters,int,std::vector<double>)> method;
    check=choice(method,parameters);
    //check on wether the user has typed an allowed value
    if(!check) {
        std::cerr << "You have not input one of the three allowed values, please compile again";
        return Nesterov(parameters);
    }
    //finite differences choice
    std::function<std::vector<double> (std::function<double (double, double)>,std::vector<double>)> fd_method;
    bool fd=choice_fd(fd_method,parameters);
    if(!fd){
        for(auto f:parameters.grad)
            grad.push_back(f(x_old[0],x_old[1]));
    }
    else{
        grad=fd_method(parameters.f,x_old);
    }
    std::vector<double> y(x.size());
    double eta=parameters.learning_rate;
    for(int k=0;k<parameters.maxiter;++k){
        alpha= method(parameters,k,x);
        y=sum(x,mul(sub(x,x_old),eta));
        grad={};
        if(!fd){
        for(auto f:parameters.grad)
            grad.push_back(f(x[0],x[1]));
        }
        else{
            grad=fd_method(parameters.f,x);
        }
        x=sub(x,mul(grad,alpha));
        x_old=x;
        if(norm(sub(x,x_old))<parameters.es)
            break;
        std::vector<double> gradx;
        for(auto f:parameters.grad)
            gradx.push_back(f(x[0],x[1]));
        if(norm(gradx)<parameters.er)
            break;
    }
    return x;
}

//gradient descent with momentum
std::vector<double> momentum(Parameters const & parameters){
    std::vector<double> x_old;
    std::vector<double> x=parameters.x0;
    double alpha=parameters.alpha;
    //strategy choice for alpha
    bool check;
    std::function<double (Parameters,int,std::vector<double>)> method;
    check=choice(method,parameters);
    //check wether the user has typed an allowed value
    if(!check) {
        std::cerr << "You have not input one of the three allowed values, please compile again";
        return momentum(parameters);
    }
    //finite differences choice
    std::function<std::vector<double> (std::function<double (double, double)>,std::vector<double>)> fd_method;
    bool fd=choice_fd(fd_method,parameters);
    std::vector<double> grad0;
    if(fd){
        grad0=fd_method(parameters.f,x);
    }
    else{
    for(auto f:parameters.grad)
            grad0.push_back(f(x[0],x[1]));
    }
    auto d=mul(grad0,-1*parameters.alpha);
    double eta=parameters.learning_rate;
    for(int k=0;k<parameters.maxiter;++k){
        x_old=x;
        alpha= method(parameters,k,x);
        x=sum(x,d);
        std::vector<double> grad;
        if(fd){
        grad=fd_method(parameters.f,x);
        }
        else{
        for(auto f:parameters.grad)
            grad.push_back(f(x[0],x[1]));
        }
        d=sub(mul(d,eta),mul(grad,alpha));
        if(norm(sub(x,x_old))<parameters.es)
            break;
        if(norm(grad)<parameters.er)
            break;
    }
    return x;
}

//FINITE DIFFERENCES
// you will see the definition of centred differences, forward differences and backward differences (in this order). I chose the parameter h (which should ideally tend to 0) as a double equal to 1e-6
std::vector<double> centred_differences(std::function<double (double, double)> f,std::vector<double> x){
    std::vector<double> grad;
    double h=1e-6;
    grad.push_back((f(x[0]+h,x[1])-f(x[0]-h,x[1]))/(2*h));
    grad.push_back((f(x[0],x[1]+h)-f(x[0],x[1]-h))/(2*h));
    return grad;
}

std::vector<double> forward_differences(std::function<double (double, double)> f, std::vector<double> x){
    std::vector<double> grad;
    double h=1e-6;
    grad.push_back((f(x[0]+h,x[1])-f(x[0],x[1]))/(2*h));
    grad.push_back((f(x[0],x[1]+h)-f(x[0],x[1]))/(2*h));
    return grad;
}

std::vector<double> bkd_differences(std::function<double (double, double)> f, std::vector<double> x){
    std::vector<double> grad;
    double h=1e-6;
    grad.push_back((f(x[0],x[1])-f(x[0]-h,x[1]))/(2*h));
    grad.push_back((f(x[0],x[1])-f(x[0],x[1]-h))/(2*h));
    return grad;
}

//user choices

bool choice(std::function<double (Parameters,int,std::vector<double>)>  & method, Parameters parameters){
    bool res=false;
    int choice=-1;
    //I ask the user to choose their preferred minimizing strategy and initialize a variable that keeps track of the choice and use it to call one auxiliary function
    std::cout<<"What strategy would you like to follow?\nType 1 for exponential decay,\nType 2 for inverse decay\nType 3 for Armijo"<<std::endl;
    std::cout<<"It is advisable not to choose Armijo with the momentum method"<<std::endl;
    std::cin>>choice;
    if(choice==1){
        method=exponential_decay;
        res=true;
    }
    else if(choice==2){
        method=inverse_decay;
        res=true;
    }
    else if(choice==3){
        method=Armijo;
        res=true;
    }
    return res;
}

bool choice_fd(std::function<std::vector<double> (std::function<double (double, double)>,std::vector<double>)>  & method, Parameters parameters){
    bool res=false;
    int choice=-1;
    char fd;
    //I ask the user to choose whether or not they want to choose a finite difference method
    while(fd!='y' &&fd!='n'){
        std::cout<<"Would you like to use a finite differences method instead of the exact gradient? (y|n)"<<std::endl;
        std::cin>>fd;
    }
    if(fd=='n')
        return false;
    //if they chose to here I ask them their preferred strategy
    std::cout<<"What strategy would you like to follow?\nType 1 for forward differences,\nType 2 for backward differences\nType 3 for centred differences"<<std::endl;
    std::cin>>choice;
    if(choice==1){
        method=forward_differences;
        res=true;
    }
    else if(choice==2){
        method=bkd_differences;
        res=true;
    }
    else if(choice==3){
        method=centred_differences;
        res=true;
    }
    return res;
}

//principal function
std::vector<double> minimize(Parameters const & parameters){
    int choice=-1;
    //I ask the user to choose their preferred minimizing strategy and initialize a variable that keeps track of the choice and use it to call one auxiliary function
    std::cout<<"Choose the preferred strategy for minimizing:\ntype 1 for gradient descent,\n2 for momentum\n or 3 for Nesterov"<<std::endl;
    std::cin>>choice;
    std::vector<double> x_min (parameters.x0.size());
    if(choice==1)
        x_min=gradient_descent(parameters);
    else if(choice==2)
        x_min=momentum(parameters);
    else if(choice==3)
        x_min=Nesterov(parameters);
    //If the value of choice does not correspond to anything I call this function again
    else{
        std::cerr<<"You have not input one of the three allowed values, please compile again"<<std::endl;
        return minimize(parameters);
        }
    return x_min;
}




