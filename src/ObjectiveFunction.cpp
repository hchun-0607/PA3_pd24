#include "ObjectiveFunction.h"

#include "cstdio"

ExampleFunction::ExampleFunction(Placement &placement) : BaseFunction(1), placement_(placement)
{
    printf("Fetch the information you need from placement database.\n");
    printf("For example:\n");
    printf("    Placement boundary: (%.f,%.f)-(%.f,%.f)\n", placement_.boundryLeft(), placement_.boundryBottom(),
           placement_.boundryRight(), placement_.boundryTop());
}

const double &ExampleFunction::operator()(const std::vector<Point2<double>> &input)
{
    // Compute the value of the function
    value_ = 3. * input[0].x * input[0].x + 2. * input[0].x * input[0].y +
             2. * input[0].y * input[0].y + 7.;
    input_ = input;
    return value_;
}

const std::vector<Point2<double>> &ExampleFunction::Backward()
{
    // Compute the gradient of the function
    grad_[0].x = 6. * input_[0].x + 2. * input_[0].y;
    grad_[0].y = 2. * input_[0].x + 4. * input_[0].y;
    return grad_;
}


const double &Wirelength::operator()(const std::vector<Point2<double>> &input){
    value_ = 0;
    input_size = input.size();
    for(int i = 0; i < input_size; i++){ 
        exp_xi_gamma[i] = exp(input[i].x / gamma);
        exp_yi_gamma[i] = exp(input[i].y / gamma);
    }

    for(int i = 0; i < placement.numNets(); i++){
        double sum1 = 0, sum2 = 0;
        double sum3 = 0, sum4 = 0;
        for(int j = 0; j < placement.net(i).numPins(); j++){
            sum1 += exp_xi_gamma[placement.net(i).pin(j).moduleId()];
            sum2 += 1/exp_xi_gamma[placement.net(i).pin(j).moduleId()];
            sum3 += exp_yi_gamma[placement.net(i).pin(j).moduleId()];
            sum4 += 1/exp_yi_gamma[placement.net(i).pin(j).moduleId()];
        }
        value_ += log(sum1);
        value_ += log(sum2);
        value_ += log(sum3);
        value_ += log(sum4);
    }
    value_ = gamma*value_;
    input_ = input;
    return value_;
}

const std::vector<Point2<double>> &Wirelength::Backward(){
    for(int i = 0; i < input_size; i++){
        grad_[0].x = 
    }
}






