#include "ObjectiveFunction.h"

#include "cstdio"
#include <cassert>

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

Wirelength::Wirelength(Placement &placement) : BaseFunction(placement.numModules()), placement(placement) {
    exp_xi_gamma.resize(placement.numModules());
    exp_yi_gamma.resize(placement.numModules());
    for(int i = 0; i < placement.numModules(); i++){
        exp_xi_gamma[i] = 0;
        exp_yi_gamma[i] = 0;
    }
    pin_netlist.clear();
    pin_netlist.resize(placement.numPins());
    for(int i = 0; i < placement.numNets(); i++){ // i is the net id
        for(int j = 0; j < placement.net(i).numPins(); j++){
            pin_netlist[placement.net(i).pin(j).pinId()].push_back(i);
        }
    }
}

const double &Wirelength::operator()(const std::vector<Point2<double>> &input){
    value_ = 0;
    exp_xi_gamma.clear();
    exp_yi_gamma.clear();
    exp_xi_gamma.resize(placement.numModules());
    exp_yi_gamma.resize(placement.numModules());
    for(int i = 0; i < placement.numModules(); i++){
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
    for(int i = 0; i < placement.numModules(); i++){
        grad_[i].x = 0;
        grad_[i].y = 0;
    }
    exp_xi_gamma.clear();
    exp_yi_gamma.clear();
    exp_xi_gamma.resize(placement.numModules());
    exp_yi_gamma.resize(placement.numModules());
    for(int i = 0; i < placement.numModules(); ++i){
        exp_xi_gamma[i] = exp(input_[i].x / gamma);
        exp_yi_gamma[i] = exp(input_[i].y / gamma);
    }
    for(int i = 0 ; i < placement.numModules(); i++){
        for(int j = 0; j < placement.module(i).numPins(); j++){
            int net_id = placement.module(i).pin(j).netId();
            double sum1 = 0,sum2 = 0;
            double sum3 = 0,sum4 = 0;
            for(int k = 0; k < placement.net(net_id).numPins(); k++){
                sum1 += exp_xi_gamma[placement.net(net_id).pin(k).moduleId()];
                sum2 += 1./exp_xi_gamma[placement.net(net_id).pin(k).moduleId()];
                sum3 += exp_yi_gamma[placement.net(net_id).pin(k).moduleId()];
                sum4 += 1./exp_yi_gamma[placement.net(net_id).pin(k).moduleId()];

            }
            
            
            grad_[i].x += (exp_xi_gamma[i]/sum1 -1/(sum2 * exp_xi_gamma[i]));
            grad_[i].y += (exp_yi_gamma[i]/sum3 - 1/(sum4 * exp_yi_gamma[i]));
            
        }
    }

    return grad_;
}

Density::Density(Placement &placement, double max_density, int num_bins_x, int num_bins_y,vector<Point2<double>> input) : BaseFunction(placement.numModules()), placement(placement), max_density(max_density), num_bins_x(num_bins_x), num_bins_y(num_bins_y), input_(input){
    bin_size_x = (placement.boundryRight() - placement.boundryLeft()) / num_bins_x;
    bin_size_y = (placement.boundryTop() - placement.boundryBottom()) / num_bins_y;
    density_map.resize(num_bins_x);
    center.resize(num_bins_x);
    for(int i = 0; i < num_bins_x; i++){
        density_map[i].resize(num_bins_y);
        center[i].resize(num_bins_y);
        for(int j = 0; j < num_bins_y; j++){
            density_map[i][j] = 0;
            center[i][j].first = placement.boundryLeft() + i * bin_size_x + bin_size_x / 2;
            center[i][j].second = placement.boundryBottom() + j * bin_size_y + bin_size_y / 2;
        }
    }
}
const double &Density::operator()(const std::vector<Point2<double>> &input){
    value_ = 0;
    for(int i = 0; i < num_bins_x; i++){
        for(int j = 0; j < num_bins_y; j++){
            density_map[i][j]= 0;
        }
    }
    for(int i = 0; i < num_bins_x; i++){
        for(int j = 0; j < num_bins_y; j++){
            for(int l = 0; l < placement.numModules(); l++){
                double Px = 0, Py = 0;
                double normalizer = 0;
                double center_x = input[l].x - center[i][j].first  ;
                double center_y = input[l].y - center[i][j].second;
                double alpha_x = 4/((placement.module(l).width()+2*bin_size_x)*(placement.module(l).width()+4*bin_size_x));
                double alpha_y = 4/((placement.module(l).height()+2*bin_size_y)*(placement.module(l).height()+4*bin_size_y));
                double beta_x = 2/((placement.module(l).width()+4*bin_size_x) * bin_size_x);
                double beta_y = 2/((placement.module(l).height()+4*bin_size_y) * bin_size_y);
                if(abs(center_x) <= placement.module(l).width()/2 + bin_size_x){
                    Px = 1 - alpha_x * center_x * center_x;
                }
                else if(abs(center_x) <= placement.module(l).width()/2 + 2*bin_size_x){
                    Px = beta_x * (center_x - placement.module(l).width()/2 - 2*bin_size_x) * (center_x - placement.module(l).width()/2 - 2*bin_size_x);
                }
                else{
                    Px = 0;
                }
                if(abs(center_y) <= placement.module(l).height()/2 + bin_size_y){
                    Py = 1 - alpha_y * center_y * center_y;
                }
                else if(abs(center_y) <= placement.module(l).height()/2 + 2*bin_size_y){
                    Py = beta_y * (center_y - placement.module(l).height()/2 - 2*bin_size_y) * (center_y - placement.module(l).height()/2 - 2*bin_size_y);
                }
                else{
                    Py = 0;
                }
                normalizer = placement.module(l).width() * placement.module(l).height() / (Px * Py);
                density_map[i][j] += Px * Py * normalizer;
            }
        }
    }
    for(int i = 0; i < bin_size_x; i++){
        for(int j = 0; j < bin_size_y; j++){
            value_ += (density_map[i][j] - max_density) * (density_map[i][j] - max_density);
        }
    }
    input_ = input;
    return value_;
}

const std::vector<Point2<double>> &Density::Backward(){
    for(int i = 0; i < placement.numModules(); i++){
        grad_[i].x = 0;
        grad_[i].y = 0;
    }
    density_diff_map.resize(num_bins_x);
    for(int i = 0; i < num_bins_x; i++){
        density_diff_map[i].resize(num_bins_y);
        for(int j = 0; j < num_bins_y; j++){
            density_diff_map[i][j] = 0;
        }
    }
    grad_.resize(placement.numModules());
    for(int i = 0 ;i < placement.numModules(); ++i){
        for(int j = 0; j < num_bins_x; j++){
            for(int k = 0; k < num_bins_y; k++){
                double cen2cen_x = input_[i].x - center[j][k].first;
                double alpha = 4/((placement.module(i).width()+2*bin_size_x)*(placement.module(i).width()+4*bin_size_x));
                double beta = 2/((placement.module(i).width()+4*bin_size_x) * bin_size_x);
                double Px_diff = 0;
                double cen2cen_y = input_[i].y - center[j][k].second;
                double alpha_y = 4/((placement.module(i).height()+2*bin_size_y)*(placement.module(i).height()+4*bin_size_y));
                double beta_y = 2/((placement.module(i).height()+4*bin_size_y) * bin_size_y);
                double Py_diff = 0;
                double Px = 0, Py = 0;
                if(abs(cen2cen_x) >= placement.module(i).width()/2 + 2*bin_size_x){
                    Px_diff = 0;
                    Px = 0;
                }
                else if(abs(cen2cen_x) >= placement.module(i).width()/2 + bin_size_x){
                    Px_diff = 2 * beta * (cen2cen_x - placement.module(i).width()/2 - 2*bin_size_x);
                    Px = beta * (cen2cen_x - placement.module(i).width()/2 - 2*bin_size_x) * (cen2cen_x - placement.module(i).width()/2 - 2*bin_size_x);
                }
                else{
                    Px_diff = -2 * alpha * cen2cen_x;
                    Px = 1 - alpha * cen2cen_x * cen2cen_x;
                }
                if(abs(cen2cen_y) >= placement.module(i).height()/2 + 2*bin_size_y){
                    Py_diff = 0;
                    Py = 0;
                }
                else if(abs(cen2cen_y) >= placement.module(i).height()/2 + bin_size_y){
                    Py_diff = 2 * beta_y * (cen2cen_y - placement.module(i).height()/2 - 2*bin_size_y);
                    Py = beta_y * (cen2cen_y - placement.module(i).height()/2 - 2*bin_size_y) * (cen2cen_y - placement.module(i).height()/2 - 2*bin_size_y);
                }
                else{
                    Py_diff = -2 * alpha_y * cen2cen_y;
                    Py = 1 - alpha_y * cen2cen_y * cen2cen_y;
                }
                double normalizer = placement.module(i).width() * placement.module(i).height() / (Px * Py);
                grad_[i].x += 2 * (density_map[j][k] - max_density) * Px_diff * Py * normalizer;
                grad_[i].y += 2 * (density_map[j][k] - max_density) * Px * Py_diff * normalizer;
                
            }
        }
    }
    return grad_;
}
ObjectiveFunction::ObjectiveFunction(Placement &placement, double lambda)
    :  BaseFunction(placement.numModules()),placement_(placement), lambda_(lambda), wirelength_(placement), 
      density_(placement, 1, 10, 10, std::vector<Point2<double>>(placement.numModules())) {
    grad_.resize(placement.numModules());
    for(int i = 0; i < placement.numModules(); i++){
        grad_[i].x = 0;
        grad_[i].y = 0;
    }
}

const double& ObjectiveFunction::operator()(const std::vector<Point2<double>> &input){
    double value = 0;
    value += wirelength_(input);
    value += lambda_ * density_(input);
    return value;
}
const std::vector<Point2<double>>& ObjectiveFunction::Backward(){
    for(int i = 0; i < placement_.numModules(); i++){
        grad_[i].x = 0;
        grad_[i].y = 0;
    }
    for(int i = 0; i < placement_.numModules(); i++){
        grad_[i].x += wirelength_.Backward()[i].x;
        grad_[i].y += wirelength_.Backward()[i].y;
        grad_[i].x += lambda_ * density_.Backward()[i].x;
        grad_[i].y += lambda_ * density_.Backward()[i].y;
    }
    return grad_;
}



