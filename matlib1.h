#pragma once
#include <vector>
#include <cstdlib>
#include<vector>
using namespace std;
const double PI=3.14159265358979;
int Fib(int N);
int sum(int N);
void print(int a, int b);
double normcdf(double x);
double norminv(double x);
double blackScholesCallPrice(double K, double T, double S0,double sigma, double r);
double rectanglerule(double a, double b, int n);
int factorial(int N);
int fib1(int n);
double normpdf(double x);
void print();
void solveQuadratic(double a, double b, double c);
double mean(std::vector<double> v);
double sd(std::vector<double> v);
double max(std::vector<double> v);
std::vector<double> randuniform(int n);
std::vector<double> randnormal(int n);
void sort(std::vector<double> &v);
class CartesianPoint {
public:
        double x;
        double y;
        double distanceTo(CartesianPoint &p1);
};

class PolarPoint {
public:
        double r;
        double theta;
};


class Circle {
public:
        double radius;
        double area();
        double circumference();
};

class BlackScholesModel {
public:
        BlackScholesModel();
        double drift;

        double stockPrice;

        double volatility;
        double riskFreeRate;
        double date;
        std::vector<double> generatePricePath(double toDate, int nSteps);
        std::vector<double> generateRiskNeutralPricePath(double toDate, int nSteps);
private:
        std::vector<double> generatePricePath(double toDate,int nSteps, double drift);
};
class PathIndependentOption {
public:
        virtual ~PathIndependentOption() {}

        virtual double payoff( double finalStockPrice) =0;

        virtual double getMaturity() =0;
};



class DigitalCall: public PathIndependentOption {
public:
        double strike;
        double maturity;
        double payoff(double stockAtMaturity);
        double getMaturity();
};


double perimeter(CartesianPoint &p1, CartesianPoint &p2, CartesianPoint &p3);

class polynomial {
public:
        vector<double> v;
        double evaluate(double x);
};

class UpAndOutOption {
public:
        double maturity;
        double strike;
        double barrier;
        double payoff(double stockAtMaturity,BlackScholesModel& model);
        };



class NormCDF {
public:
        double a;
        double b;
        int N;
        double gaussian(double x);
        double integrategaussian(double a1,double b1,int N1);
        };

class RealFunction {
public:
        virtual ~RealFunction() {};
        virtual double evaluate (double x) = 0;
        double differentiatenumerically(double x);
};

class SinFunction : public RealFunction {
        double evaluate( double x);
};

class Pair{
public:
        double x;
        double y;
        Pair();
        Pair(double x, double y);
};

int sumUsingPointer(int* toSum, int length);

class ContinuousTimeOption {
public:
        virtual ~ContinuousTimeOption() {};
        virtual double getMaturity() =0;
        virtual double payoff(std::vector<double>& stockPrices)=0;
        virtual bool isPathDependent()=0;
};

class ContinuousTimeOptionBase: public ContinuousTimeOption {
public:
        virtual ~ContinuousTimeOptionBase() {};
        double getMaturity(){
        return maturity;
        }

        void setMaturity(double maturity) {
        this->maturity=maturity;
        }
        double getStrike() {
        return strike;
        }

        void setStrike(double strike){
        this->strike=strike;
        }
        double price(BlackScholesModel& model);
private:
        double maturity;
        double strike;
};
class PutOption: public ContinuousTimeOptionBase {
public:

        double payoff(vector<double>& stockPrices);
        virtual double price(BlackScholesModel& bsm);

        bool isPathDependent() {
        return false;
        }
};


class CallOption: public ContinuousTimeOptionBase {
public:

        double payoff(vector<double>& stockPrices);
        double price(BlackScholesModel& bsm);

        bool isPathDependent() {
        return false;
        }
};

class MonteCarloPricer {
public:
        MonteCarloPricer();
        int nScenarios;
        double price(CallOption& option, BlackScholesModel& model);
        double price(PutOption& putOption, BlackScholesModel& model);
        double price(PathIndependentOption& option, BlackScholesModel& model);
        double Delta(CallOption&callOption,BlackScholesModel& model);
};
