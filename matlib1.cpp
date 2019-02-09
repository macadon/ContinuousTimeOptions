include<iostream>
#include<cmath>
#include<vector>
#include "matlib1.h"
#include<cstdlib>
#include<algorithm>
using namespace std;

void print(int a, int b) {
for(int i=a;i<=b;i++){
cout<<i<<endl;
}
}

int Fib(int N) {
if(N==0){
return 0;
}
else if(N==1){
return 1;
}
else if(N>1) {
return Fib(N-1)+Fib(N-2);
}

}




double normcdf(double x) {
double k=1/(1+0.2316419*x);
double PI=3.1415926538979;
double m=sqrt(2*PI);
if(x>=0) {
return
1-(1/m)*exp(-x*x/2)*k*(0.319381530+k*(-0.356563782+k*(1.781477937+k*(-1.821255978+1.330274429*k))));
}
else{
return 1-normcdf(-x);
}

}

int sum(int N) {
int total=0;
for(int i=1;i<N+1;i++) {
total+=i;
}
return total;
}

double hornerfunction0(double x, double a0) {
return a0;
}
double hornerfunction1(double x,double a0, double a1) {
return a0+x*hornerfunction0(x,a1);
}

double hornerfunction2(double x, double a0, double a1, double a2) {
return a0+x*hornerfunction1(x,a1,a2);
}

double hornerfunction3(double x, double a0, double a1, double a2, double a3){
return a0+x*hornerfunction2(x,a1,a2,a3);
}

double hornerfunction4(double x, double a0, double a1, double a2, double a3, double a4) {
return a0+x*hornerfunction3(x,a1,a2,a3,a4);
}
double hornerfunction5(double x, double a0, double a1, double a2, double a3, double a4, double a5) {
return a0+x*hornerfunction4(x,a1,a2,a3,a4,a5);
}

double hornerfunction6(double x, double a0, double a1, double a2,double a3, double a4, double a5, double a6){
return a0+x*hornerfunction5(x,a1,a2,a3,a4,a5,a6);
}

double hornerfunction7(double x, double a0, double a1, double a2, double a3, double a4, double a5, double a6, double a7) {
return a0+x*hornerfunction6(x,a1,a2,a3,a4,a5,a6,a7);
}

double hornerfunction8(double x, double a0, double a1, double a2, double a3, double a4, double a5, double a6, double a7, double a8) {
return a0+x*hornerfunction7(x,a1,a2,a3,a4,a5,a6,a7,a8);
}
double norminv(double x) {
double a0 = 2.50662823884;
 double a1 = -18.61500062529;
 double a2 = 41.39119773534;
 double a3 = -25.44106049637;
 double b1 = -8.47351093090;
 double b2 = 23.08336743743;
 double b3 = -21.06224101826;
 double b4 = 3.13082909833;
 double c0 = 0.3374754822726147;
 double c1 = 0.9761690190917186;
 double c2 = 0.1607979714918209;
 double c3 = 0.0276438810333863;
 double c4 = 0.0038405729373609;
double c5 = 0.0003951896511919;
 double c6 = 0.0000321767881768;
 double c7 = 0.0000002888167364;
 double c8 = 0.0000003960315187;
double r;
double y=x-0.5;
if(y<0.42 && y>-0.42) {
r=y*y;
return y*hornerfunction3(r,a0,a1,a2,a3)/hornerfunction4(r,1.0,b1,b2,b3,b4);
}
else{
if(y<0){
 r=x;
}
else {
 r=1-x;
}
double s=log(-log(r));
double t=hornerfunction8(s,c0,c1,c2,c3,c4,c5,c6,c7,c8);
if(x>0.5){
return t;
}
else{
return -t;
}

}

}
double blackScholesCallPrice(double K, double T, double S0,double sigma, double r) {
double d1=(1/(sigma*sqrt(T)))*(log(S0/K)+((r+(sigma*sigma/2))*sqrt(T)));
double d2=(1/(sigma*sqrt(T)))*(log(S0/K)+((r-(sigma*sigma/2))*sqrt(T)));
return normcdf(d1)*S0-normcdf(d2)*K*exp(-r*T);
}

int factorial(int N) {
int total;
if(N==0){
return 1;
}
else{
for(int i=1;i<=N;i++){
total*=i;
}
return total;
}
}
double rectanglerule(double a, double b, int n) {
double h=(b-a)/n;
double total;
for(int i=0;i<=n-1;i++) {
total+=pow(a+i*h+0.5*h,2);
}
return total/n;
}
double normpdf(double x) {
double total;
double a=0;
double b=1;
int N=1000;
double h=(b-a)/N;
for(int i=0;i<=N;i++){
double s=a+i*h+0.5*h;
double t=x+1-(1/s);
double f=pow(s,-2)*exp(-t*t/2);
total+=f;
}
return total/N;
}
int fib1(int n) {
int total=1;
int newtotal=1;
int Fib;
if(n==1||n==2){
return 1;
}
else{
for(int i=3;i<=n;i++) {
int temp=newtotal;
newtotal+=total;
total=temp;
}
return newtotal;
}




}

void print() {
for(int i=1;i<=10;i++) {
cout<<i<<endl;
}
}
void solveQuadratic(double a, double b, double c) {
if(b*b<4*a*c) {
cout<<"No real roots!" <<endl;
}
else if(b*b==4*a*c) {
cout<<"Exactly one real root given by: "<<-b/(2*a)<<endl;
}
else {
cout<<"Two real roots which are: "<<(-b+sqrt(b*b-4*a*c))/(2*a)<<"and "<<(-b-sqrt(b*b-4*a*c))/(2*a)<<endl;
}
}

double mean(vector<double> v) {
int n=100;
double total;
for(int i=0;i<n;i++){
total+=v[i];
}
return total/n;




}
double sd(vector<double> v){
double X=mean(v);
int n=100;
double sum;
for(int i=0;i<=n-1;i++){
sum+=(X-v[i])*(X-v[i]);
}
return sqrt(sum/(n-1));

}

double max(vector<double> v) {
int n=v.size();
double MAX=0.0;
for(int i=0;i<n;i++) {
if(MAX<=v[i]){
MAX=v[i];
}
else {
MAX+=0;
}



}
return MAX;

}

vector<double> randuniform(int n) {
vector<double> v;
for(int i=0;i<n;i++) {
double uniform=rand()/(RAND_MAX+1.0);
v.push_back(uniform);
}
return v;
}

vector<double> randnormal(int n) {
vector<double> v;
for(int i=0;i<n;i++) {
double uniform1=rand()/(RAND_MAX+1.0);
double uniform2=rand()/(RAND_MAX+1.0);
double n1=sqrt(-2*log(uniform1))*cos(2*PI*uniform2);
v.push_back(n1);
}
return v;
}
void sort(vector<double> &v){
sort(v.begin(),v.end());
}

CartesianPoint polarToCartesian(PolarPoint & p) {
CartesianPoint c;
c.x=p.r*cos(p.theta);
c.y=p.r*sin(p.theta);
return c;
}

double CartesianPoint::distanceTo(CartesianPoint &p1) {
return sqrt((p1.x-x)*(p1.x-x)+(p1.y-y)*(p1.y-y));


}

double Circle::area() {
        return PI*radius*radius;
        }

double Circle::circumference() {
        return PI*radius*2;
        }

double CallOption::payoff(std::vector<double>& stockPrices) {
double stockAtMaturity=stockPrices.back();
if(stockAtMaturity>getStrike()){
return stockAtMaturity-getStrike();
}
else{
return 0.0;
}
}
double CallOption::price(BlackScholesModel &bsm) {
double S=bsm.stockPrice;
double K=getStrike();
double sigma=bsm.volatility;
double r=bsm.riskFreeRate;
double T=getMaturity()-bsm.date;

double numerator = log(S/K)+(r+sigma*sigma*0.5)*T;
double denominator=sigma*sqrt(T);
double d1=numerator/denominator;
double d2=d1-denominator;
return S*normcdf(d1)-exp(-r*T)*K*normcdf(d2);
}


double DigitalCall::payoff(double stockAtMaturity) {
if(stockAtMaturity>strike) {
return 1;
}
else {
return 0;
}
}

double DigitalCall::getMaturity() {
return maturity;
}


double PutOption::payoff( std::vector<double>& stockPrices)  {
    double stockAtMaturity = stockPrices.back();
    if (stockAtMaturity<getStrike()) {
        return getStrike()-stockAtMaturity;
    } else {
        return 0.0;
    }
}
double PutOption::price(
         BlackScholesModel& bsm )  {
    double S = bsm.stockPrice;
    double K = getStrike();
    double sigma = bsm.volatility;
    double r = bsm.riskFreeRate;
    double T = getMaturity() - bsm.date;

    double numerator = log( S/K ) + ( r + sigma*sigma*0.5)*T;
    double denominator = sigma * sqrt(T );
    double d1 = numerator/denominator;
    double d2 = d1 - denominator;
    return -S*normcdf(-d1) + exp(-r*T)*K*normcdf(-d2);
}


double perimeter(CartesianPoint &p1, CartesianPoint &p2, CartesianPoint &p3){
return p1.distanceTo(p2)+p2.distanceTo(p3)+p3.distanceTo(p1);
}

double polynomial::evaluate(double x) {
double total;
for(int i=0;i<10;i++){
total+=v[i]*pow(x,i);
}
return total;
}

vector<double> BlackScholesModel::generatePricePath(double toDate, int nSteps, double drift) {
vector<double> path(nSteps,0.0);
vector<double> epsilon=randnormal(nSteps);
double dt=(toDate-date)/nSteps;
double a=(drift-volatility*volatility*0.5)*dt;
double b=volatility*sqrt(dt);
double currentLogS = log(stockPrice);
for(int i=0; i<nSteps;i++) {
        double dLogS=a+b*epsilon[i];
        double logS=currentLogS+dLogS;
        path[i]=exp(logS);
        currentLogS=logS;
        }
        return path;
}

vector<double> BlackScholesModel::generatePricePath(double toDate, int nSteps) {
return generatePricePath(toDate, nSteps, drift);
}
vector<double> BlackScholesModel::generateRiskNeutralPricePath(double toDate, int nSteps) {
return generatePricePath(toDate, nSteps, riskFreeRate);
}
BlackScholesModel::BlackScholesModel(){
drift=0.0;
stockPrice=0.0;
volatility=0.0;
riskFreeRate=0.0;
date=0.0;
}

MonteCarloPricer::MonteCarloPricer() :
        nScenarios(100000) {
}

double MonteCarloPricer::price(CallOption& callOption, BlackScholesModel& model) {
        double total=0.0;
        for(int i=0;i<nScenarios;i++){
        vector<double> path=model.generateRiskNeutralPricePath(callOption.getMaturity(),1);
        double stockPrice=path.back();
        double payoff=callOption.payoff(path);
        total+=payoff;
        }
        double mean=total/nScenarios;
        double r=model.riskFreeRate;
        double T=callOption.getMaturity()-model.date;
        return exp(-r*T)*mean;



}
double MonteCarloPricer::price(PutOption& putOption, BlackScholesModel& model) {
         double total=0.0;
        for(int i=0;i<nScenarios;i++) {
        vector<double> path=model.generateRiskNeutralPricePath(putOption.getMaturity(),1);
        double stockPrice=path.back();
        double payoff=putOption.payoff(path);
        total+=payoff;
}
        double mean=total/nScenarios;
        double r=model.riskFreeRate;
        double T=putOption.getMaturity()-model.date;
        return exp(-r*T)*mean;
}

double UpAndOutOption::payoff(double stockAtMaturity, BlackScholesModel& model) {
vector<double> path=model.generateRiskNeutralPricePath(maturity,1);
double total;
for(int i=0;i<path.size();i++) {
if(path[i]<barrier){
total+=1;
}
else{
total-=10000;
}
}
if(total>0 && stockAtMaturity<barrier) {
return stockAtMaturity-strike;
}
else{
return 0.0;
}

}
double MonteCarloPricer::price(PathIndependentOption& option, BlackScholesModel& model) {
         double total=0.0;
        for(int i=0;i<nScenarios;i++) {
        vector<double> path=model.generateRiskNeutralPricePath(option.getMaturity(),1);
        double stockPrice=path.back();
        double payoff=option.payoff(stockPrice);
        total+=payoff;
}
        double mean=total/nScenarios;
        double r=model.riskFreeRate;
        double T=option.getMaturity()-model.date;
        return exp(-r*T)*mean;
}

double MonteCarloPricer::Delta(CallOption& callOption,BlackScholesModel& model) {
double total1=0.0;
double total2=0.0;
double h=(model.stockPrice)/(1000000);
double T=callOption.getMaturity()-model.date;
for(int i=0;i<nScenarios;i++) {
        double uniform1=rand()/(RAND_MAX+1.0);
        double uniform2=rand()/(RAND_MAX+1.0);
        double n1=sqrt(-2*log(uniform1))*cos(2*PI*uniform2);
        double STOCK1=(model.stockPrice+h)*exp((model.riskFreeRate-model.volatility*model.volatility*0.5)*T+(model.volatility*n1*sqrt(T)));
        double STOCK2=(model.stockPrice-h)*exp((model.riskFreeRate-model.volatility*model.volatility*0.5)*T+(model.volatility*n1*sqrt(T)));
        vector<double> v1;
        vector<double> v2;
        for(int i=0;i<100;i++){
        v1.push_back(STOCK1);
        v2.push_back(STOCK2);
        }
double difference1=callOption.payoff(v1);
double difference2=callOption.payoff(v2);
total1+=difference1;
total2+=difference2;
}

double mean1=total1/nScenarios;
double mean2=total2/nScenarios;
double r=model.riskFreeRate;

return exp(-r*T)*(mean1-mean2)/(2*h);




}

double NormCDF::gaussian(double x) {
return (1/sqrt(2*PI))*exp(-x*x*0.5);
}
double NormCDF::integrategaussian(double a1,double b1, int N1) {
double h=(b1-a1)/N1;
double total=0.0;
for(int i=0;i<N1;i++) {
total+=gaussian(a1+i*h);
}
return h*total;
}

double SinFunction::evaluate(double x) {
        return sin(x);
}



double RealFunction::differentiatenumerically(double x) {
double h=0.0001;
return (evaluate(x+h)-evaluate(x-h))/(2*h);
}

int sumUsingPointer(int* toSum, int length) {
int sum=0;
for(int i=0;i<length;i++) {
sum+=toSum[i];
}
return sum;
}

double ContinuousTimeOptionBase::price(BlackScholesModel& model) {
MonteCarloPricer pricer;
return pricer.price(*this,model);
}
