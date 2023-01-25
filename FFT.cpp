#include<bits/stdc++.h>
#include<iostream>
using namespace std;

/*--------------Golbal variables and vectors----------*/
int N;
vector<double> sum(2);
vector<double> sub(2);
vector<double> mult(2);
vector<double> poly_1;
vector<double> poly_2;
vector<vector<double>> roots;

/*--------------declarations or prototypes of all the fuctions which are used in program----------*/
void find_complex_roots(int N);
int find_N(int deg_1,int deg_2);
void I_find_complex_roots(int N);
void print_polynomial(vector<double> temp, int degree);
vector<vector<double>> Eval(vector<double> poly,int N);
vector<vector<double>> I_Eval(vector<double> poly,int N);
vector<double> add_complex_num(double a,double b,double c,double d);
vector<double> sub_complex_num(double a,double b,double c,double d);
vector<double> multiply_complex_num(double a,double b,double c,double d);
vector<double> naive_polynomial_multiplication(vector<double> &poly1,vector<double> &poly2);


//sum of the complex numbers a+ib and c+id
vector<double> add_complex_num(double a,double b,double c,double d){
    sum[0]=a+c;
    sum[1]=b+d;
    return sum;
}

//subtraction of the complex numbers a+ib and c+id
vector<double> sub_complex_num(double a,double b,double c,double d){
    sub[0]=a-c;
    sub[1]=b-d;
    return sub;
    
}

//product of the complex numbers a+ib and c+id
vector<double> multiply_complex_num(double a,double b,double c,double d){
    mult[0]=((a*c)-(b*d));
    mult[1]=((a*d)+(b*c));
    return mult;
    
}

//polynomial_multiplication using naive method which takes O(N^2)
vector<double> naive_polynomial_multiplication(vector<double> &poly1,vector<double> &poly2){
    vector<double> naive_prod((poly1.size())+(poly2.size())-1);
    for(int i=0;i<poly1.size();i++)
    {
        for(int j=0;j<poly2.size();j++){
            naive_prod[i+j]=naive_prod[i+j]+poly1[i]*poly2[j];
        }
    }
    return naive_prod;
}

//fuction to print polynomial
void print_polynomial(vector<double> temp, int degree){
    for(auto it=temp.rbegin();it!=temp.rend();it++){
        if(it==temp.rbegin())
        {
            cout<<"("<<*it<<")"<<"x^"<<degree;
            degree--;
        }
        else if((*it)==0){
            degree--;
            continue;
        }
        else if(it==temp.rend()-1){
           cout<<*it; 
           break;
        }
        else{
        cout<<"+"<<"("<<*it<<")"<<"x^"<<degree;
        degree--;
        }
    }
    cout<<endl;
}

//computeing the smallest power of 2 greater than or equal to deg1 + deg2 + 1 in O(1)
int find_N(int deg_1,int deg_2){
      int var=0;
      int count=0;
      int temp_1=(deg_1+deg_2+1);
      int temp_2=temp_1-1;
      int temp_3= !((temp_1)&(temp_2));
      //if deg_1+deg_2+1 is power of two
      if((temp_1) && (temp_3))
        return temp_1;      
      //deg_1+deg_2+1 is NOT power of two
      else{                  
        while(temp_1!=0){
            temp_1=temp_1>>1;
            count++;
        }
        var=(1<<count);
      }
    return(var);//returing power of 2 which is greater than deg_1+deg_2+1
}

//computes the Nth root of unity
void find_complex_roots(int N)
{
    roots.resize(N,vector<double>(2));
    double angle = (M_PI*2)/N;
    for(int i=0; i<N; i++)
    {
        // calculate the real and imaginary part of root
        double real_part = cos(i*angle);
        double imginary_part = sin(i*angle);
        roots[i][0]=real_part;
        roots[i][1]=imginary_part;
    }
}

//This fuction is used in inverse FFT while converting point peresentation back to cofficient peresentation
void I_find_complex_roots(int N)
{
    roots.resize(N,vector<double>(2));
    double angle = (M_PI*2)/N;
    for(int i=0; i<N; i++)
    {
        double real_part = cos(i*angle);
        double imginary_part = (-1)*sin(i*angle);
        roots[i][0]=real_part;
        roots[i][1]=imginary_part;
    }
}

//evaluation part
//This is FFT fuction that convert cofficient form to point-value form, all value is calculated at Nth root of unity
vector<vector<double>> Eval(vector<double> poly,int N){ 
    //base
    if(N==1){
        vector<vector<double>> temp( 1 , vector<double> (1,poly[0]));
        return temp;
    }
    vector<double> even;
    vector<double> odd;
    for (int i = 0; i < N; i++) {  //change below also
        if(i%2==0){
        even.push_back(poly[i]);
        }
        else{
        odd.push_back(poly[i]);
        }
    }
   
    vector<vector<double>> yE=Eval(even,N/2);
    vector<vector<double>> yO=Eval(odd,N/2);
    find_complex_roots(N);
    
    vector<vector<double>> poly_evaluations( N , vector<double> (2));

    for (int j = 0; j <N / 2; j++) {
        vector<double> mult_temp= multiply_complex_num(roots[j][0],roots[j][1],yO[j][0],yO[j][1]);
        vector<double> sum_temp= add_complex_num(yE[j][0],yE[j][1],mult_temp[0],mult_temp[1]);
        poly_evaluations[j][0]=sum_temp[0];
        poly_evaluations[j][1]=sum_temp[1];
        
        vector<double> sub_temp = sub_complex_num(yE[j][0],yE[j][1],mult_temp[0],mult_temp[1]);
        poly_evaluations[j+(N/2)][0]=sub_temp[0]; 
        poly_evaluations[j+(N/2)][1]=sub_temp[1];
    }
    return poly_evaluations;
}

//interpolation part
//This is inverse FFT fuction that convert point-vaule form to cofficient form
vector<vector<double>> I_Eval(vector<vector<double>> temp_poly_product,int N){
    if(N==1){
        vector<vector<double>> temp( 1 , vector<double> (2));
        temp[0][0] = temp_poly_product[0][0];
        temp[0][1] = temp_poly_product[0][1];
        return temp;
    }

    vector<vector<double>> even;
    vector<vector<double>> odd;
    
    for (int i = 0; i < N; i++) {
        if(i%2==0){
        even.push_back(temp_poly_product[i]);
        }
        else{
        odd.push_back(temp_poly_product[i]);
        }
    }

    vector<vector<double>> yE=I_Eval(even,N/2);
    vector<vector<double>> yO=I_Eval(odd,N/2);
    I_find_complex_roots(N);
    vector<vector<double>> I_poly_evaluations( N , vector<double> (2));
    
    for (int j = 0; j < N / 2; j++) {

        vector<double> mult_temp= multiply_complex_num(roots[j][0],roots[j][1],yO[j][0],yO[j][1]);
        vector<double> sum_temp= add_complex_num(yE[j][0],yE[j][1],mult_temp[0],mult_temp[1]);
        I_poly_evaluations[j][0]=sum_temp[0];
        I_poly_evaluations[j][1]=sum_temp[1];
        
        vector<double> sub_temp = sub_complex_num(yE[j][0],yE[j][1],mult_temp[0],mult_temp[1]);
        I_poly_evaluations[j+(N/2)][0]=sub_temp[0]; 
        I_poly_evaluations[j+(N/2)][1]=sub_temp[1];
    }
    return I_poly_evaluations;
}

//This function performs the product of 2 polynomial by multipying their point-value form in O(n)
vector<vector<double>> product_polynomial_evaluations(vector<vector<double>> poly_1_evaluations, vector<vector<double>> poly_2_evaluations){
    vector<vector<double>> temp(N, vector<double>(2));
    for(int i=0; i<N; i++){
        vector<double> mul = multiply_complex_num(poly_1_evaluations[i][0], poly_1_evaluations[i][1], poly_2_evaluations[i][0], poly_2_evaluations[i][1]);
        temp[i][0]=mul[0];
        temp[i][1]=mul[1];
    }
    return temp;
}


int main()
{

  int temp;
  int deg_1;                                  //this variable will store the degree of polynomial 1
  int deg_2;                                  //this variable will store the degree of polynomial 2
  vector<vector<double>> poly_1_evaluations;  //this vector will store point-value representation of polynomial 1 at Nth root of unity
  vector<vector<double>> poly_2_evaluations;  //this vector will store point-value representation of polynomial 2 at Nth root of unity
  
  /*--------------------polynomial 1---------------------------------*/
  cout<<"enter the degree of first polynomial:";   //taking the 1st polynomial Degree as input
  cin>>deg_1;
  if(deg_1>15){
   cout<<"polynomial degree cannot be greater than 15";
   exit(0);
  }
  cout<<"Enter the "<<deg_1+1<<" coefficients of the 1st polynomial in the increasing order of the degree of the monomials:"<<endl; //taking the 1st polynomial coeffiencts as input and storing them in poly_1 vector
  for(int i=1;i<=deg_1+1;i++){
       cin>>temp;
       poly_1.push_back(temp);
    }
  print_polynomial(poly_1,deg_1);              //printing the polynomial 1
  

  /*------------------polynomial 2-------------------------------------*/
  cout<<endl<<"enter the degree of second polynomial:";  // taking the 2nd polynomial Degree as input
  cin>>deg_2;
  if(deg_2>15){
   cout<<"polynomial degree cannot be greater than 15";
   exit(0);
  }
  cout<<"Enter the "<<deg_2+1<<" coefficients of the 2nd polynomial in the increasing order of the degree of the monomials:"<<endl; //taking the 2nd polynomial coeffiencts as input and storing them in poly_2 vector
  for(int i=1;i<=deg_2+1;i++){
       cin>>temp;
       poly_2.push_back(temp);
       }
  print_polynomial(poly_2,deg_2);               //printing the polynomial 2
  
  
  /*------------------multipiying  two polynomial NAIVE METHOD O(N^2) --*/
  vector<double> naive_prod=naive_polynomial_multiplication(poly_1,poly_2);
  cout<<endl;
  cout<<"The product of the two polynomials obtained via naive polynomial multiplication is:"<<endl;
  print_polynomial(naive_prod,deg_1+deg_2);              //printing resultant polynomial which got multiplied using naive method

  
  
  /*------------------find_N computes the smallest power of 2 greater than or equal to deg1 + deg2 + 1.---------------*/
  N=find_N(deg_1,deg_2);
  
  //cout<<"value of N is:"<<N<<endl; 
  poly_1.resize(N);
  poly_2.resize(N);
  /*------------------padding zeroes at the end of Poly-1 and poly-2 according the N-----------------------------------*/
  for(int i=poly_1.size();i<N;i++)
  {
    poly_1.push_back(0);
  }
  for(int i=poly_2.size();i<N;i++)
  {
    poly_2.push_back(0);
  }


  /*-------------------performing resize operation on all related vectors-----------------------------------------------*/
  roots.resize(N, vector<double>(2,0));
  poly_1_evaluations.resize(N,vector<double>(2,0));
  poly_2_evaluations.resize(N,vector<double>(2,0));

   /*--------- Converting coefficient representation to point-value representation using divide and conquer FFT---------------------------------------*/
  /*---------Calculating point-value representation of polynomial 1 and 2 at Nth root of unity---------------*/
  poly_1_evaluations=Eval(poly_1,N);
  poly_2_evaluations=Eval(poly_2,N);
  
  /*-----------performing product on point-value representation of polynomial 1 and polynomial 2-----------*/
  vector<vector<double>> prod_poly_evaluations = product_polynomial_evaluations(poly_1_evaluations, poly_2_evaluations);
  
  /*------converting point-vaule representation of product polynomial to coefficient form-------------------- */
  vector<vector<double>> point_to_coefficient = I_Eval(prod_poly_evaluations, N);

  vector<double> result(deg_1+deg_2+1);
    for(int i=0; i<=deg_1+deg_2; i++){
        result[i]=round(point_to_coefficient[i][0]/N);
    }
    cout<<endl; 
    cout<<"The product of the two polynomials obtained via polynomial multiplication using FFT is:"<<endl;
    print_polynomial(result, deg_1+deg_2);
}