#include <bits/stdc++.h>
using namespace std;

vector<string> words;
vector<vector<double>> dp;
vector<vector<int>> roots;
double OBST(vector<double> probabilities);
double sumprob(vector<double> probabilities, int i, int j);

//function to get sum of array elements probabilities[i] to probabilities[j] which we will use while calculating cost
double sumprob(vector<double> probabilities, int i, int j) {
        double sum = 0;
        for (int k = i; i < j; i++)
            sum += probabilities[i];
        return sum;
    }

//fuction that performs preoder travesal on roots array and also prints preorder travesal
void Preorder(vector<vector<int>> roots,vector<string> words, int i, int j) {
        if (i == j)
            return;
        int temp = roots[i][j];
        cout<<words[temp - 1]<< " ";
        Preorder(roots, words, i, temp - 1);
        Preorder(roots, words, temp, j);
    }


// A recursive function to calculate cost of optimal binary search tree using dynamic programming.
double OBST(vector<double> probabilities,int n) {
        
        //dp[i][i + 1] considers all chains of length = 1 
        //so here we just have to assign probability directly.
        for (int i = 0; i < n; i++)
            dp[i][i + 1] = probabilities[i];
        

        // One by one consider all elements as root and recursively find cost
        // of the BST, compare the cost with min cost and update min if needed,
        //dp[i][j] considers all chains of length = j - i One by one consider all elements
        // as root and recursively find cost of the BST, compare the cost with min and update min if needed
        for (int len = 2; len <= n + 1; len++) {
            for (int start = 0; start <= n - len + 1; start++) {

                int end = start + len - 1;
                dp[start][end] = DBL_MAX;

                //The DP formula c[i][j] = min(i < k <= j){c[i][k-1] + c[k][j] + wt[i,j]}
                for (int k = start + 1; k <= end; k++) {
                    double temp = dp[start][k - 1] + dp[k][end] + sumprob(probabilities, start, end);
                    if (temp < dp[start][end]) {
                        dp[start][end] = temp;
                        roots[start][end] = k;
                    }
                }
            }
        }
        return dp[0][n]; //returning the minimum expected total access time
    }


int main()
{
    int n;                                     //Number of Strings that you want to insert
    int flag=0;                                //This boolean variable will help us to identify all probabilities are distinct or not                         
    double sum=0;                              //This variable will store the sum of all probabilities to check is it is equal to one or not
    bool is_sort;                              //This variable will see if string entered are sorted dictionary order or not
    string temp_string;                        
    double temp_double;
    vector<double> probabilities;              //storing frequency of strings
    
    
    cout<<"How many strings do you want to insert in the BST?";
    cin>>n;
    cout<<endl;
    
    //resizing vectors according to n
    words.resize(n);
    probabilities.resize(n);

    //taking the strings and their probabilities as input from user
    cout << "Enter " <<n<< " strings in sorted dictionary order along with their probabilities:";
    for(int i=0;i<n;i++){
        cin>>temp_string;
        words[i]=temp_string;
        cin>>temp_double;
        sum=sum+temp_double;
        probabilities[i]=(temp_double);

    }

    //checking if strings are in sorted dictionary order or not?
    is_sort=is_sorted(begin(words),end(words));
    if(is_sort==0){
        cout<<"The strings entered are not in sorted order.";
    }

    //checking if any probabilities are not distinct?
    for (int i = 1; i < probabilities.size(); ++ i) {
            for (int j = 0; j < i; ++ j) {
                if (probabilities[i] == probabilities[j]) {
                    flag=1;
                }
            }
    }
    if(flag==1){
        cout<<"probabilities are not distinct.";
    }

    //checking if The probabilities add up to 1?
    if(int(sum)!=1){
        cout<<"The probabilities donnot add up to 1.";
        return 0;
    }
    
    //return 0 is there is any error.
    if(is_sort==0 || flag==1)
        return 0;

    //resizing 2D vector dp(stores the cost) and roots(stores the root node)
    dp.resize(n+1, vector<double>(n+1));
    roots.resize(n+1, vector<int>(n+1));
    
    //mincost is the minimum expected total access time
    double mincost=OBST(probabilities,probabilities.size());
    cout<<"minimum expected total access time is ";
    cout<<mincost;
    cout<<endl;

    cout<<"Preorder traversal of the BST that provides minimum expected total access time is:"<<endl;
    //performing preorder traversal
    Preorder(roots, words, 0, n);
    return 0;
}